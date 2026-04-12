#!/usr/bin/env python3
# =============================================================
# E. COLI — SCRIPT 03: K-MER EXTRACTION (tBLASTn)
# Input:  blast_db/query_genome
#         $CARD_PROTEIN_FILE
#         $TEMPLATES_DIR/features_*.txt  (8 per-antibiotic templates)
# Output: kmer_production.csv
#         logs/03_extract_kmers.log
# NOTE:   E. coli has 8 separate per-antibiotic feature templates.
#         We take the UNION of k-mers across all templates so a
#         single tBLASTn pass covers every model's k-mer needs.
# =============================================================

import os
import re
import sys
import subprocess
import time
from collections import Counter
from datetime import datetime
from pathlib import Path

import pandas as pd
from Bio import SeqIO

# ── Paths ──────────────────────────────────────────────────────
WORK_DIR      = Path(os.environ.get("WORK_DIR",     "."))
TEMPLATES_DIR = Path(os.environ.get("TEMPLATES_DIR", str(WORK_DIR)))
CARD_FILE     = Path(os.environ.get("CARD_PROTEIN_FILE", WORK_DIR / "card_all_proteins.fasta"))
BLAST_DB      = WORK_DIR / "blast_db" / "query_genome"
LOG_FILE      = WORK_DIR / "logs" / "03_extract_kmers.log"

os.chdir(WORK_DIR)

K_SIZE              = 10
E_VALUE             = 1e-5
MIN_IDENTITY        = 80
MIN_LENGTH          = 50
MAX_PROTEINS_PER_GENE = 3

# ── Logger ─────────────────────────────────────────────────────
def log(msg: str) -> None:
    line = f"[{datetime.now():%H:%M:%S}] {msg}"
    print(line)
    with open(LOG_FILE, "a") as fh:
        fh.write(line + "\n")

(WORK_DIR / "logs").mkdir(parents=True, exist_ok=True)
with open(LOG_FILE, "w") as fh:
    fh.write("=" * 60 + "\n")
    fh.write("E. COLI SCRIPT 03: K-MER EXTRACTION\n")
    fh.write(f"Started: {datetime.now()}\n")
    fh.write("=" * 60 + "\n")

# ── Validate required files ────────────────────────────────────
for required in [CARD_FILE, Path(f"{BLAST_DB}.nhr")]:
    if not required.exists():
        log(f"✗ ERROR: Required file not found: {required}")
        sys.exit(1)

# ── 1. Collect training k-mers from ALL 8 per-antibiotic templates
log("Loading training k-mer feature list from all antibiotic templates...")
training_kmers: set = set()

template_files = sorted(TEMPLATES_DIR.glob("features_*.txt"))
if not template_files:
    log(f"⚠ WARNING: No feature templates found in {TEMPLATES_DIR}")

for tf in template_files:
    before = len(training_kmers)
    with open(tf) as fh:
        for line in fh:
            feat = line.strip()
            # K-mers: exactly 10 uppercase alphabetic characters, no digits
            if len(feat) == 10 and feat.isalpha() and feat.isupper():
                training_kmers.add(feat)
    added = len(training_kmers) - before
    log(f"  ✓ {tf.name}  (+{added} k-mers)")

log(f"\n  Total unique k-mers across all templates: {len(training_kmers):,}")

if len(training_kmers) == 0:
    log("✗ ERROR: No k-mers found in templates.")
    sys.exit(1)

# ── 2. Get resistance gene list from gene presence file ────────
log("\nLoading resistance gene list from gene_presence_production.csv...")
try:
    gene_df = pd.read_csv(WORK_DIR / "gene_presence_production.csv")
    resistance_genes = [c for c in gene_df.columns if c != "Genome_ID"]
    log(f"  Resistance genes: {len(resistance_genes)}")
except Exception as exc:
    log(f"✗ ERROR loading gene presence file: {exc}")
    sys.exit(1)

# ── 3. Filter CARD proteins to resistance genes ─────────────────
log("\nFiltering CARD proteins to resistance genes...")

def normalise_gene(gene: str) -> str:
    name = re.sub(r"^(bla|BLA)", "", gene)
    name = re.sub(r"_\d+$", "", name)
    return name.upper()

gene_to_proteins: dict = {g: [] for g in resistance_genes}
filtered_proteins = []
matched_genes: set = set()

for record in SeqIO.parse(str(CARD_FILE), "fasta"):
    desc = record.description
    for gene in resistance_genes:
        gene_norm = normalise_gene(gene)
        patterns = [
            rf'\b{re.escape(gene)}\b',
            rf'\b{re.escape(gene_norm)}\b',
        ]
        if "(" in gene:
            clean = gene.replace("(", "").replace(")", "").replace("'", "")
            patterns.append(rf'\b{re.escape(clean)}\b')
        for pat in patterns:
            if re.search(pat, desc, re.IGNORECASE):
                if len(gene_to_proteins[gene]) < MAX_PROTEINS_PER_GENE:
                    filtered_proteins.append(record)
                    matched_genes.add(gene)
                    gene_to_proteins[gene].append(record.id)
                break

filtered_fasta = WORK_DIR / "resistance_proteins_production.faa"
SeqIO.write(filtered_proteins, str(filtered_fasta), "fasta")
log(f"  Genes matched  : {len(matched_genes)} / {len(resistance_genes)}")
log(f"  Proteins written: {len(filtered_proteins)}")

if len(filtered_proteins) == 0:
    log("✗ ERROR: No proteins matched resistance genes. Check CARD file.")
    sys.exit(1)

# ── 4. Run tBLASTn ─────────────────────────────────────────────
log("\nRunning tBLASTn (may take several minutes)...")
blast_out = WORK_DIR / "blast_production.tsv"
t0 = time.time()

blast_cmd = [
    "tblastn",
    "-query",           str(filtered_fasta),
    "-db",              str(BLAST_DB),
    "-out",             str(blast_out),
    "-outfmt",          "6 qseqid sseqid pident length qseq sseq",
    "-evalue",          str(E_VALUE),
    "-max_target_seqs", "5",
    "-num_threads",     "4",
]

result = subprocess.run(blast_cmd, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
elapsed = int(time.time() - t0)

if result.returncode != 0:
    log(f"✗ BLAST ERROR: {result.stderr.decode(errors='replace')}")
    sys.exit(1)
log(f"  tBLASTn completed in {elapsed}s")

# ── 5. Handle empty output ─────────────────────────────────────
if not blast_out.exists() or blast_out.stat().st_size == 0:
    log("⚠ WARNING: No BLAST hits — all k-mers will be 0")
    pd.DataFrame([{"Genome_ID": "query_genome"}]).to_csv(
        WORK_DIR / "kmer_production.csv", index=False
    )
    log("  Saved empty kmer_production.csv")
    sys.exit(0)

# ── 6. Parse BLAST + extract k-mers ───────────────────────────
blast_df = pd.read_csv(
    blast_out, sep="\t",
    names=["query_id", "subject_id", "pident", "length", "query_seq", "subject_seq"]
)
blast_df = blast_df[
    (blast_df["pident"] >= MIN_IDENTITY) &
    (blast_df["length"] >= MIN_LENGTH)
]
log(f"  Passing BLAST hits: {len(blast_df):,}")

log("Extracting k-mers from aligned sequences...")
found_kmers: Counter = Counter()
for _, hit in blast_df.iterrows():
    seq = hit["subject_seq"].replace("-", "")
    for i in range(len(seq) - K_SIZE + 1):
        kmer = seq[i: i + K_SIZE]
        if "*" not in kmer and "X" not in kmer:
            found_kmers[kmer] += 1

log(f"  Unique k-mers extracted        : {len(found_kmers):,}")

matched_kmers = {k: v for k, v in found_kmers.items() if k in training_kmers}
log(f"  K-mers matching training set   : {len(matched_kmers):,} / {len(training_kmers):,}")
missing_pct = (1 - len(matched_kmers) / max(len(training_kmers), 1)) * 100
log(f"  Missing                        : {missing_pct:.1f}%  (filled with 0 in alignment)")

# ── 7. Save ───────────────────────────────────────────────────
kmer_row = {"Genome_ID": "query_genome"}
kmer_row.update(matched_kmers)
pd.DataFrame([kmer_row]).to_csv(WORK_DIR / "kmer_production.csv", index=False)

# Cleanup
blast_out.unlink(missing_ok=True)
filtered_fasta.unlink(missing_ok=True)

log(f"\n✓ Saved: kmer_production.csv")
log("  Next step: bash 04_extract_snps.sh")
log("=" * 60)
log(f"SCRIPT 03 COMPLETE: {datetime.now()}")
log("=" * 60)
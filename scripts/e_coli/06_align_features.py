#!/usr/bin/env python3
# =============================================================
# E. COLI — SCRIPT 06: FEATURE ALIGNMENT (PER ANTIBIOTIC)
# Input:  gene_presence_production.csv
#         kmer_production.csv
#         snp_production.csv
#         $TEMPLATES_DIR/features_ecoli.txt  (single shared template)
# Output: aligned_full.csv  — the ONE file all 8 E.COLI models use
#         logs/06_align_features.log
# NOTE:   E.COLI uses a SINGLE shared feature template.
#         Missing features are filled with 0 (standard practice for
#         genomic data where ~65% of training features are zero)
# =============================================================

import os
import sys
import pandas as pd
from datetime import datetime
from pathlib import Path

# ── Paths ──────────────────────────────────────────────────────
WORK_DIR      = Path(os.environ.get("WORK_DIR",     "."))
TEMPLATES_DIR = Path(os.environ.get("TEMPLATES_DIR", str(WORK_DIR)))
LOG_FILE      = WORK_DIR / "logs" / "06_align_features.log"

os.chdir(WORK_DIR)

# ── Logger ─────────────────────────────────────────────────────
def log(msg: str) -> None:
    line = f"[{datetime.now():%H:%M:%S}] {msg}"
    print(line)
    with open(LOG_FILE, "a") as fh:
        fh.write(line + "\n")

(WORK_DIR / "logs").mkdir(parents=True, exist_ok=True)
with open(LOG_FILE, "w") as fh:
    fh.write("=" * 60 + "\n")
    fh.write("E. COLI SCRIPT 06: FEATURE ALIGNMENT — PER ANTIBIOTIC\n")
    fh.write(f"Started: {datetime.now()}\n")
    fh.write("=" * 60 + "\n")


# ── 1. Load extracted production features ──────────────────────
log("\nLoading extracted feature files...")

for fname in ["gene_presence_production.csv", "kmer_production.csv", "snp_production.csv"]:
    if not (WORK_DIR / fname).exists():
        log(f"   ERROR: {fname} not found — ensure scripts 01b, 03, 04b completed")
        sys.exit(1)

genes_df = pd.read_csv(WORK_DIR / "gene_presence_production.csv")
kmers_df = pd.read_csv(WORK_DIR / "kmer_production.csv")
snps_df  = pd.read_csv(WORK_DIR / "snp_production.csv")

log(f"  Gene features  : {genes_df.shape[1] - 1:,}")
log(f"  K-mer features : {kmers_df.shape[1] - 1:,}")
log(f"  SNP features   : {snps_df.shape[1] - 1:,}")

snps_available = snps_df.shape[1] > 1
if not snps_available:
    log("   No SNP features available — models will predict from genes + k-mers only")

# ── 2. Merge into one raw feature row ──────────────────────────
log("\nMerging extracted features into one row...")

genome_id = genes_df["Genome_ID"].iloc[0]
merged = genes_df.copy()

if kmers_df.shape[1] > 1:
    merged = merged.merge(kmers_df, on="Genome_ID", how="left")
else:
    log("   No k-mer features — all k-mer columns will be 0")

if snps_available:
    merged = merged.merge(snps_df, on="Genome_ID", how="left")

merged.fillna(0, inplace=True)
for col in [c for c in merged.columns if c != "Genome_ID"]:
    merged[col] = merged[col].astype(int)

available = set(merged.columns) - {"Genome_ID"}
log(f"  Total raw features available: {len(available):,}")

# ── 3. Per-antibiotic alignment function ───────────────────────
def align_to_template(template_path: Path, output_name: str, label: str) -> float:
    """
    Align merged production features to a training template.
    Missing features → 0.  Extra features → dropped.
    Returns the % of template features that were missing.
    """
    log(f"\n  Aligning to: {label}")
    log(f"  Template   : {template_path.name}")

    if not template_path.exists():
        log(f"   ERROR: Template not found: {template_path}")
        sys.exit(1)

    with open(template_path) as fh:
        template_features = [line.strip() for line in fh if line.strip()]

    found   = [f for f in template_features if f in available]
    missing = [f for f in template_features if f not in available]
    missing_pct = len(missing) / max(len(template_features), 1) * 100

    log(f"  Template features : {len(template_features):,}")
    log(f"  Matched           : {len(found):,}")
    log(f"  Missing → 0       : {len(missing):,}  ({missing_pct:.1f}%)")
    log("  Note: ~65% zero-rate is normal in genomic feature data")

    # Build aligned row: matched columns + missing filled with 0
    row_matched = merged[["Genome_ID"] + found].copy()
    if missing:
        missing_df = pd.DataFrame(0, index=row_matched.index, columns=missing)
        row_matched = pd.concat([row_matched, missing_df], axis=1)

    # Enforce exact template column order
    aligned = row_matched[["Genome_ID"] + template_features]
    out_path = WORK_DIR / output_name
    aligned.to_csv(out_path, index=False)

    if missing_pct > 95:
        log(f"   WARNING: {missing_pct:.1f}% missing exceeds 95% threshold!")
        log("   Check upstream extraction scripts for failures")
    else:
        log(f"   Alignment within expected range")

    log(f"   Saved: {output_name}  (shape {aligned.shape})")
    return missing_pct


# ── 5. Run alignment for all 8 antibiotics ─────────────────────
log("\n" + "=" * 60)
log("ALIGNING TO FEATURE TEMPLATE")
log("=" * 60)

# Primary shared template
template_file = TEMPLATES_DIR / "features_ecoli.txt"
pct_full = align_to_template(template_file, "aligned_full.csv", "E. coli feature template (Genes + K-mers + SNPs)")

# ── 5. Summary ─────────────────────────────────────────────────
log("\n" + "=" * 60)
log("ALIGNMENT SUMMARY")
log("=" * 60)
log(f"  Genome ID       : {genome_id}")
log(f"  Full dataset    : {pct_full:.1f}% missing → aligned_full.csv")

if pct_full > 95:
    log("")
    log("   ALIGNMENT EXCEEDED 95% MISSING THRESHOLD")
    log("  Check scripts 01b, 03, 04b before running prediction")
    sys.exit(1)
else:
    log("")
    log("  ✓ Alignment within expected range")
    log("  ✓ Ready for prediction")

log("")
log("  Next step: python3 07_predict.py  (handled by API)")
log("=" * 60)
log(f"SCRIPT 06 COMPLETE: {datetime.now()}")
log("=" * 60)


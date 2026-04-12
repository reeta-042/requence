#!/usr/bin/env python3
# =============================================================
# E. COLI — SCRIPT 06: FEATURE ALIGNMENT (PER ANTIBIOTIC)
# Input:  gene_presence_production.csv
#         kmer_production.csv
#         snp_production.csv
#         $TEMPLATES_DIR/features_amoxillin_ecoli.txt
#         $TEMPLATES_DIR/features_coxime_ecoli.txt
#         $TEMPLATES_DIR/features_doxy_ecoli.txt
#         $TEMPLATES_DIR/features_levo_ecoli.txt
#         $TEMPLATES_DIR/features_nali_ecoli.txt
#         $TEMPLATES_DIR/features_norflo_ecoli.txt
#         $TEMPLATES_DIR/features_strep_ecoli.txt
#         $TEMPLATES_DIR/features_tetra_ecoli.txt
# Output: aligned_amoxicillin.csv
#         aligned_coxime.csv
#         aligned_doxy.csv
#         aligned_levo.csv
#         aligned_nali.csv
#         aligned_norflo.csv
#         aligned_strep.csv
#         aligned_tetra.csv
#         logs/06_align_features.log
# NOTE:   E. coli has ONE aligned CSV per antibiotic.
#         Each model loads only its own aligned file so there
#         are zero spurious zero-fills from unrelated features.
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

# ── Antibiotic → template + output file mapping ────────────────
# Keys here match the "aligned" field in app.py PATHOGEN_CONFIG.
ANTIBIOTIC_CONFIG = {
    "amoxicillin": {
        "full_name": "Amoxicillin/Clavulanic Acid",
        "template":  "features_amoxillin_ecoli.txt",
        "output":    "aligned_amoxicillin.csv",
    },
    "coxime": {
        "full_name": "Cefuroxime",
        "template":  "features_coxime_ecoli.txt",
        "output":    "aligned_coxime.csv",
    },
    "doxy": {
        "full_name": "Doxycycline",
        "template":  "features_doxy_ecoli.txt",
        "output":    "aligned_doxy.csv",
    },
    "levo": {
        "full_name": "Levofloxacin",
        "template":  "features_levo_ecoli.txt",
        "output":    "aligned_levo.csv",
    },
    "nali": {
        "full_name": "Nalidixic Acid",
        "template":  "features_nali_ecoli.txt",
        "output":    "aligned_nali.csv",
    },
    "norflo": {
        "full_name": "Norfloxacin",
        "template":  "features_norflo_ecoli.txt",
        "output":    "aligned_norflo.csv",
    },
    "strep": {
        "full_name": "Streptomycin",
        "template":  "features_strep_ecoli.txt",
        "output":    "aligned_strep.csv",
    },
    "tetra": {
        "full_name": "Tetracycline",
        "template":  "features_tetra_ecoli.txt",
        "output":    "aligned_tetra.csv",
    },
}

# ── 1. Validate template files ─────────────────────────────────
log("Validating feature template files...")
missing_templates = []
for ab, cfg in ANTIBIOTIC_CONFIG.items():
    tf = TEMPLATES_DIR / cfg["template"]
    if not tf.exists():
        log(f"  ✗ MISSING: {cfg['template']}  ({cfg['full_name']})")
        missing_templates.append(cfg["template"])
    else:
        log(f"  ✓ {cfg['template']}")

if missing_templates:
    log(f"\n✗ ERROR: {len(missing_templates)} template file(s) missing. Cannot continue.")
    sys.exit(1)

# ── 2. Load extracted production features ──────────────────────
log("\nLoading extracted feature files...")

for fname in ["gene_presence_production.csv", "kmer_production.csv", "snp_production.csv"]:
    if not (WORK_DIR / fname).exists():
        log(f"  ✗ ERROR: {fname} not found — ensure scripts 01b, 03, 04b completed")
        sys.exit(1)

genes_df = pd.read_csv(WORK_DIR / "gene_presence_production.csv")
kmers_df = pd.read_csv(WORK_DIR / "kmer_production.csv")
snps_df  = pd.read_csv(WORK_DIR / "snp_production.csv")

log(f"  Gene features  : {genes_df.shape[1] - 1:,}")
log(f"  K-mer features : {kmers_df.shape[1] - 1:,}")
log(f"  SNP features   : {snps_df.shape[1] - 1:,}")

snps_available = snps_df.shape[1] > 1
if not snps_available:
    log("  ⚠ No SNP features available — models will predict from genes + k-mers only")

# ── 3. Merge into one raw feature row ──────────────────────────
log("\nMerging extracted features into one row...")

genome_id = genes_df["Genome_ID"].iloc[0]
merged = genes_df.copy()

if kmers_df.shape[1] > 1:
    merged = merged.merge(kmers_df, on="Genome_ID", how="left")
else:
    log("  ⚠ No k-mer features — all k-mer columns will be 0")

if snps_available:
    merged = merged.merge(snps_df, on="Genome_ID", how="left")

merged.fillna(0, inplace=True)
for col in [c for c in merged.columns if c != "Genome_ID"]:
    merged[col] = merged[col].astype(int)

available = set(merged.columns) - {"Genome_ID"}
log(f"  Total raw features available: {len(available):,}")

# ── 4. Per-antibiotic alignment function ───────────────────────
def align_to_template(ab_key: str, cfg: dict) -> float:
    """
    Build one aligned CSV for a single antibiotic using its own
    feature template.  Missing features → 0.  Extra features dropped.
    Uses pd.concat to add missing columns in one operation (avoids
    DataFrame fragmentation performance warnings).
    Returns % missing.
    """
    template_path = TEMPLATES_DIR / cfg["template"]
    output_path   = WORK_DIR / cfg["output"]

    log(f"\n  {'─' * 55}")
    log(f"  {cfg['full_name']}")
    log(f"  Template : {cfg['template']}")

    with open(template_path) as fh:
        template_features = [line.strip() for line in fh if line.strip()]

    found   = [f for f in template_features if f in available]
    missing = [f for f in template_features if f not in available]
    missing_pct = len(missing) / max(len(template_features), 1) * 100

    log(f"  Template features  : {len(template_features):,}")
    log(f"  Matched            : {len(found):,}")
    log(f"  Missing → 0        : {len(missing):,}  ({missing_pct:.1f}%)")
    log("  Note: ~65–77% zeros in training is normal for genomic features")

    # Build the aligned row: matched features first
    row = merged[["Genome_ID"] + found].copy()

    # Bulk-add all missing features at once (avoids fragmentation)
    if missing:
        zeros_df = pd.DataFrame(0, index=row.index, columns=missing)
        row = pd.concat([row, zeros_df], axis=1)

    # Reorder to exact template order
    row = row[["Genome_ID"] + template_features]
    row.to_csv(output_path, index=False)

    if missing_pct > 95:
        log(f"  ⚠ WARNING: {missing_pct:.1f}% missing exceeds 95% threshold — check extraction logs")
    else:
        log(f"  ✓ Within expected range")
    log(f"  ✓ Saved: {cfg['output']}  (shape {row.shape})")

    return missing_pct


# ── 5. Run alignment for all 8 antibiotics ─────────────────────
log("\n" + "=" * 60)
log("ALIGNING TO PER-ANTIBIOTIC FEATURE TEMPLATES")
log("=" * 60)

summary_rows = []
any_failed   = False

for ab_key, cfg in ANTIBIOTIC_CONFIG.items():
    pct = align_to_template(ab_key, cfg)
    summary_rows.append((cfg["full_name"], cfg["output"], pct))
    if pct > 95:
        any_failed = True

# ── 6. Print summary table ─────────────────────────────────────
log("\n" + "=" * 60)
log("ALIGNMENT SUMMARY")
log("=" * 60)
log(f"  Genome ID      : {genome_id}")
log(f"  SNPs available : {'Yes' if snps_available else 'No — genes + k-mers only'}")
log("")
log(f"  {'Antibiotic':<35} {'Output file':<28} {'Missing':>8}")
log(f"  {'─' * 73}")
for full_name, out_file, pct in summary_rows:
    flag = "  ⚠ CHECK" if pct > 95 else ""
    log(f"  {full_name:<35} {out_file:<28} {pct:>6.1f}%{flag}")

log("")
if any_failed:
    log("  ⚠ One or more alignments exceeded 95% missing — review extraction logs")
    log("  ⚠ Pipeline will continue but prediction quality may be degraded")
else:
    log("  ✓ All alignments within expected range")
    log("  ✓ Ready for prediction")

log("")
log("  Next step: prediction handled by API (app.py)")
log("=" * 60)
log(f"SCRIPT 06 COMPLETE: {datetime.now()}")
log("=" * 60)
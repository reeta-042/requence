#!/usr/bin/env python3
# =============================================================
# E. COLI — SCRIPT 04b: PROCESS SNP EXTRACTION OUTPUT
# Input:  snippy_production_out/snps.tab
#         $TEMPLATES_DIR/features_*.txt  (JUST ONE E. COLI antibiotic template)
# Output: snp_production.csv
#         logs/04b_process_snps.log
# NOTE:   We take the UNION of SNP feature names across all 
#         templates so one pass covers every model's SNP needs.
#         Vectorised matching used for speed (vs row-by-row loop).
# =============================================================

import os
import sys
import pandas as pd
from datetime import datetime
from pathlib import Path

# ── Paths ──────────────────────────────────────────────────────
WORK_DIR      = Path(os.environ.get("WORK_DIR",     "."))
TEMPLATES_DIR = Path(os.environ.get("TEMPLATES_DIR", str(WORK_DIR)))
SNIPPY_TAB    = WORK_DIR / "snippy_production_out" / "snps.tab"
LOG_FILE      = WORK_DIR / "logs" / "04b_process_snps.log"

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
    fh.write("E. COLI SCRIPT 04b: PROCESS SNP EXTRACTION\n")
    fh.write(f"Started: {datetime.now()}\n")
    fh.write("=" * 60 + "\n")

if not SNIPPY_TAB.exists():
    log(f"✗ ERROR: {SNIPPY_TAB} not found")
    log("  Ensure 04_extract_snps.sh completed successfully")
    sys.exit(1)

# ── 1. Collect training SNP names from ALL 8 templates ─────────
log("Loading SNP feature names from all antibiotic templates...")
training_snps: set = set()

template_files = sorted(TEMPLATES_DIR.glob("features_*.txt"))
if not template_files:
    log(f"⚠ WARNING: No feature templates found in {TEMPLATES_DIR}")

for tf in template_files:
    before = len(training_snps)
    with open(tf) as fh:
        for line in fh:
            feat = line.strip()
            # E. coli SNP format: NC_XXXXXX_<POS>_<REF>><ALT>
            if ">" in feat and feat.startswith("NC_"):
                training_snps.add(feat)
    added = len(training_snps) - before
    log(f"  ✓ {tf.name}  (+{added} SNPs)")

log(f"\n  Total unique training SNPs across all templates: {len(training_snps):,}")

if len(training_snps) == 0:
    log("⚠ WARNING: No SNPs found in templates — SNP features will all be 0.")

# ── 2. Parse Snippy snps.tab ────────────────────────────────────
log("\nParsing Snippy snps.tab...")
try:
    snp_tab = pd.read_csv(SNIPPY_TAB, sep="\t", dtype=str)
    log(f"  Raw variants : {len(snp_tab)}")
    log(f"  Columns      : {list(snp_tab.columns)}")
except Exception as exc:
    log(f"✗ ERROR reading snps.tab: {exc}")
    sys.exit(1)

required_cols = {"CHROM", "POS", "REF", "ALT"}
missing_cols  = required_cols - set(snp_tab.columns)
if missing_cols:
    log(f"✗ ERROR: Missing columns in snps.tab: {missing_cols}")
    sys.exit(1)

# ── 3. Build feature names (vectorised) ────────────────────────
log("Building SNP feature names (vectorised)...")

# Strip chromosome version suffix: NC_000913.3 → NC_000913
snp_tab["CHROM_CLEAN"]  = snp_tab["CHROM"].str.split(".").str[0]
snp_tab["Feature_Name"] = (
    snp_tab["CHROM_CLEAN"]      + "_" +
    snp_tab["POS"].astype(str)  + "_" +
    snp_tab["REF"]              + ">" +
    snp_tab["ALT"]
)

# Drop rows where any key field was NaN
snp_tab.dropna(subset=["CHROM_CLEAN", "POS", "REF", "ALT"], inplace=True)

# ── 4. Match against training SNPs ─────────────────────────────
matched_df     = snp_tab[snp_tab["Feature_Name"].isin(training_snps)]
found_features = matched_df["Feature_Name"].unique()
matched_snps   = len(found_features)
total_variants = len(snp_tab)
missing_pct    = (1 - matched_snps / max(len(training_snps), 1)) * 100

log(f"  Total variants processed         : {total_variants:,}")
log(f"  SNPs matching training features  : {matched_snps:,}")
log(f"  Coverage                         : {matched_snps}/{len(training_snps):,} "
    f"({100 - missing_pct:.1f}%)")
log("  Note: Training rows average ~71% zero SNPs — high zeros are expected")

# Debug output when suspiciously few SNPs matched
if matched_snps < 10:
    log("  ⚠ Very few SNPs matched — debug info:")
    log(f"    Sample training SNP names : {list(training_snps)[:3]}")
    log("    Sample Snippy feature names (first 3):")
    for name in snp_tab["Feature_Name"].head(3).tolist():
        log(f"      Built: {name}")
    log("    Check that REFERENCE genome chromosome IDs match training data.")

# ── 5. Save ───────────────────────────────────────────────────
snp_row: dict = {feat: 1 for feat in found_features}
snp_row["Genome_ID"] = "query_genome"

pd.DataFrame([snp_row]).to_csv(WORK_DIR / "snp_production.csv", index=False)

log(f"\n✓ Saved: snp_production.csv")
log(f"  Features present (=1) : {matched_snps:,}")
log(f"  Features absent  (=0) : {len(training_snps) - matched_snps:,}")
log("  Next step: python3 06_align_features.py")
log("=" * 60)
log(f"SCRIPT 04b COMPLETE: {datetime.now()}")
log("=" * 60)
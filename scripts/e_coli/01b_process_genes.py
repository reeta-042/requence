#!/usr/bin/env python3
# =============================================================
# E. COLI — SCRIPT 01b: PROCESS GENE EXTRACTION OUTPUT
# Input:  card_production.tsv
#         resfinder_production.tsv
#         gene_summary_production.tsv
# Output: gene_presence_production.csv
#         logs/01b_process_genes.log
# =============================================================

import os
import sys
import pandas as pd
from datetime import datetime
from pathlib import Path

# ── Paths (injected by API, fallback to cwd for standalone use) ─
WORK_DIR = Path(os.environ.get("WORK_DIR", "."))
LOG_FILE = WORK_DIR / "logs" / "01b_process_genes.log"
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
    fh.write("E. COLI SCRIPT 01b: PROCESS GENE EXTRACTION\n")
    fh.write(f"Started: {datetime.now()}\n")
    fh.write("=" * 60 + "\n")

# ── 1. Resolve genome ID from the query file directly ──────────
# Using the filename avoids the ABRicate summary path-as-ID bug.
if not (WORK_DIR / "query_genome.fna").exists():
    log("✗ ERROR: query_genome.fna not found")
    sys.exit(1)
GENOME_ID = "query_genome"
log(f"Genome ID: {GENOME_ID}")

# ── 2. Load ABRicate outputs ────────────────────────────────────
log("Loading ABRicate outputs...")
try:
    card = pd.read_csv("card_production.tsv", sep="\t")
    res  = pd.read_csv("resfinder_production.tsv", sep="\t")
    log(f"  CARD hits    : {len(card)}")
    log(f"  ResFinder hits: {len(res)}")
except Exception as exc:
    log(f"✗ ERROR loading TSV files: {exc}")
    sys.exit(1)

try:
    summary = pd.read_csv("gene_summary_production.tsv", sep="\t")
    log(f"  Summary shape: {summary.shape}")
except Exception as exc:
    log(f"✗ ERROR loading summary matrix: {exc}")
    sys.exit(1)

# ── 3. Clean and binarise ───────────────────────────────────────
log("Processing gene matrix...")

summary.rename(columns={"#FILE": "Genome_ID"}, inplace=True)
summary["Genome_ID"] = GENOME_ID

# Drop bookkeeping columns produced by ABRicate
drop_cols = [c for c in summary.columns if "NUM_FOUND" in c.upper()]
if drop_cols:
    summary.drop(columns=drop_cols, inplace=True)
    log(f"  Dropped non-feature columns: {drop_cols}")

def _make_binary(val) -> int:
    return 0 if str(val).strip() in {"0", "", "0.0", ".", "nan"} else 1

feature_cols = [c for c in summary.columns if c != "Genome_ID"]
for col in feature_cols:
    summary[col] = summary[col].apply(_make_binary)

# ── 4. Collapse multi-row summary into one genome row ──────────
# ABRicate summary has one row per input TSV file (2 here).
# A gene is PRESENT if detected in EITHER database.
if len(summary) > 1:
    feat_vals = summary[feature_cols].max(axis=0)
    final_row = pd.DataFrame([feat_vals])
    final_row.insert(0, "Genome_ID", GENOME_ID)
    summary = final_row
    log(f"  Collapsed {len(summary)+1} database rows → 1 genome row using union rule")

# ── 5. Save ─────────────────────────────────────────────────────
out_path = WORK_DIR / "gene_presence_production.csv"
summary.to_csv(out_path, index=False)

present = int((summary.iloc[0, 1:] == 1).sum())
absent  = int((summary.iloc[0, 1:] == 0).sum())
log(f"\n  Total gene features : {summary.shape[1] - 1}")
log(f"  Genes present  (=1) : {present}")
log(f"  Genes absent   (=0) : {absent}")
log(f"\n✓ Saved: gene_presence_production.csv")
log("  Next step: bash 02_create_blast_db.sh")
log("=" * 60)
log(f"SCRIPT 01b COMPLETE: {datetime.now()}")
log("=" * 60)
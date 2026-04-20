#!/bin/bash
# =============================================================
# E. COLI — SCRIPT 01: GENE EXTRACTION (ABRicate)
# Input:  $WORK_DIR/query_genome.fna
# Output: card_production.tsv
#         resfinder_production.tsv
#         gene_summary_production.tsv
#         logs/01_extract_genes.log
# =============================================================

set -euo pipefail

# ── Resolve paths ─────────────────────────────────────────────
# WORK_DIR is injected by the API at runtime via the env.
# Fallback lets the script still run standalone for debugging.
WORK_DIR="${WORK_DIR:-$(pwd)}"
LOG_DIR="$WORK_DIR/logs"
LOG_FILE="$LOG_DIR/01_extract_genes.log"

cd "$WORK_DIR"
mkdir -p "$LOG_DIR"

# ── Header ────────────────────────────────────────────────────
{
echo "============================================================"
echo "E. COLI SCRIPT 01: GENE EXTRACTION"
echo "Work dir : $WORK_DIR"
echo "Started  : $(date)"
echo "============================================================"
} | tee "$LOG_FILE"

# ── Input validation ──────────────────────────────────────────
if [ ! -f "query_genome.fna" ]; then
    echo "✗ ERROR: query_genome.fna not found in $WORK_DIR" | tee -a "$LOG_FILE"
    exit 1
fi

GENOME_SIZE=$(wc -c < "query_genome.fna")
echo "✓ Input genome: query_genome.fna ($GENOME_SIZE bytes)" | tee -a "$LOG_FILE"
echo "" | tee -a "$LOG_FILE"

# ── Step 1/3 : ABRicate — CARD database ───────────────────────
echo "[1/3] ABRicate (CARD database)..." | tee -a "$LOG_FILE"
START=$(date +%s)

abricate --db card query_genome.fna \
    > card_production.tsv \
    2>> "$LOG_FILE"

if [ $? -ne 0 ]; then
    echo "✗ ERROR: ABRicate CARD failed. See $LOG_FILE" | tee -a "$LOG_FILE"
    exit 1
fi
END=$(date +%s)
CARD_HITS=$(tail -n +2 card_production.tsv | wc -l)
echo "  ✓ CARD done in $((END-START))s — $CARD_HITS hits" | tee -a "$LOG_FILE"

# ── Step 2/3 : ABRicate — ResFinder database ──────────────────
echo "[2/3] ABRicate (ResFinder database)..." | tee -a "$LOG_FILE"
START=$(date +%s)

abricate --db resfinder query_genome.fna \
    > resfinder_production.tsv \
    2>> "$LOG_FILE"

if [ $? -ne 0 ]; then
    echo "✗ ERROR: ABRicate ResFinder failed. See $LOG_FILE" | tee -a "$LOG_FILE"
    exit 1
fi
END=$(date +%s)
RES_HITS=$(tail -n +2 resfinder_production.tsv | wc -l)
echo "  ✓ ResFinder done in $((END-START))s — $RES_HITS hits" | tee -a "$LOG_FILE"

# ── Step 3/3 : Summary matrix ─────────────────────────────────
echo "[3/3] Creating ABRicate summary matrix..." | tee -a "$LOG_FILE"

abricate --summary card_production.tsv resfinder_production.tsv \
    > gene_summary_production.tsv \
    2>> "$LOG_FILE"

if [ $? -ne 0 ]; then
    echo "✗ ERROR: ABRicate summary failed. See $LOG_FILE" | tee -a "$LOG_FILE"
    exit 1
fi

# ── Summary ───────────────────────────────────────────────────
echo "" | tee -a "$LOG_FILE"
echo "  Output files:" | tee -a "$LOG_FILE"
echo "    card_production.tsv        ($CARD_HITS hits)" | tee -a "$LOG_FILE"
echo "    resfinder_production.tsv   ($RES_HITS hits)" | tee -a "$LOG_FILE"
echo "    gene_summary_production.tsv" | tee -a "$LOG_FILE"
echo "" | tee -a "$LOG_FILE"
echo "  Next step: python3 01b_process_genes.py" | tee -a "$LOG_FILE"
echo "============================================================" | tee -a "$LOG_FILE"
echo "SCRIPT 01 COMPLETE: $(date)" | tee -a "$LOG_FILE"
echo "============================================================" | tee -a "$LOG_FILE"
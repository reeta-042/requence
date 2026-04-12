#!/bin/bash
# =============================================================
# SALMONELLA — SCRIPT 04: SNP EXTRACTION (Snippy)
# Input:  query_genome.fna
#         $REFERENCE  (salmonella_LT2.gbff — injected by API)
# Output: snippy_production_out/snps.tab
#         logs/04_extract_snps.log
# NOTE:   This is typically the slowest step (~10-20 min).
#         $REFERENCE falls back to /app/reference/salmonella_LT2.gbff
# =============================================================

set -euo pipefail

WORK_DIR="${WORK_DIR:-$(pwd)}"
REFERENCE="${REFERENCE:-/app/reference/salmonella_LT2.gbff}"
LOG_DIR="$WORK_DIR/logs"
LOG_FILE="$LOG_DIR/04_extract_snps.log"
SNIPPY_OUTDIR="$WORK_DIR/snippy_production_out"

cd "$WORK_DIR"
mkdir -p "$LOG_DIR"

{
echo "============================================================"
echo "SALMONELLA SCRIPT 04: SNP EXTRACTION (Snippy)"
echo "Work dir  : $WORK_DIR"
echo "Reference : $REFERENCE"
echo "Started   : $(date)"
echo "============================================================"
} | tee "$LOG_FILE"

# ── Validate inputs ───────────────────────────────────────────
if [ ! -f "query_genome.fna" ]; then
    echo "✗ ERROR: query_genome.fna not found" | tee -a "$LOG_FILE"
    exit 1
fi

if [ ! -f "$REFERENCE" ]; then
    echo "✗ ERROR: Reference genome not found: $REFERENCE" | tee -a "$LOG_FILE"
    echo "  Set REFERENCE env variable or place salmonella_LT2.gbff in /app/reference/" | tee -a "$LOG_FILE"
    exit 1
fi

echo "✓ Genome   : query_genome.fna" | tee -a "$LOG_FILE"
echo "✓ Reference: $REFERENCE" | tee -a "$LOG_FILE"
echo "" | tee -a "$LOG_FILE"

# ── Clean previous Snippy run (ensures idempotency) ──────────
if [ -d "$SNIPPY_OUTDIR" ]; then
    echo "  Removing previous Snippy output dir..." | tee -a "$LOG_FILE"
    rm -rf "$SNIPPY_OUTDIR"
fi

# ── Run Snippy ────────────────────────────────────────────────
echo "[1/2] Running Snippy..." | tee -a "$LOG_FILE"
START=$(date +%s)

snippy \
    --outdir "$SNIPPY_OUTDIR" \
    --ref    "$REFERENCE" \
    --ctgs   "query_genome.fna" \
    --cpus   4 \
    --force \
    >> "$LOG_FILE" 2>&1

SNIPPY_EXIT=$?
END=$(date +%s)

if [ $SNIPPY_EXIT -ne 0 ]; then
    echo "✗ ERROR: Snippy exited with code $SNIPPY_EXIT. See $LOG_FILE" | tee -a "$LOG_FILE"
    exit 1
fi

echo "  ✓ Snippy completed in $((END-START))s" | tee -a "$LOG_FILE"

# ── Validate output ───────────────────────────────────────────
echo "[2/2] Validating Snippy output..." | tee -a "$LOG_FILE"

if [ ! -f "$SNIPPY_OUTDIR/snps.tab" ]; then
    echo "✗ ERROR: snps.tab not found — Snippy may have produced no variants" | tee -a "$LOG_FILE"
    # Create empty placeholder so downstream scripts degrade gracefully
    mkdir -p "$SNIPPY_OUTDIR"
    echo -e "CHROM\tPOS\tTYPE\tREF\tALT\tEVIDENCE\tFTYPE\tSTRAND\tNT_POS\tAA_POS\tEFFECT\tLOCUS_TAG\tGENE\tPRODUCT" \
        > "$SNIPPY_OUTDIR/snps.tab"
    echo "  ⚠ Created empty snps.tab — genome may lack detectable variants vs reference" | tee -a "$LOG_FILE"
else
    SNP_COUNT=$(tail -n +2 "$SNIPPY_OUTDIR/snps.tab" | wc -l)
    echo "  ✓ Raw variants found: $SNP_COUNT" | tee -a "$LOG_FILE"
fi

echo "" | tee -a "$LOG_FILE"
echo "  Output: $SNIPPY_OUTDIR/snps.tab" | tee -a "$LOG_FILE"
echo "  Next step: python3 04b_process_snps.py" | tee -a "$LOG_FILE"
echo "============================================================" | tee -a "$LOG_FILE"
echo "SCRIPT 04 COMPLETE: $(date)" | tee -a "$LOG_FILE"
echo "============================================================" | tee -a "$LOG_FILE"
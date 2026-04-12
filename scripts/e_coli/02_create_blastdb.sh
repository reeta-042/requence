#!/bin/bash
# =============================================================
# E. COLI — SCRIPT 02: CREATE BLAST DATABASE
# Input:  query_genome.fna
# Output: blast_db/query_genome.*  (BLAST index files)
#         logs/02_create_blast_db.log
# =============================================================

set -euo pipefail

WORK_DIR="${WORK_DIR:-$(pwd)}"
LOG_DIR="$WORK_DIR/logs"
LOG_FILE="$LOG_DIR/02_create_blast_db.log"
DB_DIR="$WORK_DIR/blast_db"

cd "$WORK_DIR"
mkdir -p "$LOG_DIR" "$DB_DIR"

{
echo "============================================================"
echo "E. COLI SCRIPT 02: CREATE BLAST DATABASE"
echo "Work dir : $WORK_DIR"
echo "Started  : $(date)"
echo "============================================================"
} | tee "$LOG_FILE"

if [ ! -f "query_genome.fna" ]; then
    echo "✗ ERROR: query_genome.fna not found" | tee -a "$LOG_FILE"
    exit 1
fi
echo "✓ Input: query_genome.fna" | tee -a "$LOG_FILE"
echo "" | tee -a "$LOG_FILE"

echo "Building nucleotide BLAST database..." | tee -a "$LOG_FILE"
START=$(date +%s)

makeblastdb \
    -in  "query_genome.fna" \
    -dbtype nucl \
    -out  "$DB_DIR/query_genome" \
    -parse_seqids \
    >> "$LOG_FILE" 2>&1

if [ $? -ne 0 ]; then
    echo "✗ ERROR: makeblastdb failed. See $LOG_FILE" | tee -a "$LOG_FILE"
    exit 1
fi

END=$(date +%s)
echo "  ✓ BLAST DB created in $((END-START))s" | tee -a "$LOG_FILE"
echo "  ✓ Location: $DB_DIR/" | tee -a "$LOG_FILE"
echo "" | tee -a "$LOG_FILE"
echo "  Next step: python3 03_extract_kmers.py" | tee -a "$LOG_FILE"
echo "============================================================" | tee -a "$LOG_FILE"
echo "SCRIPT 02 COMPLETE: $(date)" | tee -a "$LOG_FILE"
echo "============================================================" | tee -a "$LOG_FILE"
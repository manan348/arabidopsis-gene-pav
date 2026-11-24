#!/usr/bin/env bash
# Usage: bash scripts/trim_subsample.sh SAMPLE [MIN_LEN]
# Trims full reads with fastp, then subsamples to ~15x (1.8e9 bases).
set -euo pipefail
S=$1
LEN=${2:-36}
GENOME=119146348
TARGET=1800000000
RAW=data/reads/$S
TRIM=data/reads_trimmed_full/$S
FIN=data/reads_final/$S
mkdir -p "$TRIM" "$FIN" results/fastp

[ -f "$RAW/${S}_R1.fastq.gz" ] || { echo "missing $RAW/${S}_R1.fastq.gz"; exit 1; }

fastp -i "$RAW/${S}_R1.fastq.gz" -I "$RAW/${S}_R2.fastq.gz" \
  -o "$TRIM/${S}_R1.trim.fastq.gz" -O "$TRIM/${S}_R2.trim.fastq.gz" \
  --detect_adapter_for_pe --cut_right --cut_window_size 4 --cut_mean_quality 20 \
  --length_required "$LEN" -w 3 \
  -j "results/fastp/${S}_full.json" -h "results/fastp/${S}_full.html" \
  2> "results/fastp/${S}_full.log"

BASES=$(grep -m1 -A3 '"after_filtering"' "results/fastp/${S}_full.json" | grep total_bases | tr -dc '0-9')
FRAC=$(awk -v b="$BASES" -v t="$TARGET" 'BEGIN{f=t/b; if(f>1)f=1; printf "%.4f", f}')
DEPTH=$(awk -v b="$BASES" -v g="$GENOME" 'BEGIN{printf "%.1f", b/g}')
echo "$S trimmed_bases=$BASES depth=${DEPTH}x fraction=$FRAC"

if [ "$FRAC" = "1.0000" ]; then
  echo "trimmed depth is below 15x: keeping all reads"
  cp "$TRIM/${S}_R1.trim.fastq.gz" "$FIN/${S}_R1.fastq.gz"
  cp "$TRIM/${S}_R2.trim.fastq.gz" "$FIN/${S}_R2.fastq.gz"
else
  seqtk sample -s 100 "$TRIM/${S}_R1.trim.fastq.gz" "$FRAC" | gzip > "$FIN/${S}_R1.fastq.gz"
  seqtk sample -s 100 "$TRIM/${S}_R2.trim.fastq.gz" "$FRAC" | gzip > "$FIN/${S}_R2.fastq.gz"
fi
seqkit stats -j 3 "$FIN"/*.fastq.gz
echo "DONE $S"

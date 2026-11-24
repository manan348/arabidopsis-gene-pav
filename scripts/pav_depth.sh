#!/usr/bin/env bash
# Usage: bash scripts/pav_depth.sh SAMPLE
# Per-gene read depth for PAV calling (Phase 4). Two mosdepth runs on the markdup BAM:
#   q10 = only reads with MAPQ >= 10 (the PAV evidence)
#   q0  = all reads (used to flag genes whose reads mostly have low MAPQ = "uncertain", not "absent")
# Duplicate-marked reads are skipped by mosdepth by default.
# Output in results/pav/SAMPLE/:
#   SAMPLE.q10.regions.bed.gz      chrom start end gene mean_depth
#   SAMPLE.q10.thresholds.bed.gz   chrom start end gene bases>=1x bases>=5x
#   (same two files for q0)
set -euo pipefail

S=$1
BAM=data/bam/$S/$S.markdup.bam
BED=data/reference/genes.bed
OUT=results/pav/$S

[ -f "$BAM" ] || { echo "missing $BAM"; exit 1; }
[ -f "$BED" ] || { echo "missing $BED"; exit 1; }
[ ! -e "$OUT/$S.q10.regions.bed.gz" ] || { echo "$OUT already has results, not overwriting"; exit 1; }

mkdir -p "$OUT"
echo "$(date) start pav_depth $S"

for Q in 10 0; do
  # -n = no per-base file (not needed, saves time and disk), -t 1 = one thread,
  # --by = per-gene windows, -T = breadth thresholds (bases covered at >=1x and >=5x), -Q = min MAPQ
  mosdepth -n -t 1 --by "$BED" -T 1,5 -Q "$Q" "$OUT/$S.q$Q" "$BAM"
done

echo "DONE pav_depth $S"

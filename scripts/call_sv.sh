#!/usr/bin/env bash
# Usage: bash scripts/call_sv.sh SAMPLE
# Deletion calling with Delly for ONE accession (Phase 5), run in env `sv` (Delly 2.8.2; in Delly 2.x `call` is named `sr`).
# Delly looks at read pairs that map too far apart and at split reads, and writes candidate deletions.
# Settings: -t DEL (deletions only); Delly's default quality thresholds are used.
# Output: results/sv/SAMPLE.delly.bcf (+ index). Refuses to overwrite an existing file.
set -euo pipefail

S=$1
REF=data/reference/TAIR10.1_Col-0/GCF_000001735.4/GCF_000001735.4_TAIR10.1_genomic.fna
BAM=data/bam/$S/$S.markdup.bam
OUTDIR=results/sv
OUT=$OUTDIR/$S.delly.bcf
TMP=$OUTDIR/$S.tmp.delly.bcf

[ -f "$BAM" ] || { echo "missing $BAM"; exit 1; }
[ -f "$REF.fai" ] || { echo "missing $REF.fai"; exit 1; }
[ ! -e "$OUT" ] || { echo "$OUT already exists, not overwriting"; exit 1; }

mkdir -p "$OUTDIR"
echo "$(date) start delly $S"
delly sr -t DEL -g "$REF" -o "$TMP" "$BAM"
mv "$TMP" "$OUT"
bcftools index "$OUT"
echo "DONE $S: $(bcftools view -H "$OUT" | wc -l) deletion candidates, $(bcftools view -H -f PASS "$OUT" | wc -l) PASS"

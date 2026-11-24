#!/usr/bin/env bash
# Stream the first n_pairs read pairs of each ENA run (~15x coverage) to save disk.
# Usage: bash scripts/download_reads.sh config/samples.tsv data/meta/selected_runs.tsv data/raw
set -euo pipefail
samples=$1; runs=$2; outdir=$3
mkdir -p "$outdir"
tail -n +2 "$samples" | while IFS=$'\t' read -r sample _ _ _ run npairs _; do
  urls=$(awk -F'\t' -v r="$run" '$1==r{print $4}' "$runs")
  i=1
  for u in ${urls//;/ }; do
    out="$outdir/${sample}_R${i}.fastq.gz"
    if [ -s "$out.done" ]; then echo "skip $out"; i=$((i+1)); continue; fi
    echo "$(date +%T) $sample R$i <- $u ($npairs reads)"
    (curl -sS --retry 5 "https://$u" | gzip -dc 2>/dev/null || true) \
      | head -n $((npairs*4)) | gzip -1 > "$out"
    lines=$(gzip -dc "$out" | wc -l)
    [ "$lines" -eq $((npairs*4)) ] || { echo "ERROR $out has $lines lines"; exit 1; }
    echo "$lines" > "$out.done"
    i=$((i+1))
  done
done
echo "$(date +%T) all downloads finished"

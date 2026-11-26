#!/usr/bin/env bash
# Usage: bash scripts/run_plotsr.sh
# Figure 24: synteny and rearrangements across 5 genomes with plotsr, env `syri`.
# Needs the 4 chained SyRI runs first (each genome against the next one in the order below):
#   bash scripts/run_syri.sh Col-0 An-1; bash scripts/run_syri.sh An-1 Kyoto
#   bash scripts/run_syri.sh Kyoto Ler-1; bash scripts/run_syri.sh Ler-1 Cvi-0
# Order follows the SNP tree (fig 32), Cvi-0 (most distant) last. Each band compares neighbouring genomes,
# so a rearrangement between Ler-1 and Cvi-0 is shown there, not against Col-0.
# Output: results/figures/24_validation_rearrangements__synteny_riparian_plot.png
set -euo pipefail

ORDER=(Col-0 An-1 Kyoto Ler-1 Cvi-0)
G=results/syri/genomes.txt
printf "#file\tname\ttags\n" > "$G"
for s in "${ORDER[@]}"; do printf "data/syri/%s.fa\t%s\tlw:1.5\n" "$s" "$s" >> "$G"; done
SR=()
for i in 0 1 2 3; do
  f=results/syri/${ORDER[$i]}_${ORDER[$((i + 1))]}.syri.out
  [ -f "$f" ] || { echo "missing $f"; exit 1; }
  SR+=(--sr "$f")
done
plotsr "${SR[@]}" --genomes "$G" -o results/figures/24_validation_rearrangements__synteny_riparian_plot.png -H 8 -W 10 -f 9 \
       -s 10000 > results/syri/plotsr.log 2>&1
echo "wrote results/figures/24_validation_rearrangements__synteny_riparian_plot.png"

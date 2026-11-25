#!/usr/bin/env bash
# Usage: bash scripts/annotate_snpeff.sh
# Phase 5: predict the effect of every filtered variant on genes with SnpEff, run in env `variants`.
# Database: SnpEff `Arabidopsis_thaliana` (v5_2 build, TAIR10 gene models; downloaded 2026-11-25).
# The database names chromosomes 1-5, the VCF names them NC_...; they are renamed for SnpEff and
# renamed back afterwards, so the output uses the same names as every other file in the project.
# Output: results/variants/annotated.vcf.gz (+ .tbi), snpeff_summary.html, snpeff_summary.genes.txt
# Refuses to overwrite an existing output.
set -euo pipefail

IN=results/variants/filtered.vcf.gz
OUT=results/variants/annotated.vcf.gz
TMP=results/variants/annotated.tmp.vcf.gz
MAP=results/variants/snpeff_chrom_map.txt
BACK=results/variants/snpeff_chrom_map_back.txt

[ -f "$IN" ] || { echo "missing $IN"; exit 1; }
[ ! -e "$OUT" ] || { echo "$OUT already exists, not overwriting"; exit 1; }

echo "$(date) start snpEff"
bcftools annotate --rename-chrs "$MAP" "$IN" \
  | snpEff -Xmx2g -stats results/variants/snpeff_summary.html Arabidopsis_thaliana \
  | bcftools annotate --rename-chrs "$BACK" -Oz -o "$TMP"
mv "$TMP" "$OUT"
bcftools index -t "$OUT"
echo "$(date) DONE: $(bcftools view -H "$OUT" | wc -l) variants annotated," \
     "$(bcftools view -H "$OUT" | grep -c '|HIGH|') with at least one HIGH-impact effect"

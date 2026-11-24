#!/usr/bin/env bash
# Usage: bash scripts/call_variants.sh CHROM
# Joint SNP/indel calling of all 8 accessions on ONE chromosome (Phase 5), run in env `variants`.
#   bcftools mpileup : stacks the reads of the 8 BAMs at each position
#   bcftools call -m : multiallelic caller, keeps variant sites only (-v)
# Parameters (scientific choices, change them only on purpose):
#   -q 10   ignore reads with MAPQ < 10         (same cut-off as the PAV calling)
#   -Q 20   ignore bases with base quality < 20
#   --ploidy 1  the accessions are inbred lines, so one allele per site is assumed
# Duplicate-marked reads are skipped by mpileup by default. Organelles are not called (chr1-5 only).
# Output: results/variants/raw/CHROM.vcf.gz (+ .tbi). Refuses to overwrite an existing file.
set -euo pipefail

CHR=$1
REGION=${2:-$CHR}   # optional 2nd argument, e.g. NC_003075.7:1-500000, only for a quick test
REF=data/reference/TAIR10.1_Col-0/GCF_000001735.4/GCF_000001735.4_TAIR10.1_genomic.fna
if [ "$REGION" = "$CHR" ]; then
  OUTDIR=results/variants/raw
  OUT=$OUTDIR/$CHR.vcf.gz
else
  OUTDIR=results/variants/test   # test runs never touch the real output
  OUT=$OUTDIR/$CHR.test.vcf.gz
fi
BAMLIST=results/variants/bam_list.txt

case "$CHR" in
  NC_003070.9|NC_003071.7|NC_003074.8|NC_003075.7|NC_003076.8) ;;
  *) echo "CHROM must be one of the five nuclear chromosomes, got: $CHR"; exit 1;;
esac
[ -f "$REF.fai" ] || { echo "missing $REF.fai"; exit 1; }
[ ! -e "$OUT" ] || { echo "$OUT already exists, not overwriting"; exit 1; }

mkdir -p "$OUTDIR"
if [ ! -f "$BAMLIST" ]; then
  for S in Col-0 An-1 Kyoto Kondara Cvi-0 Ler-1 Se-0 Ale-Stenar-44-4; do
    [ -f "data/bam/$S/$S.markdup.bam" ] || { echo "missing BAM of $S"; exit 1; }
    echo "data/bam/$S/$S.markdup.bam"
  done > "$BAMLIST"
fi

echo "$(date) start calling $CHR"
bcftools mpileup -f "$REF" -r "$REGION" -b "$BAMLIST" -q 10 -Q 20 -a FORMAT/AD,FORMAT/DP -Ou \
  | bcftools call -m -v --ploidy 1 -Oz -o "$OUT.tmp"
mv "$OUT.tmp" "$OUT"
bcftools index -t "$OUT"
echo "DONE $CHR: $(bcftools view -H "$OUT" | wc -l) variant sites"

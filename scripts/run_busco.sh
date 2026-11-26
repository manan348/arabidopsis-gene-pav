#!/usr/bin/env bash
# Usage: bash scripts/run_busco.sh SAMPLE      (SAMPLE = Col-0, An-1, Cvi-0, Kyoto or Ler-1)
# Phase 6: gene-set completeness of one genome with BUSCO 6.1.0, run in env `asmqc`.
# Lineage brassicales_odb12 (7,083 conserved single-copy genes), downloaded to data/busco_downloads (--offline).
# Genome mode (-m genome); BUSCO 6 uses miniprot for eukaryote genomes. Col-0 = the TAIR10.1 reference itself,
# used as the baseline: the other assemblies are judged against what Col-0 scores.
# Output: results/assembly/busco/SAMPLE/ (short_summary*.txt). Refuses to overwrite an existing result.
set -euo pipefail

S=$1
case "$S" in
  Col-0) FA=data/reference/TAIR10.1_Col-0/GCF_000001735.4/GCF_000001735.4_TAIR10.1_genomic.fna ;;
  An-1)  FA=data/assemblies/An-1/An-1.fna ;;
  Cvi-0) FA=data/assemblies/Cvi-0/Cvi-0.fna ;;
  Kyoto) FA=data/assemblies/Kyo/Kyoto.fna ;;
  Ler-1) FA=data/assemblies/Ler/Ler_GCA_902460285.1.fna ;;
  *) echo "unknown sample $S (use Col-0, An-1, Cvi-0, Kyoto or Ler-1)"; exit 1 ;;
esac
LIN=data/busco_downloads/lineages/brassicales_odb12
OUTDIR=results/assembly/busco

[ -f "$FA" ] || { echo "missing $FA"; exit 1; }
[ -d "$LIN" ] || { echo "missing $LIN (busco --download brassicales_odb12 --download_path data/busco_downloads)"; exit 1; }
[ ! -e "$OUTDIR/$S" ] || { echo "$OUTDIR/$S already exists, not overwriting"; exit 1; }

mkdir -p "$OUTDIR"
echo "$(date) start busco $S"
busco -i "$FA" -m genome -l "$LIN" -c 3 -o "$S" --out_path "$OUTDIR" --offline \
      --download_path data/busco_downloads
echo "$(date) DONE $S: $(grep -h 'C:' "$OUTDIR/$S"/short_summary*.txt | head -n 1 | tr -s ' \t' ' ')"

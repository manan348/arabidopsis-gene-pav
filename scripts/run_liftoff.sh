#!/usr/bin/env bash
# Usage: bash scripts/run_liftoff.sh SAMPLE      (SAMPLE = An-1, Cvi-0, Kyoto or Ler-1)
# Phase 6: place the Col-0 genes on one assembly with Liftoff, run in env `liftoff` (Liftoff 1.6.3, minimap2 2.24).
# A gene that Liftoff cannot place is absent from the assembly = independent check of the read-based absent calls.
# Settings: -a 0.1 (a gene counts as found if >= 10% of it aligns, the same cut-off as the read-based absent call,
#           breadth < 0.1), -s 0.5 (default identity), -p 3 threads. Coverage is written for every found gene,
#           so stricter cut-offs can be applied afterwards.
# Assembly of Ler-1: the NCBI "Ler" assembly (GCA_902460285.1); NCBI does not say Ler-0 or Ler-1.
# Output: results/assembly/liftoff/SAMPLE.gff, SAMPLE.unmapped.txt, SAMPLE_run.log; intermediate files in SAMPLE_tmp/.
# Compressed assemblies are unpacked once to data/assemblies/<dir>/SAMPLE.fna (Liftoff needs plain FASTA).
# Refuses to overwrite an existing output.
set -euo pipefail

S=$1
REF=data/reference/TAIR10.1_Col-0/GCF_000001735.4/GCF_000001735.4_TAIR10.1_genomic.fna
GFF=data/reference/genes_only.gff
case "$S" in
  An-1)  DIR=An-1; SRC=data/assemblies/An-1/An-1_GCA_902460265.3.fna.gz ;;
  Cvi-0) DIR=Cvi-0; SRC=data/assemblies/Cvi-0/Cvi-0_GCA_902460275.1.fna.gz ;;
  Kyoto) DIR=Kyo; SRC=data/assemblies/Kyo/Kyo_GCA_902460305.1.fna.gz ;;
  Ler-1) DIR=Ler; SRC=data/assemblies/Ler/Ler_GCA_902460285.1.fna ;;
  *) echo "unknown sample $S (use An-1, Cvi-0, Kyoto or Ler-1)"; exit 1 ;;
esac
OUTDIR=results/assembly/liftoff
OUT=$OUTDIR/$S.gff
TMPDIR=$OUTDIR/${S}_tmp

[ -f "$SRC" ] || { echo "missing $SRC"; exit 1; }
[ -f "$GFF" ] || { echo "missing $GFF (run python scripts/make_liftoff_gff.py)"; exit 1; }
[ ! -e "$OUT" ] || { echo "$OUT already exists, not overwriting"; exit 1; }

TARGET=$SRC
if [[ "$SRC" == *.gz ]]; then
  TARGET=data/assemblies/$DIR/$S.fna
  [ -f "$TARGET" ] || zcat "$SRC" > "$TARGET"
fi

mkdir -p "$OUTDIR"
echo "$(date) start liftoff $S"
liftoff -g "$GFF" -o "$OUTDIR/$S.tmp.gff" -u "$OUTDIR/$S.unmapped.txt" -dir "$TMPDIR" \
        -a 0.1 -s 0.5 -p 3 "$TARGET" "$REF"
mv "$OUTDIR/$S.tmp.gff" "$OUT"
echo "$(date) DONE $S: $(awk -F'\t' '$3=="gene"' "$OUT" | wc -l) genes placed, $(wc -l < "$OUTDIR/$S.unmapped.txt") unmapped"

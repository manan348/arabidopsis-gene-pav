#!/usr/bin/env bash
# Usage: bash scripts/run_syri.sh REF QRY      (names of files in data/syri/, e.g. Col-0 An-1)
# Phase 7: whole-genome alignment (minimap2 -ax asm5 --eqx) and structural rearrangements (SyRI), env `syri`.
# Inputs: data/syri/REF.fa, data/syri/QRY.fa (Chr1-Chr5 only, made by scripts/prep_syri.py)
# Outputs: results/syri/REF_QRY.bam (alignment, not tracked by git), results/syri/REF_QRY.syri.out (+ .summary)
# Memory: minimap2 needs about 2-3 GB here; do not run while BUSCO or another large job is running.
set -euo pipefail

R=$1; Q=$2
RF=data/syri/$R.fa; QF=data/syri/$Q.fa
OUT=results/syri; P=${R}_${Q}
[ -f "$RF" ] && [ -f "$QF" ] || { echo "missing $RF or $QF (run scripts/prep_syri.py)"; exit 1; }
[ ! -e "$OUT/$P.syri.out" ] || { echo "$OUT/$P.syri.out exists, not overwriting"; exit 1; }
mkdir -p "$OUT"

echo "$(date) start minimap2 $R vs $Q"
if [ ! -s "$OUT/$P.bam" ]; then
  minimap2 -ax asm5 --eqx -t 3 "$RF" "$QF" 2> "$OUT/$P.minimap2.log" | samtools view -b -o "$OUT/$P.bam.tmp" -
  mv "$OUT/$P.bam.tmp" "$OUT/$P.bam"
fi
echo "$(date) start syri $R vs $Q"
syri -c "$OUT/$P.bam" -r "$RF" -q "$QF" -F B --dir "$OUT" --prefix "$P." --nc 1 > "$OUT/$P.syri.log" 2>&1
echo "$(date) DONE $P: $(cut -f11 "$OUT/$P.syri.out" | sort | uniq -c | awk '$2 ~ /^(SYN|INV|TRANS|INVTR|DUP|INVDP)$/ {printf "%s %s  ", $2, $1}')"

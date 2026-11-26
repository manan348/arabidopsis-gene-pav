#!/usr/bin/env bash
# Usage: bash scripts/run_quast.sh
# Phase 6: assembly statistics for the 4 assemblies with QUAST 5.3.0, run in env `asmqc`.
# Col-0 (TAIR10.1) is given as reference (-r), so QUAST also reports the share of Col-0 covered (genome fraction)
# and "misassemblies". Careful: these accessions really differ from Col-0 (inversions, translocations), so
# many "misassemblies" are probably real structural differences, not assembly errors.
# Uses -t 3 threads. Output: results/assembly/quast/ (report.txt, report.tsv, report.html). Refuses to overwrite.
set -euo pipefail

REF=data/reference/TAIR10.1_Col-0/GCF_000001735.4/GCF_000001735.4_TAIR10.1_genomic.fna
OUT=results/assembly/quast
[ ! -e "$OUT" ] || { echo "$OUT already exists, not overwriting"; exit 1; }
for f in data/assemblies/An-1/An-1.fna data/assemblies/Cvi-0/Cvi-0.fna data/assemblies/Kyo/Kyoto.fna \
         data/assemblies/Ler/Ler_GCA_902460285.1.fna; do
  [ -f "$f" ] || { echo "missing $f (run scripts/run_liftoff.sh first, it unpacks the assemblies)"; exit 1; }
done

echo "$(date) start quast"
quast.py -t 3 -r "$REF" -o "$OUT" -l "An-1,Cvi-0,Kyoto,Ler-1" \
  data/assemblies/An-1/An-1.fna data/assemblies/Cvi-0/Cvi-0.fna data/assemblies/Kyo/Kyoto.fna \
  data/assemblies/Ler/Ler_GCA_902460285.1.fna
echo "$(date) DONE quast"

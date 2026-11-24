#!/usr/bin/env bash
# Usage: bash scripts/map_sample.sh SAMPLE
# Maps data/reads_final/SAMPLE to TAIR10 with bowtie2, then fixmate -> sort -> markdup -> index,
# and writes flagstat, samtools stats and mosdepth results.
# Output BAM: data/bam/SAMPLE/SAMPLE.markdup.bam   QC files: results/mapping/SAMPLE/
set -euo pipefail

S=$1
THREADS=3
IDX=data/reference/index/TAIR10
FIN=data/reads_final/$S
BAMDIR=data/bam/$S
QC=results/mapping/$S
BAM=$BAMDIR/$S.markdup.bam

[ -f "$FIN/${S}_R1.fastq.gz" ] || { echo "missing $FIN/${S}_R1.fastq.gz"; exit 1; }
[ -f "$FIN/${S}_R2.fastq.gz" ] || { echo "missing $FIN/${S}_R2.fastq.gz"; exit 1; }
[ -f "${IDX}.1.bt2" ] || [ -f "${IDX}.1.bt2l" ] || { echo "bowtie2 index not found: $IDX"; exit 1; }
# never overwrite an existing result: delete it yourself if you really want to redo it
[ ! -e "$BAM" ] || { echo "$BAM already exists, not overwriting"; exit 1; }

mkdir -p "$BAMDIR" "$QC"
TMPBAM=$BAMDIR/$S.tmp.bam
SORTTMP=$BAMDIR/sorttmp
mkdir -p "$SORTTMP"

echo "$(date) start mapping $S"

# 1) bowtie2 (end-to-end, very sensitive) -> 2) fixmate -m (adds mate score tags needed by markdup)
# -> 3) sort by coordinate -> 4) markdup (marks duplicates, does not remove them)
# "-u" = uncompressed BAM between steps (faster, nothing written to disk in between).
# Memory: bowtie2 index about 0.3 GB + sort 2 threads x 600 MB; fits in 3.8 GB.
bowtie2 --very-sensitive -p "$THREADS" -x "$IDX" \
    -1 "$FIN/${S}_R1.fastq.gz" -2 "$FIN/${S}_R2.fastq.gz" \
    --rg-id "$S" --rg "SM:$S" --rg "PL:ILLUMINA" \
    2> "$QC/$S.bowtie2.log" \
  | samtools fixmate -m -u -@ 1 - - \
  | samtools sort -u -@ 1 -m 600M -T "$SORTTMP/$S" - \
  | samtools markdup -@ 1 -s -f "$QC/$S.markdup.txt" - "$TMPBAM"

mv "$TMPBAM" "$BAM"
rmdir "$SORTTMP" 2>/dev/null || true
samtools index -@ 2 "$BAM"
echo "$(date) BAM written and indexed"

# QC tables (MultiQC reads these)
samtools flagstat -@ 2 "$BAM" > "$QC/$S.flagstat.txt"
samtools stats -@ 2 "$BAM" > "$QC/$S.stats.txt"

# genome-wide depth per chromosome (-n = no per-base file; summary + distribution only)
mosdepth -t "$THREADS" -n "$QC/$S" "$BAM"

# a few plain numbers for the log
# flags: 4 = unmapped, 256 = secondary, 2048 = supplementary
PRIMARY=$(samtools view -c -F 2304 "$BAM")   # all reads (mapped or not), one line per read
MAPPED=$(samtools view -c -F 2308 "$BAM")    # mapped reads only
MAPQ10=$(samtools view -c -F 2308 -q 10 "$BAM")  # mapped with MAPQ >= 10
echo "$S primary_reads=$PRIMARY mapped=$MAPPED mapped_MAPQ>=10=$MAPQ10"
echo "--- bowtie2 summary ---"
tail -n 8 "$QC/$S.bowtie2.log"
echo "--- mosdepth summary ---"
cat "$QC/$S.mosdepth.summary.txt"
echo "DONE $S"

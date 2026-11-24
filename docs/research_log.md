# Research log: Arabidopsis thaliana gene presence/absence (PAV) pan-genome

One entry per working day: what was done, what was found, what is next.
Numbers come from the result files in `results/`; the full state is in `PROJECT_PLAN.md`.
Entries before 2026-11-25 were written afterwards from the project notes.

## 2026-11-22
- Decided to trim first and subsample afterwards (trim with fastp on the full reads, then subsample to about 15x of trimmed bases).

## 2026-11-23
- Phase 2 finished: trimming, subsampling and FastQC/MultiQC for all 8 accessions. Final depth is 11.2x to 15.1x.
- Phase 3 started: bowtie2 mapping to TAIR10.1 as one background queue.

## 2026-11-24
- Phase 3 finished: all 8 accessions mapped. Mapping rate 94.0-96.9%, MAPQ >= 10 share 68-78%, duplicates 2.7-18.2%, median depth of chr1-5 7-12x. No sample failed.
- MultiQC on the mapping results. Breadth at >= 1x is 100% for Col-0 and 93-95% for the other seven.
- Phase 4 first pass: per-gene depth and breadth (mosdepth) and PAV calls for 27,444 genes.
  - Col-0 control: 0 absent genes.
  - Other accessions: 173-452 absent genes each (0.6-1.7%).
  - 24,259 genes are present in all 8; 911 genes are absent in at least one accession; 18 genes are absent in all 7 other accessions.
  - Absent counts change by about 5-15% when the absent cut-off moves from 5% to 20% covered, and the order of accessions does not change.
- Phase 5 started: installed bcftools 1.25, SnpEff 5.4c and Delly 0.7.6 in the conda env `variants`. Test of SNP calling on 500 kb of chr4 gave 8,483 variant sites. Full run for chr1-5 started at 23:40.

## 2026-11-25
- Started DVC with DagsHub storage for the small results (reads and BAM files stay local).

## Open points
- Some absent genes may be genes that are too divergent from Col-0 for reads to map, not real deletions. To be checked with SNP density (Phase 5), Delly deletions (Phase 5) and the assemblies (Phase 6).

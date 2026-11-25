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
- Put the project on GitHub (`arabidopsis-gene-pav`) with the scripts, plan and this log.
- Set up DVC with DagsHub storage. The small results (fastp, figures, mapping QC, MultiQC, PAV tables, FastQC) are tracked with DVC and pushed; reads and BAM files stay local.
- SNP calling with bcftools finished for chr1-5 (01:32): 2,046,654 raw sites (1,794,047 SNPs, 252,607 indels), Ts/Tv 1.21.
- Low-quality sites have almost the same Ts/Tv (1.14) as the rest (1.22), so they are not mostly errors. Chose a mild filter (QUAL >= 20, total depth 20-130): 1,850,079 sites kept (90%).
- Delly 0.7.6 crashed (segmentation fault); installed Delly 2.8.2 in a separate env for the deletion calling. In Delly 2.x the command `call` is named `sr`; fixed in `scripts/call_sv.sh`.
- Variant figures (fig4, fig5): SNPs per accession 472k (An-1) to 634k (Cvi-0), Col-0 control has 832 non-reference sites. SNP density around absent genes is not higher than around present genes (median 3.85 vs 3.65 per kb), so there is no strong sign that the absent calls are mostly reads failing to map in divergent regions.
- Delly deletions finished for all 8 (5-12 min each). PASS deletions: Col-0 5, An-1 300, Kyoto 47, Kondara 532, Cvi-0 556, Ler-1 1,370, Se-0 226, Ale-Stenar-44-4 1,172.
- Compared the deletions (up to 100 kb; 97 longer ones dropped) with the gene calls (fig6): absent genes lie inside a deletion 8 to 76 times more often than present genes (Se-0 8x, Kyoto 12x, An-1 19x, Kondara 39x, Ale-Stenar-44-4 40x, Cvi-0 53x, Ler-1 76x). The share of absent genes confirmed rises with the number of deletions Delly found, so a low share mostly reflects Delly finding few deletions. Two independent methods agree.
- SnpEff: downloaded the `Arabidopsis_thaliana` database (v5_2 build) after asking, and annotated all 1,850,079 filtered variants (7 min). 10,215 have a HIGH-impact effect (stop gained, frameshift, splice, start or stop lost).
- Present but broken genes (fig7): 2,070 (An-1) to 2,731 (Cvi-0) genes per accession are called present but carry a HIGH-impact variant, against 76 in the Col-0 control. That is 6 to 13 times more than the number of absent genes (173-452). HIGH effects are predictions; frameshifts from short-read indel calls may include errors. 46 genes with HIGH variants are not in the 27,444-gene PAV list because SnpEff and RefSeq gene models differ.
- Phase 5 is finished.
- Pan-genome figures from the existing calls (Phase 7, early): of the 911 variable genes, 389 are absent in only one accession (145 of them in Cvi-0) and 18 in all 7 (fig8). The core-gene curve is still falling at 8 accessions (24,259 genes called present in all; 26,533 not called absent in any), so the core of this small set is an upper estimate (fig9). A SNP distance map puts Cvi-0 clearly apart from the other accessions (fig10).
- Next: Phase 6 (assemblies: QUAST, BUSCO, Liftoff; ask before downloads). First confirm accession names (Ler-1 vs Ler, Kyoto vs Kyo).

## Open points
- Some absent genes may be genes that are too divergent from Col-0 for reads to map, not real deletions. To be checked with SNP density (Phase 5), Delly deletions (Phase 5) and the assemblies (Phase 6).

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

## 2026-11-26
- Pushed the Phase 5 results: DVC (244 files to DagsHub) and git. PROJECT_PLAN.md is no longer tracked on GitHub (kept locally as working notes).
- Phase 6 started. Checked the 4 assemblies on NCBI: An-1, Cvi-0 and Kyo (= Kyoto) match the read samples; the Ler assembly is only labelled "Ler", so it is used as a stand-in for Ler-1.
- Installed Liftoff 1.6.3 and minimap2 2.24. A first Liftoff run with all transcripts, exons and CDS (about 668,000 features) ran out of memory on this machine, so Liftoff was run on gene spans only (27,444 genes), which matches how the read-based calls were made. Each assembly took 5-11 minutes.
- Liftoff could not place 311-442 genes per assembly (An-1 356, Cvi-0 442, Kyoto 311, Ler-1 318 read-absent vs 323 not placed).
- Comparison with the read-based calls (fig11): 89-94% of the genes the reads call absent are also not found intact in the assembly (missing, less than half aligned, or only a look-alike elsewhere); 22-43 per accession are found intact and are the likely errors. Of the genes missing from an assembly, only 2-7% were called present by the reads; most others were called absent, uncertain or ambiguous. Ler-1 behaves like the other three, so there is no sign that the "Ler" assembly is a different line. The cut-offs (coverage 0.5, 2 Mb window) are choices, and "placed elsewhere" can also be a real rearrangement.
- Reads vs assembly for every gene (fig12): about 26,500-26,700 genes per accession are covered by reads and found in the assembly; 306-441 have neither. Three methods together (fig13): of the read-absent genes, only 21-41 per accession are supported by neither the assembly nor Delly (likely read errors). Genes missing from an assembly are 10-47 times more often inside a Delly deletion than genes the assembly has; for Ler-1 the agreement is highest (19.8% vs 0.42%), which supports using the "Ler" assembly for the Ler-1 reads.
- Position of absent genes (fig15): 21-27% of the confirmed absent genes lie within 2 Mb of a centromere, against 4.8% of all genes (about 5 times more; centromere positions approximate, not tested). The likely read errors often recur at the same positions in several accessions, so a few specific genes are hard for the read method.
- Example gene (fig16): FAR4 (AT3G44540, chromosome 3) is missing only in Cvi-0. Read depth drops to zero exactly over a 7.7 kb Delly deletion, the other 7 accessions have normal depth over the gene, and Liftoff cannot place it in the Cvi-0 assembly.
- QUAST on the 4 assemblies: 120-122 Mb, N50 23-24 Mb, 90% of each genome in 5 pieces, about 90% of Col-0 covered (89.1% for Cvi-0, the most divergent). The 3,500-4,600 "misassemblies" against Col-0 are mostly real structural differences between accessions, not assembly errors.
- Downloaded the BUSCO lineage brassicales_odb12 (7,083 genes) after asking.
- Installed QUAST 5.3.0 and BUSCO 6.1.0 (env `asmqc`; QUAST needed setuptools because of Python 3.13).
- Orientation check of the 4 assemblies against Col-0: 0-4% of each chromosome is reversed (Chr4 4% in all four), so no chromosome needs flipping; the small shares are local inversions.
- Genes lost, broken or gained (fig17): lost 173-452 and broken 2,070-2,731 per accession are solid. Read-depth "gains" (21-452 per accession at 2x depth) rise with each sample's depth noise (Ler-1 lowest, Ale-Stenar-44-4 highest), so they are only candidates (results/pav/gene_gains.tsv).
- GO enrichment (R clusterProfiler 4.18.4, org.At.tair.db 3.22.0; background 27,444 genes, BH FDR < 0.05) (fig18, fig19, fig25). Absent genes (911): ADP binding, i.e. NLR disease-resistance genes (52 genes, FDR 3e-36), TIR-domain NLRs (10), defense response (29), RNase H (17, probably transposon remnants) and sesquiterpene synthases (5). A 14-gene "fertilization / cell surface" term comes from small unnamed genes, 8 of them clustered on Chr5 near 14.3 Mb; treated with caution. Broken genes (3,658 broken in >= 2 accessions): only ADP binding (FDR 5e-7) and RNase H. No cellular-component term.
- BUSCO run for Col-0 was killed by a laptop restart; the restart option failed (half-written HMMER folder), so the folder was deleted and the run started again. Col-0 (TAIR10.1): complete 100.0% (single 98.8%, duplicated 1.2%), fragmented 0%, missing 0% of 4,311 genes. This is the baseline for the 4 assemblies (running).
- Pan-genome categories (fig23): core (all 8) 26,533 genes (96.7%), soft-core (7) 389, shell (2-6) 463, private (Col-0 only) 18. Genes absent from Col-0 cannot be seen with this method, so shell and private are lower bounds.
- Circular cluster plot of the 911 variable genes by absence pattern, with GO groups (fig24); workflow diagram (fig21).
- Map and SNP tree (fig22): coordinates of the 8 accessions from the 1001 Genomes AraPheno API (data/meta/accession_coords.tsv); world outline Natural Earth 1:110m (downloaded after asking, not tracked). In the UPGMA tree Cvi-0 splits first, then Kondara; Kyoto groups with Ler-1.
- Synteny with SyRI 1.8.2 and plotsr 1.2.0 (fig20): minimap2 asm5 alignments chained Col-0 -> An-1 -> Kyoto -> Ler-1 -> Cvi-0, chromosomes 1-5 only. 92-95% of each genome is syntenic with its neighbour; inversions 0.15-1.75 Mb per pair. The largest is a Chr4 inversion at about 1.6-2.8 Mb between Col-0 and An-1 (fits the known Col-0 inversion on the short arm of Chr4). Rearrangements and duplications cluster near the centromeres, like the absent genes in fig15.
- Interactive PAV explorer (results/interactive/pav_explorer.html): 911 genes x 8 accessions with read depth, Delly and Liftoff on hover, filters by chromosome, GO group and number of accessions.
- Next: BUSCO for An-1, Cvi-0, Kyoto and Ler-1 (running), then the assembly quality figure (fig14).

## 2026-11-27
- Figure numbers below are the new ones. All figures were renumbered 01-39 in pipeline order, with the plot type in the file name and a caption band on each image; the old numbers in earlier entries no longer match. The list is in docs/figures.md.
- BUSCO finished for all 4 assemblies (brassicales_odb12, 4,311 genes; one run was lost when WSL closed and was restarted). An-1, Cvi-0 and Kyoto are 100% complete, Ler-1 99.9% (4 missing), the same as Col-0. Together with QUAST this is fig 17: the assemblies are complete, so a gene Liftoff cannot find is really missing.
- Length bias (fig 21): short genes are called absent more often, but they are confirmed by the assembly or Delly just as often (91-95%, Wilson 95% CI), so the bias is not mainly false calls.
- Gene loss vs SNP relatedness (fig 33): the PAV distances and the SNP distances agree, Mantel r = 0.79 (p = 0.02, exact permutations) for all 8 and r = 0.90 (p = 0.007) without Col-0.
- NLR genes (fig 37): 151 NLRs found from the RefSeq product names (results/pav/nlr_genes.tsv). 50 are variable (33%) against 3.2% of other genes, odds ratio 15.2, Fisher's exact p = 4e-36; 82% of the variable NLRs sit in clusters.
- The Chr5 "fertilization" GO cluster (fig 39) is a multi-copy region with uneven read depth, so that GO term stays flagged as uncertain.
- Nucleotide diversity (fig 15; results/variants/pi_windows.tsv, pi_genes.tsv): genome-wide pi about 0.0051 per bp, a lower bound because only filtered SNPs in covered sites are counted. Diversity is highest around the centromeres. NLR genes are 2.1 times more diverse than other genes (Mann-Whitney p = 2.2e-12).
- pN/pS with Nei-Gojobori site counts (fig 38; results/variants/pnps_genes.tsv): core genes 0.249 (95% bootstrap CI 0.245-0.252), NLR genes 0.494 (0.451-0.543), variable non-NLR genes 0.601 (0.564-0.642). NLR vs core permutation p = 1e-4. All values are below 1, so this shows weaker constraint or diversifying selection, not proof of positive selection.
- Circos overview of the genome (fig 02): gene density, SNPs, pi, absent genes, deletions and NLRs on one plot.
- Next: update the scripts to the new figure names, write the README.

## Open points
- Some absent genes may be genes that are too divergent from Col-0 for reads to map, not real deletions. Checked: SNP density, Delly and the assemblies support 89-95% of the calls; 21-41 per accession remain likely read errors.
- Gene gains and genes absent from Col-0 cannot be seen with the read-depth method (Liftoff with copies could test gains).

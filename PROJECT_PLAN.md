# Plan: Arabidopsis thaliana Gene Presence/Absence (PAV) Pan-genome Project

## CURRENT STATUS (updated 2026-11-24): read this first
This section overrides any path or step below that disagrees with it. The sections below are the original plan.

### How we work (important for a new chat)
- The USER runs every command themselves in WSL Ubuntu (`/mnt/d/genom_project`) and pastes the output. Claude gives ONE small block of simple, explained commands at a time and interprets the output in plain language. Claude does not install or run tools itself.
- Do not download or delete anything without asking. Prefer `grep`/`head` over `cat` on big files (a past `cat id_names.tsv` printed about 1,000 lines).
- Be honest about uncertainty; correct earlier mistakes openly.
- Lessons from Phase 2: always `cd /mnt/d/genom_project` first (a past check failed because the user was inside `data/`); the user mistyped folder names twice (`result/` vs `results/`, `fill` vs `full`), so give copy-paste blocks and use `mkdir -p` with the exact names; a `nohup ... &` line must be run ONCE (it was run twice by mistake; harmless here because the script is deterministic) and it returns the prompt at once with `[1] PID`, which is normal; logs written by the script stay empty until the step ends.
- Terminology: call the 8 lines "accessions" (natural inbred lines of A. thaliana, from the 1001 Genomes Project), not species or cultivars.
- STYLE (user concern about long answers, 2026-11-24): keep replies SHORT. Default = one command block, 2-4 lines on what it does, then the verdict on the numbers. No big tables, "what to look for" lists or next-step previews unless the user asks; give long explanations only when asked "why" or "explain". Read only the needed lines of big files (e.g. `samtools stats` files are mostly per-cycle tables).
- The user asks many concept questions (what is mapping, BAM, MAPQ, BUSCO, Liftoff). Answer them in plain language, briefly, with an analogy if it helps.
- Lessons from Phase 3: when giving a shell block, check that it does what it says before sending (two earlier progress blocks failed: one picked the `sh -c` wrapper PID, one assumed the wrong order of the mosdepth dist file). Do not claim a time estimate is reliable: mapping estimates were off by 10-30%.

### Machine
WSL Ubuntu, miniconda3, 4 cores (use 3 threads), 3.8 GB RAM, 1 GB swap. D: has about 177 GB free. conda env in use: `map` (`conda activate map`). Env `fingerprint-panel` belongs to another project; ignore.

### Tools installed (all in env `map`)
bowtie2 2.5.5, bwa, samtools 1.24, bedtools 2.31.1, seqkit, mosdepth 0.3.14, fastqc 0.12.1, fastp 1.0.1, multiqc 1.35, seqtk 1.4.
NOT installed yet: snakemake and dvc (Linux side), bcftools, snpEff, delly, minimap2, liftoff, quast, busco, python/R analysis packages (pandas, matplotlib, upsetplot, scikit-learn, R clusterProfiler). Install each set when its phase starts; give the user the `conda create/install -c conda-forge -c bioconda --override-channels` command and check `--version`. No git repo or DVC yet.

### Done
- Downloads: 8 accessions of full reads (about 24 GB, md5 verified), TAIR10.1 reference, 4 assemblies (An-1, Cvi-0, Kyo, Ler), metadata. See `config/samples.tsv` (column n_pairs is the old subsample target).
- **Phase 1 (reference): DONE.** md5 OK; 7 sequences (chr1-5 = 119,146,348 bp, chrM NC_037304.1, chrC NC_000932.1; total 119,668,634). Bowtie2 index: `data/reference/index/TAIR10` (prefix `TAIR10`). `data/reference/genes.bed` = 27,444 protein-coding genes on chr1-5 only, named by locus_tag (BED is 0-based). 125 organelle genes excluded (27,569 total). `data/reference/gene_names.tsv` = locus_tag and gene name. BWA index NOT built yet (only needed for the Cvi-0 comparison).
- **Phase 2 (reads and QC): DONE (2026-11-23). Trim + subsample + FastQC + MultiQC are done and interpreted for all 8 accessions. See "Phase 2 results (2026-11-23)" below the Decision block. The bullets under this heading are the older trial notes.**
  - FastQC on raw Col-0 and Ale-Stenar-44-4 only: `results/qc_raw/`.
    - Col-0 (100 bp, about 17.3M pairs, about 29x): quality good to about base 45, then decays; R2 worse (mean about 17 at base 100). Old HiSeq run. GC 38%.
    - Ale-Stenar-44-4 (76 bp, about 24.2M pairs, about 31x): passes per-base quality; adapter in about 0.3% of reads; Q2 tails in about 10% of reads from base about 52; GC 40%, duplication and overrepresented warnings (adapters, not contamination).
    - The other 6 accessions have NOT been FastQC'd.
  - Col-0 trial subsample: `data/reads_15x/Col-0/` (seqtk sample -s 100, fraction 0.519, 9,000,710 pairs, names match R1/R2).
  - Col-0 trial trim of that subset (fastp, `--detect_adapter_for_pe --cut_right --cut_window_size 4 --cut_mean_quality 20 --length_required 36 -w 3`): `data/reads_trimmed/Col-0/`, report `results/fastp/Col-0.*`. Result: only 7,033,706 pairs (78%) kept, mean length about 78/80, Q30 rose from 69% to 90%, but depth fell to about 9.4x (not 15x) because we subsampled BEFORE trimming.

### Decision made (APPLIED 2026-11-22: trim first, then subsample)
Trim first, then subsample: run fastp on the FULL reads, then subsample the trimmed pairs to about 15x of trimmed bases (about 1.8 Gb per sample; fraction = 1.8e9 / trimmed bases, same seed for R1 and R2). The trial folders `data/reads_15x/` and `data/reads_trimmed/` can stay; do not delete without asking.

### Phase 2 results (2026-11-23)
- Script `scripts/trim_subsample.sh SAMPLE [MIN_LEN]` (fastp with the settings above, then `seqtk sample -s 100` to 1.8e9 bases = about 15x; if trimmed depth is below 15x it copies all reads instead, because seqtk treats a fraction >= 1 as a read count). It was run for all 8 accessions (Ale-Stenar-44-4 with MIN_LEN 30). All finished, no failures.
- Outputs: trimmed full reads `data/reads_trimmed_full/<S>/<S>_R{1,2}.trim.fastq.gz` (kept, about 20 GB; do not delete without asking); FINAL reads for mapping `data/reads_final/<S>/<S>_R{1,2}.fastq.gz`; fastp reports `results/fastp/<S>_full.json/.html`, per-sample logs `results/fastp/<S>_run.log` (contain the depth and fraction lines). Depth = trimmed bases / 119,146,348 (chr1-5 only, so true nuclear depth is a little lower).
- Final pairs and depth (R1 and R2 counts identical for all):
  - Col-0: 11,336,064 pairs, 15.1x (trimmed 18.1x, fraction 0.8366)
  - An-1: 7,197,684 pairs, 11.2x (all reads kept; only 8.7M pairs existed)
  - Kyoto: 8,249,774 pairs, 13.3x (all kept)
  - Kondara: 9,850,439 pairs, 14.9x (all kept)
  - Cvi-0: 9,796,614 pairs, 15.0x (trimmed 29.2x, fraction 0.5167)
  - Ler-1: 9,313,005 pairs, 15.0x (trimmed 37.4x, fraction 0.4040)
  - Se-0: 9,456,036 pairs, 15.1x (trimmed 18.4x, fraction 0.8214)
  - Ale-Stenar-44-4: 12,333,578 pairs, 15.1x (trimmed 23.2x, fraction 0.6521; reads about 73 bp after trimming)
- Depth is therefore not perfectly equal (11.2x to 15.1x). Acceptable because PAV (Phase 4) normalises by each sample's own median depth; An-1 will have the least reliable low-coverage calls. Optionally equalise later.
- FastQC on the 16 final read files: DONE, 16 HTML reports in `results/qc_final/` (interpreted via MultiQC).
- MultiQC: the user was running `multiqc results/qc_final results/fastp --ignore "Col-0.json" --ignore "Col-0.html" -o results/multiqc -n multiqc_final` (ignores the old Col-0 trial fastp reports; raw FastQC in `results/qc_raw` deliberately left out because file names clash). Output: `results/multiqc/multiqc_final.html` (2.7 MB, exists, interpreted with the user; the reads are fine for mapping). MultiQC watch list: Cvi-0 (only 51% of reads kept after trimming; high duplication warning; cause unproven, possibly early quality collapse), An-1 (lowest depth, 11.2x), Ale-Stenar-44-4 (76 bp reads), Ler-1 R1 (bad cycle at base 3, 52% N, harmless), GC about 40% for Ale/Se-0.
- A stray EMPTY folder `result/` (typo) exists in the project root; offer to remove it with `rmdir result`, only after the user agrees.

### Phase 3 results (mapping, 2026-11-23/24): DONE (all 8 accessions mapped and checked 2026-11-24)
- Script `scripts/map_sample.sh SAMPLE`: bowtie2 `--very-sensitive -p 3` to `data/reference/index/TAIR10` (read group `--rg-id S --rg SM:S --rg PL:ILLUMINA`) | `samtools fixmate -m -u` | `sort -u -m 600M` | `markdup -s -f` (marks duplicates, does not remove them). Output: `data/bam/<S>/<S>.markdup.bam` + `.bai` (written as `<S>.tmp.bam`, then `mv`; the script REFUSES to overwrite an existing final BAM). QC files in `results/mapping/<S>/`: `<S>.bowtie2.log`, `.markdup.txt`, `.flagstat.txt`, `.stats.txt`, `.mosdepth.summary.txt`, `.mosdepth.global.dist.txt`. Run log: `results/mapping/<S>_run.log` (ends with `DONE <S>`).
- Run history: Col-0 (test, about 75 min) and An-1 (about 65 min) were started one by one with nohup. Then ONE background queue was started (run once): `nohup bash -c 'for S in Kyoto Kondara Cvi-0 Ler-1 Se-0 Ale-Stenar-44-4; do ...; done' > results/mapping/queue.log`; `queue.log` has `queue start/done/FAILED <S>` lines. State when this was written (2026-11-24): Kyoto done (95 min), Kondara done (71 min), Cvi-0 mapping (about 32% of R1 at last check, started 22:35 on 11-23), Ler-1, Se-0 and Ale-Stenar-44-4 still waiting. CHECK the real state first (`cat results/mapping/queue.log`; `ps aux | grep -E 'bowtie2-align|map_sample|samtools'`).
- WARNING: do NOT rerun the whole queue line. It loops over ALL listed accessions and would overwrite the `_run.log` of finished ones with an "already exists" message (the BAMs are protected, the logs are not). If the queue died (e.g. WSL was shut down), start a new loop with ONLY the unfinished accessions; an unfinished sample restarts from scratch (leftover `sorttmp/` and `*.tmp.bam` are overwritten; check them first and ask before deleting).
- Laptop sleep PAUSES the job (it resumes on wake; it only dies if WSL is shut down or the PC restarts). A restart killed the first Col-0 run once, and sleep happened during the queue. The user was given `powercfg` commands (standby/hibernate timeout 0 on AC); not confirmed applied. Lid close can also sleep the laptop.
- Progress check that works for any sample (the aligner is NOT the `sh -c` wrapper PID; do not use `pgrep -f`): find the process holding R1 open with `find /proc/[0-9]*/fd -lname '*<S>_R1.fastq.gz'`, read the `pos:` line of the matching `/proc/PID/fdinfo/N` twice 60 s apart, compare with the file size to get % and ETA (mapping step only; sort/markdup/QC add 10-15 min and show no percentage; look at the file listing instead).
- Results so far (reads = 2 x pairs in all four, nothing lost):
  - Col-0: 22,672,128 reads; mapped 96.3%; MAPQ>=10 78.4% of mapped; duplicates 3.3% (estimated library about 160 M fragments); properly paired 95.6%; insert 323+-56; error rate 0.53%; nuclear median depth 11x (mean 12.9x); chloroplast median 513x, mito 153x; organelles 8.4% of mapped reads.
  - An-1: 14,395,368 reads; mapped 94.0%; MAPQ>=10 72.1%; pairs not aligned 7.4% (probably divergence; contamination not excluded); duplicates 3.9%; properly paired 92.4%; nuclear median 7x (weakest sample); chloroplast 625x, mito 40x.
  - Kyoto: 16,499,548 reads; mapped 96.4%; MAPQ>=10 73.3%; duplicates 4.4%; properly paired 95.1%; nuclear median 8x; chloroplast 695x, mito 70x. Its slow run was NOT due to hard reads (mapping rate is like Col-0); most likely the laptop slept or was busy.
  - Kondara: 19,700,878 reads; mapped 96.6%; MAPQ>=10 69.4%; 33% of pairs concordant >1 time; duplicates 4.0%; properly paired 94.8%; nuclear median 9x; chloroplast 1350x (about 13% of its bases), which probably contributes to the low MAPQ share (chloroplast has an inverted repeat; inference, not verified).
- FINAL STATE (2026-11-24): the queue finished at 06:31 with no FAILED line; all 8 `<S>_run.log` end with `DONE <S>`, all 8 BAMs (+.bai), flagstat and mosdepth summary files exist. The sleep watcher fired at 06:31:39 (laptop went to sleep). Times: Kyoto 95 min, Kondara 71, Cvi-0 189 (laptop slept), Ler-1 117, Se-0 92, Ale-Stenar-44-4 78. The "WARNING do not rerun the queue" and the progress-check notes above are now historical.
- Results of the last four (reads = 2 x pairs, nothing lost):
  - Cvi-0: 19,593,228 reads; mapped 96.89%; properly paired 95.3%; 4.68% pairs aligned 0 times, 62.9% once, 32.4% >1 times; MAPQ>=10 69.3%; duplicates 13.1%; nuclear median 7x on all five; chloroplast 3371x, mito 58x. Usable; low depth and high duplicates need care (cause inferred: large chloroplast share plus duplicates, not tested).
  - Ler-1: 18,626,010 reads; mapped 95.25%; properly paired 91.7% (lowest, probably real divergence/rearrangements vs Col-0, not tested); 8.28% pairs aligned 0 times; MAPQ>=10 77.5%; duplicates 2.7%; nuclear median 12x on all five (best of the 7 non-Col-0); chloroplast 206x, mito 39x; error rate 1.57% (highest, most divergent).
  - Se-0: 18,912,072 reads; mapped 96.33%; properly paired 94.8%; 5.20% pairs 0 times, 33.7% >1 times; MAPQ>=10 68.4%; duplicates 4.9%; nuclear median 9x; chloroplast 784x, mito 157x.
  - Ale-Stenar-44-4: 24,667,156 reads; mapped 96.43%; properly paired 94.3%; 5.74% pairs 0 times, 33.2% >1 times; MAPQ>=10 68.8%; duplicates 18.2%; nuclear median 8x on all five; chloroplast 309x, mito 154x. Short reads (about 73 bp) plus 18% duplicates: treat its PAV calls with extra care.
- Nuclear median depth (chr1-5) per sample, the PAV normaliser: Col-0 11, An-1 7, Kyoto 8, Kondara 9, Cvi-0 7, Ler-1 12, Se-0 9, Ale-Stenar-44-4 8.
- MultiQC on the mapping QC: `results/multiqc/multiqc_mapping.html` (+ `multiqc_mapping_data/`). It found 8 samples for mosdepth, samtools stats, flagstat and bowtie2 (no name clashes). A second run made `multiqc_mapping_1.html` and `multiqc_mapping_data_1/` (redundant duplicates; remove only after asking). Side by side (`multiqc_general_stats.txt`): breadth >=1x is 100% for Col-0 but 93-95% for the other seven (so 5-7% of the reference has no reads in them: first PAV signal, but also unmappable regions); breadth >=5x: Col-0 96, Ler-1 91, Kondara 80, Kyoto 80, Se-0 78, Ale 75, Cvi-0 75, An-1 73; samtools error rate Col-0 0.53% vs 1.06-1.14% for the others and 1.57% for Ler-1 (divergence from Col-0); MQ0 1.6% (Col-0) to 3.7% (Ale); insert size average 168-367 (different libraries: expected, not a fault). Verdict: no failed sample; An-1, Cvi-0 and Ale-Stenar-44-4 are the weakest for PAV.
- Interpretation notes: (1) Use the MEDIAN depth of chr1-5, not the mean, as each sample's normaliser: in Col-0 the medians of chr1-5 are all 11 although the means are 11.3-14.7, so the means are inflated by collapsed repeats. The mosdepth `total` mean includes organelles; do not use it. (2) `*.mosdepth.global.dist.txt` is ordered from HIGH to low depth and fractions are rounded to 2 decimals; median = largest depth with fraction >= 0.5: `awk '$3>=0.5 { if (!($1 in m) || $2+0 > m[$1]+0) m[$1]=$2 } END{for(c in m) print c, m[c]}' FILE | sort`. (3) In Col-0 low MAPQ is concentrated on chr2 (69% have MAPQ>=10), chr3 (73%), chr4 (78%), chloroplast (66%), mito (42%); chr1 93%, chr5 87%. 27.8% of Col-0 pairs aligned concordantly >1 time. This is a property of repeats/paralogs/organelle-like insertions (a large mitochondrial insertion on chr2 is known from the literature, not verified here), not a mapper fault. Chr1 has 0.54% N (others <=0.04%), too little to explain anything. (4) Phase 4 must treat genes whose reads are mostly low-MAPQ as "uncertain", not "absent". (5) Reads average 79 bp after trimming (max 100). (6) Earlier I said organelles take about 4% of reads; the correct Col-0 figure is 8.4%.
- Reference sequence names: NC_003070.9, NC_003071.7, NC_003074.8, NC_003075.7, NC_003076.8 = chr1-5; NC_037304.1 = mitochondrion (367,808 bp); NC_000932.1 = chloroplast (154,478 bp). Reference FASTA: `data/reference/TAIR10.1_Col-0/GCF_000001735.4/GCF_000001735.4_TAIR10.1_genomic.fna`.
- Time estimates given for later phases (rough, not measured): Phase 4 minutes (mosdepth per gene, 15-30 min for 8 samples); Phase 5 hours (bcftools 1-3 h, Delly 2-4 h; SnpEff needs a database download, ask first); Phase 6 about 1-3 h per assembly (BUSCO and Liftoff; Liftoff RAM tight on 3.8 GB; downloads need asking); Phase 7 minutes.

### Phase 4 results (PAV calling, 2026-11-24): FIRST PASS DONE, validation pending
- Python packages installed in env `map` (conda-forge): pandas 3.1.1, numpy 2.5.5, matplotlib 3.11.2.
- `scripts/pav_depth.sh SAMPLE`: mosdepth `-n -t 1 --by data/reference/genes.bed -T 1,5 -Q {10,0}` -> `results/pav/<S>/<S>.q10.*` and `<S>.q0.*` (regions.bed.gz = mean depth per gene, thresholds.bed.gz = bases >=1x and >=5x; refuses to overwrite). Col-0 tested first (27,444 genes, depth about 10-12x, breadth 100% for AT1G01010-30), then the other 7 ran as one background queue (20:35-21:00 on 11-24, about 4 min each, no FAILED; `results/pav_queue.log`). Do not rerun the queue.
- `scripts/build_pav.py` (median depth per accession hard-coded in the script, thresholds are CLI options; defaults present = breadth >= 0.8 and depth >= 0.2 x median; absent = breadth < 0.1; "uncertain" = breadth < 0.1 at MAPQ >= 10 but all-reads breadth >= 0.5, i.e. covered by low-MAPQ reads, NOT called absent; else ambiguous). Outputs: `results/pav/pav_long.tsv.gz` (raw metrics + call), `results/pav/pav_matrix.tsv` (1 present, 0 absent, 0.5 ambiguous or uncertain), `results/pav/pav_summary.tsv`.
- Calls (27,444 genes): Col-0 present 27,041 / ambiguous 165 / uncertain 238 / absent 0 (control passes); An-1 25,948 / 848 / 282 / 366 (1.33%); Kyoto 25,930 / 918 / 316 / 280 (1.02%); Kondara 25,873 / 911 / 301 / 359 (1.31%); Cvi-0 25,607 / 1,084 / 301 / 452 (1.65%); Ler-1 26,118 / 691 / 317 / 318 (1.16%); Se-0 26,014 / 969 / 288 / 173 (0.63%); Ale-Stenar-44-4 26,124 / 752 / 269 / 299 (1.09%).
- Reading: Col-0 has 0 absences (good) but about 1.5% of its genes are still ambiguous/uncertain (repeat-rich, low-mappability genes), which is the noise floor. Absent genes are 173-452 per accession; Cvi-0 (most divergent, lowest depth) has the most, Se-0 the fewest. CAUTION (not tested): some "absent" genes may be very divergent genes whose reads fail to map to Col-0, not real deletions; validate with Delly/SNP density (Phase 5) and the assemblies (Phase 6). Thresholds have not been tuned yet (no sensitivity test).

- Matrix split (27,444 genes): present in all 8 = 24,259 (88.4%); no absence but some ambiguous/uncertain = 2,274; absent in exactly 1 accession = 389; absent in 2-6 = 504; absent in all 7 non-Col-0 = 18 (list: `results/pav/absent_in_all_others.txt`). The 18 form about 11 loci (neighbours AT1G15885/15890, AT1G80970/80980, AT2G11225/11270/11271, AT1G34810/34815/34840). 15 have normal Col-0 depth and ~0 reads in the others even with low MAPQ (real Col-0-specific candidates); AT1G34810/34815/34840 are short (~270 bp) with Col-0 depth ~0 at MAPQ>=10 (probable repeat-mappability artefacts). `gene_names.tsv` has no real names for them. Check them against the assemblies in Phase 6.
- Threshold sensitivity (absent counts An-1 / Kyoto / Kondara / Cvi-0 / Ler-1 / Se-0 / Ale; Col-0): absent-breadth < 0.05: 340/246/337/409/301/157/274, Col-0 0; < 0.10 (default): 366/280/359/452/318/173/299, Col-0 0; < 0.20: 417/328/403/500/365/204/338, Col-0 1. Order of accessions is the same under all three. Changing the PRESENT thresholds (0.9 and 0.3 x median) does not change absent counts (it only moves genes between present and ambiguous; those counts were not printed). Output folders `results/pav_sens/{abs05,abs20,pres90}` (can be removed later, ask first). `build_pav.py` now has `--indir` and `--outdir`.
- Figures: `scripts/plot_pav.py` writes `results/figures/pav_calls_per_accession.png`, `depth_hist_per_accession.png`, `pav_heatmap_variable.png` (no clustering yet; clustering needs scipy, Phase 7).

### Next steps (updated 2026-11-24)
1. DONE (2026-11-24): all 8 mapped and checked.
2. DONE (2026-11-24): MultiQC on results/mapping, 8 samples, interpreted (see Phase 3 results).
3. Optional: bowtie2 vs BWA-MEM on Cvi-0 (needs the BWA index, not built; ask first). Optional: look at a sample of unmapped reads of An-1 (7.4% unaligned) for contamination.
4. UPDATE 2026-11-24: Phase 4 first pass is DONE (see "Phase 4 results"). Remaining in Phase 4: quick matrix summary (core/dispensable/private counts), threshold sensitivity test, first plots (`scripts/plot_pav.py` -> `results/figures/`, offered, not yet agreed). The text below is the original plan for Phase 4. Phase 3 is CLOSED. Original item: write `scripts/pav_depth.sh` (mosdepth `--by data/reference/genes.bed --thresholds 1,5`, only reads with MAPQ >= 10, plus an all-MAPQ run to flag low-mappability genes) and a small script that builds the PAV matrix, normalised by each sample's chr1-5 median depth (Col-0 11, An-1 7, Kyoto 8, Kondara 9, Cvi-0 7, Ler-1 12, Se-0 9, Ale-Stenar-44-4 8). Plan thresholds: present = breadth >= 0.8 and depth >= 0.2 x median; absent = breadth < 0.1; else ambiguous. Col-0 is the control (about 0 absences expected). Test the per-gene mosdepth on Col-0 first (cheap, about a couple of minutes, one thread). Python packages (pandas etc.) are not installed: give the `conda install -c conda-forge --override-channels` command when needed and check versions.
5. Later Phases 5-8 as in the original plan. Ask before every download (SnpEff database, BUSCO lineage `brassicales_odb10`, assemblies) and before raising the `.wslconfig` memory.

### Older next-steps notes (historical, written before Phase 3)
1. Finish MultiQC (see above): confirm `results/multiqc/multiqc_final.html`, interpret it in plain language, note any accession with problems. This closes Phase 2.
2. Phase 3 mapping (no new tool installs needed: bowtie2, samtools, mosdepth are installed). Write a script like `trim_subsample.sh`, e.g. `scripts/map_sample.sh SAMPLE`, run ONE sample at a time with `nohup`, test on Col-0 first. Details: bowtie2 `--very-sensitive -p 3` to `data/reference/index/TAIR10` using `data/reads_final/<S>/`, samtools sort (`-m 1G`, few threads, RAM is 3.8 GB)/markdup/index, flagstat, mosdepth, outputs under `results/mapping/` or `data/bam/`. Expect mapping above 90% (lower for Cvi-0). Compare bowtie2 vs bwa on Cvi-0 (needs BWA index, not built yet). Add flagstat/mosdepth to MultiQC.
3. Phase 4 PAV (mosdepth per gene on genes.bed, MAPQ >= 10, thresholds in config), Col-0 control. Then Phases 5-8 as in the original plan.

### Open questions
1. Ler: reads are Ler-1 (1001G 6932); 1001G also lists Ler-0 (7213). Assembly GCA_902460285.1 is only labelled "Ler"; headers do not name the line. Kyoto vs Kyo likewise unconfirmed. An-1/Cvi-0 probable. Resolve before Phase 6 (needs NCBI metadata; ask before downloading).
2. Ale-Stenar-44-4 has 76 bp reads; treat its PAV calls with extra care.
3. RAM is 3.8 GB: BUSCO, Liftoff and snpEff may need `.wslconfig` memory raised (ask the user first).

### Current folder layout
genom_project/ PROJECT_PLAN.md, NGS-Course-MI26.pdf, config/samples.tsv, scripts/, notebooks/, tools/ (sratoolkit, unused),
data/{reads (full), reads_15x (trial), reads_trimmed (trial), reads_trimmed_full (all 8, kept), reads_final (all 8, USE THESE for mapping), reference/{TAIR10.1_Col-0, index, genes.bed, gene_names.tsv}, assemblies, meta},
results/{qc_raw, qc_final, fastp, multiqc, mapping (per-accession QC folders, `<S>_run.log`, `queue.log`)}, data/bam/<S>/<S>.markdup.bam(+.bai); scripts/trim_subsample.sh, scripts/map_sample.sh; stray empty `result/` in root (typo, ask before removing)

### Ideas to improve results (decided earlier)
- Count only reads with MAPQ >= 10 for PAV; flag low-mappability genes (NLR clusters, TE genes, paralogs).
- Cvi-0, Kondara and Se-0 are divergent: compare Bowtie2 with BWA-MEM on one accession.
- Eight accessions is small for a core/dispensable split; consider 20-30 later.
- Use all 4 assemblies (Liftoff) as ground truth to tune thresholds; add a simulated-deletion control in Col-0.
- Equalise depth (about 15x after trimming) for the main analysis; optionally map full data and compare.

---

## Context (original plan; paths below may be outdated, see CURRENT STATUS)
- `d:\genom_project` originally contained `ncbi_dataset_2.zip` (TAIR10.1 Col-0: GCF_000001735.4 RefSeq **and** GCA_000001735.2 GenBank copies of FASTA/GFF/proteins) and `NGS-Course-MI26.pdf` (course outline: QC/trimming, BWT read mapping with Bowtie2, SAM/BAM, variant calling, VCF, variant annotation, CNV/SV, assembly + QC, DGE + functional enrichment).
- **Key point:** one reference genome alone cannot show gene presence/absence. You need other accessions to compare against Col-0. The plan adds 6–8 natural accessions and asks one question: which Col-0 genes are present, absent or partly lost in each accession?
- Hardware: 8 GB RAM, 2 cores, ~15–20 GB free. The genome is only ~120 Mb, so read mapping fits easily. Full de novo pan-genome graphs (PGGB, minigraph-cactus) and assembling every accession would not fit, so they are left out.
- Environment found: WSL2 Ubuntu is installed (currently stopped), plus Anaconda, git and Python on Windows. All bioinformatics tools will run inside WSL2.
- Goal: a reproducible pipeline, published on GitHub with data versioned on DagsHub, plus a LinkedIn write-up.

## Approach (reference-based PAV, mostly using the course tools)
**Main track:** short-read WGS of each accession → QC → Bowtie2 mapping to TAIR10 → gene-level breadth and depth of coverage → PAV matrix → core vs dispensable genes → functional enrichment.
**Supporting evidence:** bcftools variants + SnpEff (loss-of-function = genes that are present but broken), Delly deletions (structural evidence for absences).
**Validation:** for 1–2 accessions that also have a chromosome-level assembly on NCBI (e.g. Ler-0), lift Col-0 genes over with Liftoff and check that the read-based "absent" calls agree.

### Accessions (6–8, chosen for diversity, ideally ones that also have published assemblies)
Candidates: Ler-0, Cvi-0, C24, Kyo, An-1, Sha, Eri-1 (the Jiao & Schneeberger 2020 set), plus Col-0 reads as a negative control (it should show ~0 absences).
Source: 1001 Genomes Project reads on SRA/ENA (BioProject PRJNA273563). Look up the actual SRR run IDs in the ENA browser and write them into `config/samples.tsv`. Don't hard-code guessed IDs.
Disk: subsample each accession to ~15× (~1.8 Gb of bases, ~1 GB fastq.gz) with `seqtk`/`rasusa`, and delete the trimmed FASTQ once the BAM is made. Peak disk use stays around 8–10 GB.

## Phases

### Phase 0 — Environment (WSL2)
1. Start Ubuntu in WSL2. Put the project at `/mnt/d/genom_project`, or under `~/` in WSL for faster I/O. Add `.wslconfig` with `memory=6GB`, `processors=2`, `swap=4GB`.
2. Install Miniforge (mamba) inside WSL. Create `envs/*.yaml` (one env per stage, pinned versions):
   - `qc.yaml`: fastqc, fastp (or trimmomatic), multiqc, seqtk, sra-tools / ena-file-downloader
   - `map.yaml`: bowtie2, samtools, mosdepth, bedtools
   - `variants.yaml`: bcftools, snpeff (db `Arabidopsis_thaliana`), delly
   - `validate.yaml`: liftoff, minimap2, quast, busco (`brassicales_odb10`, run with `-c 2`)
   - `analysis.yaml`: python (pandas, numpy, matplotlib, seaborn, scikit-learn, upsetplot), R + clusterProfiler/org.At.tair.db (or the gprofiler-official python client)
   - `snakemake` + `dvc` in a base `workflow.yaml`
3. `git init`, `dvc init`.

### Phase 1 — Reference preparation (`workflow/rules/reference.smk`)
- Extract **only the GCF_000001735.4** files from the zip into `data/reference/` (RefSeq has the curated annotation; ignore the GCA copy). Check them against `md5sum.txt`.
- `samtools faidx` and `bowtie2-build --threads 2` (~10–20 min, <2 GB RAM).
- From `genomic.gff`, write a BED of protein-coding genes (gene_biotype=protein_coding) and a gene → locus_tag (ATxGxxxxx) map. Leave out organelle chromosomes (NC_000932 chloroplast, NC_037304 mitochondrion) from the PAV calls.
- Reference stats with `seqkit stats` and QUAST on the reference (course: assembly QC). BUSCO on the proteins gives a completeness baseline.

### Phase 2 — Read download and QC (course lectures 9–15)
- Download each accession's runs (ENA FTP is faster than `fasterq-dump`), then subsample to 15×.
- Run `fastqc` on the raw reads, `fastp` (adapter + quality trimming; report HTML/JSON), `fastqc` on the trimmed reads, and combine everything with `multiqc` into one QC report.

### Phase 3 — Mapping (lectures 16–25)
- `bowtie2 --very-sensitive -p 2 -x ref` | `samtools sort -@1 -m 1G` → `samtools markdup` → index.
- `samtools flagstat`/`stats` and `mosdepth` for genome-wide depth; feed these into MultiQC.

### Phase 4 — Gene PAV calling (core of the project)
- `mosdepth --by genes.bed --thresholds 1,5` per sample → per-gene mean depth and breadth (fraction of bases ≥1×).
- Classification, normalised to each sample's median depth (put the thresholds in `config.yaml` so they are easy to change):
  - **Present:** breadth ≥ 0.8 and depth ≥ 0.2 × sample median
  - **Absent:** breadth < 0.1
  - **Partial/ambiguous:** everything in between (report it, don't force a call)
- Output `results/pav/pav_matrix.tsv` (genes × accessions, values 1/0/0.5) and a long-format table with the raw metrics.
- Col-0 control: absences should be close to zero. If not, tighten the thresholds.

### Phase 5 — Variants, LoF and SVs (lectures 26–31)
- `bcftools mpileup | bcftools call -mv` (multi-sample, run per chromosome to keep RAM low) → `bcftools filter` on QUAL/DP → `bcftools stats`.
- `snpEff Arabidopsis_thaliana` annotation; pull out HIGH-impact variants (stop-gained, frameshift, splice) → "present but disrupted" genes.
- `delly call -t DEL` → deletions that overlap genes called absent (structural confirmation).

### Phase 6 — Assembly-based validation (lectures 51–55)
- Download 1–2 NCBI chromosome-level assemblies (e.g. Ler-0) as `.fna` only.
- Run QUAST and BUSCO on them for assembly QC.
- Run `liftoff -p 2 -g genomic.gff`, keeping the protein-coding subset to save RAM. Genes in the "unmapped" list = absent in the assembly.
- Build a confusion matrix of read-based vs assembly-based calls (precision/recall) and put it in the README.

### Phase 7 — Pan-genome analysis and figures (lectures 40–43)
`notebooks/` or `workflow/scripts/analysis.py`:
- Core (present in all), shell/dispensable (in some), private-absence (missing from one accession only) gene counts, shown as an UpSet plot.
- Gene accumulation/rarefaction curve: core genes decreasing as accessions are added (permutations).
- PAV heatmap with hierarchical clustering, plus PCA/Jaccard distance tree of accessions.
- Chromosome ideogram of where absent genes sit (clusters near pericentromeres, NLR clusters).
- GO/KEGG enrichment of dispensable genes vs all genes (clusterProfiler `enrichGO` with org.At.tair.db, BH FDR — course lectures 41–43). Expected results: disease resistance (NLR), F-box and TE-related genes.
- Optional: `snakemake --report` HTML.

### Phase 8 — Reproducibility and publishing
**Repository layout**
```
genom_project/
├── README.md              # summary, figures, how to run, results, citations
├── LICENSE (MIT)  CITATION.cff  .gitignore  .dvcignore
├── config/config.yaml  config/samples.tsv
├── envs/*.yaml
├── workflow/Snakefile  workflow/rules/*.smk  workflow/scripts/*.py|R
├── notebooks/01_pav_exploration.ipynb
├── data/        (DVC tracked: reference, small test data)
├── results/     (small tables + figures in git; BAM/VCF via DVC or ignored)
├── docs/figures/  docs/methods.md
└── .github/workflows/ci.yml   # snakemake --lint + dry-run on a tiny test dataset
```
- **Test dataset:** chromosome 4 subset with 0.5× reads, so CI and reviewers can run the whole workflow in minutes.
- **GitHub:** clean commit history (one commit per phase), README badges, topics (`bioinformatics`, `pangenome`, `arabidopsis`, `snakemake`).
- **DagsHub:** connect the GitHub repo (mirror). Run `dvc remote add origin https://dagshub.com/<user>/<repo>.dvc` and `dvc push` the reference indices, PAV matrix and MultiQC report (skip raw FASTQs, which you can re-download from ENA). Optionally log per-sample QC/mapping metrics with DagsHub's MLflow tracking.
- **LinkedIn:** a 150–250 word post with 2–3 figures (UpSet, heatmap, accumulation curve), the key numbers (e.g. "X% of Col-0 genes are dispensable across N accessions; enriched for NLR defence genes"), the tool stack and links to GitHub and DagsHub. Optionally a carousel PDF showing the workflow DAG (`snakemake --dag`).

## RAM/CPU safety rules
- Use `-p/--threads 2` everywhere and `samtools sort -m 1G`. Process samples one at a time (`snakemake -j 1` for mapping and Liftoff).
- Use bowtie2 rather than BWA-MEM2 (BWA-MEM2 indexing needs more RAM). Avoid PGGB/minigraph-cactus/de novo assembly of reads.
- Delete intermediate FASTQ files once the BAM exists (`temp()` in Snakemake).

## Verification
1. `snakemake -n` (dry run) and `snakemake --lint` pass; CI runs the test dataset successfully.
2. Col-0 control: fewer than ~0.5% of genes called absent.
3. Mapping rate above 90% for most accessions (lower for Cvi-0 is expected), checked in MultiQC.
4. Read-based vs Liftoff agreement for Ler-0 (target ≥ 85% concordance on absent calls).
5. Sanity check: known Ler-0/Col-0 differences, e.g. **FLC/FRI** alleles (Ler carries an FLC transposon insertion and a FRI deletion) and RPS5/RPP loci, should show up as variants or absences.
6. A fresh clone plus `dvc pull` plus `snakemake --use-conda` reproduces `pav_matrix.tsv`.

## Execution order once approved
Phase 0 → 1 → (select SRR IDs together) → 2 → 3 → 4 → 5 → 6 → 7 → 8. Phases 2–3 take the longest (a few hours per accession on 2 cores), so they run in the background.

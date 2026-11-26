"""Known-gene check: do the PAV calls agree with well-studied genes from the literature? Env `map`.

Looks up read-depth calls (results/pav/pav_long.tsv.gz), assembly checks (results/assembly/liftoff_vs_pav_genes.tsv.gz)
and HIGH-impact SnpEff effects (results/variants/lof_genes.tsv.gz) for a fixed list of genes whose natural variation is
described in published work, and writes one row per gene and accession.

Expectations from the literature (checked on PubMed):
  RPS5  whole-locus presence/absence polymorphism kept by balancing selection (Tian et al. 2002, PNAS)
  RPM1  susceptible accessions lack the whole gene (rpm1-null deletion) (Grant et al. 1998, PNAS)
  FRI   Col-0 carries a non-functional allele; Cvi has an in-frame stop codon in exon 1
        (Johanson et al. 2000, Science; Gazzani et al. 2003, Plant Physiol)
  FLC   variation is by expression (e.g. a transposon in the Ler intron), not by deletion (Gazzani et al. 2003)
  ZAR1, ADR1, ADR1-L1, ADR1-L2, PBS1  conserved immune genes: expected present in all accessions
  RPS4 / RRS1  a head-to-head NLR pair that works together: expected to be lost or kept together

Output: results/analysis/known_genes_check.tsv
"""
import pandas as pd

GENES = {
    "AT1G12220": ("RPS5", "NLR, presence/absence polymorphism"),
    "AT3G07040": ("RPM1", "NLR, deleted in susceptible accessions"),
    "AT5G45250": ("RPS4", "NLR pair with RRS1"),
    "AT5G45260": ("RRS1", "NLR pair with RPS4"),
    "AT4G26090": ("RPS2", "NLR"),
    "AT3G46530": ("RPP13", "NLR, highly variable"),
    "AT3G50950": ("ZAR1", "conserved NLR"),
    "AT1G33560": ("ADR1", "conserved helper NLR"),
    "AT4G33300": ("ADR1-L1", "conserved helper NLR"),
    "AT5G04720": ("ADR1-L2", "conserved helper NLR"),
    "AT5G13160": ("PBS1", "conserved kinase guarded by RPS5"),
    "AT4G00650": ("FRI", "flowering; Col-0 allele non-functional"),
    "AT5G10140": ("FLC", "flowering; varies by expression"),
    "AT3G44540": ("FAR4", "example deletion in Cvi-0 (fig 23)"),
}

pav = pd.read_csv("results/pav/pav_long.tsv.gz", sep="\t", usecols=["gene", "accession", "norm_depth", "call"])
asm = pd.read_csv("results/assembly/liftoff_vs_pav_genes.tsv.gz", sep="\t", usecols=["gene", "accession", "asm_class", "coverage"])
lof = pd.read_csv("results/variants/lof_genes.tsv.gz", sep="\t", usecols=["gene", "accession", "effects"])

d = pav[pav.gene.isin(GENES)].merge(asm, on=["gene", "accession"], how="left").merge(lof, on=["gene", "accession"], how="left")
d.insert(1, "name", d.gene.map(lambda g: GENES[g][0]))
d.insert(2, "note", d.gene.map(lambda g: GENES[g][1]))
d["asm_class"] = d.asm_class.fillna("no assembly")
d["effects"] = d.effects.fillna("")
d = d.sort_values(["gene", "accession"])
d.to_csv("results/analysis/known_genes_check.tsv", sep="\t", index=False, float_format="%.3f")

wide = d.pivot(index="name", columns="accession", values="call").reindex([v[0] for v in GENES.values()])
print(wide.to_string())
print("\nwrote results/analysis/known_genes_check.tsv")

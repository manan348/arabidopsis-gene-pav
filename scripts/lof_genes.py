#!/usr/bin/env python3
"""Genes hit by HIGH-impact variants (stop gained, frameshift, splice, start lost ...), per accession (Phase 5).

Run in env `map` after exporting the HIGH sites from the SnpEff-annotated VCF:
    python scripts/lof_genes.py

Inputs:  results/variants/high_impact.tsv.gz  (CHROM POS REF ALT ANN GT x 8; one line per HIGH site)
         results/variants/samples.txt         (sample names, same order as the GT columns)
         results/pav/pav_long.tsv.gz          (gene, accession, call)
Outputs: results/variants/lof_genes.tsv.gz    (gene, accession, effects, n_sites, pav_call)
         results/variants/lof_summary.tsv     (per accession: HIGH sites carried, genes hit, genes hit by PAV call)

An accession "carries" a HIGH effect when its genotype is the alternative allele that SnpEff scored.
Inbred lines were called as haploid, so the genotype is a single allele.
"Present but broken" = the gene has a HIGH variant AND the read-depth call is `present`.
Descriptive only: no statistical test. HIGH effects are predictions, not proof that the protein is lost
(for example, a stop near the end of a gene may leave the protein working).
"""
import argparse
import gzip
import os
from collections import defaultdict

import pandas as pd

ORDER = ["Col-0", "An-1", "Kyoto", "Kondara", "Cvi-0", "Ler-1", "Se-0", "Ale-Stenar-44-4"]


def parse_ann(ann):
    """SnpEff ANN string -> {alt allele: set of (gene_id, effect)} for HIGH-impact, protein-coding entries."""
    out = defaultdict(set)
    for entry in ann.split(","):
        f = entry.split("|")
        if len(f) > 7 and f[2] == "HIGH" and f[7] == "protein_coding":
            out[f[0]].add((f[4], f[1]))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--high", default="results/variants/high_impact.tsv.gz")
    ap.add_argument("--samples", default="results/variants/samples.txt")
    ap.add_argument("--pav", default="results/pav/pav_long.tsv.gz")
    ap.add_argument("--outdir", default="results/variants")
    a = ap.parse_args()

    names = [x.strip() for x in open(a.samples) if x.strip()]
    col = {}
    for acc in ORDER:
        hit = [i for i, n in enumerate(names) if acc in n]
        # "Kyoto" must not match "Kyo..." twice and Col-0 must match exactly one name
        assert len(hit) == 1, f"cannot match {acc} to exactly one sample name: {names}"
        col[acc] = hit[0]
    assert len(set(col.values())) == len(ORDER), "two accessions matched the same sample"

    rows, n_sites, n_unmatched = [], 0, 0
    with gzip.open(a.high, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            chrom, pos, ref, alts, ann = f[0], f[1], f[2], f[3].split(","), f[4]
            gts = f[5:5 + len(names)]
            assert len(gts) == len(names), "genotype columns do not match the sample list"
            hi = parse_ann(ann)
            if not hi:
                continue  # HIGH only in a non-coding transcript
            n_sites += 1
            matched = False
            for acc in ORDER:
                g = gts[col[acc]].replace("|", "/").split("/")
                for allele in set(g):
                    if allele in (".", "0", ""):
                        continue
                    alt = alts[int(allele) - 1]
                    for gene, effect in hi.get(alt, ()):
                        rows.append((gene, acc, effect, f"{chrom}:{pos}"))
                        matched = True
            if not matched and not any(a_ in hi for a_ in alts):
                n_unmatched += 1  # SnpEff allele string does not equal any VCF ALT (should be rare)
    print(f"HIGH sites in protein-coding transcripts: {n_sites}; allele string not matched to an ALT: {n_unmatched}")

    d = pd.DataFrame(rows, columns=["gene", "accession", "effect", "site"])
    d["effect"] = d["effect"].str.split("&")  # combined effects such as stop_gained&splice_region_variant
    d["effect"] = d["effect"].str[0]
    g = (d.groupby(["gene", "accession"])
           .agg(effects=("effect", lambda x: ";".join(sorted(set(x)))), n_sites=("site", "nunique"))
           .reset_index())

    calls = pd.read_csv(a.pav, sep="\t", usecols=["gene", "accession", "call"]).rename(columns={"call": "pav_call"})
    g = g.merge(calls, on=["gene", "accession"], how="left")
    # SnpEff's gene models and the RefSeq-based gene list (genes.bed) are not identical: a few genes
    # with HIGH variants were never tested for presence/absence. Keep them, but label them.
    n_miss = g.loc[g["pav_call"].isna(), "gene"].nunique()
    print(f"genes with a HIGH variant that are not in the PAV gene list: {n_miss} (labelled not_in_pav)")
    g["pav_call"] = g["pav_call"].fillna("not_in_pav")

    os.makedirs(a.outdir, exist_ok=True)
    g.to_csv(os.path.join(a.outdir, "lof_genes.tsv.gz"), sep="\t", index=False)

    s = []
    for acc in ORDER:
        da, ga = d[d["accession"] == acc], g[g["accession"] == acc]
        r = {"accession": acc, "high_sites": da["site"].nunique(), "genes_hit": len(ga)}
        for c in ["present", "ambiguous", "uncertain", "absent", "not_in_pav"]:
            r[f"genes_{c}"] = int((ga["pav_call"] == c).sum())
        s.append(r)
    s = pd.DataFrame(s)
    s.to_csv(os.path.join(a.outdir, "lof_summary.tsv"), sep="\t", index=False)
    print(s.to_string(index=False))
    print("\ntop effect types (gene x accession):")
    print(d.drop_duplicates(["gene", "accession", "effect"])["effect"].value_counts().head(8).to_string())
    print("\nwrote", os.path.join(a.outdir, "lof_genes.tsv.gz"), "and lof_summary.tsv")


if __name__ == "__main__":
    main()

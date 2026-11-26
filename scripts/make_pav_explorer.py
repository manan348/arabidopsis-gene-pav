#!/usr/bin/env python3
"""Interactive PAV explorer: the 911 variable genes x 8 accessions, three methods side by side. Env `map`.

Run:
    python scripts/make_pav_explorer.py

Inputs (read only): results/pav/pav_long.tsv.gz, data/reference/genes.bed, data/reference/gene_names.tsv,
    results/sv/sv_vs_pav_genes.tsv.gz, results/assembly/liftoff_vs_pav_genes.tsv.gz,
    results/variants/lof_genes.tsv.gz, results/go/go_enrichment.tsv
Output: results/interactive/pav_explorer.html (one file; plotly.js is loaded from a CDN, so it needs internet)

Each cell = one gene in one accession, coloured by the read call. Hover shows read depth, Delly and Liftoff;
click a gene for a table of all 8 accessions and a link to TAIR. Filters: chromosome, GO group, minimum number of
accessions missing the gene, and gene order. GO groups come from the significant terms in go_enrichment.tsv.
"""
import argparse
import json
import os

import pandas as pd

ORDER = ["Col-0", "An-1", "Kyoto", "Kondara", "Cvi-0", "Ler-1", "Se-0", "Ale-Stenar-44-4"]
CHR = {"NC_003070.9": 1, "NC_003071.7": 2, "NC_003074.8": 3, "NC_003075.7": 4, "NC_003076.8": 5}
CALLS = ["present", "ambiguous", "uncertain", "absent"]
LIFT = ["no assembly", "found", "partial", "elsewhere", "missing"]
# first matching term wins (order = priority)
GROUPS = [("ADP binding", "NLR resistance gene"),
          ("NAD+ nucleosidase activity, cyclic ADP-ribose generating", "NLR resistance gene"),
          ("defense response to other organism", "other defense gene"),
          ("sesquiterpene biosynthetic process", "terpene synthase"),
          ("RNA-DNA hybrid ribonuclease activity", "RNase H (TE-like)"),
          ("protein localization to cell surface", "small unnamed cluster (fertilization GO)")]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="results/interactive/pav_explorer.html")
    a = ap.parse_args()

    p = pd.read_csv("results/pav/pav_long.tsv.gz", sep="\t",
                    usecols=["gene", "accession", "norm_depth", "breadth_q10", "call"])
    var = sorted(set(p.loc[p["call"] == "absent", "gene"]))
    assert len(var) == 911, len(var)
    p = p[p["gene"].isin(var)]

    bed = pd.read_csv("data/reference/genes.bed", sep="\t", header=None, usecols=[0, 1, 2, 3],
                      names=["chrom", "s", "e", "gene"]).drop_duplicates("gene").set_index("gene")
    names = pd.read_csv("data/reference/gene_names.tsv", sep="\t", header=None, names=["gene", "sym"])
    names = names.drop_duplicates("gene").set_index("gene")["sym"]
    sv = pd.read_csv("results/sv/sv_vs_pav_genes.tsv.gz", sep="\t").set_index(["gene", "accession"])["deleted_frac"]
    lo = pd.read_csv("results/assembly/liftoff_vs_pav_genes.tsv.gz", sep="\t").set_index(["gene", "accession"])
    lof = pd.read_csv("results/variants/lof_genes.tsv.gz", sep="\t").set_index(["gene", "accession"])["effects"]
    go = pd.read_csv("results/go/go_enrichment.tsv", sep="\t")
    go = go[go["set"] == "absent"]
    group = {}
    for term, label in GROUPS:
        row = go[go["Description"] == term]
        if len(row):
            for g in row.iloc[0]["geneID"].split("/"):
                group.setdefault(g, label)
    labels = list(dict.fromkeys(l for _, l in GROUPS)) + ["other"]

    pc = p.set_index(["gene", "accession"])
    genes, cells = [], []
    for g in var:
        sym = names.get(g, g)
        b = bed.loc[g]
        genes.append([g, "" if sym == g else sym, CHR[b["chrom"]], int(b["s"]) + 1, int(b["e"]),
                      labels.index(group.get(g, "other"))])
        row = []
        for acc in ORDER:
            r = pc.loc[(g, acc)]
            k = (g, acc)
            if k in lo.index:
                lr = lo.loc[k]
                lift, cov = LIFT.index(lr["asm_class"]), round(float(lr["coverage"]), 2)
            else:
                lift, cov = 0, None
            row.append([CALLS.index(r["call"]), round(float(r["norm_depth"]), 2), round(float(r["breadth_q10"]), 2),
                        round(float(sv.get(k, 0.0)), 2), lift, cov, str(lof.get(k, ""))])
        cells.append(row)
    nabs = [sum(c[0] == 3 for c in row) for row in cells]
    print(f"{len(genes)} genes; groups: " +
          ", ".join(f"{l} {sum(1 for x in genes if labels[x[5]] == l)}" for l in labels))

    data = {"acc": ORDER, "calls": CALLS, "lift": LIFT, "groups": labels, "genes": genes, "cells": cells,
            "nabs": nabs}
    with open(os.path.join(os.path.dirname(__file__), "pav_explorer_template.html"), encoding="utf-8") as f:
        html = f.read().replace("__DATA__", json.dumps(data, separators=(",", ":")))
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    with open(a.out, "w", encoding="utf-8") as f:
        f.write(html)
    print(f"wrote {a.out} ({os.path.getsize(a.out) / 1e6:.1f} MB)")


if __name__ == "__main__":
    main()

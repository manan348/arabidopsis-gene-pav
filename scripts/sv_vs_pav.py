#!/usr/bin/env python3
"""Compare Delly PASS deletions with the PAV calls (Phase 5).

For every accession and gene, compute the fraction of the gene covered by that accession's
PASS deletions (deletions longer than --max-len are dropped: very long calls are usually artefacts
near centromeres/repeats). A gene counts as "inside a deletion" if >= --min-frac of it is covered.

If absent genes are real deletions, many should be inside a Delly deletion, and present genes
should almost never be. Delly needs read pairs spanning the deletion, so a gene with no Delly
support is NOT proof of a mapping failure (low depth, repeats and very long deletions are missed).

Usage (env map):  python scripts/sv_vs_pav.py
Inputs:  results/sv/pass_deletions.tsv (accession, chrom, pos, end), results/pav/pav_long.tsv.gz,
         data/reference/genes.bed
Outputs: results/sv/sv_vs_pav_summary.tsv (counts per accession and call),
         results/sv/sv_vs_pav_genes.tsv.gz (covered fraction per gene and accession)
"""
import argparse
import os

import numpy as np
import pandas as pd

ORDER = ["Col-0", "An-1", "Kyoto", "Kondara", "Cvi-0", "Ler-1", "Se-0", "Ale-Stenar-44-4"]


def merge(starts, ends):
    """Merge overlapping intervals; returns sorted, non-overlapping starts and ends."""
    o = np.argsort(starts)
    s, e = starts[o], ends[o]
    ms, me = [], []
    for a, b in zip(s, e):
        if ms and a <= me[-1]:
            me[-1] = max(me[-1], b)
        else:
            ms.append(a)
            me.append(b)
    return np.array(ms, dtype=np.int64), np.array(me, dtype=np.int64)


def covered_fraction(gs, ge, ms, me):
    """Fraction of each gene [gs, ge) covered by merged intervals [ms, me)."""
    out = np.zeros(len(gs))
    if len(ms) == 0:
        return out
    for k, (a, b) in enumerate(zip(gs, ge)):
        i0 = np.searchsorted(me, a, side="right")  # first interval ending after the gene start
        i1 = np.searchsorted(ms, b, side="left")   # intervals starting before the gene end
        if i1 > i0:
            ov = np.minimum(me[i0:i1], b) - np.maximum(ms[i0:i1], a)
            out[k] = ov[ov > 0].sum() / (b - a)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dels", default="results/sv/pass_deletions.tsv")
    ap.add_argument("--pav", default="results/pav/pav_long.tsv.gz")
    ap.add_argument("--bed", default="data/reference/genes.bed")
    ap.add_argument("--outdir", default="results/sv")
    ap.add_argument("--max-len", type=int, default=100_000, help="drop deletions longer than this (bp)")
    ap.add_argument("--min-frac", type=float, default=0.5, help="gene fraction that must be deleted")
    a = ap.parse_args()

    d = pd.read_csv(a.dels, sep="\t", header=None, names=["accession", "chrom", "pos", "end"])
    assert set(d["accession"]) <= set(ORDER), "unknown accession in the deletion table"
    d["len"] = d["end"] - d["pos"]
    assert (d["len"] > 0).all()
    long_ = d["len"] > a.max_len
    print(f"PASS deletions: {len(d)}; dropped as longer than {a.max_len:,} bp: {long_.sum()}")
    d = d[~long_]

    bed = pd.read_csv(a.bed, sep="\t", header=None, usecols=[0, 1, 2, 3], names=["chrom", "start", "end", "gene"])
    calls = pd.read_csv(a.pav, sep="\t", usecols=["gene", "accession", "call"])
    assert set(calls["gene"]) == set(bed["gene"]), "gene lists of the BED and the PAV table differ"

    rows = []
    for acc in ORDER:
        frac = np.zeros(len(bed))
        da = d[d["accession"] == acc]
        for c in bed["chrom"].unique():
            gi = (bed["chrom"] == c).to_numpy()
            dc = da[da["chrom"] == c]
            # Delly POS is 1-based; pos-1 makes it 0-based like the BED. Breakpoints are approximate anyway.
            ms, me = merge(dc["pos"].to_numpy() - 1, dc["end"].to_numpy())
            frac[gi] = covered_fraction(bed.loc[gi, "start"].to_numpy(), bed.loc[gi, "end"].to_numpy(), ms, me)
        rows.append(pd.DataFrame({"gene": bed["gene"], "accession": acc, "deleted_frac": frac}))
    g = pd.concat(rows).merge(calls, on=["gene", "accession"], how="left")
    assert g["call"].notna().all()
    g["in_deletion"] = g["deleted_frac"] >= a.min_frac

    s = (g.groupby(["accession", "call"])
           .agg(genes=("gene", "size"), in_deletion=("in_deletion", "sum"),
                any_overlap=("deleted_frac", lambda x: int((x > 0).sum())))
           .reset_index())
    s["pct_in_deletion"] = (100 * s["in_deletion"] / s["genes"]).round(2)
    s["accession"] = pd.Categorical(s["accession"], ORDER, ordered=True)
    s = s.sort_values(["accession", "call"])

    os.makedirs(a.outdir, exist_ok=True)
    s.to_csv(os.path.join(a.outdir, "sv_vs_pav_summary.tsv"), sep="\t", index=False)
    g.drop(columns="in_deletion").to_csv(os.path.join(a.outdir, "sv_vs_pav_genes.tsv.gz"), sep="\t", index=False,
                                         float_format="%.3f")

    # short printout: absent vs present genes inside deletions, per accession
    p = s.pivot(index="accession", columns="call", values=["genes", "in_deletion", "pct_in_deletion"])
    print(f"\ngenes with >= {a.min_frac:.0%} inside a PASS deletion (<= {a.max_len:,} bp)")
    print(f"{'accession':<16}{'absent':>18}{'present':>18}")
    for acc in ORDER:
        ab = f"{int(p.loc[acc, ('in_deletion', 'absent')])}/{int(p.loc[acc, ('genes', 'absent')])}" \
             f" ({p.loc[acc, ('pct_in_deletion', 'absent')]:.1f}%)" if ("genes", "absent") in p.columns and \
            not np.isnan(p.loc[acc, ("genes", "absent")]) else "no absent"
        pr = f"{int(p.loc[acc, ('in_deletion', 'present')])}/{int(p.loc[acc, ('genes', 'present')])}" \
             f" ({p.loc[acc, ('pct_in_deletion', 'present')]:.2f}%)"
        print(f"{acc:<16}{ab:>18}{pr:>18}")
    print("\nwrote", os.path.join(a.outdir, "sv_vs_pav_summary.tsv"), "and sv_vs_pav_genes.tsv.gz")


if __name__ == "__main__":
    main()

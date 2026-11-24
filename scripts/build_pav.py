#!/usr/bin/env python3
"""Build the gene PAV table from the per-gene mosdepth results (Phase 4).

Usage (from the project root):
    python scripts/build_pav.py [--present-breadth 0.8] [--present-depth 0.2]
                                [--absent-breadth 0.1] [--outdir results/pav]

Input : results/pav/<S>/<S>.q10.{regions,thresholds}.bed.gz  (MAPQ >= 10)
        results/pav/<S>/<S>.q0.{regions,thresholds}.bed.gz   (all MAPQ)
Output: results/pav/pav_long.tsv.gz   one row per gene x accession with the raw metrics and the call
        results/pav/pav_matrix.tsv    genes x accessions, 1 present, 0 absent, 0.5 ambiguous/uncertain
        results/pav/pav_summary.tsv   counts of each call per accession

Rules (breadth = fraction of the gene's bases with >= 1 read, MAPQ >= 10):
    present   : breadth >= present_breadth AND depth >= present_depth x sample median
    absent    : breadth <  absent_breadth  AND the gene also has little coverage with all reads
    uncertain : breadth <  absent_breadth  BUT the gene is covered by low-MAPQ reads
                (repeat / paralog: reads exist but do not map uniquely) -> not called absent
    ambiguous : everything in between
"""
import argparse
import os
import sys

import pandas as pd

# Median depth of chr1-5 per accession, from results/mapping/<S>/<S>.mosdepth.global.dist.txt
# (largest depth with cumulative fraction >= 0.5; the same on all five chromosomes for each sample).
MEDIAN = {
    "Col-0": 11, "An-1": 7, "Kyoto": 8, "Kondara": 9,
    "Cvi-0": 7, "Ler-1": 12, "Se-0": 9, "Ale-Stenar-44-4": 8,
}
N_GENES = 27444  # lines of data/reference/genes.bed


def read_sample(s, q, pav_dir):
    """Return a DataFrame indexed by gene with length, mean depth and breadth for min-MAPQ q."""
    base = os.path.join(pav_dir, s, f"{s}.q{q}")
    reg = pd.read_csv(f"{base}.regions.bed.gz", sep="\t", header=None,
                      names=["chrom", "start", "end", "gene", "depth"])
    thr = pd.read_csv(f"{base}.thresholds.bed.gz", sep="\t", skiprows=1, header=None,
                      names=["chrom", "start", "end", "gene", "b1", "b5"])
    if len(reg) != N_GENES or len(thr) != N_GENES:
        sys.exit(f"{s} q{q}: expected {N_GENES} genes, got {len(reg)} / {len(thr)}")
    if not (reg["gene"].values == thr["gene"].values).all():
        sys.exit(f"{s} q{q}: gene order differs between regions and thresholds")
    out = pd.DataFrame({
        "gene": reg["gene"].values,
        "length": (reg["end"] - reg["start"]).values,
        f"depth_q{q}": reg["depth"].values,
        f"breadth_q{q}": (thr["b1"] / (thr["end"] - thr["start"])).values,
    })
    return out.set_index("gene")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--present-breadth", type=float, default=0.8)
    ap.add_argument("--present-depth", type=float, default=0.2, help="fraction of the sample median")
    ap.add_argument("--absent-breadth", type=float, default=0.1)
    ap.add_argument("--lowmap-breadth", type=float, default=0.5,
                    help="all-reads breadth at or above this marks a q10-absent gene as uncertain")
    ap.add_argument("--indir", default="results/pav", help="folder with the per-sample mosdepth results")
    ap.add_argument("--outdir", default="results/pav", help="where the PAV tables are written")
    a = ap.parse_args()

    frames = []
    for s, med in MEDIAN.items():
        q10 = read_sample(s, 10, a.indir)
        q0 = read_sample(s, 0, a.indir)
        d = q10.join(q0[["depth_q0", "breadth_q0"]])
        d["accession"] = s
        d["median_depth"] = med
        d["norm_depth"] = d["depth_q10"] / med

        present = (d["breadth_q10"] >= a.present_breadth) & (d["norm_depth"] >= a.present_depth)
        low_q10 = d["breadth_q10"] < a.absent_breadth
        lowmap = d["breadth_q0"] >= a.lowmap_breadth
        d["call"] = "ambiguous"
        d.loc[present, "call"] = "present"
        d.loc[low_q10 & ~lowmap, "call"] = "absent"
        d.loc[low_q10 & lowmap, "call"] = "uncertain"
        frames.append(d.reset_index())

    long = pd.concat(frames, ignore_index=True)
    code = {"present": 1.0, "absent": 0.0, "ambiguous": 0.5, "uncertain": 0.5}
    long["pav"] = long["call"].map(code)

    os.makedirs(a.outdir, exist_ok=True)
    long.to_csv(os.path.join(a.outdir, "pav_long.tsv.gz"), sep="\t", index=False, float_format="%.4f")
    matrix = long.pivot(index="gene", columns="accession", values="pav")[list(MEDIAN)]
    matrix.to_csv(os.path.join(a.outdir, "pav_matrix.tsv"), sep="\t")
    summary = long.groupby(["accession", "call"]).size().unstack(fill_value=0)
    summary = summary.reindex(list(MEDIAN))[["present", "ambiguous", "uncertain", "absent"]]
    summary["absent_pct"] = (100 * summary["absent"] / N_GENES).round(2)
    summary.to_csv(os.path.join(a.outdir, "pav_summary.tsv"), sep="\t")
    print(summary.to_string())


if __name__ == "__main__":
    main()

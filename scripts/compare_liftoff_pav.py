#!/usr/bin/env python3
"""Phase 6: compare the read-based PAV calls with Liftoff on 4 assemblies (env `map`).

Run after scripts/run_liftoff.sh for An-1, Cvi-0, Kyoto and Ler-1:
    python scripts/compare_liftoff_pav.py

Inputs (read only): results/assembly/liftoff/<S>.gff and <S>.unmapped.txt, results/pav/pav_long.tsv.gz,
                    data/reference/genes.bed, the assembly FASTA files (headers only, to name the chromosomes)
Outputs: results/assembly/liftoff_vs_pav_genes.tsv.gz (one row per gene and accession)
         results/assembly/liftoff_vs_pav_summary.tsv (counts per accession, read call and assembly class)

Assembly class of each gene (Liftoff run with -a 0.1 -s 0.5, gene spans only):
  missing    Liftoff could not place >= 10% of the gene
  partial    placed, but < --min-cov of the gene aligns (part of the gene may be deleted)
  elsewhere  placed with >= --min-cov, but NOT on the same chromosome within --window bp of its Col-0 position
             (often a similar gene / paralog, not the gene itself; can also be a real rearrangement)
  found      placed with >= --min-cov, on the same chromosome near the Col-0 position
"Absent in the assembly" (strict) = missing; (broad) = missing + partial + elsewhere.
Precision = share of read-absent genes that are absent in the assembly; recall = share of genes absent in the
assembly that the reads call absent. The Ler assembly is "Ler" on NCBI (Ler-0 or Ler-1 not stated).
Descriptive only: no statistical test.
"""
import argparse
import os
import re

import numpy as np
import pandas as pd

ACC = {"An-1": "data/assemblies/An-1/An-1.fna", "Cvi-0": "data/assemblies/Cvi-0/Cvi-0.fna",
       "Kyoto": "data/assemblies/Kyo/Kyoto.fna", "Ler-1": "data/assemblies/Ler/Ler_GCA_902460285.1.fna"}
REF_CHR = {"NC_003070.9": "1", "NC_003071.7": "2", "NC_003074.8": "3", "NC_003075.7": "4", "NC_003076.8": "5"}
CLASSES = ["found", "elsewhere", "partial", "missing"]
CALLS = ["present", "ambiguous", "uncertain", "absent"]


def chrom_names(fasta):
    """Assembly sequence name -> chromosome number ('1'..'5'), or None for unplaced contigs."""
    out = {}
    with open(fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split()[0]
                m = re.search(r"chromosome: (\d)", line)
                out[name] = m.group(1) if m else None
    assert sum(v is not None for v in out.values()) == 5, f"expected 5 chromosomes in {fasta}"
    return out


def read_liftoff(gff, unmapped):
    rows = []
    for line in open(gff):
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "gene":
            continue
        g = re.search(r"locus_tag=([^;]+)", f[8]).group(1)
        rows.append((g, f[0], int(f[3]), float(re.search(r"coverage=([^;]+)", f[8]).group(1)),
                     float(re.search(r"sequence_ID=([^;]+)", f[8]).group(1))))
    d = pd.DataFrame(rows, columns=["gene", "seq", "start", "coverage", "identity"])
    d = d.sort_values("coverage", ascending=False).drop_duplicates("gene")  # keep the best placement
    un = pd.DataFrame({"gene": [l.strip().replace("gene-", "") for l in open(unmapped) if l.strip()]})
    assert not set(un["gene"]) & set(d["gene"]), "a gene is both placed and unmapped"
    return pd.concat([d, un], ignore_index=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--liftoff-dir", default="results/assembly/liftoff")
    ap.add_argument("--pav", default="results/pav/pav_long.tsv.gz")
    ap.add_argument("--bed", default="data/reference/genes.bed")
    ap.add_argument("--min-cov", type=float, default=0.5)
    ap.add_argument("--window", type=float, default=2e6, help="bp; how far from the Col-0 position still counts")
    ap.add_argument("--outdir", default="results/assembly")
    a = ap.parse_args()

    bed = pd.read_csv(a.bed, sep="\t", header=None, usecols=[0, 1, 3], names=["ref_chr", "ref_start", "gene"])
    bed["ref_chr"] = bed["ref_chr"].map(REF_CHR)
    calls = pd.read_csv(a.pav, sep="\t", usecols=["gene", "accession", "call"])

    parts = []
    for acc, fasta in ACC.items():
        names = chrom_names(fasta)
        lo = read_liftoff(os.path.join(a.liftoff_dir, f"{acc}.gff"), os.path.join(a.liftoff_dir, f"{acc}.unmapped.txt"))
        assert len(lo) == len(bed) and set(lo["gene"]) == set(bed["gene"]), f"{acc}: Liftoff genes differ from genes.bed"
        lo["asm_chr"] = lo["seq"].map(names)
        d = lo.merge(bed, on="gene").merge(calls[calls["accession"] == acc], on="gene")
        assert len(d) == len(bed)
        placed = d["seq"].notna()
        near = (d["asm_chr"] == d["ref_chr"]) & ((d["start"] - d["ref_start"]).abs() <= a.window)
        d["asm_class"] = np.select(
            [~placed, d["coverage"] < a.min_cov, ~near],
            ["missing", "partial", "elsewhere"], default="found")
        parts.append(d)
    g = pd.concat(parts, ignore_index=True)

    os.makedirs(a.outdir, exist_ok=True)
    g[["gene", "accession", "call", "asm_class", "coverage", "identity", "seq", "start", "ref_chr", "ref_start"]] \
        .to_csv(os.path.join(a.outdir, "liftoff_vs_pav_genes.tsv.gz"), sep="\t", index=False, float_format="%.3f")
    s = g.groupby(["accession", "call", "asm_class"]).size().rename("genes").reset_index()
    s.to_csv(os.path.join(a.outdir, "liftoff_vs_pav_summary.tsv"), sep="\t", index=False)

    print(f"assembly classes: missing = not placed; partial = coverage < {a.min_cov}; "
          f"elsewhere = not within {a.window / 1e6:.0f} Mb of the Col-0 position\n")
    for acc in ACC:
        d = g[g["accession"] == acc]
        t = pd.crosstab(d["call"], d["asm_class"]).reindex(index=CALLS, columns=CLASSES, fill_value=0)
        print(f"== {acc}")
        print(t.to_string())
        ra = d["call"] == "absent"
        strict = d["asm_class"] == "missing"
        broad = d["asm_class"] != "found"
        print(f"   read-absent {ra.sum()}: precision strict {(ra & strict).sum() / ra.sum():.2f}, "
              f"broad {(ra & broad).sum() / ra.sum():.2f}")
        print(f"   recall: of {strict.sum()} missing, {(ra & strict).sum() / strict.sum():.2f} called absent, "
              f"{((d['call'] == 'present') & strict).sum()} called present; "
              f"of {broad.sum()} broad, {(ra & broad).sum() / broad.sum():.2f} called absent, "
              f"{((d['call'] == 'present') & broad).sum()} called present\n")
    print("wrote", os.path.join(a.outdir, "liftoff_vs_pav_genes.tsv.gz"), "and liftoff_vs_pav_summary.tsv")


if __name__ == "__main__":
    main()

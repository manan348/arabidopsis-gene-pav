#!/usr/bin/env python3
"""pN/pS per gene: are NLR resistance genes under different selection than other genes? Env `map`.

Run:
    python scripts/plot_pnps.py

Inputs (read only): data/reference/TAIR10.1_Col-0/GCF_000001735.4/genomic.gff and *_genomic.fna (CDS),
    results/variants/annotated.vcf.gz (SnpEff), results/pav/pav_long.tsv.gz, results/analysis/nlr_genes.tsv
Outputs: results/variants/pnps_genes.tsv, results/figures/38_function_nlr_selection__pnps_dot_plot_ci_scatter_plot.png

Sites: synonymous (S) and non-synonymous (N) sites of the longest RefSeq CDS of each gene, Nei-Gojobori counting
(each codon position: share of the 3 possible changes that keep the amino acid; changes to a stop count as N).
Polymorphisms: biallelic SNPs where at least one accession carries the alternative allele, classed by SnpEff on
the gene's .1 transcript (synonymous_variant = S; missense, stop_gained, stop_lost, start_lost = N).
pN = N SNPs / N sites, pS = S SNPs / S sites. Group values are ratios of sums (sum pN-SNPs / sum N-sites over
genes, etc.), with 95% bootstrap intervals over genes (2,000 resamples). Test NLR vs other core genes: permutation
of group labels (10,000), two-sided, on the difference in group pN/pS.
Caveats: polymorphism within one species (pN/pS), not divergence (dN/dS); 8 accessions give few SNPs per gene,
so single-gene values are noisy and only group values are interpreted; SnpEff and RefSeq transcripts can differ.
"""
import argparse
import gzip
import os
import re
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ACC = ["Col-0", "An-1", "Kyoto", "Kondara", "Cvi-0", "Ler-1", "Se-0", "Ale-Stenar-44-4"]
CHROMS = {"NC_003070.9", "NC_003071.7", "NC_003074.8", "NC_003075.7", "NC_003076.8"}
N_EFF = ("missense_variant", "stop_gained", "stop_lost", "start_lost")
BASES = "TCAG"
AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODE = {a + b + c: AA[16 * i + 4 * j + k] for i, a in enumerate(BASES) for j, b in enumerate(BASES)
        for k, c in enumerate(BASES)}
COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
BLUE, ORANGE, LIGHT, INK, INK2, GRID, SURFACE = "#2a78d6", "#eb6834", "#c3c2b7", "#0b0b0b", "#52514e", "#e1e0d9", "#fcfcfb"
plt.rcParams.update({
    "font.size": 9, "axes.edgecolor": LIGHT, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
    "text.color": INK, "axes.titlesize": 10, "axes.titleweight": "bold", "axes.titlelocation": "left",
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
})


def syn_sites(codon):
    aa = CODE[codon]
    s = 0.0
    for p in range(3):
        for b in "ACGT":
            if b != codon[p]:
                alt = CODE[codon[:p] + b + codon[p + 1:]]
                if alt == aa:
                    s += 1 / 3
    return s


SYN = {c: syn_sites(c) for c in CODE if CODE[c] != "*"}


def read_fasta(path, keep):
    seqs, name, buf = {}, None, []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if name in keep:
                    seqs[name] = "".join(buf).upper()
                name, buf = line[1:].split()[0], []
            elif name in keep:
                buf.append(line.strip())
    if name in keep:
        seqs[name] = "".join(buf).upper()
    return seqs


def cds_sites(gff, fna):
    tx = defaultdict(list)
    gene_of, strand_of = {}, {}
    with open(gff) as f:
        for line in f:
            if "\tCDS\t" not in line:
                continue
            c = line.split("\t")
            if c[0] not in CHROMS:
                continue
            m = re.search(r"TAIR:(AT[1-5]G\d{5})", c[8])
            par = re.search(r"Parent=([^;\n]+)", c[8])
            if not (m and par):
                continue
            tx[par.group(1)].append((int(c[3]), int(c[4])))
            gene_of[par.group(1)] = (m.group(1), c[0])
            strand_of[par.group(1)] = c[6]
    best = {}
    for t, segs in tx.items():
        g, ch = gene_of[t]
        L = sum(e - s + 1 for s, e in segs)
        if g not in best or L > best[g][1]:
            best[g] = (t, L, ch)
    genome = read_fasta(fna, CHROMS)
    rows = []
    for g, (t, L, ch) in best.items():
        segs = sorted(tx[t])
        seq = "".join(genome[ch][s - 1:e] for s, e in segs)
        if strand_of[t] == "-":
            seq = seq.translate(COMP)[::-1]
        if len(seq) % 3:
            continue
        S = 0.0
        ok = 0
        for i in range(0, len(seq) - 3, 3):  # last codon = stop
            cd = seq[i:i + 3]
            if cd in SYN:
                S += SYN[cd]
                ok += 3
        rows.append((g, S, ok - S, len(seq)))
    return pd.DataFrame(rows, columns=["gene", "S_sites", "N_sites", "cds_len"])


def count_snps(vcf):
    cnt = defaultdict(lambda: [0, 0])  # gene -> [N, S]
    with gzip.open(vcf, "rt") as f:
        for line in f:
            if line[0] == "#" or ("synonymous_variant" not in line and "missense_variant" not in line
                                  and "stop_" not in line and "start_lost" not in line):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c[3]) != 1 or len(c[4]) != 1:
                continue
            if not any(x.split(":")[0] == "1" for x in c[9:]):
                continue
            ann = re.search(r"ANN=([^;\t]+)", c[7])
            if not ann:
                continue
            seen = {}
            for e in ann.group(1).split(","):
                f2 = e.split("|")
                if len(f2) < 7 or not f2[6].endswith(".1"):
                    continue
                eff, gid = f2[1], f2[4]
                if "synonymous_variant" in eff and not any(x in eff for x in N_EFF):
                    seen.setdefault(gid, 1)
                elif any(x in eff for x in N_EFF):
                    seen[gid] = 0
            for gid, k in seen.items():
                cnt[gid][k] += 1
    return pd.DataFrame([(g, v[0], v[1]) for g, v in cnt.items()], columns=["gene", "N_snps", "S_snps"])


def ratio(d):
    return (d["N_snps"].sum() / d["N_sites"].sum()) / (d["S_snps"].sum() / d["S_sites"].sum())


def boot(d, rng, n=2000):
    a = d[["N_snps", "N_sites", "S_snps", "S_sites"]].to_numpy()
    idx = rng.integers(0, len(a), size=(n, len(a)))
    s = a[idx].sum(axis=1)
    r = (s[:, 0] / s[:, 1]) / (s[:, 2] / s[:, 3])
    return np.percentile(r, [2.5, 97.5])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gff", default="data/reference/TAIR10.1_Col-0/GCF_000001735.4/genomic.gff")
    ap.add_argument("--fna", default="data/reference/TAIR10.1_Col-0/GCF_000001735.4/GCF_000001735.4_TAIR10.1_genomic.fna")
    ap.add_argument("--vcf", default="results/variants/annotated.vcf.gz")
    ap.add_argument("--pav", default="results/pav/pav_long.tsv.gz")
    ap.add_argument("--nlr", default="results/analysis/nlr_genes.tsv")
    ap.add_argument("--table", default="results/variants/pnps_genes.tsv")
    ap.add_argument("--out", default="results/figures/38_function_nlr_selection__pnps_dot_plot_ci_scatter_plot.png")
    a = ap.parse_args()
    rng = np.random.default_rng(1)

    sites = cds_sites(a.gff, a.fna)
    print(f"CDS sites for {len(sites):,} genes")
    snps = count_snps(a.vcf)
    print(f"coding SNPs: N {snps['N_snps'].sum():,}, S {snps['S_snps'].sum():,} in {len(snps):,} genes")
    d = sites.merge(snps, on="gene", how="left").fillna({"N_snps": 0, "S_snps": 0})
    p = pd.read_csv(a.pav, sep="\t", usecols=["gene", "call"])
    nmiss = (p["call"] == "absent").groupby(p["gene"]).sum()
    d = d[d["gene"].isin(nmiss.index)].copy()
    d["n_missing"] = d["gene"].map(nmiss).astype(int)
    d["nlr"] = d["gene"].isin(set(pd.read_csv(a.nlr, sep="\t")["gene"]))
    d["pN"] = d["N_snps"] / d["N_sites"]
    d["pS"] = d["S_snps"] / d["S_sites"]
    d.to_csv(a.table, sep="\t", index=False, float_format="%.6g")

    groups = {"core genes (not NLR)": d[(~d["nlr"]) & (d["n_missing"] == 0)],
              "variable genes (not NLR)": d[(~d["nlr"]) & (d["n_missing"] > 0)],
              "NLR genes": d[d["nlr"]]}
    res = {}
    for k, v in groups.items():
        r, (lo, hi) = ratio(v), boot(v, rng)
        res[k] = (r, lo, hi, len(v), int(v["N_snps"].sum()), int(v["S_snps"].sum()))
        print(f"  {k:<26} pN/pS {r:.3f} (95% CI {lo:.3f}-{hi:.3f}), {len(v):,} genes, N {res[k][4]:,}, S {res[k][5]:,}")
    A, B = groups["NLR genes"], groups["core genes (not NLR)"]
    obs = ratio(A) - ratio(B)
    pool = pd.concat([A, B]).reset_index(drop=True)
    arr = pool[["N_snps", "N_sites", "S_snps", "S_sites"]].to_numpy()
    nA, diffs = len(A), []
    for _ in range(10000):
        perm = rng.permutation(len(arr))
        x, y = arr[perm[:nA]].sum(0), arr[perm[nA:]].sum(0)
        diffs.append((x[0] / x[1]) / (x[2] / x[3]) - (y[0] / y[1]) / (y[2] / y[3]))
    pval = (np.sum(np.abs(diffs) >= abs(obs)) + 1) / (len(diffs) + 1)
    print(f"NLR vs core: difference {obs:.3f}, permutation p = {pval:.1e}")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.2), gridspec_kw={"width_ratios": [1, 1.25], "wspace": 0.3})
    cols = [LIGHT, ORANGE, BLUE]
    for i, ((k, (r, lo, hi, n, nN, nS)), c) in enumerate(zip(res.items(), cols)):
        ax1.errorbar(r, i, xerr=[[r - lo], [hi - r]], fmt="o", ms=9, color=c, ecolor=INK2, elinewidth=1.4, capsize=4,
                     markeredgecolor=INK2)
        ax1.text(r, i + 0.18, f"{r:.2f}", ha="center", va="bottom", fontsize=9, fontweight="bold", color=INK)
    ax1.axvline(1, color=INK2, ls=":", lw=1)
    ax1.text(1, 2.45, " pN/pS = 1 (neutral)", fontsize=7.5, color=INK2, va="bottom")
    ax1.set_yticks(range(3))
    ax1.set_yticklabels([f"{k}\n{v[3]:,} genes; {v[4]:,} N / {v[5]:,} S SNPs" for k, v in res.items()],
                        fontsize=8.5)
    ax1.set_ylim(-0.6, 2.7)
    ax1.set_xlim(0, 1.2)
    ax1.set_xlabel("pN/pS (group value, 95% bootstrap CI)")
    ax1.set_title(f"(a) NLR genes: about twice the core-gene pN/pS\n(NLR vs core: permutation p = {pval:.0e})")
    ax1.grid(axis="x", color=GRID, lw=0.5)
    ax1.set_axisbelow(True)

    sel = d[(d["N_snps"] + d["S_snps"]) > 0]
    o = sel[~sel["nlr"]]
    ax2.scatter(o["pS"] * 1e3, o["pN"] * 1e3, s=5, color=LIGHT, alpha=0.5, lw=0, label="other genes")
    nl = sel[sel["nlr"]]
    ax2.scatter(nl["pS"] * 1e3, nl["pN"] * 1e3, s=26, color=BLUE, edgecolor=INK2, lw=0.4, label="NLR genes", zorder=3)
    lim = np.percentile(np.r_[sel["pS"], sel["pN"]] * 1e3, 99.8)
    xs = np.linspace(0, lim, 50)
    for r_, ls in [(1, "-"), (0.3, "--")]:
        ax2.plot(xs, r_ * xs, color=INK2, lw=0.8, ls=ls)
        ax2.text(lim * 0.98, r_ * lim * 0.98, f"pN/pS = {r_}", fontsize=7.5, color=INK2, ha="right", va="bottom")
    ax2.set_xlim(0, lim)
    ax2.set_ylim(0, lim)
    ax2.set_xlabel("pS (synonymous SNPs per 1,000 synonymous sites)")
    ax2.set_ylabel("pN (non-synonymous SNPs per 1,000 sites)")
    ax2.set_title("(b) Each gene (genes with >= 1 coding SNP)")
    ax2.legend(frameon=False, loc="upper left", fontsize=8)
    fig.suptitle("Selection on resistance genes: pN/pS of NLR genes vs other genes (8 accessions)", x=0.01,
                 ha="left", fontsize=12, fontweight="bold", y=1.02)
    fig.text(0.01, 0.955, "pN/pS < 1: protein changes removed (purifying selection). Higher values = weaker constraint "
             "or diversifying selection; all groups stay below 1, so this is not proof of positive selection.",
             fontsize=8.5, color=INK2)
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    fig.savefig(a.out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("wrote", a.out, "and", a.table)


if __name__ == "__main__":
    main()

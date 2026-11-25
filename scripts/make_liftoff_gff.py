#!/usr/bin/env python3
"""Cut the Col-0 RefSeq annotation down to the genes used in the PAV analysis (Phase 6), for Liftoff.

Run (any env with Python 3):  python scripts/make_liftoff_gff.py
Input:  data/reference/TAIR10.1_Col-0/GCF_000001735.4/genomic.gff, data/reference/genes.bed
Output: data/reference/genes_only.gff (default: one line per gene, 27,444 PAV genes, gene span only)
        data/reference/genes_protein_coding.gff (with --full: also mRNA, exon and CDS; about 668,000 lines)
The full file made Liftoff run out of memory on this 3.8 GB machine (it loads every feature into Python),
so the gene-only file is the default. The read-based calls also use the whole gene span.

Keeping only these genes makes Liftoff faster and lighter on RAM (3.8 GB machine) and makes its answer
directly comparable with the read-based calls: same genes, same names (locus_tag, e.g. AT1G01010).
Gene -> mRNA -> exon/CDS are linked through the Parent attribute. Only ID, Parent and locus_tag are kept
in column 9, which makes the file small and Liftoff's database fast to build.
"""
import argparse
import os
import re

GFF = "data/reference/TAIR10.1_Col-0/GCF_000001735.4/genomic.gff"
BED = "data/reference/genes.bed"


def attr(fields9, key):
    m = re.search(r"(?:^|;)" + key + r"=([^;]*)", fields9)
    return m.group(1) if m else None


def slim(line):
    """Keep only ID, Parent and locus_tag in column 9 (the long notes are not needed by Liftoff)."""
    f = line.rstrip("\n").split("\t")
    keep = [f"{k}={attr(f[8], k)}" for k in ("ID", "Parent", "locus_tag") if attr(f[8], k) is not None]
    f[8] = ";".join(keep)
    return "\t".join(f) + "\n"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gff", default=GFF)
    ap.add_argument("--bed", default=BED)
    ap.add_argument("--full", action="store_true",
                    help="also keep mRNA, exon and CDS lines (about 668,000 features: needs far more RAM than 3.8 GB)")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    if a.out is None:
        a.out = "data/reference/genes_protein_coding.gff" if a.full else "data/reference/genes_only.gff"

    wanted = {line.split("\t")[3] for line in open(a.bed) if line.strip()}
    gene_ids, tx_ids, kept = set(), set(), []
    counts = {}
    with open(a.gff) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            typ, at = f[2], f[8]
            if typ == "gene":
                if attr(at, "locus_tag") in wanted and attr(at, "gene_biotype") == "protein_coding":
                    gene_ids.add(attr(at, "ID"))
                    kept.append(line)
                    counts[typ] = counts.get(typ, 0) + 1
            elif typ == "mRNA":
                if attr(at, "Parent") in gene_ids:
                    tx_ids.add(attr(at, "ID"))
                    kept.append(line)
                    counts[typ] = counts.get(typ, 0) + 1
            elif typ in ("exon", "CDS"):
                if attr(at, "Parent") in tx_ids:
                    kept.append(line)
                    counts[typ] = counts.get(typ, 0) + 1
    assert len(gene_ids) == len(wanted), f"gene count {len(gene_ids)} differs from the PAV gene list {len(wanted)}"
    if not a.full:  # gene spans only: all that is needed to ask "is this gene present in the assembly?"
        kept = [x for x in kept if x.split("\t")[2] == "gene"]
        counts = {"gene": counts["gene"]}
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    with open(a.out, "w") as out:
        out.write("##gff-version 3\n")
        out.writelines(slim(x) for x in kept)
    print("kept:", counts)
    print("wrote", a.out)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Prepare genomes for SyRI/plotsr: keep chromosomes 1-5 only and rename them Chr1-Chr5 (Phase 7). Env `map`.

Run:
    python scripts/prep_syri.py

Inputs (read only): the TAIR10.1 reference and the 4 assemblies (paths below).
Outputs: data/syri/{Col-0,An-1,Kyoto,Ler-1,Cvi-0}.fa (about 120 MB each; data/syri/ is not tracked by git).
SyRI needs the same chromosome names in both genomes; unplaced contigs, mitochondrion and chloroplast are dropped.
No chromosome is reversed: the orientation check showed 0-4% of each chromosome flipped (local inversions only).
"""
import os
import re

GENOMES = {
    "Col-0": "data/reference/TAIR10.1_Col-0/GCF_000001735.4/GCF_000001735.4_TAIR10.1_genomic.fna",
    "An-1": "data/assemblies/An-1/An-1.fna",
    "Kyoto": "data/assemblies/Kyo/Kyoto.fna",
    "Ler-1": "data/assemblies/Ler/Ler_GCA_902460285.1.fna",
    "Cvi-0": "data/assemblies/Cvi-0/Cvi-0.fna",
}
PAT = re.compile(r"chromosome:? ([1-5])(?:\D|$)")
OUT = "data/syri"


def main():
    os.makedirs(OUT, exist_ok=True)
    for name, path in GENOMES.items():
        out = f"{OUT}/{name}.fa"
        if os.path.exists(out):
            print(f"{out} exists, skipped")
            continue
        seen, lengths, keep, cur = [], {}, False, None
        with open(path) as fi, open(out + ".tmp", "w") as fo:
            for line in fi:
                if line.startswith(">"):
                    m = PAT.search(line)
                    keep = bool(m) and "mitochondri" not in line and "chloroplast" not in line
                    if keep:
                        cur = f"Chr{m.group(1)}"
                        assert cur not in seen, f"{name}: {cur} found twice"
                        seen.append(cur)
                        lengths[cur] = 0
                        fo.write(f">{cur}\n")
                elif keep:
                    fo.write(line)
                    lengths[cur] += len(line.strip())
        assert sorted(seen) == [f"Chr{i}" for i in range(1, 6)], f"{name}: got {seen}"
        os.replace(out + ".tmp", out)
        print(f"{name:<6} " + ", ".join(f"{c} {lengths[c] / 1e6:.1f} Mb" for c in sorted(lengths))
              + f" | total {sum(lengths.values()) / 1e6:.1f} Mb -> {out}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Turn GenMap --raw -fs output into a BED of bisulfite-mappable bases.

Input: the converted FASTA from bisulfite_convert.py (record order defines the
layout) and GenMap's single `<prefix>.freq8` file, one uint8 k-mer frequency per
sequence position: position i of NAME is the k-mer starting at i on the
forward strand, position j of NAME__rc the k-mer starting at j on the reverse
complement.  The last K-1 positions of each record are 0 (no full k-mer), and a
0 elsewhere means the k-mer contains N.  Only frequency == 1 counts as unique.

A read-sized window [s, s+K) is *uniquely mappable* if its top-strand k-mer
(NAME[s]) or its bottom-strand k-mer (NAME__rc[L-K-s]) has frequency 1.  A base
is *mappable* if at least one uniquely mappable window covers it.  The BED
holds the maximal runs of mappable bases (0-based, half-open).

    python mappability/genmap_to_bed.py ct.fa map.freq8 OUT.bed --K 100
"""
import argparse
import sys

import numpy as np


def fasta_lengths(path):
    names, lens = [], []
    with open(path, "rb") as fh:
        for line in fh:
            if line.startswith(b">"):
                names.append(line[1:].split()[0].decode())
                lens.append(0)
            else:
                lens[-1] += len(line.strip())
    return names, lens


def runs(mask):
    """Start/end (half-open) of each True run in a boolean array."""
    d = np.diff(np.concatenate(([False], mask, [False])).astype(np.int8))
    return np.flatnonzero(d == 1), np.flatnonzero(d == -1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("fasta")
    ap.add_argument("freq8")
    ap.add_argument("out")
    ap.add_argument("--K", type=int, default=100)
    args = ap.parse_args()
    K = args.K

    names, lens = fasta_lengths(args.fasta)
    freq = np.memmap(args.freq8, dtype=np.uint8, mode="r")
    if freq.size != sum(lens):
        sys.exit(f"freq8 has {freq.size} values but the FASTA has {sum(lens)} bases")
    offs = np.concatenate(([0], np.cumsum(lens)))
    index = {n: i for i, n in enumerate(names)}

    with open(args.out, "w") as bed:
        for n in names:
            if n.endswith("__rc"):
                continue
            L = lens[index[n]]
            f = freq[offs[index[n]]:offs[index[n]] + L]
            r = freq[offs[index[n + "__rc"]]:offs[index[n + "__rc"]] + L]
            nw = L - K + 1
            if nw <= 0:
                continue
            top = f[:nw] == 1
            bottom = r[:nw][::-1] == 1            # window start s <-> rc start L-K-s
            uniq = top | bottom
            c = np.concatenate(([0], np.cumsum(uniq, dtype=np.int32)))
            b = np.arange(L)
            lo = np.maximum(b - K + 1, 0)
            hi = np.minimum(b, nw - 1)
            mappable = (c[hi + 1] - c[lo]) > 0
            s, e = runs(mappable)
            for a, z in zip(s, e):
                bed.write(f"{n}\t{a}\t{z}\n")


if __name__ == "__main__":
    main()

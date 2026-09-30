#!/usr/bin/env python3
"""Relabel a haploid simulation VCF so REF is one sample's allele, like a reference genome.

The dnAging simulated VCFs write the true ancestral allele as REF. Real maize SNPs
are instead labelled against the B73 reference, so REF is B73's allele whether it
is ancestral or derived. This picks one modern haplotype at random (seeded, or
named with --sample) and, at every site where that haplotype carries ALT, swaps
REF/ALT and flips all genotypes. Allele identities and counts are unchanged; only
the labels move. The chosen sample and the number of swapped sites are written to
a JSON sidecar.
"""

from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path

import numpy as np


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--vcf", type=Path, required=True, help="haploid GT input VCF (.vcf.gz)")
    p.add_argument("--output", type=Path, required=True, help="new output .vcf.gz")
    choose = p.add_mutually_exclusive_group(required=True)
    choose.add_argument("--seed", type=int, help="seed for choosing the reference sample")
    choose.add_argument("--sample", help="name of the reference sample")
    a = p.parse_args(argv)
    if a.output.exists():
        raise SystemExit(f"output already exists: {a.output}")

    sidecar = a.output.with_name(a.output.name.replace(".vcf.gz", "") + ".reference.json")
    swapped = total = 0
    flip = str.maketrans("01", "10")
    with gzip.open(a.vcf, "rt") as src, gzip.open(a.output, "wt") as out:
        for line in src:
            if line.startswith("##"):
                out.write(line)
                continue
            fields = line.rstrip("\n").split("\t")
            if line.startswith("#"):
                samples = fields[9:]
                if a.sample is None:
                    reference = samples[int(np.random.default_rng(a.seed).integers(len(samples)))]
                else:
                    if a.sample not in samples:
                        raise SystemExit(f"sample {a.sample!r} is not in {a.vcf}")
                    reference = a.sample
                column = 9 + samples.index(reference)
                out.write(line)
                continue
            total += 1
            gt = fields[column]
            if gt not in ("0", "1"):
                raise SystemExit(f"unexpected GT {gt!r} at {fields[0]}:{fields[1]}")
            if gt == "1":
                fields[3], fields[4] = fields[4], fields[3]
                fields[9:] = [g.translate(flip) for g in fields[9:]]
                swapped += 1
            out.write("\t".join(fields) + "\n")
    sidecar.write_text(json.dumps({
        "input": str(a.vcf), "reference_sample": reference, "seed": a.seed,
        "sites": total, "swapped_sites": swapped,
    }, indent=2) + "\n")
    print(f"reference {reference}: swapped {swapped:,} of {total:,} sites")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

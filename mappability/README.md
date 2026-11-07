# Bisulfite mappability of B73 v5

A per-base mask of where 100 bp bisulfite-sequencing reads can map uniquely to the
maize B73 v5 reference (`Zm-B73-REFERENCE-NAM-5.0`). It is computed from exact
k-mer uniqueness in bisulfite-converted space, not from simulated alignments.

**Result:** 1,018,575,920 of the 2,182,075,994 reference bases (46.7%) are covered
by at least one uniquely mappable 100-mer (47.6% for chr1-10 alone), with up to 2
mismatches allowed and reads from either strand.

## Data

[`B73v5_bisulfite_k100_e2.mappable.bed.gz`](B73v5_bisulfite_k100_e2.mappable.bed.gz):
the mask, 2,013,753 intervals of mappable bases as a 3-column BED (0-based,
half-open) covering chr1-10 and all scaffolds. About 15 MB gzipped.

Reference: [B73 v5 (Zm-B73-REFERENCE-NAM-5.0)](https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0.fa.gz)
from MaizeGDB.

## Scripts

| Script | Role |
|---|---|
| [`run_mappability.sbatch`](run_mappability.sbatch) | Slurm driver for the whole pipeline: convert, index, map, make BED. |
| [`bisulfite_convert.py`](bisulfite_convert.py) | Writes the C->T converted reference (forward and reverse complement). |
| [`genmap_to_bed.py`](genmap_to_bed.py) | Turns GenMap's raw frequencies into the BED and summary table. |
| [`environment.yml`](environment.yml) | Separate conda env holding GenMap 1.3.0 (not part of the project env). |

## What we did

1. **Defined mappability for bisulfite reads.** After conversion, C and T are
   indistinguishable in a read, so a read maps ambiguously if its C->T-collapsed
   sequence occurs elsewhere in the C->T-collapsed genome. Reads come from both
   strands, so the search set is the C->T version of every sequence and of its
   reverse complement (records `NAME` and `NAME__rc`).
2. **Converted the reference** with `bisulfite_convert.py` (N and other non-ACGT
   become N).
3. **Counted k-mer occurrences genome-wide** with GenMap: `genmap index`, then
   `genmap map -K 100 -E 2 -nc`. `-K 100` is the read length. `-E 2` allows up to
   2 mismatches: a 100-mer counts as a repeat if it occurs elsewhere in the
   converted genome with at most 2 positions differing, so near-copies count as
   multi-mapping (indels are not modelled). `-nc` (no reverse complement) is
   required because both strands are already in the FASTA; GenMap's default would
   also search the reverse complement of the converted sequence, which is not a
   bisulfite read.
4. **Ran the whole genome** (Slurm job 40473019, 16 CPUs, 96 GB requested,
   2 h 11 min total). Peak memory was not recorded.
5. **Made the BED.** A 100 bp window is *uniquely mappable* if its top-strand
   k-mer, or its bottom-strand k-mer (the matching `__rc` position), occurs
   exactly once. A base is *mappable* if at least one uniquely mappable window
   covers it. The BED holds the maximal runs of mappable bases.

## Reproduce

```bash
conda env create -p .envs/mappability -f mappability/environment.yml   # once
sbatch mappability/run_mappability.sbatch                              # K=100, E=2
# other settings: sbatch --export=ALL,K=150,E=0 mappability/run_mappability.sbatch
```

## Choices and caveats

- **Mismatches (`-E 2`) and read length (`-K 100`)** are parameters, not measured
  properties of any real library. Change them to match the data.
- **Either-strand rule.** A base counts as mappable if a unique read could cover it
  from either strand. A stricter mask (both strands) would shrink the BED.
- **Not an aligner.** Frequency 1 means no other converted k-mer within the
  mismatch limit, not that Bismark or bowtie2 will place the read. Soft clipping,
  indels and MAPQ cutoffs are not modelled.
- **Directional libraries** (original top and bottom strands) are assumed; PBAT or
  non-directional libraries add the other two conversions.
- Positions whose k-mers contain N have frequency 0 and are unmappable.
- GenMap saturates frequencies at 255; only frequency 1 is used here.

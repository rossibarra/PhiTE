# Bisulfite mappability of B73 v5

Per-base mask of where 100 bp bisulfite reads can map uniquely to B73
(`Zm-B73-REFERENCE-NAM-5.0`), computed from exact k-mer uniqueness in
bisulfite-converted space rather than from simulated alignments.

## Method

1. **Convert** (`bisulfite_convert.py`). Bisulfite reads, once read-space C and
   T are collapsed, are k-mers of the genome with every C written as T. Reads
   come from both strands of the original DNA, so the search space is the
   C->T version of each sequence *and* of its reverse complement (records
   `NAME` and `NAME__rc`). N and other non-ACGT characters become N.
2. **Index and count** (GenMap 1.3.0). `genmap index` builds the FM index of the
   converted FASTA. `genmap map -K 100 -E 2 -nc` then counts, for every
   k-mer start, how many places in the whole converted set it occurs with up to
   2 mismatches. `-nc` (no reverse complement) is essential: both strands are
   already in the FASTA, and GenMap's default would also search the reverse
   complement of the converted sequence, which is not a bisulfite read.
3. **Make the BED** (`genmap_to_bed.py`). A 100 bp window is *uniquely
   mappable* if its top-strand k-mer or its bottom-strand k-mer (the matching
   `__rc` position) has frequency 1. A base is *mappable* if at least one
   uniquely mappable window covers it. The BED holds the maximal runs of
   mappable bases.

## Run

```bash
conda env create -p .envs/mappability -f mappability/environment.yml   # once
sbatch mappability/run_mappability.sbatch                              # K=100, E=2
# other settings: sbatch --export=ALL,K=150,E=0 mappability/run_mappability.sbatch
```

The job uses 16 CPUs and 96 GB (the suffix-array build needs about 10 bytes
per base for inputs over 2 GB). It must run genome-wide: uniqueness judged
within one chromosome overstates mappability.

## Output (`results/mappability/`)

- `B73v5_bisulfite_k100_e2.mappable.bed`: 3-column BED (0-based, half-open),
  all of chr1-10 and the scaffolds in the reference.
  A gzipped copy of the K=100, E=2 result is committed as
  `mappability/B73v5_bisulfite_k100_e2.mappable.bed.gz` (about 15 MB; 2,013,753
  intervals, 46.7% of the 2.18 Gb reference mappable).
- `B73v5_bisulfite_k100_e2.summary.tsv`: per sequence, count of unique top- and
  bottom-strand windows and the mappable fraction.
- `B73v5_k100_e2/`: intermediates (converted FASTA, index, raw frequencies).
  Large; safe to delete once the BED exists.

## Choices and caveats

- **Mismatches (`-E 2`)** and **read length (`-K 100`)** are parameters, not
  measured properties of any real library. Change them to match the data.
- **Either-strand rule.** A base counts as mappable if a unique read could cover
  it from either strand. A stricter mask (both strands) would shrink the BED.
- **Not an aligner.** Frequency 1 means no other converted k-mer within the
  mismatch limit, not that Bismark or bowtie2 will place the read. Soft clipping,
  indels and MAPQ cutoffs are not modelled.
- **Directional libraries** (original top and bottom strands) are assumed;
  PBAT or non-directional libraries add the other two conversions.
- Positions whose k-mers contain N have frequency 0 and are unmappable.
- GenMap saturates frequencies at 255; only frequency 1 is used here.
- Chromosome-only pilot (chr10, same settings): 40.8% of 100-mer starts unique
  and 60.8% of bases mappable. These are not genome-wide values.

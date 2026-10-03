# v0.9.0 validation artifacts

These are compact, version-controlled copies of the reports summarized in
[`../VALIDATION.md`](../VALIDATION.md). Large arrays, match bundles, VCF-derived
tables and simulation intermediates remain outside Git under `results/`.

| validation item | directory | contents |
|---|---|---|
| V1 | [`v1/`](v1/) | per-cell rejection criteria and full JSON report |
| V2 | [`v2/`](v2/) | per-replicate and pooled depletion-capacity summaries |
| V4, original run | [`v4_1001/`](v4_1001/) | criteria, block diagnostics and report for 1,001 sets |
| V4, Amendment A | [`v4_first500/`](v4_first500/) | criteria, block diagnostics and report for the first 500 sets |
| V5 | [`v5/`](v5/) | Phi-SFS summary and provenance metadata |
| V6 | [`v6/`](v6/) | criteria, per-test results, block diagnostics and report |

The V4 first-500 report was regenerated from the unchanged match bundle and
Phi-SFS output with the amended threshold recorded in
`REMAINING_VALIDATION_PROPOSAL.md`:

```bash
python -m tools.v4_depletion_report \
  --matches results/bootstrap_matches/in_gene_rc2_first500 \
  --phi results/phi_sfs/in_gene_rc2_first500 \
  --min-qc-passes 451 \
  --output results/v4/in_gene_rc2_first500
```

The tracked reports are sufficient to audit the numerical claims in the
validation summary. Reproducing them from raw inputs still requires the large
local artifacts named in each report and the commands in `../VALIDATION.md`.

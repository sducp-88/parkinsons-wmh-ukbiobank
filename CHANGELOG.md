# Changelog

## 1.0.0 — 2026-10-07

- Replaced an inherited exploratory notebook with a PD-specific aggregate-only pipeline.
- Repaired diagnosis-column aliases using ICD-10 prefixes; reports absent fields and fails when an entire comorbidity component has no source fields.
- Uses unpenalized propensity estimation, explicit categorical coding and numeric scaling, with positivity/convergence/balance checks.
- Distinguishes weighted WMH models, unweighted primary cognitive/path models and weighted cognitive sensitivity analyses.
- Documents single imputation and the transformed-outcome estimand.
- Adds reported-hypertension sensitivity, deterministic participant bootstrap and synthetic validation tests.
- Removes participant records, notebook outputs, study result files, credentials and machine-specific paths from the public release.

This release documents manuscript-preparation refinements. It is not a claim that every analysis was prospectively specified.

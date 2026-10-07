# Reproducibility and release scope

The script records input SHA-256, local package versions, additional imputation, diagnosis-code availability, propensity diagnostics, effective sample sizes, withdrawal-filter counts and the requested bootstrap seed in a **local** run report. These records are not participant exports, but the report and study aggregate tables are not distributed in this code-only release.

`--keep-legacy-cmc` is solely for a local methodological comparison against an archived prepared score. The final workflow rebuilds available-record comorbidity categories. The release does not redistribute an inherited notebook because notebook outputs can contain participant records and its legacy blocks do not constitute a reproducible PD pipeline.

Public tests generate fresh synthetic participants entirely in memory. They check exposure time order, schema and scaling errors, comorbidity-prefix recovery, balance and weight behavior, cognitive direction, path decomposition identities, bootstrap completion and aggregate-only output. Passing them does not validate a clinical diagnosis, UK Biobank consent status or the biological interpretation of an observational association.

The public repository should remain restricted to code, documentation, dependency metadata, license and synthetic tests. Do not add real-data fixtures, notebooks with outputs, source CSVs, genotypes, exact participant dates, API credentials, withdrawal IDs or manuscript drafts. `.gitignore` is only an additional protection; inspect staged paths before every future push.

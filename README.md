# Parkinson disease and white matter hyperintensities

Code for a cross-sectional UK Biobank analysis of algorithmically recorded Parkinson disease (PD), total and regional white matter hyperintensities (WMH), and cognitive performance. This repository distributes **code, method documentation and synthetic validation tests**. Study results, participant records, genotypes, withdrawal lists and manuscript drafts are retained in the authorized local workspace.

The pipeline accepts an investigator-prepared canonical CSV. It starts after cohort extraction and upstream genotype processing; it is not a complete raw-data extraction or imaging-segmentation workflow. UK Biobank data access requires a separately approved application. Publishing this code does not grant data access.

## Run in an authorized environment

Use Python 3.12 and the versions in `requirements.txt`. Install dependencies into your own analysis environment:

```bash
python -m pip install -r requirements.txt
python -m unittest discover -s tests -v
python pd_wmh.py --input /private/approved/prepared.csv --withdrawn-list /private/approved/current_withdrawals.csv --output /private/approved/aggregate_results --paths --bootstrap 2000 --seed 20261007
python plot_results.py --results /private/approved/aggregate_results --output /private/approved/figures
```

Use a current, application-specific withdrawal list and verify consent status before running. The program cannot verify the currency of a list or your access approval. The input must already satisfy the documented cohort definition. Do not upload the input or withdrawal list to this repository or an external AI service.

See [`docs/input_schema.md`](docs/input_schema.md) for exact input column names and [`docs/methods.md`](docs/methods.md) for estimands and model specifications. A locally retained flow-count JSON can be supplied to `plot_results.py --flow-json` to generate the selection diagram; it needs the aggregate keys documented in the plot function.

## What the pipeline does

- Validates exposure dates, required columns, unique identifiers when present and positive head-size scaling factors.
- Rebuilds available-record comorbidity categories by ICD-10 prefix, retaining a report of missing code fields rather than assuming complete hypertension ascertainment.
- Fits unpenalized, standardized logistic propensity models and checks balance of modelled covariates before estimating full and restricted overlap-weighted WMH contrasts.
- Fits primary unweighted cognitive regressions and weighted sensitivity models, with endpoint-specific available samples.
- Estimates exploratory same-visit statistical path products with a deterministic participant pair bootstrap.
- Writes aggregate CSV tables and a reproducibility report. It does not write participant rows, dates or identifiers. Plotting uses those local aggregates to create PNG, PDF, SVG and 1200-dpi TIFF figures.

## Interpretation

PD is identified from algorithmic report dates, without specialist adjudication of idiopathic or vascular parkinsonism. The neurologic-code restriction is a sensitivity analysis, not a clinical diagnosis. WMH contrasts concern geometric means of **1 + normalized WMH**, because regressions use `log1p`. Cognitive outcomes and WMH are measured at the same visit; statistical path products do not identify causal mediation. HC3 intervals condition on fitted weights and do not incorporate all propensity or single-imputation uncertainty. A small case group remains small despite a large control group.

The code and analysis methods were revised during manuscript preparation. This resource is not a preregistration, a journal submission, or an author-approved publication. See [`CHANGELOG.md`](CHANGELOG.md) and [`docs/reproducibility.md`](docs/reproducibility.md).

## License and attribution

MIT covers the software only. UK Biobank data and access rights are not covered by this license. Cite the software using `CITATION.cff`; a study DOI is not assigned. OpenAI Codex assisted code review/refactoring, local validation and documentation preparation. Human investigators retain responsibility for their analyses and publications.

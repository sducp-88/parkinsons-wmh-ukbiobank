# Statistical methods and interpretation

## WMH comparisons

Each raw MRI volume is multiplied by the released T1 head-size scaling factor, then transformed with `log1p`. A logistic propensity model estimates recorded-PD status using nine main-effect covariates. Categorical variables are dummy coded and design columns standardized. The likelihood is unpenalized. PD weights are `1 - e(X)` and control weights are `e(X)`; no trimming or group normalization is used. Balance is assessed on fitted covariate terms with a pooled, unweighted within-group SD denominator.

Overlap-weighted least squares regressions also include the nine core covariates and HC3 outcome-model standard errors. Within-group effective sample size is `(sum(w))**2 / sum(w**2)`. The exponentiated PD coefficient is a geometric-mean ratio of **1 + head-size-normalized WMH**. `100 * expm1(beta)` is not an exact percentage contrast in arithmetic mean lesion volume.

The full and neurologic-code-restricted propensity models are fitted independently. The primary outcome is total WMH; Benjamini–Hochberg correction is applied to the two regional tests within each model. Conventional unweighted WMH regression is a sensitivity analysis. Additional self-reported hypertension and subgroup analyses retain full-cohort weights.

The ICD I63/I64/G45/G35/I67 restriction identifies recorded diagnoses, not specialist adjudication of vascular parkinsonism. Missing hypertension-source fields, a categorized comorbidity count and residual confounding limit claims of vascular independence.

## Cognitive outcomes

Reaction time and Trail Making Test B are negative-log transformed; fluid intelligence and digit span retain their direction. Each is standardized using the full available cohort mean and sample SD, with higher values meaning better performance. Primary cognitive contrasts are unweighted ordinary least squares with nine core covariates and additive APOE ε4 count. HC3 intervals are used. A weighted cognitive model is supplementary. Each endpoint uses its available sample without outcome imputation.

## Exploratory statistical paths

For each WMH/cognitive pair, the same complete sample supplies adjusted unweighted regressions `M ~ PD + X`, `Y ~ PD + M + X`, and `Y ~ PD + X`. The product of the PD-to-WMH and WMH-to-cognition coefficients is bootstrapped by resampling participants with replacement. Multinomial row multiplicities and exact sufficient cross-products implement ordinary pair-bootstrap resampling; fixed cognitive standardization and covariate scaling are retained. Draws with singular linear systems are skipped, and excessive nonestimability stops the analysis. Each pair must have at least 95% successful draws.

Percentile intervals are exploratory and unadjusted across the path comparisons. These concurrent regressions do not establish time order, causal mediation, a mediated proportion, or intervention efficacy. A nonsignificant overall cognitive contrast does not prove cognitive preservation. Different subgroup P values do not by themselves establish interaction.

## Uncertainty and provenance

HC3 intervals condition on fitted weights; they do not propagate all propensity or single-imputation uncertainty. The analysis begins with a prepared extract, not centrally controlled raw imaging or genotype data. Current access, withdrawals and upstream cohort derivation remain investigator responsibilities. Methods were refined during manuscript preparation and are not retrospectively labelled preregistered.

Primary methodological sources:

- Li F, Morgan KL, Zaslavsky AM. Balancing Covariates via Propensity Score Weighting. JASA 2018. https://doi.org/10.1080/01621459.2016.1260466
- Maxwell SE, Cole DA. Bias in cross-sectional analyses of longitudinal mediation. Psychol Methods 2007. https://doi.org/10.1037/1082-989X.12.1.23
- UK Biobank PD report field https://biobank.ndph.ox.ac.uk/ukb/field.cgi?id=42032
- UK Biobank self-report illness field and coding https://biobank.ndph.ox.ac.uk/ukb/field.cgi?id=20002 and https://biobank.ndph.ox.ac.uk/ukb/coding.cgi?id=6

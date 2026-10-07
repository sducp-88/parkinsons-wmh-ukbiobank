# Canonical prepared input schema

Read-only input: one row per participant in an already extracted, authorized and consent-checked cohort. The CLI does not create the cohort from the full UK Biobank resource. Exact normalized column names are listed below. Study extraction and APOE genotyping remain upstream.

## Required exposure and visit fields

`treatment_var` must be 0 for controls and 1 for recorded PD. If absent, `group` must contain `Control` or `Study` and is mapped accordingly.

| Column | Meaning |
|---|---|
| `Date_of_attending_assessment_centre_Instance_2` | First imaging assessment date |
| `Date_of_parkinson_s_disease_report` | PD-specific algorithm report date |
| `Date_of_all_cause_parkinsonism_report` | All-cause parkinsonism algorithm report date |
| `Participant_ID` | Required for withdrawal filtering; optional uniqueness validation otherwise; never exported |

PD cases require both report dates, with the PD date on or before the imaging date. Controls require both report dates missing. Unknown report-date sentinel `1900-01-01` is rejected for cases. Other diagnostic/date patterns must be resolved upstream, not silently relabelled.

## Required covariates

`Sex`, `Age_at_Instance_2`, `Townsend_deprivation_index_at_recruitment`, `Body_mass_index_BMI_Instance_0`, `Genetic_ethnic_grouping`, `Smoking_Ever`, `Alcohol_intake_frequency_ordinal`, `has_degree`, `CMC_score_cat`, `e4_count`.

Sex and genetic grouping are binary indicators in the prepared study input. Smoking means previous/current versus never. Alcohol frequency increases from never (0) to most frequent (5). Degree is binary. The comorbidity category is 0, 1 or 2 (two or more). `e4_count` is additive 0–2; its upstream allele calling is not reproduced here. Age uses assessment year minus birth year. Prepared missing-value codes must already be converted to missing values. The program completes missing covariates with median/mode single imputation and refuses an entirely missing covariate.

## Required MRI and cognitive fields

- `Volumetric_scaling_from_T1_head_image_to_standard_space_Instance_2`
- `Total_volume_of_white_matter_hyperintensities_from_T1_and_T2_FLAIR_images_Instance_2`
- `Total_volume_of_peri_ventricular_white_matter_hyperintensities_Instance_2`
- `Total_volume_of_deep_white_matter_hyperintensities_Instance_2`
- `Mean_time_to_correctly_identify_matches_Instance_2`
- `Duration_to_complete_alphanumeric_path_trail_2_Instance_2`
- `Fluid_intelligence_score_Instance_2`
- `Maximum_digits_remembered_correctly_Instance_2`

Volumes must be nonnegative, and the scaling factor must be finite and positive. Region-specific and cognitive missingness is allowed. Negative cognitive sentinel values are invalid, and timed tests must be positive. No outcome imputation is performed.

## Diagnosis columns

Columns are matched by prefixes such as `Date_E11_first_reported_`; the human-readable suffix may vary. Each comorbidity component must have at least one available source field. Expected code groups: I10–I12/I13/I15; E78; I47/I48/I49; I20–I25; I50; E10–E14; I60–I64. Dates later than `1900-01-01` and no later than imaging contribute. Missing individual code fields are explicitly reported. A dataset containing only rare hypertension-code columns is not complete hypertension ascertainment.

The neurologic restriction requires prefixes for I63, I64, G45, G35 and I67. Optional reported-hypertension analysis uses `Non_cancer_illness_code_self_reported_Instance_0_Array_*` plus `Date_of_attending_assessment_centre_Instance_0`. Available baseline arrays containing 1065 or 1072 indicate a positive hypertension report. Empty arrays are not clinically confirmed disease absence.

No identifiers, exact participant dates, participant values or withdrawal files should be committed to a repository. The withdrawal input is a private, headerless first column of current withdrawn participant IDs for the appropriate application.

# Immune-metabolic plasma proteomic signatures distinguish asymptomatic from symptomatic SARS-CoV-2 infection and predict post-COVID sequelae

Analysis code for the plasma proteomic analysis of the Gutenberg COVID-19 Study (GCS) and
the Gutenberg Post-COVID Study (GPC).

## Contents

`GCS_proteomics_analysis.R` holds the whole analysis, in the order in which it was run:

| Sections | |
| --- | --- |
| 0–1 | Setup and function definitions, including the penalised regression function, the ensemble mRMR selection frequency and the stability measure |
| 2 | Data preparation: the mass spectrometry data, the cohort data, and the severity, sequelae and covariate variables |
| 3–5 | Baseline table, ensemble mRMR feature selection with its stability assessment, ridge regression and the protein score |
| 6–14 | The analyses reported in the main figures: the signature coefficients, the protein score across severity categories, pre- and post-infection protein dynamics, acute and long-term symptoms, post-COVID sequelae, cardiovascular risk factors and the per-protein comparisons |
| 15 | The analyses added during revision: serologically confirmed infections, the comparison of participants with and without proteomic measurements, sensitivity to the method of infection ascertainment, the objective severity definition based on treatment, and the extended confounder adjustment with its forest plots |

`forest_symptoms_data.csv` and `forest_pcs_data.csv` hold the estimates from the extended
confounder adjustment and are read by the forest plots in section 15f.

## Methods implemented here

- Transformation of each protein chosen to minimise skewness, then winsorisation at six median absolute deviations from the median.
- Ensemble minimum redundancy maximum relevance (mRMR) feature selection over 100 bootstrap iterations, keeping the proteins selected in more than half of them. Selection stability is quantified with the measure of Nogueira et al. and compared against a single mRMR run, and the selection frequency is traced across 5 to 1,000 iterations.
- Ridge regression of the selected proteins on asymptomatic versus symptomatic infection, with age and sex left unpenalised and the penalty parameter chosen by 10-fold cross-validation. The protein score is the weighted sum of the standardised protein levels.
- Poisson regression with heteroscedasticity-consistent (HC0) standard errors for the prevalence ratios of post-COVID sequelae.
- Linear and logistic models for the associations with symptoms, severity and cardiovascular risk factors, under three adjustment sets: age and sex; plus BMI and smoking; plus comorbidities and vaccination before infection.

## Software

R, with: `glmnet`, `mRMRe`, `pROC`, `Hmisc`, `survival`, `doParallel`, `sandwich`,
`dplyr`, `tidyr`, `stringr`, `readr`, `readxl`, `e1071`, `openxlsx`,
`ggplot2`, `ggpubr`, `ggdist`, `ggbeeswarm`, `smplot2`, `forplo`,
`ComplexHeatmap`, `circlize`, `pheatmap`.

## Notes on reading the code

- **Paths and database access.** File paths and database identifiers are given as
  placeholders in angle brackets. The analyses were run against the study database at the
  University Medical Center Mainz.
- **In-house helper functions.** Some calls are wrappers from an in-house function library
  that is not reproduced here: `E_Net1_RM()` around the penalised regression function in
  section 1a, `baseline13()` for the baseline table, `reg_tab()`, `reg_tab_list()` and
  `reg_tab_loop()` for regression tables, `TAB_miss()` for counts of missing values, and
  `plot_save()` for figure output. These are reporting wrappers around standard R models;
  what each one computes is documented where it is called, so that every model can be
  reproduced with base R, `stats` and `sandwich`.
- **The protein score is written out in section 5** and read back in from there on, so that
  the feature selection and the ridge fit do not have to be repeated.

## Data availability

The individual-level data of the Gutenberg COVID-19 Study cannot be made publicly
available because of participant confidentiality and the terms of the ethics approval.
Access is described in the data availability statement of the manuscript.

## Licence

MIT, see `LICENSE`.

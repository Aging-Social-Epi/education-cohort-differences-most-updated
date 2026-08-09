# education-cohort-differences-most-updated
# Educational attainment and cohort differences in dementia prevalence

Analytic code for the counterfactual decomposition analysis reported in:

> *Educational Attainment and Cohort Differences in Dementia Prevalence : A Counterfactual Decomposition Analysis*

Provided to support transparency and reproducibility of the published findings.

## Data availability

Analyses use the Health and Retirement Study (HRS). HRS data are available to
researchers on application: <https://hrs.isr.umich.edu>

Reproducing these results requires the same HRS data releases used here (RAND HRS
Longitudinal File 1992–2022, plus the Langa-Weir and Gianattasio-Power dementia
classification files).

## Repository structure

    01_merge_file.R              Merge raw HRS files into the analytic dataset
    02_clean_variables.R         Construct and clean analysis variables
    03_descriptive_analysis.R    Descriptive statistics (Table 1, eTables 1-4)
    tables_and_figures_all.R     All manuscript tables and figures
    function_education.R         Counterfactual decomposition functions
    function_education_sen.R     Decomposition functions for sensitivity analyses

    main_analysis/               Primary age-period-cohort models and decomposition
      04_APC_model_all.R           Overall
      04_APC_model_{male,female}.R Stratified by sex/gender
      04_APC_model_{white,black,hisp}.R  Stratified by race and ethnicity
      *_GP_only.R                  Gianattasio-Power algorithm

    sensitivity_parental/        Additionally adjusted for parental education
    sensitivity_proxy/           Excluding proxy interviews

    RR_absolute_scale.R          Estimates on the absolute scale
    RR_mc_convergence.R          Monte Carlo and bootstrap convergence checks
    RR_drift_sensitivity.R       Sensitivity to the age-period-cohort parameterization

## Reproducing the analysis

Update the file paths at the top of each script, then run in order:

1. `01_merge_file.R`
2. `02_clean_variables.R`
3. Scripts in `main_analysis/`
4. Scripts in `sensitivity_parental/` and `sensitivity_proxy/`


Dementia classification is estimated under two algorithms. Scripts without a suffix
use the Langa-Weir classification; `*_GP_only.R` scripts use the Gianattasio-Power
algorithm.

## Software

R 4.2.1. Package dependencies are declared at the top of each script.


# education-cohort-differences-most-updated
Educational Attainment and Cohort Differences in Dementia Prevalence: A Counterfactual Decomposition Analysis
This repository contains the analytic code used to produce the results reported in the manuscript:

Educational Attainment and Cohort Differences in Dementia Prevalence: A Counterfactual Decomposition Analysis

The code is provided to support transparency and reproducibility of the published findings.

Data Availability
The analyses use data from the Health and Retirement Study (HRS) and related cohorts.
HRS data are available to researchers upon application.

Information on data access is available at:
https://hrs.isr.umich.edu

Researchers seeking to reproduce the analyses must obtain the same HRS data releases used in this study.

Repository Structure
Top-level scripts

01_merge file.R
02_clean variables.R
03_descriptive analysis.R
Tables_and_figures.R
function_education.R
function_education_sen.R
Main analyses

main_analysis/
04_APC_model_education_all.R
04_APC_model_education_black.R
04_APC_model_education_female.R
04_APC_model_education_hisp.R
04_APC_model_education_male.R
04_APC_model_education_white.R
Sensitivity analyses

sensitivity_analysis/
Additional control for parental education/
using_alternative_algorithm/
Script Description
Data preparation
01_merge file.R
Merges raw HRS files and constructs the analytic dataset.

02_clean variables.R
Cleans key variables.

03_descriptive analysis.R
Produces descriptive statistics reported in the manuscript.

Main analyses
function_education.R
Core functions implementing the counterfactual decomposition framework.

main_analysis/
Scripts generating the primary age–period–cohort models and decomposition results, including analyses stratified by sex/gender and race and ethnicity.

Sensitivity analyses
function_education_sen.R
Functions used in sensitivity analyses.

sensitivity_analysis/Additional control for parental education/
Sensitivity analyses adjusting for parental educational attainment.

Reproducibility Instructions
After obtaining access to the required HRS data:

Update file paths in the scripts as needed.
Run the scripts in the following order:
01_merge file.R
02_clean variables.R
Scripts in main_analysis/
Scripts in sensitivity_analysis/
Run Tables_and_figures.R to reproduce manuscript tables and figures.
All analyses were conducted using R. Package dependencies are specified within the scripts.

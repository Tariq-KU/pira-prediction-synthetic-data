# Synthetic RRMS Cohorts for PIRA Prediction

## Purpose

This repository contains the synthetic datasets used in the report:

"Incremental Value of Serum Biomarkers for Predicting Progression
Independent of Relapse Activity: A Synthetic Cohort Analysis."

The analysis evaluated whether serum neurofilament light chain (sNfL)
and serum glial fibrillary acidic protein (sGFAP) improved prediction
beyond clinical and MRI information.

## Dataset overview

- Development cohort: 1,000 synthetic participants.
- Test cohort: 500 synthetic participants.
- Prediction landmark: study month 12.
- Prediction horizon: the following 24 months, through study month 36.
- Event labels represent presumed qualifying PIRA onsets with
  12-month persistence, with confirmation possible through study month 48.

These are generated records, not data collected from actual patients.

## Files

List each uploaded filename here and briefly explain its contents.

## Analysis overview

Four ridge-penalized Cox models were evaluated:

1. Clinical predictors.
2. Clinical and MRI predictors.
3. Clinical and blood predictors.
4. Clinical, MRI and blood predictors.

The primary comparison evaluated the addition of blood predictors
to the clinical-MRI model.

## Data provenance

The datasets were generated for this research-methods assignment.
Dataset generation and analysis code were developed with assistance
from OpenAI's ChatGPT.

The accompanying documentation describes the variables and
data-generation assumptions.

## Version

Version 1.0.0 identifies the datasets associated with the submitted report.

## Citation

Please use the citation information in CITATION.cff and identify
the version used.

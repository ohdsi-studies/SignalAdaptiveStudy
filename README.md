Signal-Adaptive Ridge and LASSO for Patient-Level Prediction
=============

<img src="https://img.shields.io/badge/Study%20Status-Started-blue.svg" alt="Study Status: Started">

- Analytics use case(s): **Patient-Level Prediction**
- Study type: **Methods Research**
- Tags: **-**
- Study lead: **Egill Fridgeirsson**
- Study lead forums tag: **-**
- Study start date: **-**
- Study end date: **-**
- Protocol: **-**
- Publications: **-**
- Results explorer: **-**

An OHDSI study package that adds signal-adaptive ridge (SAR) and signal-adaptive LASSO (SAL) logistic regression to
[PatientLevelPrediction](https://github.com/OHDSI/PatientLevelPrediction) and compares them with ridge and LASSO on
a set of prediction tasks across OMOP CDM databases. In SAR, each covariate's ridge penalty is scaled by a weight
computed from the training rows of each fit: its prevalence, its information and its marginal association with the
outcome under the intercept-only model. The weight shape and the penalty are chosen by cross-validation, with the
weights recomputed inside every fold. SAL uses the same weights with a LASSO penalty.

- **Running the study at a site:** see [README_PRODUCTION.md](README_PRODUCTION.md) and `extras/CodeToRun.R`.
- **Cohorts:** Circe definitions in `inst/cohorts` (the rendered SQL is shipped, so no Java or CirceR is needed).
- **What leaves a site:** aggregate metrics and counts only, with small counts censored (`exportResults()`).

Status: in development; not yet for production use.

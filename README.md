---
output:
  pdf_document: default
  html_document: default
---
# R Implementations

This repository contains R functions, examples, and applications for robust estimation in different regression models.

## Folder structure

### nlm

This folder contains implementations for nonlinear regression models.

The available files include:

- estimation functions based on the normal maximum likelihood estimator (MLE);
- corrected maximum Lq-likelihood estimator (MLqE);
- minimum density power divergence estimator (MDPDE);
- maximum likelihood estimation under the Student-t distribution;
- commented examples illustrating how to use the estimation functions;
- applications to real datasets.

### lmve

This folder contains implementations for linear models with errors in variables.

The available files include:

- estimation functions based on the normal MLE;
- corrected MLqE;
- MDPDE;
- maximum likelihood estimation under the Student-t distribution;
- commented simulation and usage examples;
- applications to real datasets.

The models account for measurement errors in the explanatory variables, with known observation-specific measurement-error variances.

### lmm

This folder contains implementations for linear mixed models.

The available files include:

- estimation functions based on the normal MLE;
- corrected MLqE;
- MDPDE;
- maximum likelihood estimation under the Student-t distribution;
- commented examples illustrating model specification and estimation;
- applications to longitudinal and grouped data.

The implementations allow different structures for the random-effects covariance matrix.

## File naming

Files beginning with fit_ contain the main estimation functions.

Files beginning with example_ contain commented examples showing how to use the functions.

Files beginning with application_ contain applications to real datasets.

## Estimation methods

The following estimation methods are considered:

- normal maximum likelihood estimation;
- corrected maximum Lq-likelihood estimation;
- minimum density power divergence estimation;
- Student-t maximum likelihood estimation.

The corrected MLqE and MDPDE are robust alternatives to the usual maximum likelihood estimator and reduce the influence of atypical observations.

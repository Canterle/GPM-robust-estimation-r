---
output:
  pdf_document:
    latex_engine: xelatex
  html_document: default
---
# Nonlinear Regression Estimation: MLqE, MDPDE, and Student-t MLE

This repository contains generic R implementations and reproducible examples for robust and heavy-tailed estimation in nonlinear regression models.

The implemented procedures are:

- corrected maximum Lq-likelihood estimation (corrected MLqE) under a normal model;
- maximum density power divergence estimation (MDPDE) under a normal model;
- maximum likelihood estimation under a Student-t model with fixed degrees of freedom;
- normal maximum likelihood estimation as the special case `q = 1` of the MLqE or MDPDE implementations.

The user specifies the nonlinear mean and dispersion structures through a single model function. Homoscedastic and heteroscedastic models are supported.

## Repository contents

| File name | Description |
|---|---|
| `fit_nlm_mlqe.R` | Generic corrected MLqE implementation, including fixed-`q` fitting and automatic `q` selection by SQV. |
| `fit_nlm_mdpde.R` | Generic MDPDE implementation, including fixed-`q` fitting and automatic `q` selection by SQV. |
| `fit_nlm_t.R` | Generic Student-t maximum-likelihood implementation with fixed degrees of freedom. |
| `example1_nlm_commented.R` | Simulated heteroscedastic normal nonlinear model with one contaminated observation. Compares corrected MLqE, MDPDE, and normal MLE. |
| `example2_nlm_commented.R` | Simulated homoscedastic normal nonlinear model with one contaminated observation. Compares corrected MLqE, MDPDE, and normal MLE. |
| `example_nlm_t_commented.R` | Simulated uncontaminated heteroscedastic Student-t nonlinear model with four degrees of freedom, fitted by Student-t MLE. |
| `application_radioimmunoassay.R` | Application to the radioimmunoassay data, comparing Student-t MLE, MDPDE, corrected MLqE, normal MLE, and normal MLE after removing the identified outlier. |

The example and application scripts call:

```r
source("fit_nlm_mlqe.R")
source("fit_nlm_mdpde.R")
source("fit_nlm_t.R")
```

Keep the implementation files under these names or modify the `source()` calls accordingly.

## Requirements

The scripts use base and recommended R functionality. The normal-model implementations additionally require `numDeriv`.

```r
install.packages("numDeriv")
```

Dependency details:

- `fit_nlm_mlqe.R` requires `numDeriv`, including when analytical model derivatives are supplied, because the consistency-transformation Jacobians are computed numerically.
- `fit_nlm_mdpde.R` requires `numDeriv` only when the model does not supply `D` or `V`.
- `fit_nlm_t.R` has no non-base package dependency; missing derivatives are computed internally by finite differences.

## Directory structure

```text
project/
├── README.md
├── FUNCTION_REFERENCE.md
├── fit_nlm_mlqe.R
├── fit_nlm_mdpde.R
├── fit_nlm_t.R
├── example1_nlm_commented.R
├── example2_nlm_commented.R
├── example_nlm_t_commented.R
└── application_radioimmunoassay.R
```

## Quick start

Set the working directory to the folder containing the scripts and run one example:

```r
setwd("path/to/project")

source("example1_nlm_commented.R")
```

Because each example script sources the required fitting functions, it is not necessary to source them separately when the file names above are used.

To run the Student-t example:

```r
source("example_nlm_t_commented.R")
```

To run the real-data application:

```r
source("application_radioimmunoassay.R")
```

## Common nonlinear model interface

All three fitting functions require a user-supplied model with the interface:

```r
model <- function(theta, data) {
  # Compute the conditional mean and dispersion.
  mu <- ...
  sigma2 <- ...

  # Optional analytical derivatives.
  D <- ...
  V <- ...

  list(
    mu = mu,
    sigma2 = sigma2,
    D = D,
    V = V
  )
}
```

The required components are:

- `mu`: conditional means, either a scalar or a vector of length `length(y)`;
- `sigma2`: a positive scalar or vector of length `length(y)`.

The model may return `sigma` instead of `sigma2`. In that case, the fitting function uses `sigma^2`.

The optional derivative matrices are:

- `D = d mu / d theta^T`;
- `V = d sigma2 / d theta^T`.

Both must have `length(y)` rows and `length(theta)` columns. If one matrix is omitted, only that matrix is computed numerically.

Parameter names and ordering must agree across:

- the starting-value vector;
- the parameter extraction inside `model()`;
- the columns of `D` and `V`;
- the consistency transformation `tau()` used by corrected MLqE.

## Dispersion convention

The interpretation of `sigma2` depends on the distribution.

### Normal MLqE and MDPDE

For `fit_nlm_mlqe()` and `fit_nlm_mdpde()`:

```text
sigma2_i = Var(Y_i | covariates).
```

### Student-t MLE

For `fit_nlm_t()`:

```text
sigma2_i = squared Student-t scale parameter.
```

It is not the conditional variance. When `nu > 2`:

```text
Var(Y_i | covariates) = nu * sigma2_i / (nu - 2).
```

For `nu = 4`, the conditional variance is `2 * sigma2_i`.

## Estimation procedures

### Corrected MLqE

The original MLqE estimate is denoted by `theta_star`. The user supplies a transformation family:

```r
tau <- function(theta, r) {
  ...
}
```

The corrected estimate is:

```text
theta_hat = tau(theta_star, q).
```

The inverse mapping is evaluated with `r = 1 / q`. For the log-linear variance models used in the examples:

```r
tau <- function(theta, r) {
  result <- theta
  result["gamma0"] <- result["gamma0"] - log(r)
  result
}
```

### MDPDE

The MDPDE is fitted directly under the normal model. No consistency transformation is required.

### Student-t MLE

The Student-t degrees of freedom are treated as known and fixed. The supplied value of `nu` is used in estimation and is not estimated by `fit_nlm_t()`.

### Normal MLE

For the normal-model implementations:

```r
q = 1
```

reduces the estimator to ordinary normal maximum likelihood. Only one of the two functions is required for this reference fit.

## Automatic selection of q by SQV

For corrected MLqE and MDPDE, use:

```r
q = "SQV"
```

The current implementation uses an initial grid from `1` to `0.8`. The active controls are:

```r
q_control = list(
  m0 = 21L,
  m = 3L,
  q_min = 0.5,
  L = 0.01,
  verbose = FALSE,
  keep_fits = FALSE
)
```

`compute_vcov = TRUE` is required because SQV uses standardized estimates.

The initial lower endpoint is fixed at `0.8` in the current implementation. An `initial_lower` entry supplied in `q_control` is ignored because unrecognized controls are removed.

## Script descriptions

### `example1_nlm_commented.R`

This script simulates a heteroscedastic nonlinear normal model with:

```text
mu_i = beta1 + beta2 / (1 + beta3 * X1_i^beta4)

sigma2_i = exp(gamma0 + gamma1 * w1_i + gamma2 * w2_i).
```

It then replaces one response with an atypical value and fits:

- corrected MLqE with `q` selected by SQV;
- MDPDE with `q` selected by SQV;
- normal MLE using `q = 1`.

The script demonstrates:

- analytical derivatives `D` and `V`;
- construction of data-driven starting values;
- the MLqE transformation `tau`;
- fit summaries and confidence intervals;
- selected `q` values and SQV histories;
- corrected and original MLqE parameterizations;
- sandwich covariance matrices;
- fitted means, variances, residuals, and predictions;
- density-power score weights;
- plots of fitted mean curves and weights.

### `example2_nlm_commented.R`

This script uses a homoscedastic nonlinear normal model:

```text
mu_i = beta1 + beta2 * exp(-beta3 * x_i)

sigma2_i = exp(gamma0).
```

One response is contaminated. The script fits corrected MLqE, MDPDE, and normal MLE and reports the same main diagnostic and inferential quantities as Example 1.

This example shows that a constant variance is handled by returning either a scalar `sigma2` or a replicated vector.

### `example_nlm_t_commented.R`

This script simulates an uncontaminated heteroscedastic nonlinear Student-t sample with `nu = 4`. Its mean and squared-scale structures match those used in Example 1.

It demonstrates:

- direct Student-t data generation;
- the distinction between squared scale and variance;
- Student-t maximum-likelihood fitting;
- expected-Fisher and observed-information covariance matrices;
- response, scale, and Pearson residuals;
- Student-t likelihood weights;
- MAE calculations;
- log-likelihood, AIC, and BIC;
- prediction of means, scales, and squared scales;
- plots of fitted mean curves and likelihood weights.

### `application_radioimmunoassay.R`

This script analyzes the radioimmunoassay data using:

```text
mu_i = beta1 + beta2 / (1 + beta3 * X1_i^beta4)

sigma2_i = exp(gamma0 + gamma1 * X1_i).
```

The intended comparison contains:

- Student-t MLE with four degrees of freedom;
- normal MDPDE with `q` selected by SQV;
- normal corrected MLqE with `q` selected by SQV;
- normal MLE using all observations;
- normal MLE after excluding observation 9.

The script produces summaries, confidence intervals, fitted means, MAEs, fitted-curve plots, and density-power weight plots. The PDF outputs are:

```text
XvsYapplication.pdf
score_weights_application.pdf
```

## Accessing fitted results

Typical access patterns are:

```r
coef(fit)
vcov(fit)
confint(fit)
fitted(fit)
residuals(fit)
predict(fit, newdata = new_data, type = "mean")
summary(fit)
```

Corrected MLqE also permits access to the original parameterization:

```r
coef(fit_mlqe, type = "star")
vcov(fit_mlqe, type = "star")
confint(fit_mlqe, type = "star")
fitted(fit_mlqe, type = "star")
residuals(fit_mlqe, type = "star")
```

Detailed arguments, return components, and method-specific behavior are documented in `FUNCTION_REFERENCE.md`.

## Optimization guidance

For unconstrained problems, the examples use BFGS:

```r
method = "BFGS"
control = list(
  maxit = 2000,
  reltol = 1e-12
)
```

Use `L-BFGS-B` when finite bounds are required:

```r
method = "L-BFGS-B"
lower = ...
upper = ...
control = list(
  maxit = 2000,
  factr = 1e7,
  pgtol = 0
)
```

Starting values must satisfy the bounds. Analytical derivatives are strongly recommended for nonlinear models because they reduce numerical error and execution time.

Always inspect:

```r
fit$convergence
fit$message
fit$score
```

For Student-t fits, also inspect:

```r
fit$max.abs.score
```

A convergence code of zero is expected from `optim()`.

## Current consistency issues in `application_radioimmunoassay.R`

`application_radioimmunoassay.R` contains copied elements that should be corrected before distributing or running it as a final application:

1. The Student-t application fit is initially created with `nu = 4`, but `nu` is later overwritten with `9699` for an auxiliary normal-approximation fit. Subsequent labels that use the global `nu` can therefore display the wrong value.
2. The MAE block uses indices `78` and `80`, although the radioimmunoassay data contain only 14 observations. The intended excluded observation in this application is observation 9.
3. The MAE labels refer to `nu = 5` and to removing observations 78 and 80, which are inconsistent with the application description.
4. The auxiliary objects `fit_mle2` and `fit_mle3` are implementation checks rather than distinct methods required in the final comparison.

These points are documented rather than silently changed.

## Reproducibility

The simulated examples contain commented `set.seed()` calls. Uncomment the relevant line to reproduce the same sample on every run.

## Function documentation

See:

```text
FUNCTION_REFERENCE.md
```

for complete usage signatures, argument descriptions, output components, prediction conventions, S3 methods, and minimal examples.

---
output:
  pdf_document:
    latex_engine: xelatex
  html_document: default
---
# Linear Measurement-Error Model Estimation

This repository contains R implementations and reproducible examples for estimation in a linear measurement-error model using:

- corrected maximum Lq-likelihood estimation (corrected MLqE) under a bivariate normal model;
- maximum density power divergence estimation (MDPDE) under a bivariate normal model;
- maximum likelihood estimation under a bivariate Student-t model with fixed degrees of freedom;
- normal maximum likelihood estimation as the special case `q = 1` of the MDPDE or MLqE implementation.

The implementations allow observation-specific known measurement-error variances or Student-t scale-squared quantities.

## Model

The latent structural model is

```text
Y_i* = beta0 + beta1 X_i* + epsilon_i,
```

where, under the normal formulation,

```text
X_i*       ~ Normal(mu_x, sigma2_x),
epsilon_i  ~ Normal(0, sigma2).
```

The observed variables are

```text
Y_i = Y_i* + e_yi,
X_i = X_i* + e_xi,
```

with known observation-specific measurement-error quantities `tau_yi` and `tau_xi`.

The parameter vector is fixed as

```text
theta = (beta0, beta1, mu_x, sigma2_x, sigma2).
```

Under the normal model,

```text
tau_yi  = Var(e_yi),
tau_xi  = Var(e_xi),
sigma2_x = Var(X_i*),
sigma2   = Var(epsilon_i).
```

Consequently,

```text
(Y_i, X_i)^T ~ Normal_2(mu(theta), Sigma_i(theta)),
```

with

```text
mu(theta) =
  (beta0 + beta1 * mu_x, mu_x)^T
```

and

```text
Sigma_i(theta) =

  [ beta1^2 * sigma2_x + sigma2 + tau_yi,
    beta1 * sigma2_x                              ]

  [ beta1 * sigma2_x,
    sigma2_x + tau_xi                            ].
```

The Student-t implementation uses the same location and matrix structure, but every dispersion quantity is interpreted as a component of the Student-t scale matrix rather than as a covariance quantity.

## Repository contents

| File name | Description |
|---|---|
| `fit_lmve_mlqe.R` | Generic corrected MLqE implementation for the linear normal measurement-error model, including fixed-`q` fitting and automatic `q` selection by SQV. |
| `fit_lmve_mdpde.R` | Generic MDPDE implementation for the linear normal measurement-error model, including fixed-`q` fitting and automatic `q` selection by SQV. |
| `fit_lmve_t.R` | Generic maximum-likelihood implementation for the bivariate Student-t measurement-error model with fixed degrees of freedom. |
| `example_lmve_commented.R` | Simulated bivariate normal measurement-error example with one contaminated observed pair. It compares corrected MLqE, MDPDE, and normal MLE. |
| `example_lmve_t_commented.R` | Simulated uncontaminated bivariate Student-t measurement-error example fitted by Student-t MLE. |
| `application_sbp.R` | Application to the systolic blood pressure data from `MethComp`, comparing Student-t MLE, MDPDE, corrected MLqE, normal MLE, and normal MLE after removing observations 78 and 80. |
| `FUNCTION_REFERENCE.md` | Detailed documentation of the public fitting functions, their arguments, returned objects, prediction methods, and examples. |

The example and application scripts call:

```r
source("fit_lmve_mlqe.R")
source("fit_lmve_mdpde.R")
source("fit_lmve_t.R")
```

Keep the implementation files under these names or modify the `source()` calls accordingly.

## Requirements

Install the packages used by the examples and application:

```r
install.packages(
  c(
    "mvtnorm",
    "MethComp"
  )
)
```

The package `nloptr` is required only when `fit_lmve_t()` is called with `method = "BOBYQA"`:

```r
install.packages("nloptr")
```

Package use by file:

| File | Required package |
|---|---|
| `example_lmve_commented.R` | `mvtnorm` |
| `example_lmve_t_commented.R` | `mvtnorm` |
| `application_sbp.R` | `MethComp` |
| `fit_lmve_t.R` with `method = "BOBYQA"` | `nloptr` |

The MLqE and MDPDE implementations use analytical derivatives specialized to this five-parameter measurement-error model and do not require `numDeriv`.

## Directory structure

```text
project/
├── README.md
├── FUNCTION_REFERENCE.md
├── fit_lmve_mlqe.R
├── fit_lmve_mdpde.R
├── fit_lmve_t.R
├── example_lmve_commented.R
├── example_lmve_t_commented.R
└── application_sbp.R
```

## Quick start

Set the working directory to the repository folder and run one example:

```r
setwd("path/to/project")

source("example_lmve_commented.R")
```

To run the Student-t example:

```r
source("example_lmve_t_commented.R")
```

To run the systolic blood pressure application:

```r
source("application_sbp.R")
```

Each example or application script sources the required estimation functions.

## Data arguments

The three public fitting functions use the same four data arguments:

```r
Y
X
tau_y
tau_x
```

Their meanings are:

| Argument | Normal model | Student-t model |
|---|---|---|
| `Y` | Observed response values. | Observed response values. |
| `X` | Observed explanatory-variable values. | Observed explanatory-variable values. |
| `tau_y` | Known measurement-error variance of `Y`. | Known measurement-error scale-squared quantity of `Y`. |
| `tau_x` | Known measurement-error variance of `X`. | Known measurement-error scale-squared quantity of `X`. |

`Y` and `X` must have the same nonzero length. `tau_y` and `tau_x` may be scalars or vectors; scalar values are recycled to the sample size. They must be finite and non-negative.

## Starting values

All three functions accept:

```r
start = NULL
```

When `start = NULL`, robust moment-based starting values are calculated internally.

An explicit starting vector must contain exactly:

```r
start <- c(
  beta0 = ...,
  beta1 = ...,
  mu_x = ...,
  sigma2_x = ...,
  sigma2 = ...
)
```

Both `sigma2_x` and `sigma2` must be strictly positive.

The normal MLqE and MDPDE functions use medians, MAD-based scale estimates, median measurement-error variances, and a robust covariance expression. Initial variance values below `0.05` or non-finite initial variance values are replaced by `0.1`.

The Student-t implementation uses the same robust construction. When `nu > 2`, empirical second-moment quantities are approximately converted to the Student-t scale parameterization by multiplying them by

```text
(nu - 2) / nu.
```

## Corrected MLqE parameterization

The MLqE transformation is internal to `fit_lmve_mlqe()`; the user does not supply a `tau` function.

The corrected parameter vector is

```text
theta = (beta0, beta1, mu_x, sigma2_x, sigma2).
```

The original MLqE parameterization is

```text
theta_star =
  (beta0, beta1, mu_x, q * sigma2_x, q * sigma2).
```

The known measurement-error variances are also transformed internally:

```text
tau_y_star = q * tau_y,
tau_x_star = q * tau_x.
```

The corrected estimate is returned in:

```r
fit_mlqe$coefficients
```

and the original MLqE estimate is returned in:

```r
fit_mlqe$coefficients.star
```

## MDPDE

The MDPDE is fitted directly under the bivariate normal model. It does not require a consistency transformation.

For `q = 1`, the MDPDE reduces to ordinary maximum likelihood under the normal measurement-error model:

```r
fit_mle <- fit_lmve_mdpde(
  Y = Y,
  X = X,
  tau_y = tau_y,
  tau_x = tau_x,
  start = NULL,
  q = 1
)
```

## Student-t MLE

The Student-t degrees of freedom are fixed by the user and are not estimated:

```r
fit_t <- fit_lmve_t(
  Y = Y,
  X = X,
  tau_y = tau_y,
  tau_x = tau_x,
  start = NULL,
  nu = 5
)
```

The Student-t quantities

```text
sigma2_x, sigma2, tau_y, tau_x
```

are scale-squared quantities, not variances.

When `nu > 2`:

```text
Var[(Y_i, X_i)^T] =
  nu / (nu - 2) * Sigma_i(theta).
```

The fitted object returns both the Student-t scale matrices and the corresponding covariance matrices when the covariance exists.

## Automatic selection of q by SQV

For corrected MLqE or MDPDE, use:

```r
q = "SQV"
```

A typical configuration is:

```r
q_control <- list(
  m0 = 21L,
  m = 3L,
  q_min = 0.5,
  L = 0.01,
  verbose = FALSE,
  keep_fits = TRUE
)
```

The active defaults are:

```r
list(
  m0 = 21L,
  m = 3L,
  q_min = 0.5,
  L = 0.01,
  verbose = FALSE,
  keep_fits = FALSE
)
```

The initial grid is fixed from `q = 1` to `q = 0.8`.

`compute_vcov = TRUE` is required because the SQV procedure compares standardized estimates.

## Optimization

The normal-model functions support:

```text
"BFGS"
"L-BFGS-B"
"Nelder-Mead"
```

The Student-t function additionally supports:

```text
"CG"
"BOBYQA"
```

For positive variance or scale parameters, the examples use `L-BFGS-B`:

```r
lower_bounds <- c(
  beta0 = -Inf,
  beta1 = -Inf,
  mu_x = -Inf,
  sigma2_x = 1e-8,
  sigma2 = 1e-8
)

upper_bounds <- c(
  beta0 = Inf,
  beta1 = Inf,
  mu_x = Inf,
  sigma2_x = Inf,
  sigma2 = Inf
)

fit <- fit_lmve_mdpde(
  Y = Y,
  X = X,
  tau_y = tau_y,
  tau_x = tau_x,
  start = NULL,
  q = 1,
  method = "L-BFGS-B",
  lower = lower_bounds,
  upper = upper_bounds,
  control = list(
    maxit = 2000,
    factr = 1e4,
    pgtol = 1e-12
  )
)
```

Finite bounds are used only by `L-BFGS-B` in the normal MLqE and MDPDE implementations. In `fit_lmve_t()`, finite bounds are supported by `L-BFGS-B` and `BOBYQA`.

## Script descriptions

### `example_lmve_commented.R`

This script simulates observation-specific bivariate normal measurement-error data and replaces observation 9 with the atypical pair:

```r
c(
  Y = 35,
  X = -20
)
```

It fits:

- corrected MLqE with `q` selected by SQV;
- MDPDE with `q` selected by SQV;
- normal MLE through `fit_lmve_mdpde(q = 1)`.

The script demonstrates:

- generation with observation-specific measurement-error variances;
- internal robust starting values;
- parameter estimates and standard errors;
- corrected and original MLqE parameterizations;
- confidence intervals and sandwich covariance matrices;
- analytical mean and covariance derivative matrices;
- fitted marginal distributions and residuals;
- structural and marginal predictions;
- density-power score weights;
- structural regression-line and score-weight plots.

### `example_lmve_t_commented.R`

This script generates an uncontaminated bivariate Student-t sample with observation-specific scale matrices and fits it by Student-t maximum likelihood.

It demonstrates:

- Student-t data generation using `mvtnorm::rmvt()`;
- internal robust starting values;
- expected-Fisher and observed-information covariance matrices;
- fitted locations, scale matrices, covariance matrices, and residuals;
- conditional Student-t locations, conditional scale-squared values, and conditional degrees of freedom;
- MAEs relative to the observed response and true conditional location;
- log-likelihood, AIC, and BIC;
- structural, marginal, and conditional-observed predictions;
- Student-t likelihood weights;
- structural-line and weight plots.

The current script sets:

```r
nu <- 5
```

### `application_sbp.R`

This script analyzes systolic blood pressure measurements from the `sbp` dataset in `MethComp`.

The two selected methods are:

```text
J: observer measurement
S: machine measurement
```

Replicate measurements are reduced to item-specific means. The measurement-error variance of each item mean is estimated as:

```text
sample variance of replicates / number of replicates.
```

Thus, `tau_x` and `tau_y` vary across items.

The application compares:

- Student-t MLE;
- normal MDPDE;
- normal corrected MLqE;
- normal MLE using all observations;
- normal MLE after removing observations 78 and 80.

It computes fitted conditional locations and MAEs using all observations and after excluding observations 78 and 80. It also exports:

```text
lmve_sbp_fitted_lines.pdf
lmve_sbp_score_weights.pdf
```

The current script sets:

```r
nu <- 5
```

## Predictions

### Structural regression line

For a latent explanatory value `x*`, the structural mean is

```text
E(Y_i* | X_i* = x*) = beta0 + beta1 x*.
```

Normal MLqE example:

```r
predict(
  fit_mlqe,
  newdata = list(
    latent_x = seq(-5, 5, length.out = 100)
  ),
  type = "structural_mean"
)
```

Normal MDPDE example:

```r
predict(
  fit_mdpde,
  newdata = list(
    x = seq(-5, 5, length.out = 100)
  ),
  type = "structural_mean"
)
```

Student-t example:

```r
predict(
  fit_t,
  newx = seq(-5, 5, length.out = 100),
  type = "structural"
)
```

### Marginal observed means

Normal MLqE and MDPDE:

```r
predict(
  fit_mdpde,
  type = "mean"
)
```

Student-t:

```r
predict(
  fit_t,
  type = "marginal"
)
```

### Conditional observed location of Y given X

The Student-t function provides:

```r
predict(
  fit_t,
  newx = X_new,
  tau_x = tau_x_new,
  type = "conditional_observed"
)
```

For normal MLqE and MDPDE, the conditional mean can be calculated from the fitted coefficients:

```r
conditional_mean_y_given_x <- function(
    fit,
    X_new,
    tau_x_new
) {
  theta <- fit$coefficients

  beta0 <- theta["beta0"]
  beta1 <- theta["beta1"]
  mu_x <- theta["mu_x"]
  sigma2_x <- theta["sigma2_x"]

  beta0 +
    beta1 * mu_x +
    beta1 * sigma2_x /
      (
        sigma2_x +
          tau_x_new
      ) *
      (
        X_new -
          mu_x
      )
}
```

## Score weights

For normal MLqE:

```r
weight_mlqe <- if (
  fit_mlqe$q == 1
) {
  rep(1, length(Y))
} else {
  fit_mlqe$density.star^(
    1 -
      fit_mlqe$q
  )
}
```

For normal MDPDE:

```r
weight_mdpde <- if (
  fit_mdpde$q == 1
) {
  rep(1, length(Y))
} else {
  fit_mdpde$density^(
    1 -
      fit_mdpde$q
  )
}
```

For Student-t MLE:

```r
fit_t$weights
```

The Student-t weight for observation `i` is based on its bivariate Mahalanobis distance:

```text
weight_i = (nu + 2) / (nu + delta_i),
```

where `delta_i` is the fitted quadratic form for the two-dimensional observed pair.

## Accessing fitted results

Typical commands are:

```r
coef(fit)
vcov(fit)
confint(fit)
fitted(fit)
residuals(fit)
predict(fit, ...)
summary(fit)
```

For corrected MLqE, the original parameterization is also available:

```r
coef(
  fit_mlqe,
  type = "star"
)

vcov(
  fit_mlqe,
  type = "star"
)

confint(
  fit_mlqe,
  type = "star"
)

fitted(
  fit_mlqe,
  type = "star"
)

residuals(
  fit_mlqe,
  type = "star"
)
```

For Student-t MLE:

```r
vcov(
  fit_t,
  type = "fisher"
)

vcov(
  fit_t,
  type = "observed"
)

logLik(fit_t)
nobs(fit_t)
```

Detailed argument and output documentation is provided in `FUNCTION_REFERENCE.md`.

## Recommended checks after fitting

For every method:

```r
fit$convergence
fit$message
fit$counts
```

For normal MLqE and MDPDE:

```r
fit$score
fit$vcov.warning
```

For Student-t MLE:

```r
fit$max.abs.score
fit$fisher.information
fit$observed.information
fit$vcov.warning
fit$observed.information.warning
```

For SQV fits:

```r
fit$q
fit$q.selection$reason
fit$q.selection$history
fit$q.selection$evaluations
```

A convergence code equal to zero indicates normal convergence in the returned fit objects.

## Current degrees-of-freedom setting

Both `example_lmve_t_commented.R` and `application_sbp.R` currently set:

```r
nu <- 5
```

The fitting function itself accepts any fixed positive `nu`. To standardize the public repository to four degrees of freedom, change the two script assignments and the corresponding comments and labels to:

```r
nu <- 4
```

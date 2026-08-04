---
output:
  pdf_document:
    latex_engine: xelatex
  html_document: default
---
# Linear Mixed Model Estimation

This repository contains generic R implementations and reproducible examples for estimation in linear mixed models using:

- corrected maximum Lq-likelihood estimation (corrected MLqE) under a Gaussian marginal model;
- maximum density power divergence estimation (MDPDE) under a Gaussian marginal model;
- maximum likelihood estimation under a multivariate Student-t marginal model with fixed degrees of freedom;
- normal maximum likelihood estimation as the special case `q = 1` of the corrected MLqE or MDPDE implementations.

The implementations support grouped and unbalanced data, arbitrary fixed- and random-effects design matrices, and user-defined differentiable parameterizations of the random-effects covariance or scale matrix.

## Model

For independent groups \(i = 1,\ldots,n\), the Gaussian linear mixed model is

```text
Y_i = X_i beta + Z_i b_i + epsilon_i,
```

where

```text
b_i       ~ Normal(0, Delta(gamma)),
epsilon_i ~ Normal(0, sigma2 I_{d_i}).
```

Consequently,

```text
Y_i ~ Normal(
  X_i beta,
  Z_i Delta(gamma) Z_i^T + sigma2 I_{d_i}
).
```

The complete parameter vector is ordered as

```text
theta = (beta, gamma, sigma2).
```

The Student-t implementation uses the marginal model

```text
Y_i ~ t_{d_i}(
  X_i beta,
  Sigma_i,
  nu
),
```

where

```text
Sigma_i =
  Z_i Delta(gamma) Z_i^T +
  sigma2 I_{d_i}.
```

In the Student-t model, `Delta(gamma)`, `sigma2`, and `Sigma_i` are scale quantities, not covariance quantities. When `nu > 2`,

```text
Var(Y_i) =
  nu / (nu - 2) * Sigma_i.
```

A hierarchical representation that produces this marginal Student-t model uses one common scale-mixture variable per group:

```text
lambda_i ~ Gamma(nu / 2, rate = nu / 2),

b_i | lambda_i
  ~ Normal(
      0,
      Delta(gamma) / lambda_i
    ),

epsilon_i | lambda_i
  ~ Normal(
      0,
      sigma2 I_{d_i} / lambda_i
    ).
```

This is not equivalent to assuming independent Gaussian random effects and Student-t residual errors.

## Repository contents

| File name | Description |
|---|---|
| `fit_lmm_mlqe.R` | Generic corrected MLqE implementation for Gaussian linear mixed models, including fixed-`q` fitting and automatic `q` selection by SQV. |
| `fit_lmm_mdpde.R` | Generic MDPDE implementation for Gaussian linear mixed models, including fixed-`q` fitting and automatic `q` selection by SQV. |
| `fit_lmm_t.R` | Generic maximum-likelihood implementation for linear mixed models with a multivariate Student-t marginal distribution and fixed degrees of freedom. |
| `example_lmm_commented.R` | Simulated Gaussian mixed-model example with one contaminated group. It compares corrected MLqE, MDPDE, and normal MLE. |
| `example_lmm_t_commented.R` | Simulated uncontaminated Student-t mixed-model example with four degrees of freedom, fitted by Student-t MLE. |
| `application_orthodont.R` | Application to the orthodontic growth data, comparing Student-t MLE, MDPDE, corrected MLqE, normal MLE, and normal MLE after removing subjects M09 and M13. |
| `FUNCTION_REFERENCE.md` | Detailed documentation of the public fitting functions, their arguments, returned objects, prediction methods, and examples. |

The example and application scripts call:

```r
source("fit_lmm_mlqe.R")
source("fit_lmm_mdpde.R")
source("fit_lmm_t.R")
```

Keep the implementation files under these names or modify the `source()` calls accordingly.

## Requirements

Install the packages required by the implementations and examples:

```r
install.packages(
  c(
    "MASS",
    "numDeriv",
    "nloptr",
    "nlme",
    "lme4"
  )
)
```

Package use:

| Package | Use |
|---|---|
| `MASS` | Required when `start = NULL`, because the automatic fixed-effect and residual-scale starting values are based on `MASS::rlm()`. |
| `numDeriv` | Required by corrected MLqE for transformation Jacobians. It is also required by any method when `Delta_jacobian` is omitted. |
| `nloptr` | Required when `method = "BOBYQA"`. |
| `nlme` | Provides the `Orthodont` data used in the application. |
| `lme4` | Provides reference Gaussian maximum-likelihood fits in the application. |

The Student-t simulated example currently checks for both `MASS` and `nloptr`. Its fitted model uses BFGS, so `MASS` is required by `start = NULL`; `nloptr` is required only if the example is changed to BOBYQA.

## Directory structure

```text
project/
├── README.md
├── FUNCTION_REFERENCE.md
├── fit_lmm_mlqe.R
├── fit_lmm_mdpde.R
├── fit_lmm_t.R
├── example_lmm_commented.R
├── example_lmm_t_commented.R
└── application_orthodont.R
```

## Quick start

Set the working directory to the repository folder and run one example:

```r
setwd("path/to/project")

source("example_lmm_commented.R")
```

To run the Student-t example:

```r
source("example_lmm_t_commented.R")
```

To run the orthodontic growth application:

```r
source("application_orthodont.R")
```

Each example or application script sources the required estimation functions.

## Grouped-data interface

The public fitting functions receive:

```r
Y
X
Z
```

as group-specific objects.

The recommended representation is:

```r
Y <- list(
  group1 = numeric_vector_1,
  group2 = numeric_vector_2
)

X <- list(
  group1 = fixed_effect_matrix_1,
  group2 = fixed_effect_matrix_2
)

Z <- list(
  group1 = random_effect_matrix_1,
  group2 = random_effect_matrix_2
)
```

For each group \(i\):

```text
length(Y[[i]]) =
nrow(X[[i]]) =
nrow(Z[[i]]).
```

All `X` matrices must have the same number of columns, and all `Z` matrices must have the same number of columns.

A vector or matrix representing a single group may be supplied directly; it is converted internally to a one-element list.

The group names are preserved when named lists are supplied. This is useful when associating fitted random effects and weights with subjects.

## Random-effects covariance structure

The user supplies a function:

```r
Delta <- function(gamma) {
  ...
}
```

It must return an \(m \times m\) symmetric positive-definite matrix, where \(m\) is the number of columns of each `Z[[i]]`.

Supported structures include:

- diagonal covariance matrices;
- unstructured covariance matrices;
- compound-symmetry structures;
- covariance matrices defined through Cholesky factors;
- other differentiable positive-definite parameterizations.

### Unstructured example

```r
Delta_unstructured <- function(gamma) {
  matrix(
    c(
      gamma["gamma1"],
      gamma["gamma2"],
      gamma["gamma2"],
      gamma["gamma3"]
    ),
    nrow = 2,
    byrow = TRUE
  )
}
```

This direct parameterization does not automatically enforce positive definiteness. The optimizer may evaluate invalid values, which are rejected internally. A Cholesky-based parameterization is preferable when positive definiteness must be guaranteed over the complete parameter space.

## `gamma_structure`

When `start = NULL`, the function must know the names and locations of the parameters in `Delta(gamma)`. Supply a symmetric character matrix:

```r
gamma_structure_unstructured <- matrix(
  c(
    "gamma1",
    "gamma2",
    "gamma2",
    "gamma3"
  ),
  nrow = 2,
  byrow = TRUE
)
```

Entries may be gamma parameter names or `NA`.

The automatic starting-value rule assigns:

- `1` to gamma parameters appearing on the main diagonal;
- `0` to gamma parameters appearing only outside the main diagonal.

This rule is appropriate when gamma directly parameterizes entries of `Delta(gamma)`. For log-scale, Cholesky, or other transformed parameterizations, explicit starting values should generally be supplied.

## Analytical derivative of `Delta(gamma)`

The optional function

```r
Delta_jacobian <- function(gamma) {
  ...
}
```

must return

```text
d vec(Delta(gamma)) / d gamma^T
```

as either:

- an \(m^2 \times p_\gamma\) matrix; or
- an \(m \times m \times p_\gamma\) array.

R column-major vectorization is used.

For the unstructured two-dimensional example:

```r
Delta_unstructured_jacobian <- function(gamma) {
  matrix(
    c(
      1, 0, 0,
      0, 1, 0,
      0, 1, 0,
      0, 0, 1
    ),
    nrow = 4,
    ncol = 3,
    byrow = TRUE
  )
}
```

The rows correspond to:

```text
Delta_11
Delta_21
Delta_12
Delta_22
```

If `Delta_jacobian` is omitted, the derivative is calculated with `numDeriv::jacobian()`.

## Starting values

All three functions accept:

```r
start = NULL
```

When `start = NULL`:

- the fixed effects are initialized with `MASS::rlm()` fitted to the stacked response and fixed-effect design matrices;
- `sigma2` is initialized by the squared robust residual scale;
- gamma starting values are determined by `gamma_structure`.

An explicit starting vector must be ordered as:

```r
start <- c(
  beta,
  gamma,
  sigma2 = ...
)
```

For example:

```r
start <- c(
  beta0 = 18,
  beta1 = 10,
  gamma1 = 3,
  gamma2 = -2,
  gamma3 = 4,
  sigma2 = 5
)
```

The number of beta parameters must equal `ncol(X[[1]])`. The number of gamma parameters is inferred from the intermediate entries of `start`.

## Corrected MLqE transformation

The corrected MLqE requires a user-supplied transformation family:

```r
tau <- function(theta, r) {
  ...
}
```

For the direct covariance parameterization used in the examples, the fixed effects remain unchanged and all covariance parameters are divided by `r`:

```r
tau <- function(theta, r) {
  result <- theta

  covariance_parameters <- c(
    "gamma1",
    "gamma2",
    "gamma3",
    "sigma2"
  )

  result[covariance_parameters] <-
    result[covariance_parameters] /
    r

  result
}
```

If `theta_star` is the original MLqE estimate, then:

```text
corrected estimate =
  tau(theta_star, q),

original parameterization =
  tau(corrected estimate, 1 / q).
```

For `q = 1`, `tau` may be omitted because the identity transformation is used.

## MDPDE and normal MLE

The MDPDE is fitted directly under the Gaussian marginal model and does not require a consistency transformation.

For either normal estimator:

```r
q = 1
```

gives ordinary Gaussian maximum likelihood.

For example:

```r
fit_mle <- fit_lmm_mdpde(
  Y = Y,
  X = X,
  Z = Z,
  Delta = Delta_unstructured,
  start = NULL,
  gamma_structure =
    gamma_structure_unstructured,
  Delta_jacobian =
    Delta_unstructured_jacobian,
  q = 1
)
```

## Student-t MLE

The degrees of freedom are known and fixed:

```r
fit_t <- fit_lmm_t(
  Y = Y,
  X = X,
  Z = Z,
  Delta = Delta_unstructured,
  start = NULL,
  gamma_structure =
    gamma_structure_unstructured,
  Delta_jacobian =
    Delta_unstructured_jacobian,
  nu = 4
)
```

The implementation does not estimate `nu`.

The Student-t application and simulated example are standardized to:

```r
nu <- 4
```

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

The current implementation fits each `q` independently from the same original starting vector. When a fit at `q` fails, does not converge, or produces invalid standard errors, the SQV procedure attempts `q - 0.001`. If both attempts fail, the comparison is classified as unstable.

`compute_vcov = TRUE` is required because SQV uses standardized estimates.

## Optimization methods

The Gaussian functions support:

```text
"BFGS"
"L-BFGS-B"
"Nelder-Mead"
"BOBYQA"
```

The Student-t function additionally supports:

```text
"CG"
```

### BFGS

```r
method = "BFGS",
control = list(
  maxit = 2000,
  reltol = 1e-12
)
```

### L-BFGS-B

```r
method = "L-BFGS-B",
control = list(
  maxit = 2000,
  factr = 1e7,
  pgtol = 0
)
```

Finite bounds are supported with `L-BFGS-B`.

### BOBYQA

```r
method = "BOBYQA",
control = list(
  maxeval = 100000,
  xtol_rel = 0,
  ftol_rel = 0,
  xtol_abs = 1e-8,
  ftol_abs = 1e-8
)
```

BOBYQA uses `nloptr::bobyqa()` and does not use the analytical score.

Finite bounds are supported with `BOBYQA`.

The default lower bound constrains only `sigma2` to be positive. Positive definiteness of `Delta(gamma)` remains the responsibility of its parameterization or is checked during model evaluation.

## Script descriptions

### `example_lmm_commented.R`

This script simulates a Gaussian random-intercept and random-slope model with:

```text
beta   = (18, 10),
gamma  = (3, -2, 4),
sigma2 = 5,
```

and

```text
Delta(gamma) =

  [ gamma1  gamma2 ]
  [ gamma2  gamma3 ].
```

It uses 50 independent groups with five observations per group and contaminates group 30 through:

```r
Y[[30]] <- Y[[30]] * 3 + 10
```

It fits:

- corrected MLqE with `q` selected by SQV;
- MDPDE with `q` selected by SQV;
- Gaussian MLE through `fit_lmm_mdpde(q = 1)`.

The script demonstrates:

- grouped data generation;
- a direct unstructured `Delta(gamma)`;
- `gamma_structure`;
- an analytical `Delta_jacobian`;
- automatic robust starting values;
- corrected and original MLqE parameterizations;
- standard errors and confidence intervals;
- sensitivity and variability matrices;
- marginal and conditional fitted values;
- predicted random effects;
- group-level density-power weights;
- fitted fixed-effect lines and group-specific conditional trajectories.

### `example_lmm_t_commented.R`

This script uses the same fixed effects, covariance parameters, residual scale, and design structure as the Gaussian example, but generates each complete group from a multivariate Student-t distribution with:

```r
nu <- 4
```

It uses 5,000 groups with five observations per group and contains no artificial contamination.

The script demonstrates:

- direct multivariate Student-t generation;
- automatic robust starting values;
- Student-t MLE using BFGS;
- the estimated random-effects scale matrix;
- expected-Fisher and observed-information covariance matrices;
- marginal and conditional fitted values;
- predicted random effects;
- group-level Student-t weights;
- log-likelihood, AIC, and BIC;
- fitted fixed-effect lines and group-specific conditional trajectories.

### `application_orthodont.R`

This script analyzes the orthodontic growth data from `nlme::Orthodont`.

The fitted model is:

```text
distance_ij =
  beta0 +
  beta1 age_c_ij +
  beta2 SexFemale_i +
  beta3 age_c_ij SexFemale_i +
  b0_i +
  b1_i age_c_ij +
  epsilon_ij,
```

where age is centered at 11 years.

The application compares:

- Student-t MLE with four degrees of freedom;
- normal MDPDE;
- normal corrected MLqE;
- normal MLE using all subjects;
- normal MLE after removing subjects M09 and M13;
- reference Gaussian fits from `lme4::lmer()`.

It computes conditional fitted values and MAEs using:

- all observations;
- observations excluding subjects M09 and M13.

It also produces:

- an overall fixed-effect trajectory;
- a male fixed-effect trajectory;
- a female fixed-effect trajectory;
- a 7 by 4 grid of subject-specific conditional trajectories;
- subject-level robustness weights.

The PDF outputs are:

```text
lmm_orthodont_fitted_lines_overall.pdf
lmm_orthodont_fitted_lines_male.pdf
lmm_orthodont_fitted_lines_female.pdf
lmm_orthodont_subject_specific_conditional_trajectories.pdf
lmm_orthodont_score_weights.pdf
```

## Marginal and conditional quantities

For a fitted model:

```r
fitted(
  fit,
  type = "marginal"
)
```

returns the group-specific marginal fitted means:

```text
X_i beta.
```

```r
fitted(
  fit,
  type = "conditional"
)
```

returns:

```text
X_i beta + Z_i b_hat_i.
```

Marginal and conditional residuals are obtained similarly:

```r
residuals(
  fit,
  type = "marginal"
)

residuals(
  fit,
  type = "conditional"
)
```

The predicted random effects are stored in:

```r
fit$random.effects
```

## Predictions for new groups

For marginal means and covariance or scale matrices, `newdata` must contain `X` and `Z`:

```r
new_data <- list(
  X = list(X_new),
  Z = list(Z_new)
)
```

For random-effect or conditional-mean predictions, `newdata$Y` is also required:

```r
new_data <- list(
  Y = list(Y_new),
  X = list(X_new),
  Z = list(Z_new)
)
```

Examples:

```r
predict(
  fit_mdpde,
  newdata = new_data,
  type = "mean"
)

predict(
  fit_mdpde,
  newdata = new_data,
  type = "covariance"
)

predict(
  fit_mdpde,
  newdata = new_data,
  type = "random_effects"
)

predict(
  fit_mdpde,
  newdata = new_data,
  type = "conditional_mean"
)
```

For Student-t fits, use `"scale"` or `"variance"` instead of `"covariance"`:

```r
predict(
  fit_t,
  newdata = new_data,
  type = "scale"
)

predict(
  fit_t,
  newdata = new_data,
  type = "variance"
)
```

`type = "variance"` is available only when `nu > 2`.

## Group-level weights

For corrected MLqE:

```r
weight_mlqe <- if (
  fit_mlqe$q == 1
) {
  rep(
    1,
    fit_mlqe$n.groups
  )
} else {
  fit_mlqe$density.star^(
    1 -
      fit_mlqe$q
  )
}
```

For MDPDE:

```r
weight_mdpde <- if (
  fit_mdpde$q == 1
) {
  rep(
    1,
    fit_mdpde$n.groups
  )
} else {
  fit_mdpde$density^(
    1 -
      fit_mdpde$q
  )
}
```

These are group-level weights because each marginal density is defined for the complete response vector \(Y_i\).

For Student-t MLE:

```r
fit_t$weights
```

The weight of group \(i\), with dimension \(d_i\), is:

```text
weight_i =
  (nu + d_i) /
  (nu + delta_i),
```

where `delta_i` is the fitted group-level quadratic form.

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
  parameterization = "star"
)

residuals(
  fit_mlqe,
  parameterization = "star"
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

For corrected MLqE:

```r
fit_mlqe$q
fit_mlqe$score.star
fit_mlqe$score.corrected
fit_mlqe$tau.composition.residual
fit_mlqe$vcov.warning
```

For MDPDE:

```r
fit_mdpde$q
fit_mdpde$score
fit_mdpde$vcov.warning
```

For Student-t MLE:

```r
fit_t$max.abs.score
fit_t$fisher.information
fit_t$observed.information
fit_t$vcov.warning
fit_t$observed.information.warning
```

For SQV fits:

```r
fit$q
fit$q.selection$reason
fit$q.selection$history
fit$q.selection$evaluations
```

A convergence code equal to zero indicates normal convergence.

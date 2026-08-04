---
output:
  pdf_document: default
  html_document: default
---
# Function Reference: Linear Mixed Models

This document describes the public fitting functions:

- `fit_lmm_mlqe()`;
- `fit_lmm_mdpde()`;
- `fit_lmm_t()`.

## 1. Common grouped-data structure

The functions estimate independent group-level marginal distributions.

The inputs are:

```r
Y
X
Z
```

where:

- `Y[[i]]` is the response vector for group \(i\);
- `X[[i]]` is its fixed-effects design matrix;
- `Z[[i]]` is its random-effects design matrix.

For every group:

```text
length(Y[[i]]) =
nrow(X[[i]]) =
nrow(Z[[i]]).
```

All `X[[i]]` matrices must have the same number of columns, and all `Z[[i]]` matrices must have the same number of columns.

A vector, matrix, or data frame representing a single group is automatically converted to a one-element list.

All values must be finite and non-missing.

## 2. Common covariance or scale structure

### `Delta`

`Delta` must be a function of gamma:

```r
Delta <- function(gamma) {
  ...
}
```

For the Gaussian functions, it returns the random-effects covariance matrix.

For `fit_lmm_t()`, it returns the random-effects Student-t scale matrix.

It must have dimension:

```text
ncol(Z[[1]]) by ncol(Z[[1]]).
```

The matrix must be symmetric and positive definite at every valid parameter value evaluated by the fitting procedure.

### `gamma_structure`

`gamma_structure` is required when `start = NULL`.

It is a symmetric character matrix with the same dimensions as `Delta(gamma)`. Each non-missing entry identifies the gamma parameter occupying that position.

Example:

```r
gamma_structure <- matrix(
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

The automatic starting values are:

- `1` for gamma parameters appearing on the diagonal;
- `0` for gamma parameters appearing only off the diagonal.

For transformed parameterizations such as logarithmic or Cholesky structures, explicit starting values are generally preferable.

### `Delta_jacobian`

`Delta_jacobian` is optional:

```r
Delta_jacobian <- function(gamma) {
  ...
}
```

It returns:

```text
d vec(Delta(gamma)) / d gamma^T
```

as either:

- an \(m^2 \times p_\gamma\) matrix; or
- an \(m \times m \times p_\gamma\) array.

R column-major vectorization is used.

If it is omitted, the derivative is computed numerically with `numDeriv::jacobian()` using `jacobian_method` and `jacobian_method_args`.

## 3. Common parameter ordering

The starting vector is ordered as:

```text
theta = (beta, gamma, sigma2).
```

The first `ncol(X[[1]])` entries are fixed effects, the last entry is `sigma2`, and all intermediate entries parameterize `Delta(gamma)`.

Example:

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

Parameter names are strongly recommended.

When `start = NULL`, robust starting values are computed with `MASS::rlm()` and `gamma_structure`.

# 4. `fit_lmm_mlqe()`

## Purpose

Fits a Gaussian linear mixed model by maximum Lq-likelihood estimation and returns:

- the original MLqE parameterization;
- the Fisher-consistency-corrected parameterization.

For `q = 1`, the procedure reduces to Gaussian maximum likelihood.

## Usage

```r
fit_lmm_mlqe(
  Y,
  X,
  Z,
  Delta,
  tau = NULL,
  start = NULL,
  gamma_structure = NULL,
  Delta_jacobian = NULL,
  q = 1,
  level = 0.95,
  method = c(
    "BFGS",
    "L-BFGS-B",
    "Nelder-Mead",
    "BOBYQA"
  ),
  lower = NULL,
  upper = NULL,
  control = list(),
  use_score = FALSE,
  compute_vcov = TRUE,
  jacobian_method = c(
    "Richardson",
    "simple"
  ),
  jacobian_method_args = list(),
  q_control = list()
)
```

## Arguments

| Argument | Description |
|---|---|
| `Y` | Grouped response vectors. |
| `X` | Grouped fixed-effects design matrices. |
| `Z` | Grouped random-effects design matrices. |
| `Delta` | Function returning the random-effects covariance matrix. |
| `tau` | Transformation family `tau(theta, r)`. Required when `q < 1` or `q = "SQV"`. |
| `start` | Optional vector ordered as `(beta, gamma, sigma2)`. When `NULL`, robust starting values are calculated internally. |
| `gamma_structure` | Symmetric character matrix identifying the gamma parameters in `Delta(gamma)`. Required when `start = NULL`. |
| `Delta_jacobian` | Optional analytical derivative of `vec(Delta(gamma))`. |
| `q` | Numeric value in `(0, 1]`, or `"SQV"` for automatic selection. |
| `level` | Confidence level in `(0, 1)`. |
| `method` | `"BFGS"`, `"L-BFGS-B"`, `"Nelder-Mead"`, or `"BOBYQA"`. |
| `lower`, `upper` | Bounds supported by L-BFGS-B and BOBYQA. For MLqE they refer to the original optimization parameterization. |
| `control` | Optimizer control list. |
| `use_score` | Uses the analytical estimating equation with BFGS or L-BFGS-B. BOBYQA forces score use off; Nelder-Mead does not use a gradient. |
| `compute_vcov` | Computes `J`, `K`, sandwich covariance matrices, standard errors, and confidence intervals. Required for SQV. |
| `jacobian_method` | Numerical Jacobian method used by `numDeriv`: `"Richardson"` or `"simple"`. |
| `jacobian_method_args` | Additional `method.args` passed to `numDeriv::jacobian()`. |
| `q_control` | Controls for SQV selection. |

## Transformation family

For a standard covariance scaling correction:

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

The implementation uses:

```text
corrected estimate =
  tau(original estimate, q),

original estimate =
  tau(corrected estimate, 1 / q).
```

The Jacobian of the forward transformation is evaluated numerically, so `numDeriv` is required even when `Delta_jacobian` is supplied.

For `q = 1`, `tau` may be `NULL`; the identity transformation is used.

## Optimization defaults

For BFGS and Nelder-Mead:

```r
list(
  maxit = 1000,
  reltol = 1e-9
)
```

For L-BFGS-B:

```r
list(
  maxit = 1000,
  factr = 1e7,
  pgtol = 0
)
```

For BOBYQA:

```r
list(
  maxeval = 10000,
  xtol_rel = 1e-8,
  ftol_rel = 1e-8
)
```

For bounded methods, the default lower bound is `-Inf` for every parameter except `sigma2`, whose lower bound is `1e-10`. The default upper bounds are infinite.

The starting values must lie within the bounds.

## SQV controls

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

| Entry | Description |
|---|---|
| `m0` | Number of values in the initial grid from `1` to `0.8`. |
| `m` | Number of values in each refinement grid. It cannot exceed `m0`. |
| `q_min` | Smallest permitted `q`; it must be positive and less than `0.8`. |
| `L` | Positive stability threshold. |
| `verbose` | Prints SQV progress when `TRUE`. |
| `keep_fits` | Stores valid evaluated fits in `fit$q.selection$fits`. |

Unrecognized entries in `q_control` are discarded.

## Returned object

The returned object has class:

```text
"lmm_mlqe"
```

An SQV-selected fit additionally has class:

```text
"lmm_mlqe_sqv"
```

### Estimates and model parameters

| Component | Description |
|---|---|
| `coefficients` | Corrected parameter estimate. |
| `coefficients.star` | Original MLqE estimate. |
| `beta`, `beta.star` | Corrected and original fixed-effect estimates. |
| `gamma`, `gamma.star` | Corrected and original covariance-parameter estimates. |
| `sigma2`, `sigma2.star` | Corrected and original residual variances. |
| `Delta`, `Delta.star` | Corrected and original random-effects covariance matrices. |

### Inference

| Component | Description |
|---|---|
| `standard.error` | Corrected sandwich standard errors. |
| `standard.error.star` | Original-parameterization sandwich standard errors. |
| `vcov` | Corrected sandwich covariance matrix. |
| `vcov.star` | Sandwich covariance matrix for the original estimate. |
| `vcov.direct` | Corrected covariance calculated from the direct correction Jacobian. |
| `vcov.direct.difference` | Numerical difference between the two corrected-covariance calculations. |
| `conf.int` | Wald confidence intervals for corrected parameters. |
| `conf.int.star` | Wald confidence intervals for original parameters. |
| `J`, `J.star` | Sensitivity matrix in the original parameterization. |
| `K`, `K.star` | Variability matrix in the original parameterization. |

### Derivatives and transformation diagnostics

| Component | Description |
|---|---|
| `Delta.jacobian.star` | Derivative of `vec(Delta)` at the original estimate. |
| `mean.jacobians.star` | List of group-specific derivatives of the marginal mean. |
| `covariance.jacobians.star` | List of group-specific derivatives of `vec(Sigma_i)`. |
| `F.star` | List of combined mean and covariance derivative matrices. |
| `tau.jacobian` | Jacobian of `tau(theta, 1/q)` at the corrected estimate. |
| `tau.jacobian.inverse` | Inverse forward-transformation Jacobian. |
| `tau.q.jacobian` | Jacobian of `tau(theta_star, q)` at the original estimate. |
| `tau.composition.residual` | Numerical residual from checking the forward/inverse transformation composition. |
| `score`, `score.star` | Estimating equation at the original estimate. |
| `score.corrected` | Estimating equation evaluated at the corrected estimate. |

### Corrected fitted quantities

| Component | Description |
|---|---|
| `fitted.values`, `marginal.fitted.values` | List of group-specific marginal fitted means `X_i beta`. |
| `conditional.fitted.values` | List of conditional fitted means `X_i beta + Z_i b_hat_i`. |
| `random.effects` | List of predicted random-effect vectors. |
| `residuals`, `marginal.residuals` | List of marginal residual vectors. |
| `conditional.residuals` | List of conditional residual vectors. |
| `covariance.matrices` | List of fitted group covariance matrices. |
| `covariance.inverses` | List of fitted inverse covariance matrices. |
| `log.density`, `density` | Group-level Gaussian log-densities and densities. |

The original-parameterization versions use the `.star` suffix:

```text
fitted.values.star
marginal.fitted.values.star
conditional.fitted.values.star
random.effects.star
residuals.star
marginal.residuals.star
conditional.residuals.star
covariance.matrices.star
covariance.inverses.star
log.density.star
density.star
```

### Starting values, dimensions, and metadata

| Component | Description |
|---|---|
| `starting.values` | Corrected starting vector. |
| `starting.values.star` | Starting vector transformed to the original parameterization. |
| `automatic.start` | Whether starting values were generated internally. |
| `gamma.structure` | Gamma structure used for automatic initialization. |
| `gamma.diagonal.names` | Gamma parameters identified on the diagonal. |
| `robust.start.fit` | `MASS::rlm()` object used for automatic starting values. |
| `q`, `level` | Fitted tuning parameter and confidence level. |
| `n.groups`, `nobs`, `group.sizes` | Number of groups, total observations, and group sizes. |
| `p.beta`, `p.gamma`, `m` | Numbers of fixed effects, covariance parameters, and random effects. |
| `beta.names`, `gamma.names`, `sigma2.name` | Parameter names. |
| `objective` | Optimized Lq-likelihood objective. |
| `convergence`, `message`, `counts` | Optimizer diagnostics. |
| `method`, `jacobian.method` | Optimization and differentiation methods. |
| `Y`, `X`, `Z` | Stored grouped data. |
| `Delta.function`, `Delta.jacobian.function`, `tau` | Stored model functions. |
| `vcov.warning` | Covariance-inversion warning, or `NULL`. |
| `lower`, `upper` | Bounds used by bounded optimization. |
| `optimization.control` | Effective optimizer controls. |
| `optim` | Complete optimizer result. |
| `call` | Matched function call. |

SQV fits additionally contain:

```text
selected.q
sqv.returned.fit.source
q.selection
```

The `q.selection` object contains the selection reason, pass, stage, effective configuration, grid history, SQV history, evaluated fits, and optionally cached fit objects.

## S3 methods

```r
print(fit_mlqe)
summary(fit_mlqe)

coef(fit_mlqe)
coef(
  fit_mlqe,
  type = "star"
)

vcov(fit_mlqe)
vcov(
  fit_mlqe,
  type = "star"
)

confint(fit_mlqe)
confint(
  fit_mlqe,
  parm = c(
    "beta0",
    "gamma1"
  ),
  level = 0.90,
  type = "star"
)

fitted(
  fit_mlqe,
  type = "marginal",
  parameterization = "corrected"
)

fitted(
  fit_mlqe,
  type = "conditional",
  parameterization = "star"
)

residuals(
  fit_mlqe,
  type = "marginal",
  parameterization = "corrected"
)

residuals(
  fit_mlqe,
  type = "conditional",
  parameterization = "star"
)
```

## Prediction

```r
predict(
  fit_mlqe,
  newdata = NULL,
  type = c(
    "mean",
    "covariance",
    "random_effects",
    "conditional_mean"
  ),
  parameterization = c(
    "corrected",
    "star"
  )
)
```

Without `newdata`:

| `type` | Output |
|---|---|
| `"mean"` | Stored marginal fitted means. |
| `"covariance"` | Stored covariance matrices. |
| `"random_effects"` | Stored predicted random effects. |
| `"conditional_mean"` | Stored conditional fitted means. |

For new groups, `newdata` must contain `X` and `Z`:

```r
newdata <- list(
  X = list(X_new),
  Z = list(Z_new)
)
```

`newdata$Y` is additionally required for `"random_effects"` and `"conditional_mean"`.

## Minimal call

```r
fit_mlqe <- fit_lmm_mlqe(
  Y = Y,
  X = X,
  Z = Z,
  Delta = Delta_unstructured,
  tau = tau,
  start = NULL,
  gamma_structure =
    gamma_structure_unstructured,
  Delta_jacobian =
    Delta_unstructured_jacobian,
  q = "SQV",
  level = 0.95,
  method = "BFGS",
  control = list(
    maxit = 2000,
    reltol = 1e-12
  ),
  use_score = TRUE,
  compute_vcov = TRUE,
  jacobian_method = "Richardson",
  q_control = list(
    m0 = 21L,
    m = 3L,
    q_min = 0.5,
    L = 0.01,
    keep_fits = TRUE
  )
)
```

# 5. `fit_lmm_mdpde()`

## Purpose

Fits a Gaussian linear mixed model by maximum density power divergence estimation.

For `q = 1`, it reduces to Gaussian maximum likelihood.

## Usage

```r
fit_lmm_mdpde(
  Y,
  X,
  Z,
  Delta,
  start = NULL,
  gamma_structure = NULL,
  Delta_jacobian = NULL,
  q = 1,
  level = 0.95,
  method = c(
    "BFGS",
    "L-BFGS-B",
    "Nelder-Mead",
    "BOBYQA"
  ),
  lower = NULL,
  upper = NULL,
  control = list(),
  use_score = FALSE,
  compute_vcov = TRUE,
  jacobian_method = c(
    "Richardson",
    "simple"
  ),
  jacobian_method_args = list(),
  q_control = list()
)
```

## Arguments

The arguments have the same meanings and validation rules as in `fit_lmm_mlqe()`, except that no `tau` transformation is required.

`q` may be numeric in `(0, 1]` or equal to `"SQV"`.

## Optimization defaults

For BFGS and Nelder-Mead:

```r
list(
  maxit = 1000,
  reltol = 1e-9
)
```

For L-BFGS-B:

```r
list(
  maxit = 1000,
  factr = 1e7,
  pgtol = 0
)
```

For BOBYQA:

```r
list(
  maxeval = 10000,
  xtol_rel = 1e-8,
  ftol_rel = 1e-8
)
```

The SQV defaults and restrictions are the same as for corrected MLqE.

## Returned object

The returned object has class:

```text
"lmm_mdpde"
```

An SQV-selected fit additionally has class:

```text
"lmm_mdpde_sqv"
```

### Main components

| Component | Description |
|---|---|
| `coefficients` | MDPDE parameter estimate. |
| `beta` | Fixed-effect estimate. |
| `gamma` | Random-effects covariance-parameter estimate. |
| `sigma2` | Residual variance estimate. |
| `Delta` | Estimated random-effects covariance matrix. |
| `standard.error` | Sandwich standard errors. |
| `vcov` | Sandwich covariance matrix. |
| `conf.int` | Wald confidence intervals. |
| `J` | Sensitivity matrix. |
| `K` | Variability matrix. |
| `Delta.jacobian` | Derivative of `vec(Delta)` at the estimate. |
| `mean.jacobians` | List of group-specific marginal-mean derivative matrices. |
| `covariance.jacobians` | List of group-specific covariance derivative matrices. |
| `F` | List of combined derivative matrices. |
| `fitted.values`, `marginal.fitted.values` | List of marginal fitted means. |
| `conditional.fitted.values` | List of conditional fitted means. |
| `random.effects` | List of predicted random effects. |
| `residuals`, `marginal.residuals` | List of marginal residuals. |
| `conditional.residuals` | List of conditional residuals. |
| `covariance.matrices` | List of fitted group covariance matrices. |
| `covariance.inverses` | List of inverse covariance matrices. |
| `log.density`, `density` | Group-level fitted normal log-densities and densities. |
| `objective` | Optimized MDPDE objective. |
| `score` | Estimating equation at the estimate. |
| `starting.values` | Starting vector used by the optimizer. |
| `automatic.start` | Whether automatic starting values were used. |
| `gamma.structure`, `gamma.diagonal.names` | Automatic-start structure information. |
| `robust.start.fit` | `MASS::rlm()` object used for automatic initialization. |
| `q`, `level` | Tuning parameter and confidence level. |
| `n.groups`, `nobs`, `group.sizes` | Group and sample dimensions. |
| `p.beta`, `p.gamma`, `m` | Parameter and random-effect dimensions. |
| `beta.names`, `gamma.names`, `sigma2.name` | Parameter names. |
| `convergence`, `message`, `counts` | Optimizer diagnostics. |
| `method`, `jacobian.method` | Optimization and differentiation methods. |
| `Y`, `X`, `Z` | Stored grouped data. |
| `Delta.function`, `Delta.jacobian.function` | Stored covariance functions. |
| `vcov.warning` | Covariance-inversion warning, or `NULL`. |
| `lower`, `upper` | Bounds used by bounded optimization. |
| `optimization.control` | Effective optimizer controls. |
| `optim` | Complete optimizer result. |
| `call` | Matched function call. |

SQV fits additionally contain `selected.q`, `sqv.returned.fit.source`, and `q.selection`.

## S3 methods

```r
print(fit_mdpde)
summary(fit_mdpde)

coef(fit_mdpde)
vcov(fit_mdpde)

confint(
  fit_mdpde,
  parm = c(
    "beta0",
    "gamma1"
  ),
  level = 0.90
)

fitted(
  fit_mdpde,
  type = "marginal"
)

fitted(
  fit_mdpde,
  type = "conditional"
)

residuals(
  fit_mdpde,
  type = "marginal"
)

residuals(
  fit_mdpde,
  type = "conditional"
)
```

## Prediction

```r
predict(
  fit_mdpde,
  newdata = NULL,
  type = c(
    "mean",
    "covariance",
    "random_effects",
    "conditional_mean"
  )
)
```

For new groups, `newdata$X` and `newdata$Z` are required. `newdata$Y` is additionally required for random-effect or conditional-mean predictions.

## Minimal call

```r
fit_mdpde <- fit_lmm_mdpde(
  Y = Y,
  X = X,
  Z = Z,
  Delta = Delta_unstructured,
  start = NULL,
  gamma_structure =
    gamma_structure_unstructured,
  Delta_jacobian =
    Delta_unstructured_jacobian,
  q = "SQV",
  level = 0.95,
  method = "BFGS",
  control = list(
    maxit = 2000,
    reltol = 1e-12
  ),
  use_score = TRUE,
  compute_vcov = TRUE,
  jacobian_method = "Richardson",
  q_control = list(
    m0 = 21L,
    m = 3L,
    q_min = 0.5,
    L = 0.01,
    keep_fits = TRUE
  )
)
```

## Gaussian MLE through MDPDE

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
  q = 1,
  method = "BFGS",
  use_score = TRUE,
  compute_vcov = TRUE
)
```

# 6. `fit_lmm_t()`

## Purpose

Fits a linear mixed model by maximum likelihood under a multivariate Student-t marginal distribution with known, fixed degrees of freedom.

For group \(i\):

```text
Y_i ~ t_{d_i}(
  X_i beta,
  Z_i Delta(gamma) Z_i^T +
    sigma2 I_{d_i},
  nu
).
```

The second argument of the multivariate Student-t distribution is a scale matrix.

## Usage

```r
fit_lmm_t(
  Y,
  X,
  Z,
  Delta,
  start = NULL,
  gamma_structure = NULL,
  Delta_jacobian = NULL,
  nu,
  level = 0.95,
  method = c(
    "BFGS",
    "L-BFGS-B",
    "Nelder-Mead",
    "CG",
    "BOBYQA"
  ),
  lower = NULL,
  upper = NULL,
  control = list(),
  use_score = TRUE,
  compute_vcov = TRUE,
  hessian = TRUE,
  jacobian_method = c(
    "Richardson",
    "simple"
  ),
  jacobian_method_args = list(),
  derivative_eps =
    .Machine$double.eps^(1 / 3)
)
```

## Arguments

| Argument | Description |
|---|---|
| `Y` | Grouped response vectors. |
| `X` | Grouped fixed-effects design matrices. |
| `Z` | Grouped random-effects design matrices. |
| `Delta` | Function returning the random-effects Student-t scale matrix. |
| `start` | Optional vector ordered as `(beta, gamma, sigma2)`. |
| `gamma_structure` | Symmetric character matrix used for automatic gamma starting values when `start = NULL`. |
| `Delta_jacobian` | Optional analytical derivative of `vec(Delta(gamma))`. |
| `nu` | Fixed positive degrees of freedom. It is not estimated. |
| `level` | Confidence level in `(0, 1)`. |
| `method` | `"BFGS"`, `"L-BFGS-B"`, `"Nelder-Mead"`, `"CG"`, or `"BOBYQA"`. |
| `lower`, `upper` | Bounds supported by L-BFGS-B and BOBYQA. |
| `control` | Optimizer control list. |
| `use_score` | Uses the analytical negative score for BFGS, CG, or L-BFGS-B. It is not used by Nelder-Mead or BOBYQA. |
| `compute_vcov` | Computes expected Fisher information, its inverse, standard errors, and confidence intervals. |
| `hessian` | Computes the observed information after fitting when `TRUE`. |
| `jacobian_method` | Numerical method used when `Delta_jacobian` is omitted. |
| `jacobian_method_args` | Additional arguments passed to `numDeriv::jacobian()`. |
| `derivative_eps` | Positive relative finite-difference step used by numerical fallback derivatives. |

## Optimization defaults

For BFGS, CG, and Nelder-Mead:

```r
list(
  maxit = 2000,
  reltol = 1e-9
)
```

For L-BFGS-B:

```r
list(
  maxit = 2000,
  factr = 1e7,
  pgtol = 0
)
```

For BOBYQA:

```r
list(
  maxeval = 100000,
  xtol_rel = 1e-8,
  ftol_rel = 1e-8
)
```

For bounded methods, only `sigma2` receives a positive default lower bound. The supplied `Delta(gamma)` must still produce a positive-definite scale matrix.

## Returned object

The returned object has class:

```text
"lmm_t_fit"
```

### Estimates and scale quantities

| Component | Description |
|---|---|
| `coefficients` | Student-t maximum-likelihood estimate. |
| `beta` | Fixed-effect estimate. |
| `gamma` | Random-effects scale-parameter estimate. |
| `sigma2` | Residual scale-squared estimate. |
| `Delta` | Estimated random-effects scale matrix. |
| `Delta.variance` | Implied random-effects covariance matrix when `nu > 2`; otherwise `NULL`. |
| `residual.variance` | Implied residual variance when `nu > 2`; otherwise `NA`. |

### Inference

| Component | Description |
|---|---|
| `standard.error` | Standard errors from inverse expected Fisher information. |
| `conf.int` | Wald confidence intervals. |
| `coefficient.table` | Estimates, standard errors, z statistics, and two-sided normal-approximation p-values. |
| `vcov`, `vcov.fisher` | Inverse expected Fisher information. |
| `vcov.observed` | Inverse observed information when available. |
| `fisher.information` | Expected Fisher information matrix. |
| `observed.information` | Observed information matrix when requested. |

### Derivatives and fitted quantities

| Component | Description |
|---|---|
| `Delta.jacobian` | Derivative of `vec(Delta)` at the estimate. |
| `mean.jacobians` | List of group-specific mean derivative matrices. |
| `covariance.jacobians` | List of group-specific scale-matrix derivative matrices. |
| `F` | List of combined derivative matrices. |
| `fitted.values`, `marginal.fitted.values` | Group-specific marginal fitted locations. |
| `conditional.fitted.values` | Group-specific conditional fitted values. |
| `random.effects` | Predicted random-effect vectors. |
| `residuals`, `marginal.residuals` | Marginal residual vectors. |
| `conditional.residuals` | Conditional residual vectors. |
| `standardized.residuals` | Residuals standardized by the fitted group scale matrix. |
| `scale.matrices` | Fitted group Student-t scale matrices. |
| `covariance.matrices` | Implied group covariance matrices when `nu > 2`; otherwise a list of `NULL` entries. |
| `scale.inverses` | Inverses of the fitted scale matrices. |
| `quadratic.forms` | Group-level squared Mahalanobis distances. |
| `weights` | Group-level Student-t likelihood weights. |
| `log.density`, `density` | Group-level Student-t log-densities and densities. |

### Likelihood and diagnostics

| Component | Description |
|---|---|
| `score` | Score vector at the estimate. |
| `max.abs.score` | Largest absolute score component. |
| `logLik` | Maximized log-likelihood. |
| `objective` | Minimized negative log-likelihood. |
| `AIC` | Akaike information criterion. |
| `BIC`, `BIC.groups` | BIC using the number of groups as the sample size. |
| `BIC.observations` | BIC using the total number of observations. |
| `nu` | Fixed degrees of freedom. |
| `level` | Confidence level used during fitting. |
| `n.groups`, `nobs`, `group.sizes` | Number of groups, total observations, and group sizes. |
| `p`, `p.beta`, `p.gamma`, `m` | Parameter and design dimensions. |
| `beta.names`, `gamma.names`, `sigma2.name` | Parameter names. |
| `random.effect.names` | Names assigned to predicted random effects. |
| `convergence`, `message`, `counts` | Optimizer diagnostics. |
| `method`, `use.score`, `jacobian.method` | Fitting configuration. |
| `starting.values`, `automatic.start` | Starting-value information. |
| `gamma.structure`, `gamma.diagonal.names` | Automatic-start structure. |
| `robust.start.fit` | `MASS::rlm()` object used for automatic initialization. |
| `Y`, `X`, `Z` | Stored grouped data. |
| `Delta.function`, `Delta.jacobian.function` | Stored model functions. |
| `vcov.warning` | Expected-information inversion warning, or `NULL`. |
| `observed.information.warning` | Observed-information warning, or `NULL`. |
| `lower`, `upper` | Bounds used by bounded optimization. |
| `optimization.control` | Effective optimizer controls. |
| `optim` | Complete optimizer result. |
| `call` | Matched function call. |

## Student-t group weights

For a group of dimension \(d_i\):

```text
weight_i =
  (nu + d_i) /
  (nu + delta_i),
```

where `delta_i` is the fitted group-level quadratic form.

Access the values with:

```r
fit_t$weights
```

## S3 methods

```r
print(fit_t)
summary(fit_t)

coef(fit_t)

vcov(
  fit_t,
  type = "fisher"
)

vcov(
  fit_t,
  type = "observed"
)

confint(
  fit_t,
  parm = c(
    "beta0",
    "gamma1"
  ),
  level = 0.90
)

fitted(
  fit_t,
  type = "marginal"
)

fitted(
  fit_t,
  type = "conditional"
)

residuals(
  fit_t,
  type = "marginal"
)

residuals(
  fit_t,
  type = "conditional"
)

residuals(
  fit_t,
  type = "standardized"
)

logLik(fit_t)
nobs(fit_t)
```

`nobs(fit_t)` returns the total number of scalar observations. The `logLik` object uses the number of independent groups in its `nobs` attribute.

## Prediction

```r
predict(
  fit_t,
  newdata = NULL,
  type = c(
    "mean",
    "scale",
    "variance",
    "random_effects",
    "conditional_mean"
  )
)
```

Without `newdata`:

| `type` | Output |
|---|---|
| `"mean"` | Stored marginal fitted locations. |
| `"scale"` | Stored group scale matrices. |
| `"variance"` | Implied group covariance matrices, available only for `nu > 2`. |
| `"random_effects"` | Stored predicted random effects. |
| `"conditional_mean"` | Stored conditional fitted values. |

For new groups, `newdata` must contain `X` and `Z`. `newdata$Y` is additionally required for random-effect and conditional-mean predictions.

## Minimal BFGS call

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
  nu = 4,
  level = 0.95,
  method = "BFGS",
  control = list(
    maxit = 2000,
    reltol = 1e-12
  ),
  use_score = TRUE,
  compute_vcov = TRUE,
  hessian = TRUE,
  jacobian_method = "Richardson"
)
```

## Minimal BOBYQA call

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
  nu = 4,
  level = 0.95,
  method = "BOBYQA",
  control = list(
    maxeval = 100000,
    xtol_rel = 0,
    ftol_rel = 0,
    xtol_abs = 1e-8,
    ftol_abs = 1e-8
  ),
  use_score = FALSE,
  compute_vcov = TRUE,
  hessian = TRUE
)
```

# 7. Group-level density-power weights

The Gaussian fit objects store one marginal density per group.

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

The MLqE weight is calculated in the original parameterization because that is the parameterization used by the Lq-likelihood.

# 8. Method comparison

| Feature | Corrected MLqE | MDPDE | Student-t MLE |
|---|---|---|---|
| Marginal distribution | Multivariate normal | Multivariate normal | Multivariate Student-t |
| Robustness control | `q` | `q` | Fixed `nu` |
| Automatic tuning | SQV | SQV | Not implemented |
| Consistency transformation | User supplied | Not required | Not applicable |
| `Delta` interpretation | Covariance | Covariance | Scale matrix |
| `sigma2` interpretation | Residual variance | Residual variance | Residual scale squared |
| Default covariance estimator | Corrected sandwich | Sandwich | Inverse expected Fisher information |
| Observed-information covariance | Not separately returned | Not separately returned | Available with `hessian = TRUE` |
| Normal MLE special case | `q = 1` | `q = 1` | Approached as `nu` increases, but `nu` is fixed |

# 9. Common errors

## Different numbers of groups

`Y`, `X`, and `Z` must contain the same number of groups.

## Incompatible group dimensions

For every group:

```text
length(Y[[i]]) =
nrow(X[[i]]) =
nrow(Z[[i]]).
```

## Different design-matrix column counts

All fixed-effect matrices must have the same number of columns. The same applies to all random-effect matrices.

## Missing `gamma_structure`

It is required when:

```r
start = NULL
```

## Invalid `gamma_structure`

It must be square, symmetric, and contain valid parameter names or missing entries.

## Invalid `Delta(gamma)`

The returned matrix must be finite, symmetric, and positive definite. A direct unstructured parameterization can become indefinite during optimization.

## Incorrect `Delta_jacobian`

Its dimensions must correspond to:

```text
m^2 by p_gamma
```

or:

```text
m by m by p_gamma.
```

## Automatic starting-value failure

This may occur if `MASS::rlm()` fails or returns an invalid scale. Supply `start` explicitly.

## SQV without covariance estimation

Use:

```r
compute_vcov = TRUE
```

when `q = "SQV"`.

## MLqE without `tau`

A transformation function is required when `q < 1` or `q = "SQV"`.

## BOBYQA unavailable

Install:

```r
install.packages("nloptr")
```

## Numerical derivative package unavailable

Install:

```r
install.packages("numDeriv")
```

when `Delta_jacobian` is omitted or when corrected MLqE is fitted.

# 10. Recommended post-fit checks

For corrected MLqE:

```r
fit_mlqe$convergence
fit_mlqe$message
fit_mlqe$counts
fit_mlqe$q
fit_mlqe$score.star
fit_mlqe$score.corrected
fit_mlqe$tau.composition.residual
fit_mlqe$vcov.warning
```

For MDPDE:

```r
fit_mdpde$convergence
fit_mdpde$message
fit_mdpde$counts
fit_mdpde$q
fit_mdpde$score
fit_mdpde$vcov.warning
```

For Student-t MLE:

```r
fit_t$convergence
fit_t$message
fit_t$counts
fit_t$max.abs.score
fit_t$fisher.information
fit_t$observed.information
fit_t$vcov.warning
fit_t$observed.information.warning
```

For SQV fits:

```r
fit$q.selection$reason
fit$q.selection$history
fit$q.selection$evaluations
```

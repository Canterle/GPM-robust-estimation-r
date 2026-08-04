---
output:
  pdf_document: default
  html_document: default
---
# Function Reference: Linear Measurement-Error Model

This document describes the public fitting functions:

- `fit_lmve_mlqe()`;
- `fit_lmve_mdpde()`;
- `fit_lmve_t()`.

All functions estimate the same five-parameter structural model:

```text
theta = (beta0, beta1, mu_x, sigma2_x, sigma2).
```

## 1. Common statistical structure

The latent structural relationship is:

```text
Y_i* = beta0 + beta1 X_i* + epsilon_i.
```

The observed variables are:

```text
Y_i = Y_i* + e_yi,
X_i = X_i* + e_xi.
```

For the normal model:

```text
X_i*      ~ Normal(mu_x, sigma2_x),
epsilon_i ~ Normal(0, sigma2),

Var(e_yi) = tau_yi,
Var(e_xi) = tau_xi.
```

The marginal mean of the observed pair is:

```text
mu(theta) =
  (beta0 + beta1 * mu_x, mu_x)^T.
```

The observation-specific covariance matrix is:

```text
Sigma_i(theta) =

  [ beta1^2 * sigma2_x + sigma2 + tau_yi,
    beta1 * sigma2_x                              ]

  [ beta1 * sigma2_x,
    sigma2_x + tau_xi                            ].
```

For the Student-t model, this same matrix is a scale matrix. It is not a covariance matrix.

## 2. Common data and parameter arguments

### `Y`

Observed response vector. It must be finite, non-missing, and nonempty.

### `X`

Observed explanatory-variable vector. It must be finite, non-missing, and have the same length as `Y`.

### `tau_y`

Known observation-specific measurement-error variance under the normal model or known scale-squared quantity under the Student-t model.

A scalar is recycled to `length(Y)`.

### `tau_x`

Known observation-specific measurement-error variance under the normal model or known scale-squared quantity under the Student-t model.

A scalar is recycled to `length(Y)`.

Both `tau_y` and `tau_x` must be finite and non-negative.

### `start`

Starting-value vector:

```r
c(
  beta0 = ...,
  beta1 = ...,
  mu_x = ...,
  sigma2_x = ...,
  sigma2 = ...
)
```

It must contain five finite values. `sigma2_x` and `sigma2` must be strictly positive.

When `start = NULL`, robust moment-based starting values are computed internally.

## 3. `fit_lmve_mlqe()`

### Purpose

Fits the linear normal measurement-error model by maximum Lq-likelihood estimation and returns both the original MLqE parameterization and its Fisher-consistency correction.

For `q = 1`, the procedure reduces to ordinary normal maximum likelihood.

### Usage

```r
fit_lmve_mlqe(
  Y,
  X,
  tau_y,
  tau_x,
  start = NULL,
  q = 1,
  level = 0.95,
  method = c(
    "BFGS",
    "L-BFGS-B",
    "Nelder-Mead"
  ),
  lower = NULL,
  upper = NULL,
  control = list(),
  use_score = FALSE,
  compute_vcov = TRUE,
  q_control = list()
)
```

### Arguments

| Argument | Description |
|---|---|
| `Y` | Observed response vector. |
| `X` | Observed explanatory-variable vector. |
| `tau_y` | Known measurement-error variances for `Y`. |
| `tau_x` | Known measurement-error variances for `X`. |
| `start` | Optional five-parameter corrected starting vector. When `NULL`, robust moment-based values are used. |
| `q` | Numeric value in `(0, 1]`, or `"SQV"` for automatic selection. |
| `level` | Confidence level in `(0, 1)` used for `fit$conf.int`. |
| `method` | Optimization method: `"BFGS"`, `"L-BFGS-B"`, or `"Nelder-Mead"`. |
| `lower`, `upper` | Bounds used only with `"L-BFGS-B"`. User-supplied MLqE bounds refer to the corrected parameters and are transformed internally to the original parameterization. |
| `control` | Control list passed to `optim()`. |
| `use_score` | If `TRUE`, supplies the analytical estimating-equation vector to BFGS or L-BFGS-B. It is ignored for Nelder-Mead. |
| `compute_vcov` | If `TRUE`, computes sensitivity, variability, sandwich covariance matrices, standard errors, and confidence intervals. Required for `"SQV"`. |
| `q_control` | SQV controls used when `q = "SQV"`. |

The function also supports an unambiguous positional shortcut:

```r
fit_lmve_mlqe(
  Y,
  X,
  tau_y,
  tau_x,
  "SQV"
)
```

or:

```r
fit_lmve_mlqe(
  Y,
  X,
  tau_y,
  tau_x,
  0.8
)
```

### Internal consistency transformation

The transformation is built into the implementation:

```text
tau_r(theta) =
  (beta0, beta1, mu_x, sigma2_x / r, sigma2 / r).
```

The original parameterization is:

```text
theta_star =
  tau_(1/q)(theta) =
  (beta0, beta1, mu_x, q * sigma2_x, q * sigma2).
```

The known measurement-error variances are transformed as:

```text
tau_y_star = q * tau_y,
tau_x_star = q * tau_x.
```

The Jacobian of the forward mapping is:

```text
A = diag(1, 1, 1, q, q).
```

The corrected covariance matrix is:

```text
V = A^(-1) V_star A^(-T).
```

### Optimization defaults

For `method = "L-BFGS-B"`:

```r
list(
  maxit = 1000,
  factr = 1e7,
  pgtol = 0
)
```

The default bounds are:

```r
c(
  beta0 = -Inf,
  beta1 = -Inf,
  mu_x = -Inf,
  sigma2_x = 1e-10,
  sigma2 = 1e-10
)
```

with infinite upper bounds.

For BFGS and Nelder-Mead:

```r
list(
  maxit = 1000,
  reltol = 1e-9
)
```

The objective is maximized internally by setting `control$fnscale = -1`.

### SQV controls

The defaults are:

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
| `m0` | Number of points in the initial grid from `1` to `0.8`. |
| `m` | Number of points in each refinement grid. |
| `q_min` | Smallest allowed `q`; it must be positive and less than `0.8`. |
| `L` | Positive stability threshold. |
| `verbose` | Prints SQV progress when `TRUE`. |
| `keep_fits` | Retains the evaluated fixed-`q` fits when `TRUE`. |

Only these entries are retained from `q_control`.

### Returned object

The returned object has class:

```text
"lmve_mlqe"
```

An SQV-selected fit additionally has class:

```text
"lmve_mlqe_sqv"
```

#### Estimates and inference

| Component | Description |
|---|---|
| `coefficients` | Corrected MLqE estimate. |
| `coefficients.star` | Original MLqE estimate. |
| `standard.error` | Corrected sandwich standard errors. |
| `standard.error.star` | Original-parameterization sandwich standard errors. |
| `vcov` | Corrected sandwich covariance matrix. |
| `vcov.star` | Sandwich covariance matrix for the original MLqE estimate. |
| `vcov.direct` | Corrected covariance computed through the direct correction Jacobian. |
| `vcov.direct.difference` | Numerical comparison of the two corrected covariance calculations. |
| `conf.int` | Wald confidence intervals for corrected parameters. |
| `conf.int.star` | Wald confidence intervals for original parameters. |

#### Estimating-function matrices and derivatives

| Component | Description |
|---|---|
| `J`, `J.star` | Sensitivity matrix in the original parameterization. |
| `K`, `K.star` | Variability matrix in the original parameterization. |
| `mean.jacobian.star` | Derivative of the bivariate mean with respect to the original parameter vector. |
| `covariance.jacobian.star` | Derivative of `vec(Sigma_i)` with respect to the original parameter vector. |
| `F.star` | Combined derivative matrix formed from the mean and covariance derivatives. |
| `tau.jacobian` | Jacobian `A` of the corrected-to-original transformation. |
| `tau.jacobian.inverse` | Inverse of `A`. |
| `tau.q.jacobian` | Direct correction Jacobian. |
| `tau.composition.residual` | Numerical check of the transformation composition. |
| `score`, `score.star` | Estimating-equation vector in the original parameterization. |
| `score.corrected` | Estimating-equation vector evaluated in the corrected parameterization. |

#### Corrected fitted quantities

| Component | Description |
|---|---|
| `fitted.values` | `n x 2` matrix of fitted marginal locations for `Y` and `X`. |
| `fitted.y`, `fitted.x` | Fitted marginal location vectors. |
| `covariance.array` | `2 x 2 x n` array of corrected fitted covariance matrices. |
| `variance.y`, `variance.x` | Corrected fitted marginal variances. |
| `covariance.yx` | Corrected fitted covariance between `Y` and `X`. |
| `sigma.y`, `sigma.x` | Square roots of the corrected marginal variances. |
| `residuals` | `n x 2` matrix of corrected marginal residuals. |
| `residual.y`, `residual.x` | Corrected residual vectors. |
| `log.density`, `density` | Corrected fitted bivariate normal log-densities and densities. |

The corresponding original-parameterization components use the `.star` suffix, including:

```text
fitted.values.star
fitted.y.star
fitted.x.star
covariance.array.star
variance.y.star
variance.x.star
covariance.yx.star
sigma.y.star
sigma.x.star
residuals.star
residual.y.star
residual.x.star
log.density.star
density.star
```

#### Optimization and metadata

| Component | Description |
|---|---|
| `starting.values` | Corrected starting vector. |
| `starting.values.star` | Transformed original-parameterization starting vector. |
| `automatic.start` | Whether starting values were generated internally. |
| `objective` | Maximized Lq-likelihood objective. |
| `convergence` | Optimizer convergence code. |
| `message` | Optimizer message. |
| `counts` | Function and gradient evaluation counts. |
| `q` | Fixed or selected `q`. |
| `level` | Confidence level used during fitting. |
| `nobs` | Number of observed pairs. |
| `method` | Optimization method. |
| `Y`, `X`, `tau_y`, `tau_x` | Stored data. |
| `tau_y.star`, `tau_x.star` | Transformed known measurement-error variances. |
| `tau` | Internal transformation function. |
| `vcov.warning` | Covariance-inversion warning, or `NULL`. |
| `optimization.control` | Effective optimization controls. |
| `optim` | Complete optimizer result. |
| `call` | Matched function call. |

SQV fits additionally contain:

```text
selected.q
sqv.returned.fit.source
q.selection
```

The `q.selection` object contains the selection reason, pass, stage, configuration, grid history, SQV history, evaluations, and optionally cached fits.

### S3 methods

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
    "beta1"
  ),
  level = 0.90,
  type = "star"
)

fitted(fit_mlqe)
fitted(
  fit_mlqe,
  type = "star"
)

residuals(fit_mlqe)
residuals(
  fit_mlqe,
  type = "star"
)
```

### Prediction

```r
predict(
  fit_mlqe,
  newdata = NULL,
  type = c(
    "mean",
    "mean_y",
    "mean_x",
    "structural_mean",
    "covariance",
    "variance_y",
    "variance_x",
    "covariance_yx",
    "sigma_y",
    "sigma_x"
  ),
  parameterization = c(
    "corrected",
    "star"
  )
)
```

Prediction types:

| `type` | Output |
|---|---|
| `"mean"` | Matrix with marginal fitted locations for `Y` and `X`. |
| `"mean_y"` | Marginal location of `Y`. |
| `"mean_x"` | Marginal location of `X`. |
| `"structural_mean"` | Structural line `beta0 + beta1 * latent_x`. |
| `"covariance"` | `2 x 2 x n` covariance array. |
| `"variance_y"` | Marginal variance of observed `Y`. |
| `"variance_x"` | Marginal variance of observed `X`. |
| `"covariance_yx"` | Covariance between observed `Y` and `X`. |
| `"sigma_y"` | Square root of `variance_y`. |
| `"sigma_x"` | Square root of `variance_x`. |

For structural predictions:

```r
predict(
  fit_mlqe,
  newdata = list(
    latent_x = c(
      -2,
      0,
      2
    )
  ),
  type = "structural_mean"
)
```

`newdata$x` may be used instead of `newdata$latent_x`.

For covariance predictions at new measurement-error variances:

```r
predict(
  fit_mlqe,
  newdata = list(
    tau_y = c(
      1,
      2
    ),
    tau_x = c(
      0.5,
      0.8
    )
  ),
  type = "covariance"
)
```

### Minimal call

```r
fit_mlqe <- fit_lmve_mlqe(
  Y = Y,
  X = X,
  tau_y = tau_y,
  tau_x = tau_x,
  start = NULL,
  q = "SQV",
  level = 0.95,
  method = "L-BFGS-B",
  lower = c(
    beta0 = -Inf,
    beta1 = -Inf,
    mu_x = -Inf,
    sigma2_x = 1e-8,
    sigma2 = 1e-8
  ),
  upper = rep(
    Inf,
    5
  ),
  control = list(
    maxit = 2000,
    factr = 1e4,
    pgtol = 1e-12
  ),
  use_score = TRUE,
  compute_vcov = TRUE,
  q_control = list(
    m0 = 21L,
    m = 3L,
    q_min = 0.5,
    L = 0.01,
    keep_fits = TRUE
  )
)
```

## 4. `fit_lmve_mdpde()`

### Purpose

Fits the linear normal measurement-error model by maximum density power divergence estimation.

For `q = 1`, it reduces to ordinary maximum likelihood under the bivariate normal measurement-error model.

### Usage

```r
fit_lmve_mdpde(
  Y,
  X,
  tau_y,
  tau_x,
  start = NULL,
  q = 1,
  level = 0.95,
  method = c(
    "BFGS",
    "L-BFGS-B",
    "Nelder-Mead"
  ),
  lower = NULL,
  upper = NULL,
  control = list(),
  use_score = FALSE,
  compute_vcov = TRUE,
  q_control = list()
)
```

### Arguments

The arguments have the same meanings and validation rules as in `fit_lmve_mlqe()`, except that no consistency transformation is required.

The function also supports:

```r
fit_lmve_mdpde(
  Y,
  X,
  tau_y,
  tau_x,
  "SQV"
)
```

and:

```r
fit_lmve_mdpde(
  Y,
  X,
  tau_y,
  tau_x,
  0.8
)
```

### Optimization defaults

For `L-BFGS-B`:

```r
list(
  maxit = 1000,
  factr = 1e7,
  pgtol = 0
)
```

with default lower bounds:

```r
c(
  beta0 = -Inf,
  beta1 = -Inf,
  mu_x = -Inf,
  sigma2_x = 1e-10,
  sigma2 = 1e-10
)
```

For BFGS and Nelder-Mead:

```r
list(
  maxit = 1000,
  reltol = 1e-9
)
```

The SQV controls and defaults are the same as for `fit_lmve_mlqe()`.

### Returned object

The returned object has class:

```text
"lmve_mdpde"
```

An SQV-selected fit additionally has class:

```text
"lmve_mdpde_sqv"
```

#### Main components

| Component | Description |
|---|---|
| `coefficients` | MDPDE estimate. |
| `standard.error` | Sandwich standard errors. |
| `vcov` | Sandwich covariance matrix. |
| `conf.int` | Wald confidence intervals. |
| `J` | Sensitivity matrix. |
| `K` | Variability matrix. |
| `mean.jacobian` | Mean derivative matrix. |
| `covariance.jacobian` | Covariance derivative matrix. |
| `F` | Combined mean and covariance derivative matrix. |
| `fitted.values` | `n x 2` matrix of fitted marginal locations. |
| `fitted.y`, `fitted.x` | Fitted marginal location vectors. |
| `covariance.array` | `2 x 2 x n` array of fitted covariance matrices. |
| `variance.y`, `variance.x` | Fitted marginal variances. |
| `covariance.yx` | Fitted covariance between `Y` and `X`. |
| `sigma.y`, `sigma.x` | Marginal standard deviations. |
| `residuals` | `n x 2` matrix of marginal residuals. |
| `residual.y`, `residual.x` | Marginal residual vectors. |
| `log.density`, `density` | Fitted bivariate normal log-densities and densities. |
| `objective` | Maximized MDPDE objective. |
| `score` | Estimating-equation vector at the estimate. |
| `convergence`, `message`, `counts` | Optimizer diagnostics. |
| `q` | Fixed or selected `q`. |
| `level`, `nobs`, `method` | Fit metadata. |
| `starting.values` | Starting vector used by the optimizer. |
| `automatic.start` | Whether starting values were generated internally. |
| `Y`, `X`, `tau_y`, `tau_x` | Stored data. |
| `vcov.warning` | Covariance-inversion warning, or `NULL`. |
| `optimization.control` | Effective optimization controls. |
| `optim` | Complete optimizer result. |
| `call` | Matched function call. |

SQV fits additionally contain `selected.q`, `sqv.returned.fit.source`, and `q.selection`.

### S3 methods

```r
print(fit_mdpde)
summary(fit_mdpde)

coef(fit_mdpde)
vcov(fit_mdpde)

confint(
  fit_mdpde,
  parm = c(
    "beta0",
    "beta1"
  ),
  level = 0.90
)

fitted(fit_mdpde)
residuals(fit_mdpde)
```

### Prediction

```r
predict(
  fit_mdpde,
  newdata = NULL,
  type = c(
    "mean",
    "mean_y",
    "mean_x",
    "structural_mean",
    "covariance",
    "variance_y",
    "variance_x",
    "covariance_yx",
    "sigma_y",
    "sigma_x"
  )
)
```

The prediction types have the same meanings as for the corrected MLqE.

### Minimal call

```r
fit_mdpde <- fit_lmve_mdpde(
  Y = Y,
  X = X,
  tau_y = tau_y,
  tau_x = tau_x,
  start = NULL,
  q = "SQV",
  level = 0.95,
  method = "L-BFGS-B",
  lower = c(
    beta0 = -Inf,
    beta1 = -Inf,
    mu_x = -Inf,
    sigma2_x = 1e-8,
    sigma2 = 1e-8
  ),
  upper = rep(
    Inf,
    5
  ),
  control = list(
    maxit = 2000,
    factr = 1e4,
    pgtol = 1e-12
  ),
  use_score = TRUE,
  compute_vcov = TRUE,
  q_control = list(
    m0 = 21L,
    m = 3L,
    q_min = 0.5,
    L = 0.01,
    keep_fits = TRUE
  )
)
```

### Normal MLE through MDPDE

```r
fit_mle <- fit_lmve_mdpde(
  Y = Y,
  X = X,
  tau_y = tau_y,
  tau_x = tau_x,
  start = NULL,
  q = 1,
  level = 0.95,
  method = "L-BFGS-B",
  lower = c(
    beta0 = -Inf,
    beta1 = -Inf,
    mu_x = -Inf,
    sigma2_x = 1e-8,
    sigma2 = 1e-8
  ),
  upper = rep(
    Inf,
    5
  ),
  use_score = TRUE,
  compute_vcov = TRUE
)
```

## 5. `fit_lmve_t()`

### Purpose

Fits the linear measurement-error model by maximum likelihood under a bivariate Student-t distribution with fixed degrees of freedom.

The observed pair satisfies:

```text
(Y_i, X_i)^T ~ t_2(
  mu(theta),
  Sigma_i(theta),
  nu
).
```

`Sigma_i(theta)` is the Student-t scale matrix.

When `nu > 2`:

```text
Cov[(Y_i, X_i)^T] =
  nu / (nu - 2) * Sigma_i(theta).
```

### Usage

```r
fit_lmve_t(
  Y,
  X,
  tau_y,
  tau_x,
  start = NULL,
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
  derivative_eps =
    .Machine$double.eps^(1 / 3)
)
```

### Arguments

| Argument | Description |
|---|---|
| `Y` | Observed response vector. |
| `X` | Observed explanatory-variable vector. |
| `tau_y` | Known non-negative Student-t measurement-error scale-squared quantities for `Y`. |
| `tau_x` | Known non-negative Student-t measurement-error scale-squared quantities for `X`. |
| `start` | Optional five-parameter starting vector. When `NULL`, robust moment-based starting values are calculated internally. |
| `nu` | Fixed positive degrees of freedom. It is not estimated. |
| `level` | Confidence level used for expected-information intervals stored in `fit$conf.int`. |
| `method` | `"BFGS"`, `"L-BFGS-B"`, `"Nelder-Mead"`, `"CG"`, or `"BOBYQA"`. |
| `lower`, `upper` | Parameter bounds. Supported for `"L-BFGS-B"` and `"BOBYQA"`. |
| `control` | Optimizer control list. |
| `use_score` | Uses the analytical negative score for BFGS, CG, and L-BFGS-B. Ignored for Nelder-Mead and BOBYQA. |
| `compute_vcov` | Computes the inverse expected Fisher information, standard errors, and confidence intervals. |
| `hessian` | If `TRUE`, computes the observed information using `optimHess()` after fitting. |
| `derivative_eps` | Positive relative finite-difference step used for numerical fallback gradients and Hessian-related calculations. |

### Optimization defaults

For `L-BFGS-B`:

```r
list(
  maxit = 2000,
  factr = 1e7,
  pgtol = 0
)
```

For BFGS, CG, and Nelder-Mead:

```r
list(
  maxit = 2000,
  reltol = 1e-9
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

If `control$maxit`, `control$reltol`, or `control$abstol` are supplied for BOBYQA, they are translated to `maxeval`, `ftol_rel`, or `ftol_abs` when the corresponding BOBYQA control is absent.

The default bounded lower vector is:

```r
c(
  beta0 = -Inf,
  beta1 = -Inf,
  mu_x = -Inf,
  sigma2_x = 1e-10,
  sigma2 = 1e-10
)
```

### Returned object

The returned object has class:

```text
"lmve_t_fit"
```

#### Estimates and inference

| Component | Description |
|---|---|
| `coefficients` | Student-t MLE. |
| `beta` | Vector containing `beta0` and `beta1`. |
| `mu.x` | Estimate of `mu_x`. |
| `sigma2.x`, `sigma2` | Estimated scale-squared parameters. |
| `sigma2.x.scale`, `sigma2.scale` | Explicit aliases for the estimated latent-X and equation-error scale-squared quantities. |
| `latent.x.variance` | Implied latent-X variance when `nu > 2`; otherwise `NA`. |
| `equation.error.variance` | Implied equation-error variance when `nu > 2`; otherwise `NA`. |
| `standard.error` | Standard errors from the inverse expected Fisher information. |
| `conf.int` | Wald confidence intervals based on expected Fisher information. |
| `coefficient.table` | Estimates, standard errors, z statistics, and two-sided normal-approximation p-values. |
| `vcov`, `vcov.fisher` | Inverse expected Fisher information. |
| `vcov.observed` | Inverse observed information when available. |
| `fisher.information` | Expected Fisher information matrix. |
| `observed.information` | Numerical observed information matrix when requested. |

#### Derivatives and fitted distributions

| Component | Description |
|---|---|
| `mean.jacobian` | Derivative of the bivariate location with respect to `theta`. |
| `scale.jacobian` | Derivative of `vec(Sigma_i)` with respect to `theta`. |
| `F` | Combined mean and scale derivative matrix. |
| `fitted.values` | `n x 2` matrix of fitted marginal locations. |
| `fitted.y`, `fitted.x` | Fitted marginal location vectors. |
| `residuals` | `n x 2` matrix of marginal residuals. |
| `residual.y`, `residual.x` | Marginal residual vectors. |
| `scale.y`, `scale.x`, `scale.yx` | Observation-specific fitted Student-t scale-matrix elements. |
| `scale.array` | `2 x 2 x n` array of fitted Student-t scale matrices. |
| `covariance.array` | `2 x 2 x n` covariance array when `nu > 2`; otherwise `NULL`. |
| `scale.inverse.array` | Inverse of each fitted scale matrix. |
| `quadratic.forms` | Observation-specific squared Mahalanobis distances. |
| `weights` | Student-t likelihood weights. |
| `log.density`, `density` | Fitted bivariate Student-t log-densities and densities. |

#### Conditional Student-t distribution

| Component | Description |
|---|---|
| `conditional.observed.location.y.given.x` | Conditional location of observed `Y` given observed `X`. |
| `conditional.observed.scale2.y.given.x` | Conditional Student-t scale-squared value. |
| `conditional.observed.df.y.given.x` | Conditional degrees of freedom, equal to `nu + 1`. |

#### Likelihood and diagnostics

| Component | Description |
|---|---|
| `score` | Score vector at the estimate. |
| `max.abs.score` | Largest absolute score component. |
| `logLik` | Maximized log-likelihood. |
| `objective` | Minimized negative log-likelihood. |
| `AIC`, `BIC` | Information criteria. |
| `nu` | Fixed degrees of freedom. |
| `convergence`, `message`, `counts` | Optimizer diagnostics. |
| `method` | Optimization method. |
| `use.score` | Whether score-based optimization was requested. |
| `vcov.warning` | Expected-information inversion warning, or `NULL`. |
| `observed.information.warning` | Observed-information computation or inversion warning, or `NULL`. |
| `lower`, `upper` | Bounds used by bounded optimization. |
| `optimization.control` | Effective controls. |
| `optim` | Complete optimizer result. |
| `call` | Matched function call. |

### Student-t weights

For dimension `d = 2`, the stored weight is:

```text
weight_i =
  (nu + 2) /
  (nu + quadratic_form_i).
```

Access it with:

```r
fit_t$weights
```

### S3 methods

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
    "beta1"
  ),
  level = 0.90,
  type = "fisher"
)

fitted(fit_t)
fitted(
  fit_t,
  component = "Y"
)

residuals(fit_t)
residuals(
  fit_t,
  component = "X"
)

logLik(fit_t)
nobs(fit_t)
```

The `component` argument may be:

```text
"both"
"Y"
"X"
```

### Prediction

```r
predict(
  fit_t,
  newx = NULL,
  type = c(
    "structural",
    "marginal",
    "conditional_observed"
  ),
  tau_x = NULL
)
```

Prediction types:

| `type` | Output |
|---|---|
| `"structural"` | Structural regression line `beta0 + beta1 * newx`. |
| `"marginal"` | Marginal observed-`Y` location `beta0 + beta1 * mu_x`. |
| `"conditional_observed"` | Conditional location of observed `Y` given observed `X = newx`. |

For conditional predictions at new observed values:

```r
predict(
  fit_t,
  newx = X_new,
  tau_x = tau_x_new,
  type = "conditional_observed"
)
```

If `newx` has the original sample length and `tau_x` is omitted, the stored `object$tau.x` is used.

### Minimal L-BFGS-B call

```r
fit_t <- fit_lmve_t(
  Y = Y,
  X = X,
  tau_y = tau_y,
  tau_x = tau_x,
  start = NULL,
  nu = 5,
  level = 0.95,
  method = "L-BFGS-B",
  lower = c(
    beta0 = -Inf,
    beta1 = -Inf,
    mu_x = -Inf,
    sigma2_x = 1e-8,
    sigma2 = 1e-8
  ),
  upper = rep(
    Inf,
    5
  ),
  control = list(
    maxit = 2000,
    factr = 1e4,
    pgtol = 1e-12
  ),
  use_score = TRUE,
  compute_vcov = TRUE,
  hessian = TRUE
)
```

### Minimal BOBYQA call

```r
fit_t <- fit_lmve_t(
  Y = Y,
  X = X,
  tau_y = tau_y,
  tau_x = tau_x,
  start = NULL,
  nu = 5,
  method = "BOBYQA",
  control = list(
    maxeval = 100000,
    xtol_rel = 1e-8,
    ftol_rel = 1e-8
  ),
  use_score = FALSE,
  compute_vcov = TRUE,
  hessian = TRUE
)
```

## 6. Conditional location under the normal model

The normal MLqE and MDPDE prediction methods do not directly expose a `"conditional_observed"` type.

The fitted conditional mean of observed `Y` given observed `X = x` is:

```text
beta0 + beta1 * mu_x
+ beta1 * sigma2_x / (sigma2_x + tau_x)
  * (x - mu_x).
```

Use:

```r
conditional_prediction <- function(
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

## 7. Method comparison

| Feature | Corrected MLqE | MDPDE | Student-t MLE |
|---|---|---|---|
| Observed distribution | Bivariate normal | Bivariate normal | Bivariate Student-t |
| Robustness control | `q` | `q` | fixed `nu` |
| Automatic tuning | SQV | SQV | Not implemented |
| Consistency correction | Internal | Not required | Not applicable |
| `tau_y`, `tau_x` interpretation | Variances | Variances | Scale-squared quantities |
| `sigma2_x`, `sigma2` interpretation | Variances | Variances | Scale-squared quantities |
| Default covariance estimator | Corrected sandwich | Sandwich | Inverse expected Fisher information |
| Observed-information covariance | Not separately returned | Not separately returned | Available when `hessian = TRUE` |
| Normal MLE special case | `q = 1` | `q = 1` | Approached as `nu` increases, but `nu` is fixed by the user |

## 8. Common errors

### `Y` and `X` have different lengths

They must contain one value for each observed pair.

### Negative `tau_y` or `tau_x`

Known measurement-error variances or scale-squared quantities must be non-negative.

### Invalid starting vector

The starting vector must contain exactly the five named parameters:

```text
beta0
beta1
mu_x
sigma2_x
sigma2
```

The last two must be strictly positive.

### Automatic starting values fail

This can occur when:

```text
mad(X)^2 - median(tau_x)
```

is zero or non-finite. Supply `start` explicitly.

### Bounds do not contain the starting values

All starting values must lie inside the selected lower and upper bounds.

### SQV is used without covariance estimation

Use:

```r
compute_vcov = TRUE
```

when `q = "SQV"`.

### BOBYQA is unavailable

Install:

```r
install.packages("nloptr")
```

### Observed covariance is unavailable for Student-t MLE

Use:

```r
hessian = TRUE
```

and inspect:

```r
fit_t$observed.information.warning
```

if `vcov(fit_t, type = "observed")` is unavailable.

## 9. Recommended post-fit checks

```r
fit$convergence
fit$message
fit$counts
```

For corrected MLqE:

```r
fit_mlqe$q
fit_mlqe$score
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

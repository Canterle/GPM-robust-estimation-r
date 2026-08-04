---
output:
  pdf_document: default
  html_document: default
---
# Function Reference

This document describes the public nonlinear-regression fitting functions:

- `fit_nlm_mlqe()`;
- `fit_nlm_mdpde()`;
- `fit_nlm_t()`.

The descriptions follow the current implementations.

# 1. Common model specification

The three functions accept a model function of the form:

```r
model <- function(theta, data) {
  list(
    mu = ...,
    sigma2 = ...,
    D = ...,
    V = ...
  )
}
```

## Required model outputs

| Component | Meaning |
|---|---|
| `mu` | Conditional mean. A scalar is replicated; otherwise it must have length `length(y)`. |
| `sigma2` | Positive conditional variance for the normal functions, or positive squared scale for the Student-t function. A scalar is replicated. |

The model may return `sigma` instead of `sigma2`. The fitting function then sets `sigma2 = sigma^2`.

## Optional analytical derivatives

| Component | Definition | Required dimension |
|---|---|---|
| `D` | `d mu / d theta^T` | `length(y)` by `length(theta)` |
| `V` | `d sigma2 / d theta^T` | `length(y)` by `length(theta)` |

If `D` or `V` is omitted, only the missing matrix is computed numerically.

The names and order of the columns should agree with the parameter names in the starting-value vector.

## Example model

```r
nlm_model <- function(theta, data) {
  beta1 <- theta["beta1"]
  beta2 <- theta["beta2"]
  beta3 <- theta["beta3"]
  gamma0 <- theta["gamma0"]

  x <- data$x
  n <- length(x)

  exponential_term <- exp(-beta3 * x)

  mu <- beta1 +
    beta2 * exponential_term

  sigma2 <- rep(
    exp(gamma0),
    n
  )

  D <- cbind(
    beta1 = rep(1, n),
    beta2 = exponential_term,
    beta3 = -beta2 * x * exponential_term,
    gamma0 = rep(0, n)
  )

  V <- cbind(
    beta1 = rep(0, n),
    beta2 = rep(0, n),
    beta3 = rep(0, n),
    gamma0 = sigma2
  )

  list(
    mu = mu,
    sigma2 = sigma2,
    D = D,
    V = V
  )
}
```

# 2. `fit_nlm_mlqe()`

## Purpose

Fits a nonlinear normal model by maximum Lq-likelihood estimation and returns both:

- the original MLqE parameterization, `theta_star`;
- the Fisher-consistency-corrected parameterization, `theta`.

For `q = 1`, the procedure reduces to ordinary normal maximum likelihood.

## Usage

```r
fit_nlm_mlqe(
  y,
  model,
  start,
  data = NULL,
  q = 1,
  tau = NULL,
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
| `y` | Finite numeric response vector. |
| `model` | User-supplied model function returning `mu` and `sigma2` or `sigma`; it may also return `D` and `V`. |
| `start` | Finite numeric starting-value vector. Named parameters are strongly recommended. |
| `data` | Object passed unchanged to `model(theta, data)`. Usually a list containing covariates and design matrices. |
| `q` | Numeric value in `(0, 1]`, or `"SQV"` for automatic selection. |
| `tau` | Transformation family `tau(theta, r)`. Required for `q < 1` and for `"SQV"`. It must return a finite vector with the same length and parameter ordering as `theta`. |
| `level` | Confidence level in `(0, 1)`. Used for `fit$conf.int`. |
| `method` | Optimization method passed to `optim()`: `"BFGS"`, `"L-BFGS-B"`, or `"Nelder-Mead"`. |
| `lower`, `upper` | Parameter bounds used only with `"L-BFGS-B"`. Scalars are recycled to the parameter-vector length. |
| `control` | Control list passed to `optim()`. For BFGS and Nelder-Mead, defaults are `maxit = 1000` and `reltol = 1e-9`. For L-BFGS-B, defaults are `maxit = 1000`, `factr = 1e7`, and `pgtol = 0`. |
| `use_score` | If `TRUE`, supplies the analytical estimating-equation vector to gradient-based `optim()` methods. It is ignored for Nelder-Mead. |
| `compute_vcov` | If `TRUE`, computes `J`, `K`, sandwich covariance matrices, standard errors, and confidence intervals. It must be `TRUE` for `"SQV"`. |
| `jacobian_method` | Numerical differentiation method used by `numDeriv`: `"Richardson"` or `"simple"`. |
| `jacobian_method_args` | List passed as `method.args` to `numDeriv::jacobian()`. |
| `q_control` | Controls for automatic SQV selection. Ignored for fixed numeric `q`. |

## Consistency transformation

The user supplies:

```r
tau <- function(theta, r) {
  ...
}
```

The implementation uses:

```text
corrected estimate = tau(original estimate, q)
```

and verifies the inverse composition through `r = 1 / q`.

For the log-linear variance models in the examples:

```r
tau <- function(theta, r) {
  result <- theta
  result["gamma0"] <-
    result["gamma0"] -
    log(r)
  result
}
```

For `q = 1`, `tau` may be `NULL`; the identity transformation is used internally.

## SQV controls

The active defaults are:

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

| Entry | Description |
|---|---|
| `m0` | Number of points in the initial grid from `1` to `0.8`. Must be an integer at least 2. |
| `m` | Number of points in each refinement grid. Must be an integer at least 2 and cannot exceed `m0`. |
| `q_min` | Smallest permitted `q`; it must be positive and smaller than `0.8`. |
| `L` | Positive SQV stability threshold. |
| `verbose` | Prints selection progress when `TRUE`. |
| `keep_fits` | Stores the cached fixed-`q` fits in `fit$q.selection$fits` when `TRUE`. |

The initial lower grid endpoint is fixed at `0.8`. Other entries in `q_control`, including `initial_lower`, are removed by the current implementation.

## Main return components

The returned object has class `nlm_mlqe`. An SQV-selected object additionally has class `nlm_mlqe_sqv`.

### Estimates and inference

| Component | Description |
|---|---|
| `coefficients` | Corrected MLqE estimate. |
| `coefficients.star` | Original MLqE estimate before consistency correction. |
| `standard.error` | Standard errors for the corrected estimate. |
| `standard.error.star` | Standard errors for the original estimate. |
| `vcov` | Corrected sandwich covariance matrix. |
| `vcov.star` | Sandwich covariance matrix for the original estimate. |
| `vcov.direct` | Corrected covariance obtained by the direct Jacobian transformation. |
| `vcov.direct.difference` | Numerical difference used to compare the two corrected-covariance calculations. |
| `conf.int` | Wald confidence intervals for corrected parameters at the fitted `level`. |
| `conf.int.star` | Wald confidence intervals for original parameters. |

### Estimating-function and derivative quantities

| Component | Description |
|---|---|
| `J`, `J.star` | Sensitivity matrix in the original MLqE parameterization. `J` is an alias of `J.star`. |
| `K`, `K.star` | Variability matrix in the original MLqE parameterization. `K` is an alias of `K.star`. |
| `mean.jacobian.star` | Matrix `D` evaluated at the original estimate. |
| `variance.jacobian.star` | Matrix `V` evaluated at the original estimate. |
| `tau.jacobian` | Jacobian of `tau(theta, 1/q)` evaluated at the corrected estimate. |
| `tau.jacobian.inverse` | Inverse of `tau.jacobian`. |
| `tau.q.jacobian` | Jacobian of `tau(theta_star, q)` evaluated at the original estimate. |
| `tau.composition.residual` | Maximum scaled residual from checking the forward/inverse transformation composition. |
| `score`, `score.star` | Estimating-equation vector evaluated at the original estimate. |

### Fitted quantities

| Component | Description |
|---|---|
| `fitted.values` | Corrected fitted conditional means. |
| `sigma` | Corrected fitted conditional standard deviations. |
| `sigma2` | Corrected fitted conditional variances. |
| `residuals` | Corrected response residuals. |
| `fitted.values.star` | Fitted means under the original parameterization. |
| `sigma.star` | Fitted standard deviations under the original parameterization. |
| `sigma2.star` | Fitted variances under the original parameterization. |
| `residuals.star` | Response residuals under the original parameterization. |

### Optimization and metadata

| Component | Description |
|---|---|
| `objective` | Optimized MLqE objective value. |
| `convergence` | `optim()` convergence code. Zero indicates normal convergence. |
| `message` | Optimizer message, when available. |
| `counts` | Function and gradient evaluation counts returned by `optim()`. |
| `q` | Fixed or selected value of `q`. |
| `level` | Confidence level used at fitting. |
| `nobs` | Number of observations. |
| `method` | Optimization method. |
| `jacobian.method` | Numerical Jacobian method. |
| `model`, `data`, `tau` | Stored model function, data object, and transformation. |
| `vcov.warning` | Warning generated during covariance inversion, or `NULL`. |
| `optimization.control` | Effective optimization controls. |
| `optim` | Complete object returned by `optim()`. |
| `call` | Matched function call. |
| `implementation.version` | Implementation identifier. |

### Additional SQV components

| Component | Description |
|---|---|
| `selected.q` | Selected value of `q`. |
| `sqv.returned.fit.source` | Source of the returned cached/refitted object. |
| `q.selection$method` | Selection method, `"SQV"`. |
| `q.selection$selected.q` | Selected `q`. |
| `q.selection$reason` | Reason the algorithm stopped or selected the returned value. |
| `q.selection$pass` | Selection pass. |
| `q.selection$stage` | Selection stage. |
| `q.selection$configuration` | Effective SQV settings. |
| `q.selection$history` | Stage-level SQV comparison history. |
| `q.selection$grids` | Sequence of evaluated grids. |
| `q.selection$evaluations` | Fit and stability information for evaluated `q` values. |
| `q.selection$fits` | Cached fits when `keep_fits = TRUE`; otherwise empty. |

## S3 methods

```r
print(fit_mlqe)
summary(fit_mlqe)

coef(fit_mlqe)
coef(fit_mlqe, type = "star")

vcov(fit_mlqe)
vcov(fit_mlqe, type = "star")

confint(fit_mlqe)
confint(
  fit_mlqe,
  parm = c("beta1", "gamma0"),
  level = 0.90,
  type = "star"
)

fitted(fit_mlqe)
fitted(fit_mlqe, type = "star")

residuals(fit_mlqe)
residuals(fit_mlqe, type = "star")

predict(
  fit_mlqe,
  newdata = new_data,
  type = "mean",
  parameterization = "corrected"
)
```

Prediction `type` may be `"mean"`, `"sigma"`, or `"variance"`.

## Minimal call

```r
fit_mlqe <- fit_nlm_mlqe(
  y = y,
  model = nlm_model,
  start = start_values,
  data = model_data,
  q = "SQV",
  tau = tau,
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

# 3. `fit_nlm_mdpde()`

## Purpose

Fits a nonlinear normal model by maximum density power divergence estimation. For `q = 1`, the procedure reduces to ordinary normal maximum likelihood.

## Usage

```r
fit_nlm_mdpde(
  y,
  model,
  start,
  data = NULL,
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
  jacobian_method = c(
    "Richardson",
    "simple"
  ),
  jacobian_method_args = list(),
  q_control = list()
)
```

## Arguments

The arguments have the same meanings as in `fit_nlm_mlqe()`, except that there is no `tau` argument.

`q` may be numeric in `(0, 1]` or `"SQV"`.

The same `q_control` defaults and restrictions are used:

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

`compute_vcov = TRUE` is required for SQV.

The MDPDE function requires `numDeriv` only if `D` or `V` is missing from the model output.

## Main return components

The current implementation returns class `nlhm_mdpde`. The name contains `nlhm`, although the public fitting function is named `fit_nlm_mdpde()`. An SQV-selected object additionally has class `nlhm_mdpde_sqv`.

| Component | Description |
|---|---|
| `coefficients` | MDPDE estimate. |
| `standard.error` | Sandwich standard errors. |
| `vcov` | Sandwich covariance matrix `J^{-1} K J^{-T}`. |
| `conf.int` | Wald confidence intervals at the fitted `level`. |
| `J` | Sensitivity matrix. |
| `K` | Variability matrix. |
| `mean.jacobian` | Matrix `D` at the estimate. |
| `variance.jacobian` | Matrix `V` at the estimate. |
| `fitted.values` | Fitted conditional means. |
| `sigma` | Fitted conditional standard deviations. |
| `sigma2` | Fitted conditional variances. |
| `residuals` | Response residuals. |
| `objective` | Optimized MDPDE objective value. |
| `score` | Estimating-equation vector at the estimate. |
| `convergence` | `optim()` convergence code. |
| `message` | Optimizer message, when available. |
| `counts` | Function and gradient evaluation counts. |
| `q` | Fixed or selected value of `q`. |
| `level` | Confidence level used at fitting. |
| `nobs` | Number of observations. |
| `method` | Optimization method. |
| `jacobian.method` | Numerical Jacobian method. |
| `model`, `data` | Stored model function and data. |
| `vcov.warning` | Covariance-inversion warning, or `NULL`. |
| `optimization.control` | Effective optimization controls. |
| `optim` | Complete object returned by `optim()`. |
| `call` | Matched function call. |
| `implementation.version` | Implementation identifier. |

SQV fits additionally contain `selected.q`, `sqv.returned.fit.source`, and `q.selection` with the same structure described for MLqE.

## S3 methods

```r
print(fit_mdpde)
summary(fit_mdpde)

coef(fit_mdpde)
vcov(fit_mdpde)

confint(
  fit_mdpde,
  parm = c("beta1", "gamma0"),
  level = 0.90
)

fitted(fit_mdpde)
residuals(fit_mdpde)

predict(
  fit_mdpde,
  newdata = new_data,
  type = "mean"
)
```

Prediction `type` may be `"mean"`, `"sigma"`, or `"variance"`.

## Minimal call

```r
fit_mdpde <- fit_nlm_mdpde(
  y = y,
  model = nlm_model,
  start = start_values,
  data = model_data,
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

## Normal MLE through MDPDE

```r
fit_mle <- fit_nlm_mdpde(
  y = y,
  model = nlm_model,
  start = start_values,
  data = model_data,
  q = 1,
  method = "BFGS",
  use_score = TRUE,
  compute_vcov = TRUE
)
```

# 4. `fit_nlm_t()`

## Purpose

Fits an independent nonlinear Student-t regression model by maximum likelihood with known, fixed degrees of freedom.

The model is parameterized as:

```text
Y_i = mu_i + sqrt(sigma2_i) T_i,
T_i ~ t_nu.
```

`sigma2_i` is the squared scale, not the conditional variance.

When `nu > 2`:

```text
Var(Y_i | covariates) = nu * sigma2_i / (nu - 2).
```

## Usage

```r
fit_nlm_t(
  y,
  model,
  start_theta,
  data = NULL,
  nu,
  method = c(
    "BFGS",
    "L-BFGS-B",
    "Nelder-Mead",
    "CG"
  ),
  lower = -Inf,
  upper = Inf,
  control = list(),
  use_gradient = TRUE,
  hessian = TRUE,
  derivative_eps =
    .Machine$double.eps^(1 / 3),
  conf_level = 0.95
)
```

## Arguments

| Argument | Description |
|---|---|
| `y` | Finite nonempty numeric response vector. |
| `model` | Model function returning `mu` and `sigma2` or `sigma`; it may also return `D` and `V`. Here `sigma2` is squared scale. |
| `start_theta` | Finite nonempty starting-value vector. Named parameters are strongly recommended. |
| `data` | Object passed unchanged to `model(theta, data)`. |
| `nu` | Fixed positive Student-t degrees of freedom. It is not estimated. |
| `method` | Optimization method: `"BFGS"`, `"L-BFGS-B"`, `"Nelder-Mead"`, or `"CG"`. |
| `lower`, `upper` | Bounds of length one or `length(start_theta)`. Finite bounds require `method = "L-BFGS-B"`. |
| `control` | Control list passed to `optim()`. Defaults are `maxit = 2000`, `reltol = 1e-10` for unconstrained methods, and `maxit = 2000`, `factr = 1e7`, `pgtol = 0` for L-BFGS-B. |
| `use_gradient` | If `TRUE`, supplies the negative analytical score to BFGS, CG, or L-BFGS-B. At invalid trial points, the implementation falls back to a numerical gradient. |
| `hessian` | Requests the optimizer Hessian. When available, it is stored as observed information and inverted to obtain `vcov.observed`. |
| `derivative_eps` | Positive relative finite-difference step used when `D`, `V`, or a fallback gradient must be calculated numerically. |
| `conf_level` | Confidence level in `(0, 1)` used to create `fit$conf.int`. |

## Main return components

The returned object has class `nlm_t_fit`.

### Estimates and inference

| Component | Description |
|---|---|
| `coefficients` | Student-t MLE. |
| `standard.error` | Standard errors from the inverse expected Fisher information. |
| `conf.int` | Wald confidence intervals using `conf_level`. |
| `coefficient.table` | Estimates, standard errors, z statistics, and two-sided normal-approximation p-values. |
| `vcov` | Alias of `vcov.fisher`. |
| `vcov.fisher` | Inverse expected Fisher information. |
| `vcov.observed` | Inverse observed information when `hessian = TRUE` and inversion succeeds. |
| `fisher.information` | Expected Fisher information matrix. |
| `observed.information` | Optimizer Hessian when requested and available. |

### Likelihood, score, and diagnostics

| Component | Description |
|---|---|
| `score` | Score vector at the estimate. |
| `max.abs.score` | Largest absolute score component. |
| `logLik` | Maximized log-likelihood. |
| `objective` | Minimized negative log-likelihood returned by `optim()`. |
| `AIC` | Akaike information criterion. |
| `BIC` | Bayesian information criterion. |
| `u` | Squared scale residuals, `(y - mu)^2 / sigma2`. |
| `weights` | Student-t likelihood weights, `(nu + 1) / (nu + u)`. |

### Fitted quantities

| Component | Description |
|---|---|
| `fitted.values` | Fitted conditional locations/means. |
| `sigma2` | Fitted squared Student-t scales. |
| `scale` | Fitted Student-t scales, `sqrt(sigma2)`. |
| `residuals` | Response residuals. |
| `scale.residuals` | Residuals divided by `sqrt(sigma2)`. |
| `pearson.residuals` | Residuals divided by the fitted conditional standard deviation. Returned as `NA` when `nu <= 2`, because the variance does not exist. |
| `mean.jacobian` | Matrix `D` at the estimate. |
| `scale2.jacobian` | Matrix `V` at the estimate. |

### Optimization and metadata

| Component | Description |
|---|---|
| `nu` | Fixed degrees of freedom. |
| `nobs`, `n` | Number of observations. |
| `p` | Number of estimated parameters. |
| `convergence` | `optim()` convergence code. |
| `message` | Optimizer message, when available. |
| `counts` | Function and gradient evaluation counts. |
| `method` | Optimization method. |
| `model`, `data`, `y` | Stored model, data object, and response. |
| `optimization.control` | Effective optimizer controls. |
| `optim` | Complete object returned by `optim()`. |
| `call` | Matched function call. |

## S3 methods

```r
print(fit_t)
summary(fit_t)

coef(fit_t)

vcov(fit_t)
vcov(fit_t, type = "observed")

confint(
  fit_t,
  parm = c("beta1", "gamma0"),
  level = 0.90
)

fitted(fit_t)

residuals(fit_t)
residuals(fit_t, type = "scale")
residuals(fit_t, type = "pearson")

predict(
  fit_t,
  newdata = new_data,
  type = "mean"
)

logLik(fit_t)
nobs(fit_t)
```

Residual `type` may be:

- `"response"`;
- `"scale"`;
- `"pearson"`.

Prediction `type` may be:

- `"mean"`;
- `"sigma"`;
- `"variance"`.

### Important prediction convention

In the current implementation:

```r
predict(fit_t, type = "variance")
```

returns the model's `sigma2` component, which is the squared Student-t scale. It does **not** multiply by `nu / (nu - 2)`.

For the actual conditional variance when `nu > 2`, use:

```r
conditional_variance <-
  predict(
    fit_t,
    newdata = new_data,
    type = "variance"
  ) *
  fit_t$nu /
  (
    fit_t$nu - 2
  )
```

## Minimal call

```r
fit_t <- fit_nlm_t(
  y = y,
  model = nlm_model,
  start_theta = start_values,
  data = model_data,
  nu = 4,
  method = "BFGS",
  control = list(
    maxit = 2000,
    reltol = 1e-12
  ),
  use_gradient = TRUE,
  hessian = TRUE,
  conf_level = 0.95
)
```

# 5. Computing robustness or likelihood weights

## MLqE and MDPDE density-power component

The normal-model fit objects do not store the density-power weights directly. Compute them from the fitted density:

```r
weight_mlqe <- dnorm(
  y,
  mean = fit_mlqe$fitted.values,
  sd = fit_mlqe$sigma
)^(1 - fit_mlqe$q)

weight_mdpde <- dnorm(
  y,
  mean = fit_mdpde$fitted.values,
  sd = fit_mdpde$sigma
)^(1 - fit_mdpde$q)
```

If `q = 1`, these values equal one.

## Student-t likelihood weights

The Student-t fit stores:

```r
fit_t$weights
```

where:

```text
u_i = (y_i - mu_i)^2 / sigma2_i

weight_i = (nu + 1) / (nu + u_i).
```

# 6. Method comparison

| Feature | Corrected MLqE | MDPDE | Student-t MLE |
|---|---|---|---|
| Distribution | Normal | Normal | Student-t |
| Robustness control | `q` | `q` | fixed `nu` |
| Automatic tuning | SQV for `q` | SQV for `q` | Not implemented |
| Consistency transformation | Required for `q < 1` | Not required | Not applicable |
| `sigma2` interpretation | Variance | Variance | Squared scale |
| Covariance estimator | Corrected sandwich | Sandwich | Inverse expected Fisher information by default |
| Observed Hessian covariance | No separate method | No separate method | Available with `hessian = TRUE` |
| Normal MLE special case | `q = 1` | `q = 1` | Approached as `nu` becomes large, but `nu` remains fixed by the user |

# 7. Recommended checks after fitting

```r
fit$convergence
fit$message
fit$counts
```

For MLqE and MDPDE:

```r
max(
  abs(
    fit$score
  ),
  na.rm = TRUE
)

fit$vcov.warning
```

For Student-t MLE:

```r
fit$max.abs.score
fit$fisher.information
fit$observed.information
```

For SQV fits:

```r
fit$q
fit$q.selection$reason
fit$q.selection$history
fit$q.selection$evaluations
```

# 8. Common errors

## Model returns an invalid dispersion

All values of `sigma2` must be finite and strictly positive.

A log-link is often convenient:

```r
sigma2 <- exp(linear_predictor)
```

## Derivative dimensions do not match

Both `D` and `V` must have:

```r
nrow = length(y)
ncol = length(start)
```

## Parameter names are inconsistent

Use the same names in `start`, `model()`, `D`, `V`, and `tau()`.

## Bounds are used with the wrong optimizer

Finite bounds require `"L-BFGS-B"` for `fit_nlm_t()`. For the normal implementations, bounds are only passed to `optim()` when `"L-BFGS-B"` is selected.

## SQV fails because covariance was disabled

Use:

```r
compute_vcov = TRUE
```

for `q = "SQV"`.

## Corrected MLqE transformation is not invertible

The implementation checks approximately:

```text
tau(tau(theta, q), 1/q) = theta.
```

Inspect:

```r
fit_mlqe$tau.composition.residual
```

and verify the definition of `tau()` if a warning is produced.

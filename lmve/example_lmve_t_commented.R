# ============================================================
# Complete simulated example for fit_lmve_t()
#
# Linear Student-t measurement-error model
#
# Latent structural model:
#
#   Y_i^* = beta0 + beta1 * X_i^* + epsilon_i.
#
# Observed variables:
#
#   Y_i = Y_i^* + e_yi,
#   X_i = X_i^* + e_xi.
#
# For each observation i, the joint latent vector satisfies
#
#   (X_i^*, epsilon_i, e_yi, e_xi)^T
#     ~ t_nu(
#         (mu_x, 0, 0, 0)^T,
#         diag(sigma2_x, sigma2, tau_yi, tau_xi)
#       ),
#
# independently across observations. All dispersion quantities
# are components of the Student-t scale matrix, not covariance
# quantities.
#
# Consequently,
#
#   (Y_i, X_i)^T
#     ~ t_2(mu(theta), Sigma_i(theta), nu),
#
# where
#
#   mu(theta) =
#     (beta0 + beta1 * mu_x, mu_x)^T
#
# and
#
#   Sigma_i(theta) =
#
#     [ beta1^2 * sigma2_x + sigma2 + tau_yi,
#       beta1 * sigma2_x                              ]
#     [ beta1 * sigma2_x,
#       sigma2_x + tau_xi                            ].
#
# The script illustrates how to:
#
#   1. simulate an uncontaminated bivariate Student-t sample;
#   2. calculate robust starting values internally;
#   3. fit the model by Student-t maximum likelihood;
#   4. inspect estimates, information matrices, fitted
#      distributions, conditional predictions, MAEs, likelihood
#      criteria, structural lines, and score weights.
# ============================================================

# Load the current Student-t maximum-likelihood implementation.
source("fit_lmve_t.R")

# Uncomment the following line to reproduce exactly the same
# simulated sample every time the script is run.
# set.seed(20260722)


# ============================================================
# 1. Simulation settings
# ============================================================

# Number of observations in the simulated sample.
n <- 5000

# Fixed degrees of freedom of the multivariate Student-t
# distribution.
nu <- 5


# ============================================================
# 2. True parameter values
# ============================================================

# beta0: intercept of the structural regression model;
# beta1: slope of the structural regression model;
# mu_x: location of the latent explanatory variable;
# sigma2_x: scale-squared parameter of the latent explanatory
#   variable;
# sigma2: scale-squared parameter of the structural equation
#   error.
true_theta <- c(
  beta0 = 2,
  beta1 = 1,
  mu_x = -2,
  sigma2_x = 5,
  sigma2 = 10
)


# ============================================================
# 3. Known measurement-error scale quantities
# ============================================================

# Generate the known measurement-error scale-squared quantities
# for the observed response Y.
tau_y <- runif(
  n,
  min = 0.5,
  max = 4
)^2

# Generate the known measurement-error scale-squared quantities
# for the observed covariate X.
tau_x <- runif(
  n,
  min = 0.5,
  max = 1.5
)^2


# ============================================================
# 4. Observed-data generation
# ============================================================

# Marginal location of the observed pair (Y_i, X_i).
true_mean <- c(
  Y =
    true_theta["beta0"] +
    true_theta["beta1"] *
    true_theta["mu_x"],
  X =
    true_theta["mu_x"]
)

# Generate each observed pair directly from its bivariate
# Student-t distribution. The scale matrix varies across
# observations because tau_y and tau_x are observation-specific.
observed_data <- t(
  vapply(
    seq_len(n),
    function(i) {
      Sigma_i <- matrix(
        c(
          true_theta["beta1"]^2 *
            true_theta["sigma2_x"] +
            true_theta["sigma2"] +
            tau_y[i],

          true_theta["beta1"] *
            true_theta["sigma2_x"],

          true_theta["beta1"] *
            true_theta["sigma2_x"],

          true_theta["sigma2_x"] +
            tau_x[i]
        ),
        nrow = 2
      )

      as.numeric(
        mvtnorm::rmvt(
          n = 1,
          sigma = Sigma_i,
          df = nu,
          delta = true_mean,
          type = "shifted"
        )
      )
    },
    numeric(2)
  )
)

colnames(observed_data) <- c(
  "Y",
  "X"
)

Y <- observed_data[, "Y"]
X <- observed_data[, "X"]


# ============================================================
# 5. Auxiliary functions
# ============================================================

# Conditional Student-t distribution of observed Y given
# observed X for a supplied parameter vector.
conditional_distribution <- function(
    theta,
    X_new,
    tau_y_new,
    tau_x_new,
    nu
) {
  beta0 <- unname(
    theta["beta0"]
  )

  beta1 <- unname(
    theta["beta1"]
  )

  mu_x <- unname(
    theta["mu_x"]
  )

  sigma2_x <- unname(
    theta["sigma2_x"]
  )

  sigma2 <- unname(
    theta["sigma2"]
  )

  mu_y <- beta0 +
    beta1 * mu_x

  scale_yy <- beta1^2 *
    sigma2_x +
    sigma2 +
    tau_y_new

  scale_yx <- beta1 *
    sigma2_x

  scale_xx <- sigma2_x +
    tau_x_new

  quadratic_x <- (
    X_new -
      mu_x
  )^2 /
    scale_xx

  conditional_location <- mu_y +
    scale_yx /
    scale_xx *
    (
      X_new -
        mu_x
    )

  conditional_scale2 <- (
    nu +
      quadratic_x
  ) /
    (
      nu +
        1
    ) *
    (
      scale_yy -
        scale_yx^2 /
        scale_xx
    )

  list(
    location = conditional_location,
    scale2 = conditional_scale2,
    df = rep(
      nu + 1,
      length(X_new)
    )
  )
}

# Structural mean relating the latent variables.
structural_mean <- function(
    theta,
    latent_x
) {
  unname(
    theta["beta0"] +
      theta["beta1"] *
      latent_x
  )
}


# ============================================================
# 6. Starting values
# ============================================================

# Starting values are not calculated explicitly in this example.
#
# With start = NULL, fit_lmve_t() calculates robust moment-based
# starting values internally from Y, X, tau_y, tau_x, and nu.


# ============================================================
# 7. Optimization settings
# ============================================================

# The scale parameters are estimated on their original scales.
# L-BFGS-B is used to keep them strictly positive.
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

optimization_control <- list(
  maxit = 2000,
  factr = 1e4,
  pgtol = 1e-12
)

# Confidence level used for standard errors and intervals.
conf_level <- 0.95


# ============================================================
# 8. Student-t maximum likelihood fit
# ============================================================

# The same fixed degrees of freedom are used in data generation
# and estimation.
fit_t <- fit_lmve_t(
  Y = Y,
  X = X,
  tau_y = tau_y,
  tau_x = tau_x,
  start = NULL,
  nu = nu,
  level = conf_level,
  method = "L-BFGS-B",
  lower = lower_bounds,
  upper = upper_bounds,
  control = optimization_control,
  use_score = TRUE,
  compute_vcov = TRUE,
  hessian = TRUE
)


# ============================================================
# 9. Fit summary
# ============================================================

cat(
  "\n",
  paste(rep("=", 70), collapse = ""),
  "\nSTUDENT-t MAXIMUM LIKELIHOOD FIT\n",
  paste(rep("=", 70), collapse = ""),
  "\n",
  sep = ""
)

summary(fit_t)


# ============================================================
# 10. Parameter estimates
# ============================================================

parameter_estimates <- cbind(
  True = true_theta,
  `Student-t MLE` =
    fit_t$coefficients
)

print(parameter_estimates)


# ============================================================
# 11. Standard errors and confidence intervals
# ============================================================

standard_error_results <- cbind(
  `Student-t MLE estimate` =
    fit_t$coefficients,
  `Student-t MLE SE` =
    fit_t$standard.error
)

print(standard_error_results)

fit_t$conf.int


# ============================================================
# 12. Covariance and information matrices
# ============================================================

# Covariance matrix based on the expected Fisher information.
fit_t$vcov

# Covariance matrix based on the observed information.
fit_t$vcov.observed

# Expected Fisher information matrix.
fit_t$fisher.information

# Observed information matrix.
fit_t$observed.information


# ============================================================
# 13. Analytical derivative matrices
# ============================================================

# D = d mu(theta) / d theta^T.
fit_t$mean.jacobian

# V = d vec(Sigma_i(theta)) / d theta^T.
fit_t$scale.jacobian

# F = (D^T, V^T)^T.
fit_t$F


# ============================================================
# 14. Fitted distributions and residuals
# ============================================================

fitted_results <- data.frame(
  index = seq_len(n),
  Y = Y,
  X = X,
  tau_y = tau_y,
  tau_x = tau_x,
  fitted_y = fit_t$fitted.y,
  fitted_x = fit_t$fitted.x,
  residual_y = fit_t$residual.y,
  residual_x = fit_t$residual.x,
  quadratic_form = fit_t$quadratic.forms,
  log_density = fit_t$log.density
)

print(
  head(fitted_results)
)

# The third array index identifies the observation.
fit_t$scale.array[, , 1]

# For nu > 2, the covariance matrices are obtained as
#
#   Cov(Y_i, X_i) = nu / (nu - 2) * Sigma_i.
fit_t$covariance.array[, , 1]

# Estimated latent-X and equation-error scale quantities.
fit_t$sigma2.x.scale
fit_t$sigma2.scale

# Corresponding variances when nu > 2.
fit_t$latent.x.variance
fit_t$equation.error.variance


# ============================================================
# 15. Conditional predictions and mean absolute errors
# ============================================================

# Fitted conditional Student-t distribution at the observed X
# values.
conditional_fitted <- data.frame(
  index = seq_len(n),
  X = X,
  fitted_location =
    fit_t$conditional.observed.location.y.given.x,
  fitted_scale2 =
    fit_t$conditional.observed.scale2.y.given.x,
  fitted_df =
    fit_t$conditional.observed.df.y.given.x
)

print(
  head(conditional_fitted)
)

# True conditional distribution at the observed X values.
true_conditional <- conditional_distribution(
  theta = true_theta,
  X_new = X,
  tau_y_new = tau_y,
  tau_x_new = tau_x,
  nu = nu
)

mae_observed <- mean(
  abs(
    Y -
      fit_t$conditional.observed.location.y.given.x
  )
)

mae_true_location <- mean(
  abs(
    true_conditional$location -
      fit_t$conditional.observed.location.y.given.x
  )
)

mae_results <- c(
  `MAE relative to observed Y` =
    mae_observed,
  `MAE relative to true conditional location` =
    mae_true_location
)

print(mae_results)


# ============================================================
# 16. Likelihood criteria
# ============================================================

log_likelihood <- sum(
  fit_t$log.density
)

number_of_parameters <- length(
  fit_t$coefficients
)

likelihood_criteria <- c(
  logLik = log_likelihood,
  AIC =
    -2 *
    log_likelihood +
    2 *
    number_of_parameters,
  BIC =
    -2 *
    log_likelihood +
    log(n) *
    number_of_parameters
)

print(likelihood_criteria)


# ============================================================
# 17. Predictions for new data
# ============================================================

# Evaluate the conditional Student-t distribution on a regular
# grid. The new measurement-error scale quantities are fixed at
# their sample medians.
X_new <- seq(
  min(X),
  max(X),
  length.out = 200
)

tau_y_new <- rep(
  median(tau_y),
  length(X_new)
)

tau_x_new <- rep(
  median(tau_x),
  length(X_new)
)

true_prediction <- conditional_distribution(
  theta = true_theta,
  X_new = X_new,
  tau_y_new = tau_y_new,
  tau_x_new = tau_x_new,
  nu = nu
)

fitted_prediction <- conditional_distribution(
  theta = fit_t$coefficients,
  X_new = X_new,
  tau_y_new = tau_y_new,
  tau_x_new = tau_x_new,
  nu = nu
)

prediction_results <- data.frame(
  X = X_new,
  tau_y = tau_y_new,
  tau_x = tau_x_new,
  true_location =
    true_prediction$location,
  true_scale2 =
    true_prediction$scale2,
  true_df =
    true_prediction$df,
  fitted_location =
    fitted_prediction$location,
  fitted_scale2 =
    fitted_prediction$scale2,
  fitted_df =
    fitted_prediction$df
)

print(
  head(prediction_results)
)


# ============================================================
# 18. Student-t score weights
# ============================================================

# For dimension d = 2, the multiplier of the score contribution
# of observation i is
#
#   w_i = (nu + 2) / (nu + u_i),
#
# where u_i is the fitted Mahalanobis squared distance.
score_weights <- data.frame(
  index = seq_len(n),
  weight = fit_t$weights,
  quadratic_form = fit_t$quadratic.forms
)

print(score_weights)


# ============================================================
# 19. True and fitted structural regression lines
# ============================================================

latent_x_new <- seq(
  min(X),
  max(X),
  length.out = 200
)

structural_mean_true <- structural_mean(
  true_theta,
  latent_x_new
)

structural_mean_t <- structural_mean(
  fit_t$coefficients,
  latent_x_new
)

plot(
  X,
  Y,
  pch = 16,
  cex = 0.6,
  xlab = "Observed X",
  ylab = "Observed Y",
  main = "True and fitted structural regression lines"
)

lines(
  latent_x_new,
  structural_mean_true,
  lty = 2,
  lwd = 2
)

lines(
  latent_x_new,
  structural_mean_t,
  lty = 1,
  lwd = 2
)

legend(
  "topleft",
  legend = c(
    "True structural mean",
    paste0(
      "Student-t MLE (nu = ",
      nu,
      ")"
    )
  ),
  lty = c(
    2,
    1
  ),
  lwd = 2,
  bty = "n"
)


# ============================================================
# 20. Student-t score weights
# ============================================================

plot(
  score_weights$index,
  score_weights$weight,
  pch = 1,
  xlab = "Index",
  ylab = "Student-t score weight",
  main = "Student-t score weights"
)

abline(
  h = 1,
  lty = 2
)

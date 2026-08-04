# ============================================================
# Complete simulated example for fit_lmve_mlqe() and
# fit_lmve_mdpde()
#
# Linear normal measurement-error model
#
# Latent structural model:
#
#   Y_i^* = beta0 + beta1 * X_i^* + epsilon_i,
#
# where
#
#   X_i^* ~ Normal(mu_x, sigma2_x),
#   epsilon_i ~ Normal(0, sigma2).
#
# Observed variables:
#
#   Y_i = Y_i^* + e_yi,
#   X_i = X_i^* + e_xi,
#
# with known observation-specific measurement-error variances
#
#   Var(e_yi) = tau_yi,
#   Var(e_xi) = tau_xi.
#
# Consequently,
#
#   (Y_i, X_i)^T ~ Normal_2(mu(theta), Sigma_i(theta)),
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
#   1. simulate a heteroscedastic measurement-error dataset;
#   2. contaminate one observed pair;
#   3. calculate robust starting values internally;
#   4. fit the corrected MLqE, MDPDE, and normal MLE;
#   5. compare estimates, fitted distributions, predictions,
#      structural regression lines, and score weights.
# ============================================================

# Load the current implementations of the two estimators.
source("fit_lmve_mlqe.R")
source("fit_lmve_mdpde.R")

# Uncomment the following line to reproduce exactly the same
# simulated sample every time the script is run.
# set.seed(20260722)


# ============================================================
# 1. Simulation settings
# ============================================================

# Number of observations in the simulated sample.
n <- 50


# ============================================================
# 2. True parameter values
# ============================================================

# beta0: intercept of the structural regression model;
# beta1: slope of the structural regression model;
# mu_x: mean of the latent explanatory variable;
# sigma2_x: variance of the latent explanatory variable;
# sigma2: variance of the structural regression error.
true_theta <- c(
  beta0 = 2,
  beta1 = 1,
  mu_x = -2,
  sigma2_x = 5,
  sigma2 = 10
)


# ============================================================
# 3. Known measurement-error variances
# ============================================================

# Generate the known measurement-error variances for Y.
tau_y <- runif(
  n,
  min = 0.5,
  max = 4
)^2

# Generate the known measurement-error variances for X.
tau_x <- runif(
  n,
  min = 0.5,
  max = 1.5
)^2


# ============================================================
# 4. Observed-data generation
# ============================================================

# Marginal mean of the observed pair (Y_i, X_i).
true_mean <- c(
  Y =
    true_theta["beta0"] +
    true_theta["beta1"] *
    true_theta["mu_x"],
  X =
    true_theta["mu_x"]
)

# Generate each observed pair directly from its bivariate normal
# distribution. The covariance matrix varies across observations
# because tau_y and tau_x are observation-specific.
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
        mvtnorm::rmvnorm(
          n = 1,
          mean = true_mean,
          sigma = Sigma_i
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

# Replace one observed pair by atypical values to illustrate the
# behavior of the robust estimators.
contaminated_index <- 9L

observed_data[
  contaminated_index,
] <- c(
  35,
  -20
)

Y <- observed_data[, "Y"]
X <- observed_data[, "X"]


# ============================================================
# 5. Auxiliary functions
# ============================================================

# Conditional mean of observed Y given observed X under the
# fitted bivariate normal model.
conditional_prediction <- function(
    theta,
    X_new,
    tau_x_new
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
# With start = NULL, fit_lmve_mlqe() and fit_lmve_mdpde()
# calculate the robust moment-based starting values internally
# from Y, X, tau_y, and tau_x.


# ============================================================
# 7. Common optimization settings
# ============================================================

# The variance parameters are estimated on their original
# scales. L-BFGS-B is used to keep them strictly positive.
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

# Configuration of the automatic SQV procedure.
sqv_control <- list(
  # Number of values in the initial grid.
  m0 = 21,

  # Number of q values in each subsequent refinement grid.
  m = 3,

  # Smallest q value that may be evaluated.
  q_min = 0.5,

  # Stability threshold used to compare adjacent estimates.
  L = 0.01,

  # Retain all valid evaluated fits.
  keep_fits = TRUE
)


# ============================================================
# 8. Corrected MLqE fit with q selected by SQV
# ============================================================

fit_mlqe <- fit_lmve_mlqe(
  Y = Y,
  X = X,
  tau_y = tau_y,
  tau_x = tau_x,
  start = NULL,
  q = "SQV",
  level = conf_level,
  method = "L-BFGS-B",
  lower = lower_bounds,
  upper = upper_bounds,
  control = optimization_control,
  use_score = TRUE,
  compute_vcov = TRUE,
  q_control = sqv_control
)


# ============================================================
# 9. MDPDE fit with q selected by SQV
# ============================================================

# The same internal starting-value rule, bounds, optimization
# controls, and SQV configuration are used for a direct
# comparison with MLqE.
fit_mdpde <- fit_lmve_mdpde(
  Y = Y,
  X = X,
  tau_y = tau_y,
  tau_x = tau_x,
  start = NULL,
  q = "SQV",
  level = conf_level,
  method = "L-BFGS-B",
  lower = lower_bounds,
  upper = upper_bounds,
  control = optimization_control,
  use_score = TRUE,
  compute_vcov = TRUE,
  q_control = sqv_control
)


# ============================================================
# 10. Normal maximum likelihood fit
# ============================================================

# For q = 1, the MDPDE reduces to maximum likelihood under the
# normal measurement-error model.
fit_mle <- fit_lmve_mdpde(
  Y = Y,
  X = X,
  tau_y = tau_y,
  tau_x = tau_x,
  start = NULL,
  q = 1,
  level = conf_level,
  method = "L-BFGS-B",
  lower = lower_bounds,
  upper = upper_bounds,
  control = optimization_control,
  use_score = TRUE,
  compute_vcov = TRUE
)


# ============================================================
# 11. Fit summaries
# ============================================================

cat(
  "\n",
  paste(rep("=", 70), collapse = ""),
  "\nCORRECTED MLqE FIT\n",
  paste(rep("=", 70), collapse = ""),
  "\n",
  sep = ""
)

summary(fit_mlqe)

cat(
  "\n",
  paste(rep("=", 70), collapse = ""),
  "\nMDPDE FIT\n",
  paste(rep("=", 70), collapse = ""),
  "\n",
  sep = ""
)

summary(fit_mdpde)

cat(
  "\n",
  paste(rep("=", 70), collapse = ""),
  "\nNORMAL MAXIMUM LIKELIHOOD FIT\n",
  paste(rep("=", 70), collapse = ""),
  "\n",
  sep = ""
)

summary(fit_mle)


# ============================================================
# 12. Selected q values and SQV information
# ============================================================

selected_q <- c(
  MLqE = fit_mlqe$q,
  MDPDE = fit_mdpde$q
)

print(selected_q)

fit_mlqe$q.selection$history
fit_mdpde$q.selection$history

fit_mlqe$q.selection$evaluations
fit_mdpde$q.selection$evaluations

names(
  fit_mlqe$q.selection$fits
)

names(
  fit_mdpde$q.selection$fits
)


# ============================================================
# 13. Parameter estimates
# ============================================================

parameter_estimates <- cbind(
  True = true_theta,
  `Corrected MLqE` =
    fit_mlqe$coefficients,
  MDPDE =
    fit_mdpde$coefficients,
  MLE =
    fit_mle$coefficients
)

print(parameter_estimates)

mlqe_original_and_corrected <- cbind(
  `Original MLqE` =
    fit_mlqe$coefficients.star,
  `Corrected MLqE` =
    fit_mlqe$coefficients
)

print(mlqe_original_and_corrected)


# ============================================================
# 14. Standard errors and confidence intervals
# ============================================================

standard_error_results <- cbind(
  `MLqE estimate` =
    fit_mlqe$coefficients,
  `MLqE SE` =
    fit_mlqe$standard.error,
  `MDPDE estimate` =
    fit_mdpde$coefficients,
  `MDPDE SE` =
    fit_mdpde$standard.error,
  `MLE estimate` =
    fit_mle$coefficients,
  `MLE SE` =
    fit_mle$standard.error
)

print(standard_error_results)

fit_mlqe$conf.int
fit_mdpde$conf.int
fit_mle$conf.int


# ============================================================
# 15. Covariance and estimating-function matrices
# ============================================================

fit_mlqe$vcov
fit_mlqe$vcov.star
fit_mdpde$vcov
fit_mle$vcov

fit_mlqe$J.star
fit_mlqe$K.star

fit_mdpde$J
fit_mdpde$K


# ============================================================
# 16. Analytical derivative matrices
# ============================================================

fit_mlqe$mean.jacobian.star
fit_mlqe$covariance.jacobian.star
fit_mlqe$F.star

fit_mdpde$mean.jacobian
fit_mdpde$covariance.jacobian
fit_mdpde$F


# ============================================================
# 17. Fitted distributions and residuals
# ============================================================

fitted_results <- data.frame(
  index = seq_len(n),
  Y = Y,
  X = X,
  tau_y = tau_y,
  tau_x = tau_x,
  mlqe_fitted_y =
    fit_mlqe$fitted.y,
  mlqe_fitted_x =
    fit_mlqe$fitted.x,
  mlqe_residual_y =
    fit_mlqe$residual.y,
  mlqe_residual_x =
    fit_mlqe$residual.x,
  mlqe_log_density =
    fit_mlqe$log.density,
  mdpde_fitted_y =
    fit_mdpde$fitted.y,
  mdpde_fitted_x =
    fit_mdpde$fitted.x,
  mdpde_residual_y =
    fit_mdpde$residual.y,
  mdpde_residual_x =
    fit_mdpde$residual.x,
  mdpde_log_density =
    fit_mdpde$log.density,
  mle_fitted_y =
    fit_mle$fitted.y,
  mle_fitted_x =
    fit_mle$fitted.x,
  mle_residual_y =
    fit_mle$residual.y,
  mle_residual_x =
    fit_mle$residual.x,
  mle_log_density =
    fit_mle$log.density
)

print(
  head(fitted_results)
)


# ============================================================
# 18. Predictions for observed and new covariate values
# ============================================================

# Conditional fitted means at the observed X values.
conditional_fitted_results <- data.frame(
  index = seq_len(n),
  X = X,
  tau_x = tau_x,
  true_conditional_mean =
    conditional_prediction(
      theta = true_theta,
      X_new = X,
      tau_x_new = tau_x
    ),
  mlqe_conditional_mean =
    conditional_prediction(
      theta = fit_mlqe$coefficients,
      X_new = X,
      tau_x_new = tau_x
    ),
  mdpde_conditional_mean =
    conditional_prediction(
      theta = fit_mdpde$coefficients,
      X_new = X,
      tau_x_new = tau_x
    ),
  mle_conditional_mean =
    conditional_prediction(
      theta = fit_mle$coefficients,
      X_new = X,
      tau_x_new = tau_x
    )
)

print(
  head(conditional_fitted_results)
)

# Prediction at a regular grid of new observed X values. The
# measurement-error variance for the new observations is fixed
# at the median of the simulated tau_x values.
X_new <- seq(
  min(X),
  max(X),
  length.out = 200
)

tau_x_new <- rep(
  median(tau_x),
  length(X_new)
)

prediction_results <- data.frame(
  X = X_new,
  tau_x = tau_x_new,
  true_conditional_mean =
    conditional_prediction(
      theta = true_theta,
      X_new = X_new,
      tau_x_new = tau_x_new
    ),
  mlqe_conditional_mean =
    conditional_prediction(
      theta = fit_mlqe$coefficients,
      X_new = X_new,
      tau_x_new = tau_x_new
    ),
  mdpde_conditional_mean =
    conditional_prediction(
      theta = fit_mdpde$coefficients,
      X_new = X_new,
      tau_x_new = tau_x_new
    ),
  mle_conditional_mean =
    conditional_prediction(
      theta = fit_mle$coefficients,
      X_new = X_new,
      tau_x_new = tau_x_new
    )
)

print(
  head(prediction_results)
)


# ============================================================
# 19. Score weights
# ============================================================

# The MLqE estimating function is evaluated in the original
# parameterization, so its density-power component uses the
# fitted original density.
weight_mlqe <- if (
  fit_mlqe$q == 1
) {
  rep(
    1,
    n
  )
} else {
  fit_mlqe$density.star^(
    1 -
      fit_mlqe$q
  )
}

# The MDPDE density-power component uses the fitted corrected
# bivariate normal density.
weight_mdpde <- if (
  fit_mdpde$q == 1
) {
  rep(
    1,
    n
  )
} else {
  fit_mdpde$density^(
    1 -
      fit_mdpde$q
  )
}

score_weights <- data.frame(
  index = seq_len(n),
  MLqE = weight_mlqe,
  MDPDE = weight_mdpde
)

print(score_weights)


# ============================================================
# 20. Graphical comparison of structural regression lines
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

structural_mean_mlqe <- structural_mean(
  fit_mlqe$coefficients,
  latent_x_new
)

structural_mean_mdpde <- structural_mean(
  fit_mdpde$coefficients,
  latent_x_new
)

structural_mean_mle <- structural_mean(
  fit_mle$coefficients,
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

points(
  X[contaminated_index],
  Y[contaminated_index],
  pch = 4,
  cex = 1.2,
  lwd = 2
)

lines(
  latent_x_new,
  structural_mean_true,
  lty = 2,
  lwd = 2
)

lines(
  latent_x_new,
  structural_mean_mlqe,
  lty = 1,
  lwd = 2
)

lines(
  latent_x_new,
  structural_mean_mdpde,
  lty = 3,
  lwd = 2
)

lines(
  latent_x_new,
  structural_mean_mle,
  lty = 4,
  lwd = 2
)

legend(
  "topleft",
  legend = c(
    "True structural mean",
    paste0(
      "Corrected MLqE (q = ",
      format(
        fit_mlqe$q,
        digits = 3
      ),
      ")"
    ),
    paste0(
      "MDPDE (q = ",
      format(
        fit_mdpde$q,
        digits = 3
      ),
      ")"
    ),
    "Normal MLE",
    "Contaminated observation"
  ),
  lty = c(
    2,
    1,
    3,
    4,
    NA
  ),
  lwd = c(
    2,
    2,
    2,
    2,
    NA
  ),
  pch = c(
    NA,
    NA,
    NA,
    NA,
    4
  ),
  bty = "n"
)


# ============================================================
# 21. Graphical comparison of score weights
# ============================================================

weight_range <- range(
  c(
    score_weights$MDPDE,
    score_weights$MLqE
  ),
  finite = TRUE
)

plot(
  score_weights$index,
  score_weights$MDPDE,
  pch = 1,
  ylim = weight_range,
  xlab = "Index",
  ylab = "Score weight",
  main = "MLqE and MDPDE score weights"
)

points(
  score_weights$index,
  score_weights$MLqE,
  pch = 3
)

legend(
  "topright",
  legend = c(
    "MDPDE",
    "Corrected MLqE"
  ),
  pch = c(
    1,
    3
  ),
  bty = "n"
)


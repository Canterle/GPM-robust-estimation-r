# ============================================================
# Complete simulated example for fit_nlm_t()
#
# Heteroscedastic nonlinear Student-t regression model
#
# Conditional mean:
#
#   mu_i = beta1 + beta2 /
#          (1 + beta3 * X1_i^beta4)
#
# Student-t squared scale:
#
#   sigma2_i = exp(
#     gamma0 + gamma1 * w1_i + gamma2 * w2_i
#   )
#
# Data-generating model:
#
#   Y_i = mu_i + sqrt(sigma2_i) * T_i,
#   T_i ~ t_nu,
#
# with nu = 4 fixed degrees of freedom.
#
# In this parameterization, sigma2_i is the squared scale of the
# Student-t distribution. Since nu > 2, the conditional variance is
#
#   Var(Y_i | X1_i, w1_i, w2_i)
#     = nu * sigma2_i / (nu - 2).
#
# The script illustrates how to:
#
#   1. simulate an uncontaminated heteroscedastic nonlinear
#      Student-t dataset;
#   2. define the model and its analytical derivatives;
#   3. construct data-driven starting values;
#   4. fit the model by Student-t maximum likelihood;
#   5. inspect estimates, confidence intervals, fitted values,
#      residuals, likelihood weights, and fitted curves.
# ============================================================

# Load the current implementation of the Student-t estimator.
source("fit_nlm_t.R")

# Uncomment the following line to reproduce exactly the same
# simulated sample every time the script is run.
# set.seed(20260721)


# ============================================================
# 1. Simulation settings
# ============================================================

# Number of observations in the simulated sample.
n <- 50

# Fixed degrees of freedom of the Student-t distribution.
nu <- 4


# ============================================================
# 2. True parameter values
# ============================================================

# Parameters of the nonlinear conditional mean:
#
# beta1: lower asymptote of the mean curve;
# beta2: difference between the response level near X1 = 0
#        and the lower asymptote;
# beta3: controls the rate of decay;
# beta4: controls the shape of the decay.
true_beta <- c(
  beta1 = 50,
  beta2 = 500,
  beta3 = 0.5,
  beta4 = 2
)

# Parameters of the log-squared-scale model:
#
# gamma0: intercept;
# gamma1: effect of w1;
# gamma2: effect of w2.
true_gamma <- c(
  gamma0 = 5.3,
  gamma1 = -0.5,
  gamma2 = 0.4
)

# Complete true parameter vector. The names and order must agree
# with those used by nlm_model().
true_theta <- c(
  true_beta,
  true_gamma
)


# ============================================================
# 3. Covariates and response generation
# ============================================================

# Generate the covariate used in the conditional mean.
# Sorting is not required for estimation; it only facilitates
# plotting the true and fitted curves.
X1 <- sort(
  runif(
    n,
    min = 0,
    max = 15
  )
)

# Generate the covariates used in the squared-scale model.
w1 <- runif(
  n,
  min = 0,
  max = 1
)

w2 <- runif(
  n,
  min = 0,
  max = 1
)

# Design matrix of the log-squared-scale model. The first column
# is the intercept and the other columns are observed covariates.
W1 <- cbind(
  intercept = 1,
  w1 = w1,
  w2 = w2
)

# Compute the true conditional mean for every observation.
true_mu <- true_beta["beta1"] +
  true_beta["beta2"] /
  (
    1 +
      true_beta["beta3"] *
      X1^true_beta["beta4"]
  )

# Compute the true Student-t squared scale. Because W1 contains
# an intercept column, W1 %*% true_gamma evaluates
# gamma0 + gamma1 * w1 + gamma2 * w2.
true_sigma2 <- exp(
  drop(
    W1 %*% true_gamma
  )
)

# Student-t scale parameter.
true_scale <- sqrt(
  true_sigma2
)

# Conditional variance and standard deviation. For nu = 4,
# Var(Y_i | covariates) = 2 * sigma2_i.
true_variance <-
  nu *
  true_sigma2 /
  (
    nu - 2
  )

true_sd <- sqrt(
  true_variance
)

# Generate an uncontaminated response from the heteroscedastic
# Student-t model with four degrees of freedom.
y <- true_mu +
  true_scale *
  rt(
    n,
    df = nu
  )

# Store all covariates in a single object. This object is passed
# to nlm_model() during every objective, score, and derivative
# evaluation.
model_data <- list(
  X1 = X1,
  W1 = W1
)


# ============================================================
# 4. Model specification
# ============================================================

# The model function must have the interface:
#
#   nlm_model <- function(theta, data)
#
# theta is the complete named parameter vector. The function
# must return at least mu and sigma2, where sigma2 is the squared
# scale of the Student-t distribution. It may also return:
#
#   D = d mu / d theta^T;
#   V = d sigma2 / d theta^T.
#
# Providing D and V analytically avoids numerical differentiation.
nlm_model <- function(theta, data) {
  # Extract the mean parameters by name.
  beta1 <- theta["beta1"]
  beta2 <- theta["beta2"]
  beta3 <- theta["beta3"]
  beta4 <- theta["beta4"]

  # Extract the squared-scale parameters.
  gamma <- theta[
    c(
      "gamma0",
      "gamma1",
      "gamma2"
    )
  ]

  # Recover the covariates from the supplied data object.
  X1 <- data$X1
  W1 <- data$W1
  n <- length(X1)

  # Quantities reused in the mean and its derivatives.
  X1_beta4 <- X1^beta4
  denominator <- 1 + beta3 * X1_beta4

  # Conditional mean.
  mu <- beta1 +
    beta2 /
    denominator

  # Student-t squared scale.
  sigma2 <- exp(
    drop(
      W1 %*% gamma
    )
  )

  # Derivative of X1^beta4 with respect to beta4. The value is
  # set to zero at X1 = 0 to avoid evaluating log(0).
  X1_beta4_log <- ifelse(
    X1 > 0,
    X1_beta4 * log(X1),
    0
  )

  # Matrix D = d mu / d theta^T.
  # The mean does not depend on the gamma parameters, so the
  # corresponding columns are zero.
  D <- cbind(
    beta1 = rep(1, n),
    beta2 = 1 / denominator,
    beta3 = -beta2 *
      X1_beta4 /
      denominator^2,
    beta4 = -beta2 *
      beta3 *
      X1_beta4_log /
      denominator^2,
    gamma0 = rep(0, n),
    gamma1 = rep(0, n),
    gamma2 = rep(0, n)
  )

  # Since sigma2_i = exp(W1_i gamma), its derivative with
  # respect to gamma is sigma2_i * W1_i.
  V_gamma <- sweep(
    W1,
    MARGIN = 1,
    STATS = sigma2,
    FUN = "*"
  )

  colnames(V_gamma) <- c(
    "gamma0",
    "gamma1",
    "gamma2"
  )

  # Matrix V = d sigma2 / d theta^T.
  # The squared scale does not depend on the beta parameters, so
  # the corresponding columns are zero.
  V <- cbind(
    beta1 = rep(0, n),
    beta2 = rep(0, n),
    beta3 = rep(0, n),
    beta4 = rep(0, n),
    V_gamma
  )

  # Return all quantities required by fit_nlm_t().
  list(
    mu = mu,
    sigma2 = sigma2,
    D = D,
    V = V
  )
}


# ============================================================
# 5. Starting values
# ============================================================

# Starting values based on the suggestion of
# Thiede and Pagano (1979).
X1_round <- round(X1)

beta1_start <- median(
  y[X1_round == max(X1_round)]
)

beta2_start <- median(
  y[X1_round == min(X1_round)]
) - beta1_start

beta4_start <- 1

N1_bar <- median(
  y[X1_round == unique(X1_round)[2]]
)

N2_bar <- median(
  y[X1_round == unique(X1_round)[3]]
)

beta3_start <- (
  N1_bar -
    N2_bar
) / (
  beta2_start *
    (
      unique(X1_round)[3] -
        unique(X1_round)[2]
    )
)

# Preliminary fitted means.
mu_start <- beta1_start +
  beta2_start /
  (
    1 +
      beta3_start *
      X1^beta4_start
  )

# Initial squared scale based on the robust residual scale.
sigma2_start <- mad(
  y -
    mu_start
)^2

# The squared-scale coefficients other than the intercept
# initially equal zero.
start_values <- c(
  beta1 = beta1_start,
  beta2 = beta2_start,
  beta3 = beta3_start,
  beta4 = beta4_start,
  gamma0 = log(sigma2_start),
  gamma1 = 0,
  gamma2 = 0
)


# ============================================================
# 6. Optimization settings
# ============================================================

# Controls passed to optim() through fit_nlm_t().
optimization_control <- list(
  maxit = 2000,
  reltol = 1e-12
)

# Confidence level used for standard errors and intervals.
conf_level <- 0.95


# ============================================================
# 7. Student-t maximum likelihood fit
# ============================================================

# The same fixed degrees of freedom are used in data generation
# and estimation.
fit_t <- fit_nlm_t(
  y = y,
  model = nlm_model,
  start_theta = start_values,
  data = model_data,
  nu = nu,
  method = "BFGS",
  control = optimization_control,
  use_gradient = TRUE,
  hessian = TRUE,
  conf_level = conf_level
)


# ============================================================
# 8. Fit summary
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
# 9. Parameter estimates
# ============================================================

parameter_estimates <- cbind(
  True = true_theta,
  `Student-t MLE` = fit_t$coefficients
)

print(
  parameter_estimates
)


# ============================================================
# 10. Standard errors and confidence intervals
# ============================================================

standard_error_results <- cbind(
  Estimate = fit_t$coefficients,
  `Standard error` = fit_t$standard.error
)

print(
  standard_error_results
)

fit_t$conf.int


# ============================================================
# 11. Covariance and information matrices
# ============================================================

# By default, fit_t$vcov equals the covariance matrix based on
# the expected Fisher information.
fit_t$vcov
fit_t$vcov.fisher
fit_t$vcov.observed

fit_t$fisher.information
fit_t$observed.information


# ============================================================
# 12. Analytical model derivatives
# ============================================================

fit_t$mean.jacobian
fit_t$scale2.jacobian


# ============================================================
# 13. Fitted values, scales, variances, and residuals
# ============================================================

# fit_t$sigma2 contains the fitted squared scales. The fitted
# conditional variances are nu * sigma2 / (nu - 2).
fitted_variance <-
  nu *
  fit_t$sigma2 /
  (
    nu - 2
  )

fitted_results <- data.frame(
  index = seq_len(n),
  X1 = X1,
  y = y,
  true_mu = true_mu,
  true_scale2 = true_sigma2,
  true_variance = true_variance,
  fitted_mu = fit_t$fitted.values,
  fitted_scale2 = fit_t$sigma2,
  fitted_variance = fitted_variance,
  residual = fit_t$residuals,
  scale_residual = fit_t$scale.residuals,
  pearson_residual = fit_t$pearson.residuals
)

print(
  head(fitted_results)
)


# ============================================================
# 14. Mean absolute errors
# ============================================================

# MAE relative to the observed responses measures in-sample
# predictive fit. MAE relative to true_mu measures recovery of
# the data-generating conditional mean.
mae_results <- data.frame(
  Method = paste0(
    "Student-t MLE (nu = ",
    nu,
    ")"
  ),
  MAE_observed = mean(
    abs(
      y -
        fit_t$fitted.values
    )
  ),
  MAE_true_mean = mean(
    abs(
      true_mu -
        fit_t$fitted.values
    )
  ),
  RMSE_true_mean = sqrt(
    mean(
      (
        true_mu -
          fit_t$fitted.values
      )^2
    )
  )
)

print(
  mae_results
)


# ============================================================
# 15. Likelihood criteria
# ============================================================

likelihood_results <- data.frame(
  Method = paste0(
    "Student-t MLE (nu = ",
    nu,
    ")"
  ),
  logLik = fit_t$logLik,
  AIC = fit_t$AIC,
  BIC = fit_t$BIC
)

print(
  likelihood_results
)


# ============================================================
# 16. Predictions for new data
# ============================================================

# Evaluate the fitted mean and squared-scale curves on a regular
# X1 grid while holding w1 and w2 at their sample means.
X1_new <- seq(
  min(X1),
  max(X1),
  length.out = 200
)

W1_new <- cbind(
  intercept = rep(1, length(X1_new)),
  w1 = rep(mean(w1), length(X1_new)),
  w2 = rep(mean(w2), length(X1_new))
)

new_data <- list(
  X1 = X1_new,
  W1 = W1_new
)

predicted_mean <- predict(
  fit_t,
  newdata = new_data,
  type = "mean"
)

predicted_scale2 <- predict(
  fit_t,
  newdata = new_data,
  type = "variance"
)

predicted_variance <-
  nu *
  predicted_scale2 /
  (
    nu - 2
  )

prediction_results <- data.frame(
  X1 = X1_new,
  fitted_mean = predicted_mean,
  fitted_scale2 = predicted_scale2,
  fitted_variance = predicted_variance
)

print(
  head(prediction_results)
)


# ============================================================
# 17. Student-t likelihood weights
# ============================================================

# fit_t$weights contains the likelihood weights
#
#   (nu + 1) / (nu + u_i),
#
# where
#
#   u_i = (y_i - mu_i)^2 / sigma2_i.
#
# These weights arise directly from the Student-t likelihood.
t_weight_results <- data.frame(
  index = seq_len(n),
  scale_residual = fit_t$scale.residuals,
  squared_scale_residual = fit_t$u,
  weight = fit_t$weights
)

print(
  t_weight_results
)


# ============================================================
# 18. True and fitted conditional mean curves
# ============================================================

plot(
  X1,
  y,
  pch = 16,
  cex = 0.6,
  xlab = "X1",
  ylab = "y",
  main = "True and fitted conditional means"
)

lines(
  X1,
  true_mu,
  lty = 2,
  lwd = 2
)

lines(
  X1_new,
  predicted_mean,
  lty = 1,
  lwd = 2
)

legend(
  "topright",
  legend = c(
    "True mean",
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
# 19. Student-t likelihood weights
# ============================================================

plot(
  t_weight_results$index,
  t_weight_results$weight,
  pch = 16,
  xlab = "Index",
  ylab = "Student-t likelihood weight",
  main = "Student-t likelihood weights"
)

abline(
  h = 1,
  lty = 2
)

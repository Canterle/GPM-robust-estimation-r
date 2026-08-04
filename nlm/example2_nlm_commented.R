# ============================================================
# Complete simulated example for fit_nlm_mlqe() and
# fit_nlm_mdpde()
#
# Homoscedastic nonlinear normal regression model
#
# Conditional mean:
#
#   mu_i = beta1 + beta2 * exp(-beta3 * x_i)
#
# Conditional variance:
#
#   sigma2_i = exp(gamma0)
#
# The script illustrates how to:
#
#   1. simulate a homoscedastic nonlinear dataset;
#   2. contaminate one observation;
#   3. define the model and its analytical derivatives;
#   4. construct data-driven starting values;
#   5. fit the corrected MLqE, MDPDE, and normal MLE;
#   6. compare estimates, fitted values, and score weights.
#
# To adapt the script to another dataset, the sections that
# usually require changes are the data input, model definition,
# and starting-value construction.
# ============================================================

# Load the current implementations of the two estimators.
source("fit_nlm_mlqe.R")
source("fit_nlm_mdpde.R")

# Uncomment the following line to reproduce exactly the same
# simulated sample every time the script is run.
# set.seed(20260720)


# ============================================================
# 1. Simulation settings
# ============================================================

# Number of observations in the simulated sample.
n <- 50


# ============================================================
# 2. True parameter values
# ============================================================

# Parameters of the nonlinear conditional mean:
#
# beta1: lower asymptote of the mean curve;
# beta2: difference between the response level at x = 0 and
#        the lower asymptote;
# beta3: exponential decay rate.
true_beta <- c(
  beta1 = 20,
  beta2 = 80,
  beta3 = 0.45
)

# Intercept of the homoscedastic log-variance model.
true_gamma <- c(
  gamma0 = 3.22
)

# Complete true parameter vector. The names and order must agree
# with those used by nlm_model().
true_theta <- c(
  true_beta,
  true_gamma
)


# ============================================================
# 3. Covariate and response generation
# ============================================================

# Generate the covariate used in the conditional mean.
# Sorting is not required for estimation; it only facilitates
# plotting the true and fitted curves.
x <- sort(
  runif(
    n,
    min = 0,
    max = 10
  )
)

# Compute the true conditional mean for every observation.
true_mu <- true_beta["beta1"] +
  true_beta["beta2"] *
  exp(
    -true_beta["beta3"] *
      x
  )

# Compute the true constant conditional variance.
true_sigma2 <- exp(
  true_gamma["gamma0"]
)

# rnorm() requires the standard deviation rather than the
# variance.
true_sigma <- sqrt(
  true_sigma2
)

# Generate the response from the homoscedastic normal model.
y <- rnorm(
  n,
  mean = true_mu,
  sd = true_sigma
)

# Replace one observation by an atypical value to illustrate
# the behavior of the robust estimators. Remove this line when
# a clean simulated sample is desired.
y[2] <- 5

# Store the covariate in a single object. This object is passed
# to nlm_model() during every objective, score, and derivative
# evaluation.
model_data <- list(
  x = x
)


# ============================================================
# 4. Model specification
# ============================================================

# The model function must have the interface:
#
#   nlm_model <- function(theta, data)
#
# theta is the complete named parameter vector. The function
# must return at least mu and sigma2. It may also return:
#
#   D = d mu / d theta^T;
#   V = d sigma2 / d theta^T.
#
# Providing D and V analytically avoids numerical
# differentiation.

nlm_model <- function(theta, data) {
  # Extract the parameters by name.
  beta1 <- theta["beta1"]
  beta2 <- theta["beta2"]
  beta3 <- theta["beta3"]
  gamma0 <- theta["gamma0"]

  # Recover the covariate from the supplied data object.
  x <- data$x
  n <- length(x)

  # Exponential term reused in the mean and its derivatives.
  exponential_term <- exp(
    -beta3 * x
  )

  # Conditional mean.
  mu <- beta1 +
    beta2 *
    exponential_term

  # Constant conditional variance.
  sigma2 <- rep(
    exp(gamma0),
    n
  )

  # Matrix D = d mu / d theta^T.
  # The mean does not depend on gamma0, so the corresponding
  # column is zero.
  D <- cbind(
    beta1 = rep(1, n),
    beta2 = exponential_term,
    beta3 = -beta2 *
      x *
      exponential_term,
    gamma0 = rep(0, n)
  )

  # Matrix V = d sigma2 / d theta^T.
  # The variance does not depend on the beta parameters.
  V <- cbind(
    beta1 = rep(0, n),
    beta2 = rep(0, n),
    beta3 = rep(0, n),
    gamma0 = sigma2
  )

  # Return all quantities required by the fitting functions.
  list(
    mu = mu,
    sigma2 = sigma2,
    D = D,
    V = V
  )
}


# ============================================================
# 5. Transformation tau_r for the corrected MLqE
# ============================================================

# The corrected MLqE requires a transformation family tau_r.
# For the homoscedastic log-variance model, the correction
# changes only gamma0:
#
#   tau_r(gamma0) = gamma0 - log(r).
#
# The mean parameters remain unchanged.
tau <- function(theta, r) {
  result <- theta

  result["gamma0"] <-
    result["gamma0"] -
    log(r)

  result
}


# ============================================================
# 6. Starting values
# ============================================================

x_round <- round(x)

beta1_start <- median(
  y[x_round == max(x_round)]
)

beta2_start <- median(
  y[x_round == min(x_round )]
) - beta1_start

N1_bar <- median(
  y[x_round == unique(x_round)[2]]
)

N2_bar <- median(
  y[x_round == unique(x_round)[3]]
)

# For the exponential mean model,
#
#   N(x) - beta1 = beta2 * exp(-beta3 * x).
#
# Therefore, the response levels at x = 1 and x = 2 give the
# following initial value for beta3.
beta3_start <- log(
  (
    N1_bar -
      beta1_start
  ) / (
    N2_bar -
      beta1_start
  )
) / (
  unique(x_round)[3] - unique(x_round)[2]
)

# Preliminary fitted means.
mu_start <- beta1_start +
  beta2_start *
  exp(
    -beta3_start *
      x
  )

# Initial variance based on the robust residual scale.
sigma2_start <- mad(
  y -
    mu_start
)^2

start_values <- c(
  beta1 = beta1_start,
  beta2 = beta2_start,
  beta3 = beta3_start,
  gamma0 = log(sigma2_start)
)


# ============================================================
# 7. Common optimization settings
# ============================================================

# Controls passed to optim() for MLqE, MDPDE, and MLE.
optimization_control <- list(
  maxit = 2000,
  reltol = 1e-12
)

# Method used by numDeriv when a numerical Jacobian is required.
# D and V are analytical in this example, but numerical
# derivatives may still be needed internally for transformations.
jacobian_method <- "Richardson"
jacobian_method_args <- list()

# Confidence level used for standard errors and intervals.
conf_level <- 0.95

# Configuration of the automatic SQV procedure.
sqv_control <- list(
  # Initial grid from q = 1 down to q = 0.8.
  initial_lower = 0.8,

  # Number of values in the initial grid. With the endpoints
  # above, 21 points produce spacing 0.01.
  m0 = 21,

  # Number of q values in each subsequent refinement grid.
  m = 3,

  # Smallest q value that may be evaluated.
  q_min = 0.5,

  # Stability threshold used to compare adjacent estimates.
  L = 0.01,

  # Retain all evaluated fits in fit$q.selection$fits.
  keep_fits = TRUE
)


# ============================================================
# 8. Corrected MLqE fit with q selected by SQV
# ============================================================

# The SQV procedure selects q automatically. compute_vcov must
# be TRUE because SQV compares standardized estimates.
fit_mlqe <- fit_nlm_mlqe(
  y = y,
  model = nlm_model,
  start = start_values,
  data = model_data,
  q = "SQV",
  tau = tau,
  level = conf_level,
  method = "BFGS",
  control = optimization_control,
  use_score = TRUE,
  compute_vcov = TRUE,
  jacobian_method = jacobian_method,
  jacobian_method_args = jacobian_method_args,
  q_control = sqv_control
)


# ============================================================
# 9. MDPDE fit with q selected by SQV
# ============================================================

# The same model, starting values, optimization controls, and
# SQV configuration are used to make the comparison with MLqE
# as direct as possible.
fit_mdpde <- fit_nlm_mdpde(
  y = y,
  model = nlm_model,
  start = start_values,
  data = model_data,
  q = "SQV",
  level = conf_level,
  method = "BFGS",
  control = optimization_control,
  use_score = TRUE,
  compute_vcov = TRUE,
  jacobian_method = jacobian_method,
  jacobian_method_args = jacobian_method_args,
  q_control = sqv_control
)


# ============================================================
# 10. Normal maximum likelihood fit
# ============================================================

# For q = 1, both procedures reduce to maximum likelihood
# estimation under the normal model. Only one reference fit is
# needed; here it is obtained through fit_nlm_mdpde().
fit_mle <- fit_nlm_mdpde(
  y = y,
  model = nlm_model,
  start = start_values,
  data = model_data,
  q = 1,
  level = conf_level,
  method = "BFGS",
  control = optimization_control,
  use_score = TRUE,
  compute_vcov = TRUE,
  jacobian_method = jacobian_method,
  jacobian_method_args = jacobian_method_args
)


# ============================================================
# 11. Fit summaries
#
# Print separate headings so the console output from each
# estimator can be identified easily.
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
#
# Inspect the selected robustness parameters and the sequence of
# q values evaluated by the SQV procedure.
# ============================================================

selected_q <- c(
  MLqE = fit_mlqe$q,
  MDPDE = fit_mdpde$q
)

print(
  selected_q
)

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
#
# Compare the estimates with the true values used to generate
# the simulated data.
# ============================================================

parameter_estimates <- cbind(
  True = true_theta,
  `Corrected MLqE` = fit_mlqe$coefficients,
  MDPDE = fit_mdpde$coefficients,
  MLE = fit_mle$coefficients
)

print(
  parameter_estimates
)

# Compare the original MLqE estimate with its consistency-
# corrected version.
mlqe_original_and_corrected <- cbind(
  `Original MLqE` = fit_mlqe$coefficients.star,
  `Corrected MLqE` = fit_mlqe$coefficients
)

print(
  mlqe_original_and_corrected
)


# ============================================================
# 14. Standard errors and confidence intervals
#
# Combine the estimates and standard errors from all methods and
# inspect their confidence intervals.
# ============================================================

standard_error_results <- cbind(
  `MLqE estimate` = fit_mlqe$coefficients,
  `MLqE SE` = fit_mlqe$standard.error,
  `MDPDE estimate` = fit_mdpde$coefficients,
  `MDPDE SE` = fit_mdpde$standard.error,
  `MLE estimate` = fit_mle$coefficients,
  `MLE SE` = fit_mle$standard.error
)

print(
  standard_error_results
)

fit_mlqe$conf.int
fit_mdpde$conf.int
fit_mle$conf.int


# ============================================================
# 15. Covariance and estimating-function matrices
#
# These objects are useful for checking the sandwich covariance
# calculations and for methodological analyses.
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
# 16. Analytical model derivatives
#
# Confirm the derivative matrices evaluated at the final fitted
# parameter values.
# ============================================================

fit_mlqe$mean.jacobian.star
fit_mlqe$variance.jacobian.star

fit_mdpde$mean.jacobian
fit_mdpde$variance.jacobian


# ============================================================
# 17. Fitted values, variances, and residuals
#
# Collect observation-level results in one table for comparison.
# ============================================================

fitted_results <- data.frame(
  index = seq_len(n),
  x = x,
  y = y,
  true_mu = true_mu,
  true_sigma2 = rep(true_sigma2, n),
  mlqe_mu = fit_mlqe$fitted.values,
  mlqe_sigma2 = fit_mlqe$sigma2,
  mlqe_residual = fit_mlqe$residuals,
  mdpde_mu = fit_mdpde$fitted.values,
  mdpde_sigma2 = fit_mdpde$sigma2,
  mdpde_residual = fit_mdpde$residuals,
  mle_mu = fit_mle$fitted.values,
  mle_sigma2 = fit_mle$sigma2,
  mle_residual = fit_mle$residuals
)

print(
  head(
    fitted_results
  )
)


# ============================================================
# 18. Predictions for new data
#
# Evaluate the fitted mean curves on a regular x grid.
# ============================================================

x_new <- seq(
  min(x),
  max(x),
  length.out = 200
)

new_data <- list(
  x = x_new
)

predicted_mean_mlqe <- predict(
  fit_mlqe,
  newdata = new_data,
  type = "mean"
)

predicted_sigma2_mlqe <- predict(
  fit_mlqe,
  newdata = new_data,
  type = "variance"
)

predicted_mean_mdpde <- predict(
  fit_mdpde,
  newdata = new_data,
  type = "mean"
)

predicted_sigma2_mdpde <- predict(
  fit_mdpde,
  newdata = new_data,
  type = "variance"
)

predicted_mean_mle <- predict(
  fit_mle,
  newdata = new_data,
  type = "mean"
)


# ============================================================
# 19. Score weights
#
# Compute the density-power component f(y_i; theta)^(1 - q).
# Smaller values indicate observations receiving less influence.
# ============================================================

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

score_weights <- data.frame(
  index = seq_len(n),
  MLqE = weight_mlqe,
  MDPDE = weight_mdpde
)

print(
  score_weights
)


# ============================================================
# 20. Graphical comparison of the true and fitted mean curves
#
# In a real-data application, the true mean is unavailable and
# should be omitted from this plot.
# ============================================================

plot(
  x,
  y,
  pch = 16,
  cex = 0.6,
  xlab = "x",
  ylab = "y",
  main = "True and fitted conditional means"
)

lines(
  x,
  true_mu,
  lty = 2,
  lwd = 2
)

lines(
  x_new,
  predicted_mean_mlqe,
  lty = 1,
  lwd = 2
)

lines(
  x_new,
  predicted_mean_mdpde,
  lty = 3,
  lwd = 2
)

lines(
  x_new,
  predicted_mean_mle,
  lty = 4,
  lwd = 2
)

legend(
  "topright",
  legend = c(
    "True mean",
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
    "Normal MLE"
  ),
  lty = c(
    2,
    1,
    3,
    4
  ),
  lwd = 2,
  bty = "n"
)


# ============================================================
# 21. Graphical comparison of score weights
#
# Compare how strongly MLqE and MDPDE downweight each
# observation.
# ============================================================

plot(
  score_weights$index,
  score_weights$MDPDE,
  pch = 1,
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

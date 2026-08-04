# ============================================================
# Complete simulated example for fit_nlm_mlqe() and
# fit_nlm_mdpde()
#
# Heteroscedastic nonlinear normal regression model
#
# Conditional mean:
#
#   mu_i = beta1 + beta2 /
#          (1 + beta3 * X1_i^beta4)
#
# Conditional variance:
#
#   sigma2_i = exp(
#     gamma0 + gamma1 * w1_i + gamma2 * w2_i
#   )
#
# The script illustrates how to:
#
#   1. simulate a heteroscedastic nonlinear dataset;
#   2. contaminate one observation;
#   3. define the model and its analytical derivatives;
#   4. construct data-driven starting values;
#   5. fit the corrected MLqE, MDPDE, and normal MLE;
#   6. compare estimates, fitted values, and score weights.
#
# To adapt the script to another dataset, the sections that
# usually require changes are the data input, model definition,
# variance design matrix, and starting-value construction.
# ============================================================

# Load the current implementations of the two estimators.
source("fit_nlm_mlqe.R")
source("fit_nlm_mdpde.R")

# Uncomment the following line to reproduce exactly the same
# simulated sample every time the script is run.
# set.seed(20260721)


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

# Parameters of the log-variance model:
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

# Generate the covariates used in the variance model.
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

# Design matrix of the log-variance model. The first column is
# the intercept and the other columns are observed covariates.
# For another variance structure, change W1, true_gamma,
# nlm_model(), and the corresponding starting values together.
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

# Compute the true conditional variance. Because W1 contains an
# intercept column, W1 %*% true_gamma evaluates
# gamma0 + gamma1 * w1 + gamma2 * w2.
true_sigma2 <- exp(
  drop(
    W1 %*% true_gamma
  )
)

# rnorm() requires standard deviations rather than variances.
true_sigma <- sqrt(
  true_sigma2
)

# Generate the response from the heteroscedastic normal model.
y <- rnorm(
  n,
  mean = true_mu,
  sd = true_sigma
)

# Replace one observation by an atypical value to illustrate
# the behavior of the robust estimators. Remove this line when
# a clean simulated sample is desired.
y[4] <- 800


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
# must return at least mu and sigma2. It may also return:
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

  # Extract the variance parameters.
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
    beta2 / denominator

  # Conditional variance.
  sigma2 <- exp(
    drop(
      W1 %*% gamma
    )
  )

  # Derivative of X1^beta4 with respect to beta4. The value
  # is set to zero at X1 = 0 to avoid evaluating log(0).
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
  # respect to gamma is sigma2_i * W1_i. sweep() multiplies each
  # row of W1 by the corresponding conditional variance.
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
  # The variance does not depend on the beta parameters, so the
  # corresponding columns are zero.
  V <- cbind(
    beta1 = rep(0, n),
    beta2 = rep(0, n),
    beta3 = rep(0, n),
    beta4 = rep(0, n),
    V_gamma
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
# For this log-linear variance model, the correction changes
# only the log-variance intercept:
#
#   tau_r(gamma0) = gamma0 - log(r).
#
# The mean parameters and variance slopes remain unchanged.
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
    (unique(X1_round)[3] - unique(X1_round)[2])
)

# Preliminary fitted means.
mu_start <- beta1_start +
  beta2_start /
  (
    1 +
      beta3_start *
      X1^beta4_start
  )

# Initial variance based on the robust residual scale.
sigma2_start <- mad(
  y -
    mu_start
)^2

# The variance coefficients other than the intercept initially
# equal zero.
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

print(parameter_estimates)

mlqe_original_and_corrected <- cbind(
  `Original MLqE` = fit_mlqe$coefficients.star,
  `Corrected MLqE` = fit_mlqe$coefficients
)

print(mlqe_original_and_corrected)


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

print(standard_error_results)

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
  X1 = X1,
  y = y,
  true_mu = true_mu,
  true_sigma2 = true_sigma2,
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
  head(fitted_results)
)


# ============================================================
# 18. Predictions for new data
#
# Evaluate fitted mean curves on a regular X1 grid while holding
# the variance covariates at their sample means. For another
# application, choose scientifically meaningful fixed values.
# ============================================================

X1_new <- seq(
  min(X1),
  max(X1),
  length.out = 200
)

W1_new <- cbind(
  intercept = 1,
  w1 = mean(w1),
  w2 = mean(w2)
)

W1_new <- W1_new[
  rep(
    1,
    length(X1_new)
  ),
  ,
  drop = FALSE
]

new_data <- list(
  X1 = X1_new,
  W1 = W1_new
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

print(score_weights)


# ============================================================
# 20. Graphical comparison of the true and fitted mean curves
#
# In a real-data application, the true mean is unavailable and
# should be omitted from this plot.
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
  predicted_mean_mlqe,
  lty = 1,
  lwd = 2
)

lines(
  X1_new,
  predicted_mean_mdpde,
  lty = 3,
  lwd = 2
)

lines(
  X1_new,
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
      format(fit_mlqe$q, digits = 3),
      ")"
    ),
    paste0(
      "MDPDE (q = ",
      format(fit_mdpde$q, digits = 3),
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
# Compare how strongly MLqE and MDPDE downweight each observation.
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


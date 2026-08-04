# ============================================================
# Application to a heteroscedastic nonlinear regression model
#
# The following models and estimators are considered:
#
#   1. Student-t maximum likelihood estimator;
#   2. Normal maximum density power divergence estimator
#      (MDPDE);
#   3. Normal corrected maximum Lq-likelihood estimator
#      (MLqE);
#   4. Normal maximum likelihood estimator using all data;
#   5. Normal maximum likelihood estimator after removing
#      the outlier.
#
# Mean structure:
#
#   mu_i = beta1 + beta2 /
#          (1 + beta3 * X1_i^beta4)
#
# Scale structure:
#
#   sigma2_i = exp(gamma0 + gamma1 * X1_i)
#
# Under the Student-t model, sigma2_i denotes the square of
# the scale parameter. Therefore, for nu > 2,
#
#   Var(Y_i) = nu * sigma2_i / (nu - 2).
# ============================================================


# ------------------------------------------------------------
# 1. Load the estimation functions
# ------------------------------------------------------------

source("fit_nlm_t.R")
source("fit_nlm_mlqe.R")
source("fit_nlm_mdpde.R")


# ------------------------------------------------------------
# 2. Data
# ------------------------------------------------------------

y <- c(
  7720,
  8113,
  6664,
  6804,
  4994,
  4948,
  3410,
  3208,
  4478,
  2396,
  1302,
  1377,
  1025,
  1096
)

X1 <- c(
  0.000001,
  0.000001,
  2,
  2,
  5,
  5,
  10,
  10,
  20,
  20,
  50,
  50,
  100,
  100
)

# Perform basic data checks.
stopifnot(
  length(y) == length(X1),
  all(is.finite(y)),
  all(is.finite(X1)),
  all(X1 > 0)
)

n <- length(y)

# ------------------------------------------------------------
# 3. Auxiliary function for computing fitted means
# ------------------------------------------------------------

mu_f <- function(beta, data) {
  beta1 <- beta["beta1"]
  beta2 <- beta["beta2"]
  beta3 <- beta["beta3"]
  beta4 <- beta["beta4"]
  
  X1 <- data$X1
  
  beta1 +
    beta2 /
    (
      1 +
        beta3 *
        X1^beta4
    )
}

# ------------------------------------------------------------
# 4. Starting values
# ------------------------------------------------------------
# All sets of starting values produce similar results.
# Starting values for the mean parameters suggested by 
# Thiede and Pagano (1979).
start_beta <- c(
  beta1 = 1060.5,
  beta2 = 6856,
  beta3 = 0.085716,
  beta4 = 1
)

# # The starting values for the mean parameters are the
# # estimates reported by Lemonte and Patriota (2011).
# start_beta <- c(
#   beta1 = 929.2840,
#   beta2 = 6881.7149,
#   beta3 = 0.0781,
#   beta4 = 1.3562
# )
# 
# # Starting values for the mean parameters are the Hubber
# # estimates reported by Thiede and Pagano (2011).
# start_beta <- c(
#   beta1 = 919.9,
#   beta2 = 7022.7149,
#   beta3 = 0.086,
#   beta4 = 1.33
# )

# Evaluate the mean function at the starting values 
mu_start <- mu_f(start_beta, list(X1 = X1))

# The starting value for gamma0 is obtained from a robust
# estimate of the residual scale, whereas gamma1 is initially
# set equal to zero.
sigma2_start <- mad(y - mu_start)^2

stopifnot(
  is.finite(sigma2_start),
  sigma2_start > 0
)

# Construct the complete vector of starting values.
start_values <- c(
  start_beta,
  gamma0 = log(sigma2_start),
  gamma1 = 0
)


# ------------------------------------------------------------
# 5. Nonlinear model specification
#
# The model function must have the following interface:
#
#   function(theta, data)
#
# It returns:
#
#   mu:     vector of means;
#   sigma2: vector of squared scale parameters;
#   D:      derivative of mu with respect to theta;
#   V:      derivative of sigma2 with respect to theta.
# ------------------------------------------------------------

nlm_model <- function(theta, data) {
  beta1 <- theta["beta1"]
  beta2 <- theta["beta2"]
  beta3 <- theta["beta3"]
  beta4 <- theta["beta4"]
  gamma0 <- theta["gamma0"]
  gamma1 <- theta["gamma1"]
  
  X1 <- data$X1
  n <- length(X1)
  
  # Auxiliary quantities.
  X1_beta4 <- X1^beta4
  denominator <- 1 + beta3 * X1_beta4
  
  # Mean structure.
  mu <- beta1 +
    beta2 / denominator
  
  # Squared scale structure.
  sigma2 <- exp(
    gamma0 +
      gamma1 * X1
  )
  
  # Derivative of X1^beta4 with respect to beta4.
  X1_beta4_log <- X1_beta4 * log(X1)
  
  # Matrix D = d mu / d theta^T.
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
    
    gamma1 = rep(0, n)
  )
  
  # Matrix V = d sigma2 / d theta^T.
  V <- cbind(
    beta1 = rep(0, n),
    
    beta2 = rep(0, n),
    
    beta3 = rep(0, n),
    
    beta4 = rep(0, n),
    
    gamma0 = sigma2,
    
    gamma1 = X1 * sigma2
  )
  
  list(
    mu = mu,
    sigma2 = sigma2,
    D = D,
    V = V
  )
}

# ------------------------------------------------------------
# 6. Transformation tau_r for the corrected MLqE
#
# Since
#
#   sigma2_i = exp(gamma0 + gamma1 * X1_i),
#
# the transformation modifies only gamma0. The parameter
# gamma1 remains unchanged.
# ------------------------------------------------------------

tau <- function(theta, r) {
  result <- theta
  
  result["gamma0"] <-
    result["gamma0"] -
    log(r)
  
  result
}

# ------------------------------------------------------------
# 7. Common optimization settings
# ------------------------------------------------------------

optimization_control <- list(
  maxit = 100000,
  reltol = 1e-12
)

jacobian_method <- "Richardson"
jacobian_method_args <- list()

conf_level <- 0.95

# ------------------------------------------------------------
# 8. Student-t maximum likelihood fit
#
# The degrees of freedom are assumed to be known and are
# fixed at nu = 4.
# ------------------------------------------------------------

nu <- 4

fit_t <- fit_nlm_t(
  y = y,
  model = nlm_model,
  start_theta = start_values,
  data = list(X1 = X1),
  nu = nu,
  method = "BFGS",
  control = optimization_control,
  use_gradient = TRUE,
  hessian = TRUE,
  conf_level = conf_level
)


# ------------------------------------------------------------
# 9. Normal MDPDE fit
#
# The robustness parameter is selected using the SQV
# procedure.
# ------------------------------------------------------------

fit_mdpde <- fit_nlm_mdpde(
  y = y,
  model = nlm_model,
  start = start_values,
  data = list(X1 = X1),
  q = "SQV",
  level = conf_level,
  method = "BFGS",
  control = optimization_control,
  use_score = TRUE,
  compute_vcov = TRUE,
  jacobian_method = jacobian_method,
  jacobian_method_args = jacobian_method_args,
  q_control = list(
    m0 = 21L,
    m = 3L,
    q_min = 0.5,
    L = 0.01
  )
)


# ------------------------------------------------------------
# 10. Normal corrected MLqE fit
#
# The distortion parameter q is selected using the SQV
# procedure.
# ------------------------------------------------------------

fit_mlqe <- fit_nlm_mlqe(
  y = y,
  model = nlm_model,
  start = start_values,
  data = list(X1 = X1),
  q = "SQV",
  tau = tau,
  level = conf_level,
  method = "BFGS",
  control = optimization_control,
  use_score = TRUE,
  compute_vcov = TRUE,
  jacobian_method = jacobian_method,
  jacobian_method_args = jacobian_method_args,
  q_control = list(
    m0 = 21L,
    m = 3L,
    q_min = 0.75,
    L = 0.01
  )
)


# ------------------------------------------------------------
# 11. Normal maximum likelihood fit using all observations
#
# For q = 1, the MDPDE and the MLqE reduces to the maximum 
# likelihood estimator under the normal model.
#
# With 9699 degrees of freedom, the maximum likelihood 
# estimator under the Student's t model approximates that 
# under the normal model.
# ------------------------------------------------------------

fit_mle <- fit_nlm_mdpde(
  y = y,
  model = nlm_model,
  start = start_values,
  data = list(X1 = X1),
  q = 1,
  level = conf_level,
  method = "BFGS",
  control = optimization_control,
  use_score = TRUE,
  compute_vcov = TRUE,
  jacobian_method = jacobian_method,
  jacobian_method_args = jacobian_method_args
)

fit_mle2 <- fit_nlm_mlqe(
  y = y,
  model = nlm_model,
  start = start_values,
  data = list(X1 = X1),
  q = 1,
  tau = tau,
  level = conf_level,
  method = "BFGS",
  control = optimization_control,
  use_score = TRUE,
  compute_vcov = TRUE,
  jacobian_method = jacobian_method,
  jacobian_method_args = jacobian_method_args
)

nu <- 9699

fit_mle3 <- fit_nlm_t(
  y = y,
  model = nlm_model,
  start_theta = start_values,
  data = list(X1 = X1),
  nu = nu,
  method = "BFGS",
  control = optimization_control,
  use_gradient = TRUE,
  hessian = TRUE,
  conf_level = conf_level
)

# ------------------------------------------------------------
# 12. Normal maximum likelihood fit after removing the outlier
#
# Observation 9, corresponding to
#
#   X1 = 20
#   y  = 4478,
#
# is excluded from this fit.
# ------------------------------------------------------------

outlier_index <- 9L

y_without_outlier <- y[-outlier_index]
X1_without_outlier <- X1[-outlier_index]

fit_mle_wo <- fit_nlm_mdpde(
  y = y_without_outlier,
  model = nlm_model,
  start = start_values,
  data = list(X1 = X1_without_outlier),
  q = 1,
  level = conf_level,
  method = "BFGS",
  control = optimization_control,
  use_score = TRUE,
  compute_vcov = TRUE,
  jacobian_method = jacobian_method,
  jacobian_method_args = jacobian_method_args
)


# ------------------------------------------------------------
# 13. Fit summaries and confidence intervals
# ------------------------------------------------------------

cat(
  "\n",
  paste(rep("=", 70), collapse = ""),
  "\nSTUDENT-t MAXIMUM LIKELIHOOD FIT\n",
  paste(rep("=", 70), collapse = ""),
  "\n",
  sep = ""
)

summary(fit_t)
fit_t$conf.int


cat(
  "\n",
  paste(rep("=", 70), collapse = ""),
  "\nNORMAL MDPDE FIT\n",
  paste(rep("=", 70), collapse = ""),
  "\n",
  sep = ""
)

summary(fit_mdpde)
fit_mdpde$conf.int

cat(
  "\n",
  paste(rep("=", 70), collapse = ""),
  "\nNORMAL CORRECTED MLqE FIT\n",
  paste(rep("=", 70), collapse = ""),
  "\n",
  sep = ""
)

summary(fit_mlqe)
fit_mlqe$conf.int

cat(
  "\n",
  paste(rep("=", 70), collapse = ""),
  "\nNORMAL MAXIMUM LIKELIHOOD FIT USING ALL DATA\n",
  paste(rep("=", 70), collapse = ""),
  "\n",
  sep = ""
)

summary(fit_mle)
fit_mle$conf.int

cat(
  "\n",
  paste(rep("=", 70), collapse = ""),
  "\nNORMAL MAXIMUM LIKELIHOOD FIT WITHOUT THE OUTLIER\n",
  paste(rep("=", 70), collapse = ""),
  "\n",
  sep = ""
)

summary(fit_mle_wo)
fit_mle_wo$conf.int

# ------------------------------------------------------------
# 14. Fitted means and MAE
#
# The fitted means are evaluated at all original values of X1,
# including the curve obtained from the model fitted after
# removing the outlier.
# ------------------------------------------------------------

beta_names <- c(
  "beta1",
  "beta2",
  "beta3",
  "beta4"
)

mu_t <- mu_f(
  beta = fit_t$coefficients[beta_names],
  data = list(X1 = X1)
)

mu_mdpde <- mu_f(
  beta = fit_mdpde$coefficients[beta_names],
  data = list(X1 = X1)
)

mu_mlqe <- mu_f(
  beta = fit_mlqe$coefficients[beta_names],
  data = list(X1 = X1)
)

mu_mle <- mu_f(
  beta = fit_mle$coefficients[beta_names],
  data = list(X1 = X1)
)

mu_mle_wo <- mu_f(
  beta = fit_mle_wo$coefficients[beta_names],
  data = list(X1 = X1)
)

# Indices of the observations identified as outliers
outlier_indices <- c(78, 80)

# Observations retained after removing observations 78 and 80
non_outlier_indices <- setdiff(
  seq_along(y),
  outlier_indices
)

# MAE computed using all observations

mae_t_all <- mean(
  abs(y - mu_t)
)

mae_mdpde_all <- mean(
  abs(y - mu_mdpde)
)

mae_mlqe_all <- mean(
  abs(y - mu_mlqe)
)

mae_mle_all <- mean(
  abs(y - mu_mle)
)

mae_mle_wo_all <- mean(
  abs(y - mu_mle_wo)
)

# MAE computed after removing observations 78 and 80

mae_t_wo <- mean(
  abs(
    y[non_outlier_indices] -
      mu_t[non_outlier_indices]
  )
)

mae_mdpde_wo <- mean(
  abs(
    y[non_outlier_indices] -
      mu_mdpde[non_outlier_indices]
  )
)

mae_mlqe_wo <- mean(
  abs(
    y[non_outlier_indices] -
      mu_mlqe[non_outlier_indices]
  )
)

mae_mle_wo_outliers <- mean(
  abs(
    y[non_outlier_indices] -
      mu_mle[non_outlier_indices]
  )
)

mae_mle_wo_wo <- mean(
  abs(
    y[non_outlier_indices] -
      mu_mle_wo[non_outlier_indices]
  )
)

# MAE comparison

mae_comparison <- data.frame(
  Method = c(
    "Student-t MLE (nu = 5)",
    "MDPDE",
    "Corrected MLqE",
    "Normal MLE",
    "Normal MLE without observations 78 and 80"
  ),
  MAE_all_observations = c(
    mae_t_all,
    mae_mdpde_all,
    mae_mlqe_all,
    mae_mle_all,
    mae_mle_wo_all
  ),
  MAE_without_observations_78_80 = c(
    mae_t_wo,
    mae_mdpde_wo,
    mae_mlqe_wo,
    mae_mle_wo_outliers,
    mae_mle_wo_wo
  ),
  row.names = NULL
)

print(
  mae_comparison,
  digits = 6
)


# ------------------------------------------------------------
# 15. Plot of the fitted curves
# ------------------------------------------------------------

# Order the covariate values before drawing the fitted curves.
plot_order <- order(X1)

# Value displayed in the MLqE legend. This value must agree
# with the value selected by the SQV procedure.
q_mlqe_plot <- fit_mlqe$q

# Value displayed in the MDPDE legend. This value must agree
# with the value selected by the SQV procedure.
q_mdpde_plot <- fit_mdpde$q

mlqe_legend <- bquote(
  "Normal corrected ML"[q] *
    "E (q = " *
    .(sprintf("%.2f", q_mlqe_plot)) *
    ")"
)

mdpde_legend <- bquote(
  "Normal MDPDE (q = " *
    .(sprintf("%.2f", q_mdpde_plot)) *
    ")"
)

pdf(
  file = "XvsYapplication.pdf",
  width = 6.5,
  height = 3
)

# Set the bottom, left, top, and right margins.
par(
  mar = c(2.3, 2.3, 0.05, 0.05),
  mgp = c(1.4, 0.5, 0)
)

plot(
  X1,
  y,
  xlab = "Thyrotropin dose",
  ylab = "Observed radioactivity",
  pch = 1,
  lwd = 2
)

# Normal maximum likelihood fit using all observations.
lines(
  X1[plot_order],
  mu_mle[plot_order],
  lty = 2,
  lwd = 2,
  col = "red"
)

# Normal maximum likelihood fit after removing the outlier.
lines(
  X1[plot_order],
  mu_mle_wo[plot_order],
  lty = 4,
  lwd = 2,
  col = "darkgreen"
)

# Student-t maximum likelihood fit.
lines(
  X1[plot_order],
  mu_t[plot_order],
  lty = 4,
  lwd = 2,
  col = "black"
)

# Normal MDPDE fit.
lines(
  X1[plot_order],
  mu_mdpde[plot_order],
  lty = 5,
  lwd = 2,
  col = "purple4"
)

# Normal corrected MLqE fit.
lines(
  X1[plot_order],
  mu_mlqe[plot_order],
  lty = 3,
  lwd = 2,
  col = "blue"
)

legend(
  "topright",
  
  legend = as.expression(
    list(
      "Normal MLE (all data)",
      "Normal MLE (outlier removed)",
      "Student-t MLE",
      mdpde_legend,
      mlqe_legend,
      "Data"
    )
  ),
  
  col = c(
    "red",
    "darkgreen",
    "black",
    "purple4",
    "blue",
    "black"
  ),
  
  lty = c(
    2,
    4,
    4,
    5,
    3,
    NA
  ),
  
  lwd = c(
    2,
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
    NA,
    1
  ),
  
  bty = "n"
)

dev.off()

# ------------------------------------------------------------
# 16. Plot of the score weights
# ------------------------------------------------------------

# The score weights f(y_i; theta)^(1 - q), 
# where f is the normal density.



pdf(
  file = "score_weights_application.pdf",
  width = 6.5,
  height = 3
)

weight_mdpde <- dnorm(y,fit_mdpde$fitted.values,fit_mdpde$sigma)^(1-fit_mdpde$q)
weight_mlqe <- dnorm(y,fit_mlqe$fitted.values,fit_mlqe$sigma)^(1-fit_mlqe$q)

# Set the bottom, left, top, and right margins.
par(
  mar = c(2.3, 2.3, 0.05, 0.05),
  mgp = c(1.4, 0.5, 0)
)

plot(weight_mdpde,
     xlab = "Index",
     ylab = "Score weight",
     pch = 1,
     ylim = c(0,0.7)
)

points(weight_mlqe,
       xlab = "Index",
       ylab = "Score weight",
       pch = 3,
)

# Add a horizontal reference line at zero.
abline(
  h = 0,
  lty = 2,
  lwd = 1
)

legend(
  x = 0.35,
  y = 0.175,
  legend = as.expression(
    list(
      mdpde_legend,
      mlqe_legend
    )
  ),
  

  pch = c(
    1,
    3
  ),
  
  bty = "n"

)

dev.off()


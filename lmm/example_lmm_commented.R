# ============================================================
# Complete simulated example for fit_lmm_mlqe() and
# fit_lmm_mdpde()
#
# Gaussian linear mixed model
#
# For group i:
#
#   Y_i = X_i beta + Z_i b_i + epsilon_i
#
# where:
#
#   b_i       ~ Normal(0, Delta)
#   epsilon_i ~ Normal(0, sigma2 I)
#
# The simulation uses:
#
#   beta   = (18, 10)
#   gamma  = (3, -2, 4)
#   sigma2 = 5
#
# with the unstructured random-effects covariance matrix:
#
#   Delta(gamma) =
#
#     [ gamma1  gamma2 ]
#     [ gamma2  gamma3 ].
#
# The script illustrates how to:
#
#   1. simulate grouped data from a Gaussian linear mixed model;
#   2. contaminate one complete group;
#   3. define Delta(gamma), its analytical Jacobian, and tau_r;
#   4. calculate robust starting values internally;
#   5. fit the corrected MLqE, MDPDE, and normal MLE;
#   6. inspect the SQV selection results;
#   7. compare estimates, fitted values, random effects, and
#      group-level score weights.
# ============================================================

# Load the current implementations of the two estimators.
source("fit_lmm_mlqe.R")
source("fit_lmm_mdpde.R")

# # Fix the seed so the simulated example is reproducible.
# set.seed(1234)


# ============================================================
# 1. Simulation settings
# ============================================================

# Number of independent groups.
n_groups <- 50

# Number of observations in every group.
group_size <- 5

group_sizes <- rep(
  group_size,
  n_groups
)


# ============================================================
# 2. True parameter values
# ============================================================

# Fixed-effect parameters.
true_beta <- c(
  beta0 = 18,
  beta1 = 10
)

# Parameters of the unstructured random-effects covariance
# matrix.
true_gamma <- c(
  gamma1 = 3,
  gamma2 = -2,
  gamma3 = 4
)

# Residual variance.
true_sigma2 <- c(
  sigma2 = 5
)

# Complete true parameter vector.
true_theta <- c(
  true_beta,
  true_gamma,
  true_sigma2
)


# ============================================================
# 3. Random-effects covariance structure
# ============================================================

# Construct the unstructured 2 by 2 covariance matrix
# Delta(gamma).
Delta_unstructured <- function(gamma) {
  matrix(
    c(
      gamma["gamma1"],
      gamma["gamma2"],
      gamma["gamma2"],
      gamma["gamma3"]
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      c(
        "random_intercept",
        "random_slope"
      ),
      c(
        "random_intercept",
        "random_slope"
      )
    )
  )
}


# This matrix identifies where each gamma parameter appears in
# Delta(gamma). It is required when start = NULL.
#
# Parameters appearing on the main diagonal receive automatic
# starting value 1. Parameters appearing only off the diagonal
# receive automatic starting value 0.
gamma_structure_unstructured <- matrix(
  c(
    "gamma1",
    "gamma2",
    "gamma2",
    "gamma3"
  ),
  nrow = 2,
  byrow = TRUE,
  dimnames = list(
    c(
      "random_intercept",
      "random_slope"
    ),
    c(
      "random_intercept",
      "random_slope"
    )
  )
)


# Analytical Jacobian of vec(Delta(gamma)) with respect to
# gamma. R vectorizes matrices column by column, so the rows
# correspond to Delta_11, Delta_21, Delta_12, and Delta_22.
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
    byrow = TRUE,
    dimnames = list(
      c(
        "Delta_11",
        "Delta_21",
        "Delta_12",
        "Delta_22"
      ),
      names(gamma)
    )
  )
}


# True covariance matrix of the random effects.
true_Delta <- Delta_unstructured(
  true_gamma
)

true_Delta

# Positive eigenvalues confirm that true_Delta is positive
# definite.
eigen(
  true_Delta,
  symmetric = TRUE
)$values


# ============================================================
# 4. Transformation tau_r for the corrected MLqE
# ============================================================

# The consistency correction leaves beta unchanged and divides
# all covariance parameters by r:
#
#   tau_r(theta) =
#     (beta0, beta1,
#      gamma1 / r, gamma2 / r, gamma3 / r,
#      sigma2 / r).
#
# Therefore, tau_{1/q} maps corrected parameters to the original
# MLqE parameterization, whereas tau_q performs the consistency
# correction.
tau <- function(theta,
                r) {
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


# ============================================================
# 5. Design matrices
# ============================================================

X <- vector(
  "list",
  n_groups
)

Z <- vector(
  "list",
  n_groups
)

for (i in seq_len(n_groups)) {
  x1 <- runif(
    group_sizes[i],
    min = -2,
    max = 2
  )

  # Fixed intercept and fixed slope.
  X[[i]] <- cbind(
    beta0 = 1,
    beta1 = x1
  )

  # Random intercept and random slope.
  Z[[i]] <- cbind(
    random_intercept = 1,
    random_slope = x1
  )
}


# ============================================================
# 6. Response generation
# ============================================================

Y <- vector(
  "list",
  n_groups
)

true_random_effects <- vector(
  "list",
  n_groups
)

Delta_cholesky <- t(
  chol(
    true_Delta
  )
)

for (i in seq_len(n_groups)) {
  random_effect_i <- drop(
    Delta_cholesky %*%
      rnorm(2)
  )

  residual_error_i <- rnorm(
    group_sizes[i],
    mean = 0,
    sd = sqrt(
      true_sigma2["sigma2"]
    )
  )

  Y[[i]] <-
    drop(
      X[[i]] %*%
        true_beta
    ) +
    drop(
      Z[[i]] %*%
        random_effect_i
    ) +
    residual_error_i

  true_random_effects[[i]] <-
    setNames(
      random_effect_i,
      colnames(
        Z[[i]]
      )
    )
}


# Contaminate group 30. Remove these two lines when a clean
# simulated sample is desired.
Y[[30]] <- Y[[30]]*3 + 10

# Y[[30]][3:5] <-
#   Y[[30]][3:5] +
#   20


# ============================================================
# 7. Starting values
# ============================================================

# Starting values are not calculated in this example.
#
# With start = NULL, fit_lmm_mlqe() and fit_lmm_mdpde()
# calculate them internally:
#
#   - beta is initialized by MASS::rlm() fitted to the stacked
#     response and fixed-effect design matrices;
#   - sigma2 is initialized by the squared robust residual scale;
#   - gamma parameters appearing on the diagonal of
#     gamma_structure receive 1;
#   - gamma parameters appearing only off the diagonal receive 0.
#
# The same internally calculated vector is used as the initial
# vector for every q evaluated by the SQV procedure.


# ============================================================
# 8. Common estimation settings
# ============================================================

# Controls passed to optim() for MLqE, MDPDE, and MLE.
optimization_control <- list(
  maxit = 2000,
  reltol = 1e-12
)

# Numerical differentiation settings. The analytical
# Delta_jacobian is supplied in this example, but numerical
# differentiation is still used internally for the tau
# transformation Jacobians.
jacobian_method <- "Richardson"
jacobian_method_args <- list()

# Confidence level used for standard errors and intervals.
conf_level <- 0.95

# Configuration of the current SQV procedure.
#
# Each q is fitted independently from the same original starting
# vector.
sqv_control <- list(
  # Number of values in the initial grid. With the endpoints
  # above, 21 values produce spacing 0.01.
  m0 = 21L,

  # Number of q values in each subsequent refinement grid.
  m = 3L,

  # Smallest q value that may be evaluated.
  q_min = 0.5,

  # Stability threshold.
  L = 0.01,

  # Retain all valid evaluated fits.
  keep_fits = TRUE
)


# ============================================================
# 9. Corrected MLqE fit with q selected by SQV
# ============================================================

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
  level = conf_level,
  method = "BFGS",
  control = optimization_control,
  use_score = TRUE,
  compute_vcov = TRUE,
  jacobian_method = jacobian_method,
  jacobian_method_args =
    jacobian_method_args,
  q_control = sqv_control
)


# ============================================================
# 10. MDPDE fit with q selected by SQV
# ============================================================

# The same data, covariance structure, automatic starting-value
# rule, optimization controls, and SQV configuration are used
# for a direct comparison with the corrected MLqE.
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
  level = conf_level,
  method = "BFGS",
  control = optimization_control,
  use_score = TRUE,
  compute_vcov = TRUE,
  jacobian_method = jacobian_method,
  jacobian_method_args =
    jacobian_method_args,
  q_control = sqv_control
)


# ============================================================
# 11. Normal maximum likelihood fit
# ============================================================

# For q = 1, both procedures reduce to normal maximum likelihood
# estimation. Only one reference fit is needed; here it is
# obtained through fit_lmm_mdpde().
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
  level = conf_level,
  method = "BFGS",
  control = optimization_control,
  use_score = TRUE,
  compute_vcov = TRUE,
  jacobian_method = jacobian_method,
  jacobian_method_args =
    jacobian_method_args
)


# ============================================================
# 12. Fit summaries
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
# 13. Selected q values and SQV information
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
# 14. Parameter estimates
# ============================================================

parameter_estimates <- cbind(
  True =
    true_theta,
  `Corrected MLqE` =
    fit_mlqe$coefficients,
  MDPDE =
    fit_mdpde$coefficients,
  MLE =
    fit_mle$coefficients
)

print(
  parameter_estimates
)

# Compare the original MLqE estimates with their consistency-
# corrected values.
mlqe_original_and_corrected <- cbind(
  `Original MLqE` =
    fit_mlqe$coefficients.star,
  `Corrected MLqE` =
    fit_mlqe$coefficients
)

print(
  mlqe_original_and_corrected
)


# ============================================================
# 15. Estimated random-effects covariance matrices
# ============================================================

true_Delta
fit_mlqe$Delta
fit_mlqe$Delta.star
fit_mdpde$Delta
fit_mle$Delta


# ============================================================
# 16. Standard errors and confidence intervals
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

print(
  standard_error_results
)

fit_mlqe$conf.int
fit_mdpde$conf.int
fit_mle$conf.int

# Original MLqE confidence intervals before the consistency
# correction.
fit_mlqe$conf.int.star


# ============================================================
# 17. Covariance and estimating-function matrices
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
# 18. Analytical derivative matrices
# ============================================================

fit_mlqe$Delta.jacobian.star
fit_mlqe$mean.jacobians.star[[1]]
fit_mlqe$covariance.jacobians.star[[1]]
fit_mlqe$F.star[[1]]

fit_mdpde$Delta.jacobian
fit_mdpde$mean.jacobians[[1]]
fit_mdpde$covariance.jacobians[[1]]
fit_mdpde$F[[1]]


# ============================================================
# 19. Fitted values, residuals, and random effects
# ============================================================

# Marginal fitted means.
fitted(
  fit_mlqe
)

fitted(
  fit_mdpde
)

fitted(
  fit_mle
)

# Conditional fitted means including predicted random effects.
fitted(
  fit_mlqe,
  type = "conditional"
)

fitted(
  fit_mdpde,
  type = "conditional"
)

fitted(
  fit_mle,
  type = "conditional"
)

# Marginal residuals.
residuals(
  fit_mlqe
)

residuals(
  fit_mdpde
)

residuals(
  fit_mle
)

# Conditional residuals.
residuals(
  fit_mlqe,
  type = "conditional"
)

residuals(
  fit_mdpde,
  type = "conditional"
)

residuals(
  fit_mle,
  type = "conditional"
)

# Predicted random effects for the contaminated group.
random_effect_results_group_30 <- cbind(
  True =
    true_random_effects[[30]],
  `Corrected MLqE` =
    fit_mlqe$random.effects[[30]],
  MDPDE =
    fit_mdpde$random.effects[[30]],
  MLE =
    fit_mle$random.effects[[30]]
)

print(
  random_effect_results_group_30
)


# ============================================================
# 20. Predictions
# ============================================================

# Marginal means for the original groups.
predict(
  fit_mlqe,
  type = "mean"
)

predict(
  fit_mdpde,
  type = "mean"
)

# Marginal covariance matrices.
predict(
  fit_mlqe,
  type = "covariance"
)

predict(
  fit_mdpde,
  type = "covariance"
)

# Predicted random effects.
predict(
  fit_mlqe,
  type = "random_effects"
)

predict(
  fit_mdpde,
  type = "random_effects"
)

# Conditional fitted means.
predict(
  fit_mlqe,
  type = "conditional_mean"
)

predict(
  fit_mdpde,
  type = "conditional_mean"
)


# ============================================================
# 21. Group-level score weights
# ============================================================

# The MLqE estimating function is evaluated in the original
# parameterization, so its density-power component uses the
# original fitted group densities.
weight_mlqe <-
  fit_mlqe$density.star^
  (
    1 -
      fit_mlqe$q
  )

# The MDPDE density-power component uses the fitted group
# densities in the standard parameterization.
weight_mdpde <-
  fit_mdpde$density^
  (
    1 -
      fit_mdpde$q
  )

score_weights <- data.frame(
  group = seq_len(n_groups),
  MLqE = weight_mlqe,
  MDPDE = weight_mdpde
)

print(
  score_weights
)


# ============================================================
# 22. Optimization and estimating-function information
# ============================================================

fit_mlqe$objective
fit_mlqe$score.star
fit_mlqe$score.corrected
fit_mlqe$convergence
fit_mlqe$message
fit_mlqe$counts

fit_mdpde$objective
fit_mdpde$score
fit_mdpde$convergence
fit_mdpde$message
fit_mdpde$counts

fit_mle$objective
fit_mle$score
fit_mle$convergence
fit_mle$message
fit_mle$counts


# ============================================================
# 23. Graphical comparison of fixed-effect means
# ============================================================

x_observed <- unlist(
  lapply(
    X,
    function(matrix_i) {
      matrix_i[, "beta1"]
    }
  ),
  use.names = FALSE
)

Y_observed <- unlist(
  Y,
  use.names = FALSE
)

x_new <- seq(
  min(x_observed),
  max(x_observed),
  length.out = 200
)

fixed_mean_true <-
  true_beta["beta0"] +
  true_beta["beta1"] *
  x_new

fixed_mean_mlqe <-
  fit_mlqe$beta["beta0"] +
  fit_mlqe$beta["beta1"] *
  x_new

fixed_mean_mdpde <-
  fit_mdpde$beta["beta0"] +
  fit_mdpde$beta["beta1"] *
  x_new

fixed_mean_mle <-
  fit_mle$beta["beta0"] +
  fit_mle$beta["beta1"] *
  x_new

plot(
  x_observed,
  Y_observed,
  pch = 16,
  cex = 0.45,
  xlab = "x",
  ylab = "Y",
  main = "True and fitted fixed-effect means"
)

lines(
  x_new,
  fixed_mean_true,
  lty = 2,
  lwd = 2
)

lines(
  x_new,
  fixed_mean_mlqe,
  lty = 1,
  lwd = 2
)

lines(
  x_new,
  fixed_mean_mdpde,
  lty = 3,
  lwd = 2
)

lines(
  x_new,
  fixed_mean_mle,
  lty = 4,
  lwd = 2
)

legend(
  "topleft",
  legend = c(
    "True fixed-effect mean",
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
# 24. Scatterplots and fitted conditional lines for the first
#     50 experimental groups
# ============================================================

# Number of groups displayed.
n_groups_plot <- min(
  50,
  length(Y),
  length(X),
  length(fit_mlqe$random.effects),
  length(fit_mdpde$random.effects),
  length(fit_mle$random.effects)
)

# Save the current graphical settings.
old_par <- par(
  no.readonly = TRUE
)

# Restore the graphical settings even if an error occurs while
# producing the plots.
tryCatch(
  {
    
    # Display 25 groups per page.
    par(
      mfrow = c(
        5,
        5
      ),
      mar = c(
        2.5,
        2.5,
        2,
        0.5
      ),
      mgp = c(
        1.5,
        0.5,
        0
      ),
      oma = c(
        0,
        0,
        2,
        0
      ),
      ask = TRUE
    )
    
    for (i in seq_len(n_groups_plot)) {
      
      # Covariate values for group i.
      x_i <- drop(
        X[[i]][, "beta1"]
      )
      
      # Observed responses for group i.
      y_i <- as.numeric(
        Y[[i]]
      )
      
      # Sort the covariate values to construct an ordered grid.
      x_grid <- seq(
        from = min(
          x_i,
          na.rm = TRUE
        ),
        to = max(
          x_i,
          na.rm = TRUE
        ),
        length.out = 100
      )
      
      # Extract the predicted random effects for group i.
      random_mlqe_i <-
        fit_mlqe$random.effects[[i]]
      
      random_mdpde_i <-
        fit_mdpde$random.effects[[i]]
      
      random_mle_i <-
        fit_mle$random.effects[[i]]
      
      # Corrected MLqE conditional fitted line:
      #
      # (beta0 + b0_i) + (beta1 + b1_i) * x.
      line_mlqe <-
        (
          fit_mlqe$beta["beta0"] +
            random_mlqe_i["random_intercept"]
        ) +
        (
          fit_mlqe$beta["beta1"] +
            random_mlqe_i["random_slope"]
        ) *
        x_grid
      
      # MDPDE conditional fitted line.
      line_mdpde <-
        (
          fit_mdpde$beta["beta0"] +
            random_mdpde_i["random_intercept"]
        ) +
        (
          fit_mdpde$beta["beta1"] +
            random_mdpde_i["random_slope"]
        ) *
        x_grid
      
      # Normal MLE conditional fitted line.
      line_mle <-
        (
          fit_mle$beta["beta0"] +
            random_mle_i["random_intercept"]
        ) +
        (
          fit_mle$beta["beta1"] +
            random_mle_i["random_slope"]
        ) *
        x_grid
      
      # Limits including the observations and all fitted lines.
      # extendrange() adds a small margin above and below.
      y_limits <- extendrange(
        range(
          y_i,
          line_mlqe,
          line_mdpde,
          line_mle,
          finite = TRUE
        ),
        f = 0.05
      )
      
      plot(
        x_i,
        y_i,
        pch = 16,
        cex = 0.7,
        col = "gray30",
        xlab = "x",
        ylab = "Y",
        ylim = y_limits,
        main = paste(
          "Group",
          i
        )
      )
      
      # Draw the MLE first:
      # black dot-dashed line.
      lines(
        x_grid,
        line_mle,
        col = adjustcolor(
          "black",
          alpha.f = 0.75
        ),
        lty = "dotdash",
        lwd = 1.4
      )
      
      # Draw the corrected MLqE second:
      # red long-dashed line.
      lines(
        x_grid,
        line_mlqe,
        col = adjustcolor(
          "red",
          alpha.f = 0.80
        ),
        lty = "longdash",
        lwd = 2.2
      )
      
      # Draw the MDPDE last:
      # blue dotted line, thinner than the MLqE line.
      lines(
        x_grid,
        line_mdpde,
        col = adjustcolor(
          "blue",
          alpha.f = 0.90
        ),
        lty = "dotted",
        lwd = 1.3
      )
      
      # Add the legend only to the first panel.
      if (i == 1) {
        legend(
          "topleft",
          legend = c(
            "Corrected MLqE",
            "MDPDE",
            "Normal MLE"
          ),
          col = c(
            "red",
            "blue",
            "black"
          ),
          lty = c(
            "longdash",
            "dotted",
            "dotdash"
          ),
          lwd = c(
            2.2,
            1.3,
            1.4
          ),
          cex = 0.62,
          bty = "n"
        )
      }
    }
    
    mtext(
      "Observed responses and fitted conditional lines",
      outer = TRUE,
      side = 3,
      line = 0.5
    )
  },
  finally = {
    
    # Restore the previous graphical settings.
    par(
      old_par
    )
  }
)


# ============================================================
# 25. Graphical comparison of group-level score weights
# ============================================================

plot(
  score_weights$group,
  score_weights$MDPDE,
  pch = 1,
  xlab = "Group",
  ylab = "Score weight",
  main = "MLqE and MDPDE group-level score weights"
)

points(
  score_weights$group,
  score_weights$MLqE,
  pch = 3
)

abline(
  h = 0,
  lty = 2
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


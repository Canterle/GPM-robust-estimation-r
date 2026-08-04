# ============================================================
# Complete simulated example for fit_lmm_t()
#
# Linear mixed model with a marginal multivariate Student-t
# distribution
#
# For group i:
#
#   Y_i ~ t_{d_i}(X_i beta, Sigma_i, nu),
#
# where:
#
#   Sigma_i = Z_i Delta(gamma) Z_i^T + sigma2 I_{d_i}.
#
# The simulation uses exactly the same parameter values and the
# same direct unstructured Delta(gamma) parameterization used in
# the Gaussian MLqE and MDPDE examples:
#
#   beta   = (18, 10)
#   gamma  = (3, -2, 4)
#   sigma2 = 5
#   nu     = 4
#
# with
#
#   Delta(gamma) =
#
#     [ gamma1  gamma2 ]
#     [ gamma2  gamma3 ].
#
# The only change in the data generation is that Y_i is generated
# from a multivariate Student-t distribution with four degrees of
# freedom instead of a multivariate normal distribution.
#
# The script illustrates how to:
#
#   1. simulate grouped data from the Student-t mixed model;
#   2. define Delta(gamma) and its analytical Jacobian;
#   3. calculate starting values internally with start = NULL;
#   4. fit the model by maximum likelihood using BFGS;
#   5. inspect estimates, standard errors, fitted values, random
#      effects, weights, and optimization information.
# ============================================================


# Load the Student-t maximum-likelihood implementation.
source("fit_lmm_t.R")


# Packages used in this example.
required_packages <- c(
  "MASS",
  "nloptr"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1L),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0L) {
  stop(
    "Install the required packages: ",
    paste(
      missing_packages,
      collapse = ", "
    )
  )
}


# Fix the seed so the simulated example is reproducible.
set.seed(
  1234
)


# ============================================================
# 1. Simulation settings
# ============================================================

# Number of independent groups.
n_groups <- 5000

# Number of observations in every group.
group_size <- 5

group_sizes <- rep(
  group_size,
  n_groups
)

# Fixed Student-t degrees of freedom.
nu <- 4


# ============================================================
# 2. True parameter values
# ============================================================

# Fixed-effect parameters.
true_beta <- c(
  beta0 = 18,
  beta1 = 10
)

# Parameters of the unstructured random-effects scale matrix.
true_gamma <- c(
  gamma1 = 3,
  gamma2 = -2,
  gamma3 = 4
)

# Residual squared-scale parameter.
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
# 3. Random-effects scale structure
# ============================================================

# Construct the direct unstructured 2 by 2 matrix
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
      names(
        gamma
      )
    )
  )
}


# True random-effects scale matrix.
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
# 4. Design matrices
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
# 5. Multivariate Student-t generator
# ============================================================

# Generate one vector from
#
#   t_d(mu, Sigma, nu),
#
# where Sigma is the Student-t scale matrix.
rmvt_one <- function(mu,
                     Sigma,
                     nu) {
  mu <- as.numeric(
    mu
  )

  Sigma <- as.matrix(
    Sigma
  )

  d <- length(
    mu
  )

  L <- t(
    chol(
      Sigma
    )
  )

  z <- rnorm(
    d
  )

  w <- rchisq(
    1,
    df = nu
  )

  drop(
    mu +
      L %*%
      z /
      sqrt(
        w /
          nu
      )
  )
}


# ============================================================
# 6. Response generation
# ============================================================

Y <- vector(
  "list",
  n_groups
)

for (i in seq_len(n_groups)) {
  d_i <- group_sizes[i]

  mean_i <- drop(
    X[[i]] %*%
      true_beta
  )

  Sigma_i <-
    Z[[i]] %*%
      true_Delta %*%
      t(
        Z[[i]]
      ) +
    true_sigma2["sigma2"] *
      diag(
        d_i
      )

  Y[[i]] <- rmvt_one(
    mu = mean_i,
    Sigma = Sigma_i,
    nu = nu
  )
}


# ============================================================
# 7. Starting values
# ============================================================

# Starting values are not supplied in this example.
#
# With start = NULL, fit_lmm_t() calculates them internally:
#
#   - beta is initialized by MASS::rlm() fitted to the stacked
#     response and fixed-effect design matrices;
#   - sigma2 is initialized by the squared robust residual scale;
#   - gamma1 and gamma3 receive starting value 1 because they
#     appear on the main diagonal of gamma_structure;
#   - gamma2 receives starting value 0 because it appears only
#     outside the main diagonal.
#
# Explicit starting values remain available through:
#
#   start = c(beta0, beta1, gamma1, gamma2, gamma3, sigma2).


# ============================================================
# 8. Estimation settings
# ============================================================

# Confidence level.
conf_level <- 0.95

# Controls passed to optim() for BFGS.
optimization_control_bfgs <- list(
  maxit = 2000,
  reltol = 1e-12
)

# Numerical differentiation settings. The analytical derivative
# of Delta(gamma) is supplied in this example.
jacobian_method <- "Richardson"
jacobian_method_args <- list()


# ============================================================
# 9. Student-t maximum-likelihood fit using BFGS
# ============================================================

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
  nu = nu,
  level = conf_level,
  method = "BFGS",
  control =
    optimization_control_bfgs,
  use_score = TRUE,
  compute_vcov = TRUE,
  hessian = TRUE,
  jacobian_method =
    jacobian_method,
  jacobian_method_args =
    jacobian_method_args
)


# ============================================================
# 10. Fit summary
# ============================================================

summary(
  fit_t
)


# ============================================================
# 11. Parameter estimates
# ============================================================

parameter_estimates <- cbind(
  True =
    true_theta,
  Estimate =
    fit_t$coefficients
)

print(
  parameter_estimates
)


# ============================================================
# 12. Estimated random-effects scale matrix
# ============================================================

true_Delta
fit_t$Delta


# ============================================================
# 13. Standard errors and confidence intervals
# ============================================================

standard_error_results <- cbind(
  Estimate =
    fit_t$coefficients,
  `Std. Error` =
    fit_t$standard.error
)

print(
  standard_error_results
)

fit_t$conf.int

vcov(
  fit_t,
  type = "fisher"
)

vcov(
  fit_t,
  type = "observed"
)


# ============================================================
# 14. Score and Fisher information
# ============================================================

fit_t$score
fit_t$max.abs.score
fit_t$fisher.information
fit_t$observed.information


# ============================================================
# 15. Derivative matrices
# ============================================================

fit_t$Delta.jacobian
fit_t$mean.jacobians[[1]]
fit_t$covariance.jacobians[[1]]
fit_t$F[[1]]


# ============================================================
# 16. Fitted values, residuals, and random effects
# ============================================================

# Marginal fitted means.
fitted(
  fit_t,
  type = "marginal"
)

# Conditional fitted means.
fitted(
  fit_t,
  type = "conditional"
)

# Marginal residuals.
residuals(
  fit_t,
  type = "marginal"
)

# Conditional residuals.
residuals(
  fit_t,
  type = "conditional"
)

# Predicted random effects.
fit_t$random.effects


# ============================================================
# 17. Group-level Student-t weights
# ============================================================

student_t_weights <- data.frame(
  group = seq_len(
    n_groups
  ),
  weight = fit_t$weights
)

print(
  student_t_weights
)


# ============================================================
# 18. Predictions
# ============================================================

predict(
  fit_t,
  type = "mean"
)

predict(
  fit_t,
  type = "scale"
)

predict(
  fit_t,
  type = "random_effects"
)

predict(
  fit_t,
  type = "conditional_mean"
)


# ============================================================
# 19. Optimization information
# ============================================================

fit_t$objective
fit_t$logLik
fit_t$AIC
fit_t$BIC
fit_t$convergence
fit_t$message
fit_t$counts
fit_t$method


# ============================================================
# 20. Graphical comparison of the true and fitted fixed-effect
#     means
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
  min(
    x_observed
  ),
  max(
    x_observed
  ),
  length.out = 200
)

fixed_mean_true <-
  true_beta["beta0"] +
  true_beta["beta1"] *
  x_new

fixed_mean_fitted <-
  fit_t$beta["beta0"] +
  fit_t$beta["beta1"] *
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
  lwd = 2,
  col = "red"
)

lines(
  x_new,
  fixed_mean_fitted,
  lty = 1,
  lwd = 2,
  col = "red"
)

legend(
  "topleft",
  legend = c(
    "True fixed-effect mean",
    "Student-t MLE"
  ),
  lty = c(
    2,
    1
  ),
  lwd = 2,
  col = c(
    "red",
    "red"
  ),
  bty = "n"
)


# ============================================================
# 21. Scatterplots and fitted conditional lines for the first
#     50 experimental groups
# ============================================================

# Number of groups displayed.
n_groups_plot <- min(
  n_groups,
  length(Y),
  length(X),
  length(fit_t$random.effects)
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
      
      # Construct an ordered covariate grid.
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
      random_t_i <-
        fit_t$random.effects[[i]]
      
      # Student-t MLE conditional fitted line:
      #
      # (beta0 + b0_i) + (beta1 + b1_i) * x.
      line_t <-
        (
          fit_t$beta["beta0"] +
            random_t_i["random_intercept"]
        ) +
        (
          fit_t$beta["beta1"] +
            random_t_i["random_slope"]
        ) *
        x_grid
      
      # Limits including the observations and the fitted line.
      y_limits <- extendrange(
        range(
          y_i,
          line_t,
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
      
      # Draw the Student-t MLE conditional fitted line.
      lines(
        x_grid,
        line_t,
        col = adjustcolor(
          "black",
          alpha.f = 0.85
        ),
        lty = "solid",
        lwd = 1.8
      )
      
      # Add the legend only to the first panel.
      if (i == 1) {
        legend(
          "topleft",
          legend = "Student-t MLE",
          col = "black",
          lty = "solid",
          lwd = 1.8,
          cex = 0.68,
          bty = "n"
        )
      }
    }
    
    mtext(
      "Observed responses and Student-t conditional fitted lines",
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
# 22. Graphical comparison of group-level Student-t weights
# ============================================================

plot(
  student_t_weights$group,
  student_t_weights$weight,
  pch = 16,
  xlab = "Group",
  ylab = "Student-t weight",
  main = "Group-level Student-t weights"
)

abline(
  h = 1,
  lty = 2
)


# ============================================================
# Application to a linear measurement-error model
#
# Systolic blood pressure data: observer J versus machine S
#
# The following models and estimators are considered:
#
#   1. Student-t maximum likelihood estimator with four
#      degrees of freedom;
#   2. normal maximum density power divergence estimator
#      (MDPDE);
#   3. normal corrected maximum Lq-likelihood estimator
#      (MLqE);
#   4. normal maximum likelihood estimator using all data;
#   5. normal maximum likelihood estimator after removing
#      observations 78 and 80.
#
# Replicate measurements are replaced by their item-specific
# means. For each item and method, the measurement-error variance
# of the mean is estimated by:
#
#   sample variance of the replicates / number of replicates.
#
# Thus, tau_x and tau_y vary across items.
# ============================================================


# ------------------------------------------------------------
# 1. Load the estimation functions and data
# ------------------------------------------------------------

source("fit_lmve_mlqe.R")
source("fit_lmve_mdpde.R")
source("fit_lmve_t.R")

data("sbp", package = "MethComp")


# ------------------------------------------------------------
# 2. Select the two measurement methods
# ------------------------------------------------------------

data_x <- subset(
  sbp,
  meth == "J"
)

data_y <- subset(
  sbp,
  meth == "S"
)


# ------------------------------------------------------------
# 3. Calculate item-specific means and variances
# ------------------------------------------------------------

mean_x <- aggregate(
  y ~ item,
  data = data_x,
  FUN = mean
)

mean_y <- aggregate(
  y ~ item,
  data = data_y,
  FUN = mean
)

variance_x <- aggregate(
  y ~ item,
  data = data_x,
  FUN = var
)

variance_y <- aggregate(
  y ~ item,
  data = data_y,
  FUN = var
)

number_x <- aggregate(
  y ~ item,
  data = data_x,
  FUN = length
)

number_y <- aggregate(
  y ~ item,
  data = data_y,
  FUN = length
)

names(mean_x)[2] <- "X"
names(mean_y)[2] <- "Y"

names(variance_x)[2] <- "s2_x"
names(variance_y)[2] <- "s2_y"

names(number_x)[2] <- "m_x"
names(number_y)[2] <- "m_y"

application_data <- Reduce(
  function(x, y) {
    merge(
      x,
      y,
      by = "item"
    )
  },
  list(
    mean_x,
    mean_y,
    variance_x,
    variance_y,
    number_x,
    number_y
  )
)

# At least two replicates are required to estimate each
# item-specific variance.
application_data <- subset(
  application_data,
  m_x >= 2 &
    m_y >= 2 &
    is.finite(s2_x) &
    is.finite(s2_y)
)

application_data <- application_data[
  order(application_data$item),
]

rownames(application_data) <- NULL


# ------------------------------------------------------------
# 4. Construct the observed data and error variances
# ------------------------------------------------------------

application_data$tau_x <-
  application_data$s2_x /
  application_data$m_x

application_data$tau_y <-
  application_data$s2_y /
  application_data$m_y

Y <- application_data$Y
X <- application_data$X

tau_y <- application_data$tau_y
tau_x <- application_data$tau_x

n <- length(Y)

# Observations excluded only from the reduced normal MLE fit.
excluded_indices <- c(
  78L,
  80L
)

if (any(excluded_indices > n)) {
  stop(
    "Observations 78 and 80 are not both available after data preparation."
  )
}

included_indices <- setdiff(
  seq_len(n),
  excluded_indices
)

cat(
  "\nItem-specific measurement-error variances:\n"
)

print(
  application_data[
    ,
    c(
      "item",
      "X",
      "Y",
      "m_x",
      "m_y",
      "s2_x",
      "s2_y",
      "tau_x",
      "tau_y"
    )
  ],
  row.names = FALSE,
  digits = 6
)

cat(
  "\nObservations excluded from the reduced normal MLE fit:\n"
)

print(
  data.frame(
    index = excluded_indices,
    item = application_data$item[excluded_indices],
    Y = Y[excluded_indices],
    X = X[excluded_indices],
    tau_y = tau_y[excluded_indices],
    tau_x = tau_x[excluded_indices],
    row.names = NULL
  ),
  row.names = FALSE,
  digits = 6
)

# ------------------------------------------------------------
# 5. Auxiliary functions
# ------------------------------------------------------------

# Conditional location of the observed response given the
# observed covariate. Under the normal model this is the
# conditional mean. Under the Student-t model it is also the
# conditional mean because nu > 1.
conditional_prediction <- function(
    fit,
    X_new,
    tau_x_new
) {
  theta <- fit$coefficients

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

# Structural regression line relating the latent variables.
structural_mean <- function(
    fit,
    latent_x
) {
  theta <- fit$coefficients

  unname(
    theta["beta0"] +
      theta["beta1"] *
      latent_x
  )
}

# Mean absolute error evaluated on a selected set of indices.
mae <- function(
    observed,
    predicted,
    indices
) {
  mean(
    abs(
      observed[indices] -
        predicted[indices]
    )
  )
}


# ------------------------------------------------------------
# 6. Starting values
# ------------------------------------------------------------

# Starting values are not supplied explicitly.
#
# With start = NULL, each fitting function calculates the same
# robust moment-based starting values internally from Y, X,
# tau_y, and tau_x.


# ------------------------------------------------------------
# 7. Common optimization settings
# ------------------------------------------------------------

# The variance or squared-scale parameters are estimated on
# their original scales. L-BFGS-B is used to keep them strictly
# positive.
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

sqv_control <- list(
  m0 = 21,
  m = 3,
  q_min = 0.5,
  L = 0.01,
  verbose = FALSE,
  keep_fits = TRUE
)

conf_level <- 0.95


# ------------------------------------------------------------
# 8. Student-t maximum likelihood fit
#
# The degrees of freedom are assumed to be known and are fixed
# at nu = 4.
# ------------------------------------------------------------

nu <- 4

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


# ------------------------------------------------------------
# 9. Normal MDPDE fit
#
# The robustness parameter is selected using the SQV procedure.
# ------------------------------------------------------------

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


# ------------------------------------------------------------
# 10. Normal corrected MLqE fit
#
# The distortion parameter is selected using the SQV procedure.
# ------------------------------------------------------------

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


# ------------------------------------------------------------
# 11. Normal maximum likelihood fit using all observations
#
# For q = 1, the MDPDE reduces to maximum likelihood under the
# normal measurement-error model.
# ------------------------------------------------------------

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


# ------------------------------------------------------------
# 12. Normal maximum likelihood fit without observations 78
#     and 80
# ------------------------------------------------------------

fit_mle_wo <- fit_lmve_mdpde(
  Y = Y[included_indices],
  X = X[included_indices],
  tau_y = tau_y[included_indices],
  tau_x = tau_x[included_indices],
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
  "\nNORMAL MAXIMUM LIKELIHOOD FIT WITHOUT OBSERVATIONS 78 AND 80\n",
  paste(rep("=", 70), collapse = ""),
  "\n",
  sep = ""
)

summary(fit_mle_wo)
fit_mle_wo$conf.int

selected_q <- c(
  MLqE = fit_mlqe$q,
  MDPDE = fit_mdpde$q
)

cat(
  "\nSelected q values:\n"
)

print(selected_q)

parameter_results <- data.frame(
  parameter = names(
    fit_mlqe$coefficients
  ),
  `Student-t MLE` = as.numeric(
    fit_t$coefficients
  ),
  `Student-t SE` = as.numeric(
    fit_t$standard.error
  ),
  MDPDE = as.numeric(
    fit_mdpde$coefficients
  ),
  `MDPDE SE` = as.numeric(
    fit_mdpde$standard.error
  ),
  `Corrected MLqE` = as.numeric(
    fit_mlqe$coefficients
  ),
  `MLqE SE` = as.numeric(
    fit_mlqe$standard.error
  ),
  `Normal MLE` = as.numeric(
    fit_mle$coefficients
  ),
  `Normal MLE SE` = as.numeric(
    fit_mle$standard.error
  ),
  `Normal MLE without 78 and 80` = as.numeric(
    fit_mle_wo$coefficients
  ),
  `Reduced MLE SE` = as.numeric(
    fit_mle_wo$standard.error
  ),
  row.names = NULL,
  check.names = FALSE
)

cat(
  "\nParameter estimates and standard errors:\n"
)

print(
  parameter_results,
  row.names = FALSE,
  digits = 6
)


# ------------------------------------------------------------
# 14. Conditional fitted means and MAE
#
# Predictions from every fitted model are evaluated at all
# original observed X values, including the model fitted after
# removing observations 78 and 80.
# ------------------------------------------------------------

prediction_t <- conditional_prediction(
  fit = fit_t,
  X_new = X,
  tau_x_new = tau_x
)

prediction_mdpde <- conditional_prediction(
  fit = fit_mdpde,
  X_new = X,
  tau_x_new = tau_x
)

prediction_mlqe <- conditional_prediction(
  fit = fit_mlqe,
  X_new = X,
  tau_x_new = tau_x
)

prediction_mle <- conditional_prediction(
  fit = fit_mle,
  X_new = X,
  tau_x_new = tau_x
)

prediction_mle_wo <- conditional_prediction(
  fit = fit_mle_wo,
  X_new = X,
  tau_x_new = tau_x
)

prediction_results <- data.frame(
  index = seq_len(n),
  item = application_data$item,
  Y = Y,
  X = X,
  `Student-t MLE` = prediction_t,
  MDPDE = prediction_mdpde,
  `Corrected MLqE` = prediction_mlqe,
  `Normal MLE` = prediction_mle,
  `Normal MLE without 78 and 80` =
    prediction_mle_wo,
  row.names = NULL,
  check.names = FALSE
)

mae_results <- data.frame(
  Method = c(
    paste0(
      "Student-t MLE (nu = ",
      nu,
      ")"
    ),
    "MDPDE",
    "Corrected MLqE",
    "Normal MLE",
    "Normal MLE without observations 78 and 80"
  ),
  MAE_all_observations = c(
    mae(
      Y,
      prediction_t,
      seq_len(n)
    ),
    mae(
      Y,
      prediction_mdpde,
      seq_len(n)
    ),
    mae(
      Y,
      prediction_mlqe,
      seq_len(n)
    ),
    mae(
      Y,
      prediction_mle,
      seq_len(n)
    ),
    mae(
      Y,
      prediction_mle_wo,
      seq_len(n)
    )
  ),
  MAE_without_observations_78_80 = c(
    mae(
      Y,
      prediction_t,
      included_indices
    ),
    mae(
      Y,
      prediction_mdpde,
      included_indices
    ),
    mae(
      Y,
      prediction_mlqe,
      included_indices
    ),
    mae(
      Y,
      prediction_mle,
      included_indices
    ),
    mae(
      Y,
      prediction_mle_wo,
      included_indices
    )
  ),
  row.names = NULL,
  check.names = FALSE
)

cat(
  "\nConditional fitted means:\n"
)

print(
  prediction_results,
  row.names = FALSE,
  digits = 6
)

cat(
  "\nMean absolute errors:\n"
)

print(
  mae_results,
  row.names = FALSE,
  digits = 6
)


# ------------------------------------------------------------
# 15. Plot of the fitted structural regression lines
# ------------------------------------------------------------

latent_x_grid <- seq(
  min(X, na.rm = TRUE),
  max(X, na.rm = TRUE),
  length.out = 500
)

line_t <- structural_mean(
  fit_t,
  latent_x_grid
)

line_mdpde <- structural_mean(
  fit_mdpde,
  latent_x_grid
)

line_mlqe <- structural_mean(
  fit_mlqe,
  latent_x_grid
)

line_mle <- structural_mean(
  fit_mle,
  latent_x_grid
)

line_mle_wo <- structural_mean(
  fit_mle_wo,
  latent_x_grid
)

mlqe_legend <- bquote(
  "Normal corrected ML"[q] *
    "E (q = " *
    .(sprintf("%.2f", fit_mlqe$q)) *
    ")"
)

mdpde_legend <- bquote(
  "Normal MDPDE (q = " *
    .(sprintf("%.2f", fit_mdpde$q)) *
    ")"
)

pdf(
  file = "lmve_sbp_fitted_lines.pdf",
  width = 6.5,
  height = 3
)

par(
  mar = c(2.3, 2.3, 0.05, 0.05),
  mgp = c(1.4, 0.5, 0)
)

plot(
  X,
  Y,
  xlab = "Observer J",
  ylab = "Machine S",
  pch = 1,
  lwd = 2,
  xlim = c(80,230)
)

points(
  X[excluded_indices],
  Y[excluded_indices],
  pch = 4,
  cex = 1.2,
  lwd = 2
)

lines(
  latent_x_grid,
  line_mle,
  lty = 2,
  lwd = 2,
  col = "red"
)

lines(
  latent_x_grid,
  line_mle_wo,
  lty = 4,
  lwd = 2,
  col = "darkgreen"
)

lines(
  latent_x_grid,
  line_t,
  lty = 4,
  lwd = 2,
  col = "black"
)

lines(
  latent_x_grid,
  line_mdpde,
  lty = 5,
  lwd = 2,
  col = "purple4"
)

lines(
  latent_x_grid,
  line_mlqe,
  lty = 3,
  lwd = 2,
  col = "blue"
)

legend(
  "bottomright",
  legend = as.expression(
    list(
      "Normal MLE (all data)",
      "Normal MLE (observations 78 and 80 removed)",
      paste0(
        "Student-t MLE (nu = ",
        nu,
        ")"
      ),
      mdpde_legend,
      mlqe_legend,
      "Data",
      "Observations 78 and 80"
    )
  ),
  col = c(
    "red",
    "darkgreen",
    "black",
    "purple4",
    "blue",
    "black",
    "black"
  ),
  lty = c(
    2,
    4,
    4,
    5,
    3,
    NA,
    NA
  ),
  lwd = c(
    2,
    2,
    2,
    2,
    2,
    NA,
    NA
  ),
  pch = c(
    NA,
    NA,
    NA,
    NA,
    NA,
    1,
    4
  ),
  bty = "n",
  cex = 0.75
)

dev.off()


# ------------------------------------------------------------
# 16. Plot of the score weights
# ------------------------------------------------------------

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

weight_results <- data.frame(
  index = seq_len(n),
  item = application_data$item,
  Y = Y,
  X = X,
  MLqE = weight_mlqe,
  MDPDE = weight_mdpde,
  row.names = NULL
)

cat(
  "\nDensity-power weights:\n"
)

print(
  weight_results,
  row.names = FALSE,
  digits = 6
)

pdf(
  file = "lmve_sbp_score_weights.pdf",
  width = 6.5,
  height = 3
)

par(
  mar = c(2.3, 2.3, 0.05, 0.05),
  mgp = c(1.4, 0.5, 0)
)

weight_range <- range(
  c(
    weight_mdpde,
    weight_mlqe
  ),
  finite = TRUE
)

plot(
  seq_len(n),
  weight_mdpde,
  xlab = "Index",
  ylab = "Score weight",
  pch = 1,
  ylim = weight_range
)

points(
  seq_len(n),
  weight_mlqe,
  pch = 3
)

abline(
  h = 0,
  lty = 2,
  lwd = 1
)

legend(
  "bottomleft",
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

# ============================================================
# Application to a linear mixed model
#
# Orthodontic growth data
#
# Linear mixed model:
#
#   distance_ij = beta0 + beta1 age_c_ij
#                 + beta2 SexFemale_i
#                 + beta3 age_c_ij SexFemale_i
#                 + b0_i + b1_i age_c_ij + epsilon_ij.
#
# The following estimators are considered:
#
#   1. Student-t maximum likelihood estimator with four
#      degrees of freedom;
#   2. normal maximum density power divergence estimator
#      (MDPDE);
#   3. normal corrected maximum Lq-likelihood estimator
#      (MLqE);
#   4. normal maximum likelihood estimator using all subjects;
#   5. normal maximum likelihood estimator after removing
#      subjects M09 and M13.
#
# Normal maximum-likelihood fits from lme4 are also calculated
# as reference fits.
#
# The script produces:
#
#   1. fit summaries, standard errors, and confidence intervals;
#   2. conditional fitted values and mean absolute errors;
#   3. an overall fixed-effect trajectory for all observations;
#   4. fixed-effect trajectories for males and females;
#   5. subject-specific conditional trajectories that include
#      the predicted random effects;
#   6. subject-level robustness weights.
# ============================================================


# ------------------------------------------------------------
# 1. Load packages, estimation functions, and data
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(nlme)
  library(lme4)
})

source("fit_lmm_t.R")
source("fit_lmm_mdpde.R")
source("fit_lmm_mlqe.R")

data(
  "Orthodont",
  package = "nlme"
)


# ------------------------------------------------------------
# 2. Data preparation
# ------------------------------------------------------------

orthodont <- as.data.frame(
  Orthodont
)

# Use Male as the reference category.
orthodont$Sex <- factor(
  orthodont$Sex,
  levels = c(
    "Male",
    "Female"
  )
)

# Center age at 11 years.
orthodont$age_c <- orthodont$age - 11

# Keep observations ordered within subjects.
orthodont <- orthodont[
  order(
    orthodont$Subject,
    orthodont$age
  ),
]

rownames(orthodont) <- NULL

# Subjects excluded only from the reduced normal MLE fits.
subjects_removed <- c(
  "M09",
  "M13"
)

orthodont_reduced <- droplevels(
  orthodont[
    !orthodont$Subject %in%
      subjects_removed,
  ]
)


# ------------------------------------------------------------
# 3. Auxiliary functions
# ------------------------------------------------------------

# Fixed-effect mean for a specified value of the female
# indicator. The indicator may also be a number in [0, 1], which
# gives a weighted average of the male and female trajectories.
fitted_fixed_mean <- function(
    beta,
    age_c,
    female
) {
  beta["(Intercept)"] +
    beta["age_c"] *
      age_c +
    beta["SexFemale"] *
      female +
    beta["age_c:SexFemale"] *
      age_c *
      female
}

# Standardize a two-dimensional predicted random-effect vector.
standardize_random_effect <- function(
    random_effect
) {
  result <- as.numeric(
    random_effect
  )

  if (length(result) != 2L) {
    stop(
      "Each predicted random-effect vector must have length two."
    )
  }

  names(result) <- c(
    "(Intercept)",
    "age_c"
  )

  result
}

# Obtain the predicted random effects for a subject. A zero
# vector is returned when the subject was not included in the
# fitted sample, as occurs for M09 and M13 under the reduced MLE.
subject_random_effect <- function(
    fit,
    subject
) {
  random_effects <- fit$random.effects

  if (is.null(names(random_effects))) {
    stop(
      "The fitted random effects must be named by subject."
    )
  }

  if (!subject %in% names(random_effects)) {
    return(
      c(
        `(Intercept)` = 0,
        age_c = 0
      )
    )
  }

  standardize_random_effect(
    random_effects[[subject]]
  )
}

# Conditional fitted trajectory for one subject.
fitted_conditional_mean <- function(
    beta,
    random_effect,
    age_c,
    female
) {
  fitted_fixed_mean(
    beta = beta,
    age_c = age_c,
    female = female
  ) +
    random_effect["(Intercept)"] +
    random_effect["age_c"] *
      age_c
}

# Conditional fitted values for all observations. For a subject
# not present in the fitted sample, only the fixed effects are
# used because subject_random_effect() returns zero.
conditional_predictions <- function(
    fit,
    X_list,
    Z_list
) {
  subject_names <- names(X_list)

  if (is.null(subject_names)) {
    stop(
      "The grouped design matrices must be named by subject."
    )
  }

  unlist(
    lapply(
      subject_names,
      function(subject) {
        beta <- fit$beta
        random_effect <- subject_random_effect(
          fit = fit,
          subject = subject
        )

        drop(
          X_list[[subject]] %*%
            beta +
            Z_list[[subject]] %*%
              random_effect
        )
      }
    ),
    use.names = FALSE
  )
}

# Mean absolute error evaluated on selected observations.
mae <- function(
    observed,
    predicted,
    indices = seq_along(observed)
) {
  mean(
    abs(
      observed[indices] -
        predicted[indices]
    )
  )
}


# ------------------------------------------------------------
# 4. Grouped response and design matrices
# ------------------------------------------------------------

X_complete <- model.matrix(
  ~ age_c * Sex,
  data = orthodont
)

Z_complete <- model.matrix(
  ~ age_c,
  data = orthodont
)

subject_indices <- split(
  seq_len(
    nrow(orthodont)
  ),
  orthodont$Subject
)

Y <- lapply(
  subject_indices,
  function(index) {
    as.numeric(
      orthodont$distance[index]
    )
  }
)

X <- lapply(
  subject_indices,
  function(index) {
    X_complete[
      index,
      ,
      drop = FALSE
    ]
  }
)

Z <- lapply(
  subject_indices,
  function(index) {
    Z_complete[
      index,
      ,
      drop = FALSE
    ]
  }
)

# Preserve the subject names in every grouped object.
names(Y) <- names(subject_indices)
names(X) <- names(subject_indices)
names(Z) <- names(subject_indices)

keep_subjects <- !names(Y) %in%
  subjects_removed

keep_observations <- !orthodont$Subject %in%
  subjects_removed


# ------------------------------------------------------------
# 5. Random-effects covariance structure
# ------------------------------------------------------------

# Direct unstructured covariance matrix:
#
#   Delta(gamma) =
#
#     [ gamma1  gamma2 ]
#     [ gamma2  gamma3 ].
Delta_unstructured <- function(
    gamma
) {
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
        "(Intercept)",
        "age_c"
      ),
      c(
        "(Intercept)",
        "age_c"
      )
    )
  )
}

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
      "(Intercept)",
      "age_c"
    ),
    c(
      "(Intercept)",
      "age_c"
    )
  )
)

# Analytical Jacobian of vec(Delta(gamma)) with respect to
# gamma. R vectorizes matrices column by column.
Delta_unstructured_jacobian <- function(
    gamma
) {
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


# ------------------------------------------------------------
# 6. Consistency correction for the MLqE
# ------------------------------------------------------------

# The fixed effects remain unchanged. The random-effects
# covariance parameters and residual variance are divided by r.
tau <- function(
    theta,
    r
) {
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


# ------------------------------------------------------------
# 7. Starting values
# ------------------------------------------------------------

# Starting values are not supplied explicitly.
#
# With start = NULL, the fitting functions calculate the same
# robust starting values internally from Y, X, Z, and the
# supplied random-effects covariance structure.


# ------------------------------------------------------------
# 8. Common optimization settings
# ------------------------------------------------------------

bobyqa_control <- list(
  maxeval = 100000,
  xtol_rel = 0,
  ftol_rel = 0,
  xtol_abs = 1e-8,
  ftol_abs = 1e-8
)

jacobian_method <- "Richardson"
jacobian_method_args <- list()

conf_level <- 0.95

sqv_control <- list(
  m0 = 21,
  m = 3,
  q_min = 0.5,
  L = 0.01,
  keep_fits = TRUE
)

# Fixed degrees of freedom used in every Student-t calculation.
nu <- 4


# ------------------------------------------------------------
# 9. Student-t maximum likelihood fit
# ------------------------------------------------------------

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
  method = "BOBYQA",
  control = bobyqa_control,
  use_score = FALSE,
  compute_vcov = TRUE,
  hessian = TRUE
)


# ------------------------------------------------------------
# 10. Normal MDPDE fit
# ------------------------------------------------------------

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
  method = "BOBYQA",
  control = bobyqa_control,
  use_score = FALSE,
  compute_vcov = TRUE,
  jacobian_method = jacobian_method,
  jacobian_method_args =
    jacobian_method_args,
  q_control = sqv_control
)


# ------------------------------------------------------------
# 11. Normal corrected MLqE fit
# ------------------------------------------------------------

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
  method = "BOBYQA",
  control = bobyqa_control,
  use_score = FALSE,
  compute_vcov = TRUE,
  jacobian_method = jacobian_method,
  jacobian_method_args =
    jacobian_method_args,
  q_control = sqv_control
)


# ------------------------------------------------------------
# 12. Normal maximum likelihood fit using all subjects
# ------------------------------------------------------------

# For q = 1, the corrected MLqE reduces to normal maximum
# likelihood estimation.
fit_mle <- fit_lmm_mlqe(
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
  q = 1,
  level = conf_level,
  method = "BOBYQA",
  control = bobyqa_control,
  use_score = FALSE,
  compute_vcov = TRUE,
  jacobian_method = jacobian_method,
  jacobian_method_args =
    jacobian_method_args
)


# ------------------------------------------------------------
# 13. Normal maximum likelihood fit without M09 and M13
# ------------------------------------------------------------

fit_mle_wo <- fit_lmm_mlqe(
  Y = Y[keep_subjects],
  X = X[keep_subjects],
  Z = Z[keep_subjects],
  Delta = Delta_unstructured,
  tau = tau,
  start = NULL,
  gamma_structure =
    gamma_structure_unstructured,
  Delta_jacobian =
    Delta_unstructured_jacobian,
  q = 1,
  level = conf_level,
  method = "BOBYQA",
  control = bobyqa_control,
  use_score = FALSE,
  compute_vcov = TRUE,
  jacobian_method = jacobian_method,
  jacobian_method_args =
    jacobian_method_args
)


# ------------------------------------------------------------
# 14. Reference normal maximum-likelihood fits from lme4
# ------------------------------------------------------------

fit_lmer <- lmer(
  distance ~ age_c * Sex +
    (1 + age_c | Subject),
  data = orthodont,
  REML = FALSE,
  control = lmerControl(
    optimizer = "nloptwrap",
    optCtrl = list(
      algorithm = "NLOPT_LN_BOBYQA",
      xtol_abs = 1e-8,
      ftol_abs = 1e-8,
      maxeval = 100000
    )
  )
)

fit_lmer_wo <- lmer(
  distance ~ age_c * Sex +
    (1 + age_c | Subject),
  data = orthodont_reduced,
  REML = FALSE,
  control = lmerControl(
    optimizer = "nloptwrap",
    optCtrl = list(
      algorithm = "NLOPT_LN_BOBYQA",
      xtol_abs = 1e-8,
      ftol_abs = 1e-8,
      maxeval = 100000
    )
  )
)


# ------------------------------------------------------------
# 15. Fit summaries and confidence intervals
# ------------------------------------------------------------

# Ensure that the predicted random effects are named by subject.
names(fit_t$random.effects) <- names(Y)
names(fit_mdpde$random.effects) <- names(Y)
names(fit_mlqe$random.effects) <- names(Y)
names(fit_mle$random.effects) <- names(Y)
names(fit_mle_wo$random.effects) <- names(Y[keep_subjects])

cat(
  "\n",
  paste(rep("=", 70), collapse = ""),
  "\nSTUDENT-t MAXIMUM LIKELIHOOD FIT (nu = 4)\n",
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
  "\nNORMAL MAXIMUM LIKELIHOOD FIT USING ALL SUBJECTS\n",
  paste(rep("=", 70), collapse = ""),
  "\n",
  sep = ""
)

summary(fit_mle)
fit_mle$conf.int

cat(
  "\n",
  paste(rep("=", 70), collapse = ""),
  "\nNORMAL MAXIMUM LIKELIHOOD FIT WITHOUT M09 AND M13\n",
  paste(rep("=", 70), collapse = ""),
  "\n",
  sep = ""
)

summary(fit_mle_wo)
fit_mle_wo$conf.int

cat(
  "\n",
  paste(rep("=", 70), collapse = ""),
  "\nlme4 REFERENCE FIT USING ALL SUBJECTS\n",
  paste(rep("=", 70), collapse = ""),
  "\n",
  sep = ""
)

summary(fit_lmer)

cat(
  "\n",
  paste(rep("=", 70), collapse = ""),
  "\nlme4 REFERENCE FIT WITHOUT M09 AND M13\n",
  paste(rep("=", 70), collapse = ""),
  "\n",
  sep = ""
)

summary(fit_lmer_wo)

selected_tuning_parameters <- c(
  `Student-t degrees of freedom` = nu,
  `Corrected MLqE q` = fit_mlqe$q,
  `MDPDE q` = fit_mdpde$q
)

cat(
  "\nSelected tuning parameters:\n"
)

print(
  selected_tuning_parameters
)

parameter_results <- data.frame(
  parameter = names(
    fit_mle$coefficients
  ),
  `Student-t MLE` = as.numeric(
    fit_t$coefficients
  ),
  `Student-t MLE SE` = as.numeric(
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
  `Corrected MLqE SE` = as.numeric(
    fit_mlqe$standard.error
  ),
  `Normal MLE` = as.numeric(
    fit_mle$coefficients
  ),
  `Normal MLE SE` = as.numeric(
    fit_mle$standard.error
  ),
  `Normal MLE without M09 and M13` = as.numeric(
    fit_mle_wo$coefficients
  ),
  `Normal MLE without M09 and M13 SE` = as.numeric(
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
# 16. Conditional fitted values and mean absolute errors
# ------------------------------------------------------------

pred_t <- conditional_predictions(
  fit = fit_t,
  X_list = X,
  Z_list = Z
)

pred_mdpde <- conditional_predictions(
  fit = fit_mdpde,
  X_list = X,
  Z_list = Z
)

pred_mlqe <- conditional_predictions(
  fit = fit_mlqe,
  X_list = X,
  Z_list = Z
)

pred_mle <- conditional_predictions(
  fit = fit_mle,
  X_list = X,
  Z_list = Z
)

# M09 and M13 receive fixed-effect predictions because they were
# not included in the reduced fit.
pred_mle_wo <- conditional_predictions(
  fit = fit_mle_wo,
  X_list = X,
  Z_list = Z
)

observed_distance <- orthodont$distance

mae_results <- data.frame(
  Method = c(
    "Student-t MLE (nu = 4)",
    "MDPDE",
    "Corrected MLqE",
    "Normal MLE",
    "Normal MLE without M09 and M13"
  ),
  MAE_all_observations = c(
    mae(observed_distance, pred_t),
    mae(observed_distance, pred_mdpde),
    mae(observed_distance, pred_mlqe),
    mae(observed_distance, pred_mle),
    mae(observed_distance, pred_mle_wo)
  ),
  MAE_without_M09_M13 = c(
    mae(
      observed_distance,
      pred_t,
      which(keep_observations)
    ),
    mae(
      observed_distance,
      pred_mdpde,
      which(keep_observations)
    ),
    mae(
      observed_distance,
      pred_mlqe,
      which(keep_observations)
    ),
    mae(
      observed_distance,
      pred_mle,
      which(keep_observations)
    ),
    mae(
      observed_distance,
      pred_mle_wo,
      which(keep_observations)
    )
  ),
  row.names = NULL,
  check.names = FALSE
)

cat(
  "\nMean absolute errors based on conditional fitted values:\n"
)

print(
  mae_results,
  row.names = FALSE,
  digits = 6
)

fitted_results <- data.frame(
  Subject = orthodont$Subject,
  Sex = orthodont$Sex,
  age = orthodont$age,
  age_c = orthodont$age_c,
  observed = observed_distance,
  Student_t_MLE = pred_t,
  MDPDE = pred_mdpde,
  Corrected_MLqE = pred_mlqe,
  Normal_MLE = pred_mle,
  Normal_MLE_without_M09_M13 = pred_mle_wo,
  row.names = NULL,
  check.names = FALSE
)

print(
  head(fitted_results)
)


# ------------------------------------------------------------
# 17. Graphical settings and legends
# ------------------------------------------------------------

beta_t <- fit_t$beta
beta_mdpde <- fit_mdpde$beta
beta_mlqe <- fit_mlqe$beta
beta_mle <- fit_mle$beta
beta_mle_wo <- fit_mle_wo$beta

age_grid <- seq(
  min(
    orthodont$age_c
  ),
  max(
    orthodont$age_c
  ),
  length.out = 200
)

female_proportion <- mean(
  orthodont$Sex == "Female"
)

mlqe_legend <- bquote(
  "Normal corrected ML"[q] *
    "E (q = " *
    .(
      sprintf(
        "%.2f",
        fit_mlqe$q
      )
    ) *
    ")"
)

mdpde_legend <- bquote(
  "Normal MDPDE (q = " *
    .(
      sprintf(
        "%.2f",
        fit_mdpde$q
      )
    ) *
    ")"
)

method_labels <- as.expression(
  list(
    "Normal MLE (all data)",
    "Normal MLE (M09 and M13 removed)",
    paste0(
      "Student-t MLE (nu = ",
      nu,
      ")"
    ),
    mdpde_legend,
    mlqe_legend,
    "Data"
  )
)

method_colors <- c(
  "red",
  "darkgreen",
  "black",
  "purple4",
  "blue",
  "black"
)

method_line_types <- c(
  2,
  4,
  4,
  5,
  3,
  NA
)

method_line_widths <- c(
  2,
  2,
  2,
  2,
  2,
  NA
)

method_point_types <- c(
  NA,
  NA,
  NA,
  NA,
  NA,
  1
)

add_fixed_effect_lines <- function(
    female
) {
  lines(
    age_grid,
    fitted_fixed_mean(
      beta = beta_mle,
      age_c = age_grid,
      female = female
    ),
    col = method_colors[1],
    lty = method_line_types[1],
    lwd = method_line_widths[1]
  )
  
  lines(
    age_grid,
    fitted_fixed_mean(
      beta = beta_mle_wo,
      age_c = age_grid,
      female = female
    ),
    col = method_colors[2],
    lty = method_line_types[2],
    lwd = method_line_widths[2]
  )
  
  lines(
    age_grid,
    fitted_fixed_mean(
      beta = beta_t,
      age_c = age_grid,
      female = female
    ),
    col = method_colors[3],
    lty = method_line_types[3],
    lwd = method_line_widths[3]
  )
  
  lines(
    age_grid,
    fitted_fixed_mean(
      beta = beta_mdpde,
      age_c = age_grid,
      female = female
    ),
    col = method_colors[4],
    lty = method_line_types[4],
    lwd = method_line_widths[4]
  )
  
  lines(
    age_grid,
    fitted_fixed_mean(
      beta = beta_mlqe,
      age_c = age_grid,
      female = female
    ),
    col = method_colors[5],
    lty = method_line_types[5],
    lwd = method_line_widths[5]
  )
}

add_method_legend <- function(
    location = "bottomright",
    cex = 0.58
) {
  legend(
    location,
    legend = method_labels,
    col = method_colors,
    lty = method_line_types,
    lwd = method_line_widths,
    pch = method_point_types,
    bty = "n",
    cex = cex,
    y.intersp = 0.82,
    x.intersp = 0.75,
    seg.len = 2.2,
    pt.cex = 0.9
  )
}


# ------------------------------------------------------------
# 18. Overall fixed-effect trajectories
# ------------------------------------------------------------

pdf(
  file = "lmm_orthodont_fitted_lines_overall.pdf",
  width = 6.5,
  height = 3,
  useDingbats = FALSE
)

par(
  mar = c(
    2.3,
    2.3,
    0.05,
    0.05
  ),
  mgp = c(
    1.4,
    0.5,
    0
  )
)

plot(
  orthodont$age_c,
  orthodont$distance,
  xlab = "Centered age",
  ylab = "Distance",
  pch = 1,
  lwd = 2,
  xlim = range(
    orthodont$age_c,
    finite = TRUE
  ),
  ylim = range(
    orthodont$distance,
    finite = TRUE
  )
)

add_fixed_effect_lines(
  female = female_proportion
)

add_method_legend(
  location = "topleft",
  cex = 0.55
)

dev.off()


# ------------------------------------------------------------
# 19. Fixed-effect trajectories for male subjects
# ------------------------------------------------------------

male_data <- orthodont[
  orthodont$Sex == "Male",
]

pdf(
  file = "lmm_orthodont_fitted_lines_male.pdf",
  width = 6.5,
  height = 3,
  useDingbats = FALSE
)

par(
  mar = c(
    2.3,
    2.3,
    0.05,
    0.05
  ),
  mgp = c(
    1.4,
    0.5,
    0
  )
)

plot(
  male_data$age_c,
  male_data$distance,
  xlab = "Centered age",
  ylab = "Distance",
  pch = 1,
  lwd = 2,
  xlim = range(
    orthodont$age_c,
    finite = TRUE
  ),
  ylim = range(
    orthodont$distance,
    finite = TRUE
  )
)

add_fixed_effect_lines(
  female = 0
)

add_method_legend(
  location = "topleft",
  cex = 0.55
)

dev.off()


# ------------------------------------------------------------
# 20. Fixed-effect trajectories for female subjects
# ------------------------------------------------------------

female_data <- orthodont[
  orthodont$Sex == "Female",
]

pdf(
  file = "lmm_orthodont_fitted_lines_female.pdf",
  width = 6.5,
  height = 3,
  useDingbats = FALSE
)

par(
  mar = c(
    2.3,
    2.3,
    0.05,
    0.05
  ),
  mgp = c(
    1.4,
    0.5,
    0
  )
)

plot(
  female_data$age_c,
  female_data$distance,
  xlab = "Centered age",
  ylab = "Distance",
  pch = 1,
  lwd = 2,
  xlim = range(
    orthodont$age_c,
    finite = TRUE
  ),
  ylim = range(
    orthodont$distance,
    finite = TRUE
  )
)

add_fixed_effect_lines(
  female = 1
)

add_method_legend(
  location = "topleft",
  cex = 0.55
)

dev.off()


# ------------------------------------------------------------
# 21. Subject-specific conditional trajectories
# ------------------------------------------------------------

n_subjects_plot <- length(
  subject_indices
)

pdf(
  file =
    "lmm_orthodont_subject_specific_conditional_trajectories.pdf",
  width = 15.5,
  height = 12.2,
  pointsize = 10,
  useDingbats = FALSE
)

par(
  mfrow = c(
    7,
    4
  ),
  mar = c(
    2.55,
    3.00,
    0.10,
    0.10
  ),
  mgp = c(
    1.60,
    0.68,
    0
  ),
  oma = c(
    2.80,
    3.50,
    0.10,
    0.10
  ),
  tcl = -0.28
)

for (i in seq_len(
  n_subjects_plot
)) {
  subject <- names(
    subject_indices
  )[i]
  
  subject_rows <- subject_indices[[i]]
  
  subject_sex <- as.character(
    orthodont$Sex[
      subject_rows[1]
    ]
  )
  
  female <- as.numeric(
    subject_sex == "Female"
  )
  
  age_c_i <- orthodont$age_c[
    subject_rows
  ]
  
  y_i <- orthodont$distance[
    subject_rows
  ]
  
  age_c_subject_grid <- seq(
    min(
      age_c_i,
      na.rm = TRUE
    ),
    max(
      age_c_i,
      na.rm = TRUE
    ),
    length.out = 100
  )
  
  line_mle <- fitted_conditional_mean(
    beta = beta_mle,
    random_effect = subject_random_effect(
      fit_mle,
      subject
    ),
    age_c = age_c_subject_grid,
    female = female
  )
  
  line_mle_wo <- fitted_conditional_mean(
    beta = beta_mle_wo,
    random_effect = subject_random_effect(
      fit_mle_wo,
      subject
    ),
    age_c = age_c_subject_grid,
    female = female
  )
  
  line_t <- fitted_conditional_mean(
    beta = beta_t,
    random_effect = subject_random_effect(
      fit_t,
      subject
    ),
    age_c = age_c_subject_grid,
    female = female
  )
  
  line_mdpde <- fitted_conditional_mean(
    beta = beta_mdpde,
    random_effect = subject_random_effect(
      fit_mdpde,
      subject
    ),
    age_c = age_c_subject_grid,
    female = female
  )
  
  line_mlqe <- fitted_conditional_mean(
    beta = beta_mlqe,
    random_effect = subject_random_effect(
      fit_mlqe,
      subject
    ),
    age_c = age_c_subject_grid,
    female = female
  )
  
  y_limits <- extendrange(
    range(
      y_i,
      line_mle,
      line_mle_wo,
      line_t,
      line_mdpde,
      line_mlqe,
      finite = TRUE
    ),
    f = 0.06
  )
  
  x_ticks <- pretty(
    range(
      age_c_i,
      finite = TRUE
    ),
    n = 3
  )
  
  x_ticks <- x_ticks[
    x_ticks >= min(
      age_c_i,
      na.rm = TRUE
    ) &
      x_ticks <= max(
        age_c_i,
        na.rm = TRUE
      )
  ]
  
  y_ticks <- pretty(
    y_limits,
    n = 3
  )
  
  y_ticks <- y_ticks[
    y_ticks >= y_limits[1] &
      y_ticks <= y_limits[2]
  ]
  
  plot(
    age_c_i,
    y_i,
    xlab = "",
    ylab = "",
    pch = 1,
    cex = 1.00,
    lwd = 1.4,
    xlim = range(
      age_c_i,
      finite = TRUE
    ),
    ylim = y_limits,
    axes = FALSE,
    frame.plot = FALSE
  )
  
  axis(
    side = 1,
    at = x_ticks,
    labels = x_ticks,
    cex.axis = 1.75,
    lwd = 1,
    lwd.ticks = 1
  )
  
  axis(
    side = 2,
    at = y_ticks,
    labels = y_ticks,
    las = 1,
    cex.axis = 1.75,
    lwd = 1,
    lwd.ticks = 1
  )
  
  box(
    lwd = 1
  )
  
  lines(
    age_c_subject_grid,
    line_mle,
    col = method_colors[1],
    lty = method_line_types[1],
    lwd = 2
  )
  
  lines(
    age_c_subject_grid,
    line_mle_wo,
    col = method_colors[2],
    lty = method_line_types[2],
    lwd = 2
  )
  
  lines(
    age_c_subject_grid,
    line_t,
    col = method_colors[3],
    lty = method_line_types[3],
    lwd = 2
  )
  
  lines(
    age_c_subject_grid,
    line_mdpde,
    col = method_colors[4],
    lty = method_line_types[4],
    lwd = 2
  )
  
  lines(
    age_c_subject_grid,
    line_mlqe,
    col = method_colors[5],
    lty = method_line_types[5],
    lwd = 2
  )
  
  plotting_region <- par(
    "usr"
  )
  
  text(
    x =
      plotting_region[1] +
      0.035 *
      diff(
        plotting_region[1:2]
      ),
    y =
      plotting_region[4] -
      0.055 *
      diff(
        plotting_region[3:4]
      ),
    labels = subject,
    adj = c(
      0,
      1
    ),
    cex = 1.75,
    font = 2
  )
}



plot.new()

legend(
  "center",
  legend = method_labels,
  col = method_colors,
  lty = method_line_types,
  lwd = c(
    rep(
      2,
      5
    ),
    NA
  ),
  pch = method_point_types,
  bty = "n",
  cex = 2.00,
  y.intersp = 0.95,
  x.intersp = 0.75,
  seg.len = 2.0,
  pt.cex = 1.00
)



mtext(
  "Centered age",
  side = 1,
  outer = TRUE,
  line = 1.35,
  cex = 1.75
)

mtext(
  "Distance",
  side = 2,
  outer = TRUE,
  line = 1.45,
  cex = 1.75
)

dev.off()


# ------------------------------------------------------------
# 22. Subject-level robustness weights
# ------------------------------------------------------------

weight_mlqe <- if (
  fit_mlqe$q == 1
) {
  rep(
    1,
    length(
      subject_indices
    )
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
    length(
      subject_indices
    )
  )
} else {
  fit_mdpde$density^(
    1 -
      fit_mdpde$q
  )
}

group_score_weights <- data.frame(
  Subject = names(
    subject_indices
  ),
  MLqE = weight_mlqe,
  MDPDE = weight_mdpde,
  row.names = NULL,
  check.names = FALSE
)

group_score_weights_ordered <-
  group_score_weights[
    order(
      pmin(
        group_score_weights$MLqE,
        group_score_weights$MDPDE
      )
    ),
  ]

cat(
  "\nSubject-level density-power weights:\n"
)

print(
  group_score_weights_ordered,
  row.names = FALSE,
  digits = 6
)

subject_position <- seq_len(
  nrow(
    group_score_weights
  )
)

pdf(
  file = "lmm_orthodont_score_weights.pdf",
  width = 6.5,
  height = 3,
  useDingbats = FALSE
)

par(
  mar = c(
    2.3,
    2.3,
    0.05,
    0.05
  ),
  mgp = c(
    1.4,
    0.5,
    0
  )
)

weight_range <- range(
  c(
    weight_mdpde,
    weight_mlqe
  ),
  finite = TRUE
)

plot(
  subject_position,
  weight_mdpde,
  xlab = "Subject",
  ylab = "Score weight",
  pch = 1,
  xaxt = "n",
  ylim = weight_range
)

points(
  subject_position,
  weight_mlqe,
  pch = 3
)

axis(
  side = 1,
  at = subject_position,
  labels = group_score_weights$Subject,
  las = 2,
  cex.axis = 0.7
)

abline(
  h = 0,
  lty = 2,
  lwd = 1
)

legend(
  "bottomright",
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
  bty = "n",
  cex = 0.78
)

dev.off()
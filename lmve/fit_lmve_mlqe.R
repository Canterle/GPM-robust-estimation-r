# Generic MLqE fitting for a linear normal measurement-error model
#
# The latent structural model is
#
#   Y_i^* = beta0 + beta1 X_i^* + epsilon_i,
#
# where
#
#   X_i^* ~ Normal(mu_x, sigma2_x),
#   epsilon_i ~ Normal(0, sigma2).
#
# The observed variables are
#
#   Y_i = Y_i^* + e_yi,
#   X_i = X_i^* + e_xi,
#
# with known measurement-error variances
#
#   Var(e_yi) = tau_yi,
#   Var(e_xi) = tau_xi.
#
# The corrected parameter vector is
#
#   theta = (beta0, beta1, mu_x, sigma2_x, sigma2).
#
# For MLqE, the original reparameterized vector is
#
#   theta_star = tau_{1/q}(theta)
#              = (beta0, beta1, mu_x, q * sigma2_x, q * sigma2).
#
# The known measurement-error variances are also reparameterized:
#
#   tau_yi_star = q * tau_yi,
#   tau_xi_star = q * tau_xi.
#
# Consequently, the original covariance matrix used in the
# Lq-likelihood is q times the corrected covariance matrix.
#
# The transformation family is internal:
#
#   tau_r(theta) =
#     (beta0, beta1, mu_x, sigma2_x / r, sigma2 / r).
#
# Therefore:
#
#   theta_star = tau_{1/q}(theta),
#   theta      = tau_q(theta_star).
#
# The Jacobian of tau_{1/q}(theta) is
#
#   A = diag(1, 1, 1, q, q),
#
# and the corrected covariance matrix is
#
#   V = A^{-1} V_star A^{-T}.


# ============================================================
# Internal utility functions
# ============================================================

.lmve_mlqe_parameter_names <- c(
  "beta0",
  "beta1",
  "mu_x",
  "sigma2_x",
  "sigma2"
)


.lmve_mlqe_validate_start <- function(start) {
  if (!is.numeric(start) || length(start) != 5L) {
    stop(
      "'start' must be a numeric vector of length 5 containing ",
      "beta0, beta1, mu_x, sigma2_x and sigma2."
    )
  }

  if (is.null(names(start))) {
    names(start) <- .lmve_mlqe_parameter_names
  }

  missing_names <- setdiff(
    .lmve_mlqe_parameter_names,
    names(start)
  )

  if (length(missing_names) > 0L) {
    stop(
      "'start' must contain the following named parameters: ",
      paste(.lmve_mlqe_parameter_names, collapse = ", "),
      "."
    )
  }

  start <- start[.lmve_mlqe_parameter_names]

  start <- setNames(
    as.numeric(start),
    .lmve_mlqe_parameter_names
  )

  if (anyNA(start) || any(!is.finite(start))) {
    stop(
      "'start' must contain only finite, non-missing values."
    )
  }

  if (start["sigma2_x"] <= 0 ||
      start["sigma2"] <= 0) {
    stop(
      "The starting values of 'sigma2_x' and 'sigma2' ",
      "must be strictly positive."
    )
  }

  start
}


.lmve_mlqe_prepare_data <- function(Y,
                                    X,
                                    tau_y,
                                    tau_x) {
  Y <- as.numeric(Y)
  X <- as.numeric(X)

  n <- length(Y)

  if (n == 0L) {
    stop("'Y' must not be empty.")
  }

  if (length(X) != n) {
    stop("'Y' and 'X' must have the same length.")
  }

  tau_y <- rep_len(
    as.numeric(tau_y),
    n
  )

  tau_x <- rep_len(
    as.numeric(tau_x),
    n
  )

  if (anyNA(Y) ||
      anyNA(X) ||
      anyNA(tau_y) ||
      anyNA(tau_x) ||
      any(!is.finite(Y)) ||
      any(!is.finite(X)) ||
      any(!is.finite(tau_y)) ||
      any(!is.finite(tau_x))) {
    stop(
      "'Y', 'X', 'tau_y' and 'tau_x' must contain only ",
      "finite, non-missing values."
    )
  }

  if (any(tau_y < 0) ||
      any(tau_x < 0)) {
    stop(
      "The known measurement-error variances 'tau_y' and ",
      "'tau_x' must be non-negative."
    )
  }

  list(
    Y = Y,
    X = X,
    tau_y = tau_y,
    tau_x = tau_x,
    n = n
  )
}


# Compute the corrected starting values when start = NULL.
.lmve_mlqe_default_start <- function(data) {
  n <- data$n
  X <- data$X
  Y <- data$Y
  tau_x <- data$tau_x
  tau_y <- data$tau_y

  Mxy <- sum(
    (X - median(X)) * Y
  ) / (n - 1)

  Mx <- mad(X)^2
  My <- mad(Y)^2

  Xbar <- median(X)
  Ybar <- median(Y)

  tau_xbar <- median(tau_x)
  tau_ybar <- median(tau_y)

  denominator <- Mx - tau_xbar

  if (!is.finite(denominator) ||
      denominator == 0) {
    stop(
      "The automatic starting values could not be computed because ",
      "mad(X)^2 - median(tau_x) is zero or non-finite. ",
      "Supply 'start' explicitly."
    )
  }

  b1ini <- Mxy / denominator

  betaini <- c(
    beta0 = Ybar - Xbar * b1ini,
    beta1 = b1ini
  )

  mu_xini <- Xbar
  sigma2_xini <- Mx - tau_xbar

  sigma2ini <- My -
    Mxy^2 / denominator -
    tau_ybar

  # Protect the automatic variance starting values.
  # Whenever an initial variance is smaller than 0.05,
  # replace it with 0.1 before the reparameterization.
  if (!is.finite(sigma2_xini) ||
      sigma2_xini < 0.05) {
    sigma2_xini <- 0.1
  }

  if (!is.finite(sigma2ini) ||
      sigma2ini < 0.05) {
    sigma2ini <- 0.1
  }

  start <- c(
    betaini,
    mu_x = mu_xini,
    sigma2_x = sigma2_xini,
    sigma2 = sigma2ini
  )

  if (anyNA(start) ||
      any(!is.finite(start))) {
    stop(
      "The automatic starting-value calculation produced non-finite ",
      "values. Supply 'start' explicitly."
    )
  }

  setNames(
    as.numeric(start),
    .lmve_mlqe_parameter_names
  )
}


# Internal transformation family:
#
# tau_r(theta) =
#   (beta0, beta1, mu_x, sigma2_x / r, sigma2 / r).
.lmve_mlqe_tau <- function(theta,
                           r) {
  theta <- setNames(
    as.numeric(theta),
    .lmve_mlqe_parameter_names
  )

  if (length(r) != 1L ||
      !is.numeric(r) ||
      !is.finite(r) ||
      r <= 0) {
    stop("'r' must be a strictly positive finite number.")
  }

  result <- theta

  result[c("sigma2_x", "sigma2")] <-
    result[c("sigma2_x", "sigma2")] / r

  result
}


# A = d tau_{1/q}(theta) / d theta.
.lmve_mlqe_tau_jacobian <- function(q) {
  result <- diag(
    c(
      1,
      1,
      1,
      q,
      q
    )
  )

  dimnames(result) <- list(
    .lmve_mlqe_parameter_names,
    .lmve_mlqe_parameter_names
  )

  result
}


# A^{-1} = d tau_q(theta_star) / d theta_star.
.lmve_mlqe_tau_jacobian_inverse <- function(q) {
  result <- diag(
    c(
      1,
      1,
      1,
      1 / q,
      1 / q
    )
  )

  dimnames(result) <- list(
    .lmve_mlqe_parameter_names,
    .lmve_mlqe_parameter_names
  )

  result
}


.lmve_mlqe_transform_bounds <- function(bounds,
                                        q,
                                        label) {
  p <- length(.lmve_mlqe_parameter_names)

  if (!is.null(names(bounds))) {
    missing_names <- setdiff(
      .lmve_mlqe_parameter_names,
      names(bounds)
    )

    if (length(missing_names) > 0L) {
      stop(
        "'", label, "' must contain all parameter names."
      )
    }

    bounds <- bounds[.lmve_mlqe_parameter_names]
  }

  bounds <- setNames(
    rep_len(
      as.numeric(bounds),
      p
    ),
    .lmve_mlqe_parameter_names
  )

  bounds[c("sigma2_x", "sigma2")] <-
    q * bounds[c("sigma2_x", "sigma2")]

  bounds
}


.lmve_mlqe_model_components <- function(theta,
                                        data,
                                        lambda_scale = 1,
                                        strict = FALSE) {
  theta <- setNames(
    as.numeric(theta),
    .lmve_mlqe_parameter_names
  )

  beta0 <- theta["beta0"]
  beta1 <- theta["beta1"]
  mu_x <- theta["mu_x"]
  sigma2_x <- theta["sigma2_x"]
  sigma2 <- theta["sigma2"]

  invalid_variances <-
    !is.finite(sigma2_x) ||
    !is.finite(sigma2) ||
    sigma2_x <= 0 ||
    sigma2 <= 0

  if (invalid_variances) {
    if (strict) {
      stop(
        "'sigma2_x' and 'sigma2' must be strictly positive."
      )
    }

    return(NULL)
  }

  if (length(lambda_scale) != 1L ||
      !is.numeric(lambda_scale) ||
      !is.finite(lambda_scale) ||
      lambda_scale <= 0) {
    stop("'lambda_scale' must be strictly positive.")
  }

  mean_y <- beta0 + beta1 * mu_x
  mean_x <- mu_x

  variance_y <-
    beta1^2 * sigma2_x +
    sigma2 +
    lambda_scale * data$tau_y

  covariance_yx <- rep(
    beta1 * sigma2_x,
    data$n
  )

  variance_x <-
    sigma2_x +
    lambda_scale * data$tau_x

  determinant <-
    variance_y * variance_x -
    covariance_yx^2

  valid <-
    all(is.finite(variance_y)) &&
    all(is.finite(covariance_yx)) &&
    all(is.finite(variance_x)) &&
    all(is.finite(determinant)) &&
    all(variance_y > 0) &&
    all(variance_x > 0) &&
    all(determinant > 0)

  if (!valid) {
    if (strict) {
      stop(
        "At least one fitted covariance matrix is not ",
        "positive definite."
      )
    }

    return(NULL)
  }

  residual_y <- data$Y - mean_y
  residual_x <- data$X - mean_x

  # Analytical bivariate normal calculation. This is equivalent
  # to a Cholesky-based lMvn()/dMvn() calculation, but avoids one
  # matrix decomposition for each observation.
  quadratic_form <-
    (
      variance_x * residual_y^2 -
      2 * covariance_yx * residual_y * residual_x +
      variance_y * residual_x^2
    ) / determinant

  log_density <-
    -log(2 * pi) -
    0.5 * log(determinant) -
    0.5 * quadratic_form

  if (any(!is.finite(log_density))) {
    if (strict) {
      stop(
        "The model produced non-finite bivariate normal ",
        "log-density values."
      )
    }

    return(NULL)
  }

  fitted_values <- cbind(
    Y = rep(mean_y, data$n),
    X = rep(mean_x, data$n)
  )

  residuals <- cbind(
    Y = residual_y,
    X = residual_x
  )

  covariance_array <- array(
    NA_real_,
    dim = c(2L, 2L, data$n),
    dimnames = list(
      c("Y", "X"),
      c("Y", "X"),
      NULL
    )
  )

  covariance_array[1L, 1L, ] <- variance_y
  covariance_array[1L, 2L, ] <- covariance_yx
  covariance_array[2L, 1L, ] <- covariance_yx
  covariance_array[2L, 2L, ] <- variance_x

  list(
    mean = c(
      Y = mean_y,
      X = mean_x
    ),
    fitted.values = fitted_values,
    fitted.y = fitted_values[, "Y"],
    fitted.x = fitted_values[, "X"],
    residuals = residuals,
    residual.y = residual_y,
    residual.x = residual_x,
    variance.y = variance_y,
    variance.x = variance_x,
    covariance.yx = covariance_yx,
    sigma.y = sqrt(variance_y),
    sigma.x = sqrt(variance_x),
    covariance.array = covariance_array,
    determinant = determinant,
    log_density = log_density,
    density = exp(log_density)
  )
}


# Derivatives with respect to the original parameterization theta_star.
.lmve_mlqe_derivative_matrices_star <- function(theta_star) {
  theta_star <- setNames(
    as.numeric(theta_star),
    .lmve_mlqe_parameter_names
  )

  beta1 <- theta_star["beta1"]
  mu_x <- theta_star["mu_x"]
  sigma2_x_star <- theta_star["sigma2_x"]

  mean_jacobian <- matrix(
    c(
      1,
      mu_x,
      beta1,
      0,
      0,
      0,
      0,
      1,
      0,
      0
    ),
    nrow = 2L,
    byrow = TRUE,
    dimnames = list(
      c("mean_Y", "mean_X"),
      .lmve_mlqe_parameter_names
    )
  )

  covariance_jacobian <- matrix(
    c(
      0,
      2 * sigma2_x_star * beta1,
      0,
      beta1^2,
      1,
      0,
      sigma2_x_star,
      0,
      beta1,
      0,
      0,
      sigma2_x_star,
      0,
      beta1,
      0,
      0,
      0,
      0,
      1,
      0
    ),
    nrow = 4L,
    byrow = TRUE,
    dimnames = list(
      c(
        "Sigma_11",
        "Sigma_21",
        "Sigma_12",
        "Sigma_22"
      ),
      .lmve_mlqe_parameter_names
    )
  )

  F <- rbind(
    mean_jacobian,
    covariance_jacobian
  )

  list(
    mean = mean_jacobian,
    covariance = covariance_jacobian,
    F = F
  )
}


.lmve_mlqe_block_matrix <- function(top_left,
                                    bottom_right) {
  top_dimension <- nrow(top_left)
  bottom_dimension <- nrow(bottom_right)

  result <- matrix(
    0,
    nrow = top_dimension + bottom_dimension,
    ncol = top_dimension + bottom_dimension
  )

  result[
    seq_len(top_dimension),
    seq_len(top_dimension)
  ] <- top_left

  result[
    top_dimension + seq_len(bottom_dimension),
    top_dimension + seq_len(bottom_dimension)
  ] <- bottom_right

  result
}


.lmve_mlqe_inverse_covariance <- function(variance_y,
                                          covariance_yx,
                                          variance_x,
                                          determinant) {
  matrix(
    c(
      variance_x,
      -covariance_yx,
      -covariance_yx,
      variance_y
    ),
    nrow = 2L
  ) / determinant
}


.lmve_mlqe_invert_matrix <- function(matrix_value,
                                     label) {
  inverse <- tryCatch(
    solve(matrix_value),
    error = function(e) NULL
  )

  if (is.null(inverse)) {
    inverse <- tryCatch(
      qr.solve(matrix_value),
      error = function(e) {
        stop(
          label,
          " could not be inverted: ",
          conditionMessage(e)
        )
      }
    )
  }

  inverse
}


# ============================================================
# Fixed-q MLqE fit
# ============================================================

.fit_lmve_mlqe_fixed <- function(Y,
                                 X,
                                 tau_y,
                                 tau_x,
                                 start = NULL,
                                 q = 1,
                                 level = 0.95,
                                 method = c(
                                   "BFGS",
                                   "L-BFGS-B",
                                   "Nelder-Mead"
                                 ),
                                 lower = NULL,
                                 upper = NULL,
                                 control = list(),
                                 use_score = FALSE,
                                 compute_vcov = TRUE) {
  call <- match.call()
  method <- match.arg(method)

  data <- .lmve_mlqe_prepare_data(
    Y = Y,
    X = X,
    tau_y = tau_y,
    tau_x = tau_x
  )

  automatic_start <- is.null(start)

  start_corrected <- if (automatic_start) {
    .lmve_mlqe_default_start(data)
  } else {
    .lmve_mlqe_validate_start(start)
  }

  p <- length(start_corrected)

  if (length(q) != 1L ||
      !is.numeric(q) ||
      !is.finite(q) ||
      q <= 0 ||
      q > 1) {
    stop("'q' must be in the interval (0, 1].")
  }

  if (length(level) != 1L ||
      !is.numeric(level) ||
      !is.finite(level) ||
      level <= 0 ||
      level >= 1) {
    stop("'level' must be in the interval (0, 1).")
  }

  if (!is.list(control)) {
    stop("'control' must be a list.")
  }

  # Reparameterize the starting values:
  # theta_star = tau_{1/q}(theta).
  start_star <- .lmve_mlqe_tau(
    start_corrected,
    r = 1 / q
  )

  tau_jacobian <- .lmve_mlqe_tau_jacobian(q)
  tau_jacobian_inverse <-
    .lmve_mlqe_tau_jacobian_inverse(q)

  # The direct correction Jacobian equals A^{-1}.
  tau_q_jacobian <- tau_jacobian_inverse

  objective <- function(theta_star) {
    theta_star <- setNames(
      as.numeric(theta_star),
      .lmve_mlqe_parameter_names
    )

    sigma2_x_star <- theta_star["sigma2_x"]
    sigma2_star <- theta_star["sigma2"]

    if (!is.finite(sigma2_x_star) ||
        !is.finite(sigma2_star) ||
        sigma2_x_star <= 0 ||
        sigma2_star <= 0) {
      return(NaN)
    }

    # The known measurement-error variances are multiplied by q
    # in the original reparameterized model.
    components_star <- .lmve_mlqe_model_components(
      theta = theta_star,
      data = data,
      lambda_scale = q
    )

    if (is.null(components_star)) {
      return(NaN)
    }

    if (q == 1) {
      value <- sum(
        components_star$log_density
      )
    } else {
      one_minus_q <- 1 - q

      value <- sum(
        expm1(
          one_minus_q *
          components_star$log_density
        ) /
        one_minus_q
      )
    }

    if (is.finite(value)) {
      value
    } else {
      NaN
    }
  }

  score_star <- function(theta_star,
                         strict = FALSE) {
    theta_star <- setNames(
      as.numeric(theta_star),
      .lmve_mlqe_parameter_names
    )

    components_star <- .lmve_mlqe_model_components(
      theta = theta_star,
      data = data,
      lambda_scale = q,
      strict = strict
    )

    if (is.null(components_star)) {
      return(
        setNames(
          rep(NaN, p),
          .lmve_mlqe_parameter_names
        )
      )
    }

    derivatives <- .lmve_mlqe_derivative_matrices_star(
      theta_star
    )

    F <- derivatives$F

    density_power <- if (q == 1) {
      rep(1, data$n)
    } else {
      exp(
        (1 - q) *
        components_star$log_density
      )
    }

    score_sum <- rep(0, p)

    for (i in seq_len(data$n)) {
      Sigma_inverse <- .lmve_mlqe_inverse_covariance(
        variance_y = components_star$variance.y[i],
        covariance_yx = components_star$covariance.yx[i],
        variance_x = components_star$variance.x[i],
        determinant = components_star$determinant[i]
      )

      H <- .lmve_mlqe_block_matrix(
        top_left = Sigma_inverse,
        bottom_right =
          0.5 *
          kronecker(
            Sigma_inverse,
            Sigma_inverse
          )
      )

      z <- c(
        components_star$residual.y[i],
        components_star$residual.x[i]
      )

      centered_second_moment <-
        tcrossprod(z) -
        matrix(
          c(
            components_star$variance.y[i],
            components_star$covariance.yx[i],
            components_star$covariance.yx[i],
            components_star$variance.x[i]
          ),
          nrow = 2L
        )

      s <- c(
        z,
        c(centered_second_moment)
      )

      score_sum <- score_sum +
        as.numeric(
          t(F) %*%
          H %*%
          s
        ) *
        density_power[i]
    }

    setNames(
      score_sum,
      .lmve_mlqe_parameter_names
    )
  }

  compute_JK <- function(theta_star) {
    components_star <- .lmve_mlqe_model_components(
      theta = theta_star,
      data = data,
      lambda_scale = q,
      strict = TRUE
    )

    derivatives <- .lmve_mlqe_derivative_matrices_star(
      theta_star
    )

    F <- derivatives$F
    one_minus_q <- 1 - q

    J_star <- matrix(
      0,
      nrow = p,
      ncol = p
    )

    K_star <- matrix(
      0,
      nrow = p,
      ncol = p
    )

    for (i in seq_len(data$n)) {
      Sigma_inverse <- .lmve_mlqe_inverse_covariance(
        variance_y = components_star$variance.y[i],
        covariance_yx = components_star$covariance.yx[i],
        variance_x = components_star$variance.x[i],
        determinant = components_star$determinant[i]
      )

      H <- .lmve_mlqe_block_matrix(
        top_left = Sigma_inverse,
        bottom_right =
          0.5 *
          kronecker(
            Sigma_inverse,
            Sigma_inverse
          )
      )

      c1 <- q *
        exp(
          -one_minus_q * log(2 * pi) -
          0.5 * one_minus_q *
            log(components_star$determinant[i])
        )

      J_star <- J_star +
        c1 *
        t(F) %*%
        H %*%
        F

      Sigma_inverse_vector <- c(
        Sigma_inverse
      )

      H_q <- .lmve_mlqe_block_matrix(
        top_left = Sigma_inverse,
        bottom_right =
          (
            0.5 *
            kronecker(
              Sigma_inverse,
              Sigma_inverse
            ) +
            0.25 *
            one_minus_q^2 *
            tcrossprod(
              Sigma_inverse_vector
            )
          ) /
          (2 - q)
      )

      c2 <- c1^2 /
        (
          q *
          (2 - q)
        )

      K_star <- K_star +
        c2 *
        t(F) %*%
        H_q %*%
        F
    }

    J_star <- q * J_star
    K_star <- K_star / (2 - q)

    J_star <- (
      J_star +
      t(J_star)
    ) / 2

    K_star <- (
      K_star +
      t(K_star)
    ) / 2

    dimnames(J_star) <- list(
      .lmve_mlqe_parameter_names,
      .lmve_mlqe_parameter_names
    )

    dimnames(K_star) <- list(
      .lmve_mlqe_parameter_names,
      .lmve_mlqe_parameter_names
    )

    list(
      J = J_star,
      K = K_star,
      mean.jacobian = derivatives$mean,
      covariance.jacobian =
        derivatives$covariance,
      F = F
    )
  }

  .lmve_mlqe_model_components(
    theta = start_star,
    data = data,
    lambda_scale = q,
    strict = TRUE
  )

  if (method == "L-BFGS-B") {
    default_control <- list(
      maxit = 1000,
      factr = 1e7,
      pgtol = 0
    )

    if (!is.null(control$reltol)) {
      if (is.null(control$factr)) {
        control$factr <-
          control$reltol /
          .Machine$double.eps
      }

      control$reltol <- NULL
    }

    control$abstol <- NULL
  } else {
    default_control <- list(
      maxit = 1000,
      reltol = 1e-9
    )

    control$factr <- NULL
    control$pgtol <- NULL
  }

  control <- modifyList(
    default_control,
    control
  )

  if (method == "L-BFGS-B") {
    control$reltol <- NULL
    control$abstol <- NULL
  } else {
    control$factr <- NULL
    control$pgtol <- NULL
  }

  control$fnscale <- -1

  optim_arguments <- list(
    par = start_star,
    fn = objective,
    method = method,
    control = control
  )

  if (use_score) {
    if (method == "Nelder-Mead") {
      warning(
        "'use_score' is ignored when method = 'Nelder-Mead'."
      )
    } else {
      optim_arguments$gr <- score_star
    }
  }

  if (method == "L-BFGS-B") {
    default_lower <- c(
      beta0 = -Inf,
      beta1 = -Inf,
      mu_x = -Inf,
      sigma2_x = 1e-10,
      sigma2 = 1e-10
    )

    default_upper <- setNames(
      rep(Inf, p),
      .lmve_mlqe_parameter_names
    )

    if (is.null(lower)) {
      lower <- default_lower
    }

    if (is.null(upper)) {
      upper <- default_upper
    }

    if (!is.null(names(lower))) {
      lower <- lower[.lmve_mlqe_parameter_names]
    }

    if (!is.null(names(upper))) {
      upper <- upper[.lmve_mlqe_parameter_names]
    }

    lower_corrected <- setNames(
      rep_len(as.numeric(lower), p),
      .lmve_mlqe_parameter_names
    )

    upper_corrected <- setNames(
      rep_len(as.numeric(upper), p),
      .lmve_mlqe_parameter_names
    )

    if (anyNA(lower_corrected) ||
        anyNA(upper_corrected)) {
      stop(
        "'lower' and 'upper' must not contain missing values."
      )
    }

    if (any(lower_corrected > upper_corrected)) {
      stop(
        "Each lower bound must be less than or equal ",
        "to its upper bound."
      )
    }

    if (any(start_corrected < lower_corrected) ||
        any(start_corrected > upper_corrected)) {
      stop(
        "All corrected starting values must lie within ",
        "the specified bounds."
      )
    }

    # User-supplied bounds refer to the corrected parameters.
    # They are transformed internally to the original scale.
    lower_star <- .lmve_mlqe_transform_bounds(
      lower_corrected,
      q = q,
      label = "lower"
    )

    upper_star <- .lmve_mlqe_transform_bounds(
      upper_corrected,
      q = q,
      label = "upper"
    )

    optim_arguments$lower <- lower_star
    optim_arguments$upper <- upper_star
  }

  optimization <- do.call(
    optim,
    optim_arguments
  )

  estimate_star <- setNames(
    as.numeric(optimization$par),
    .lmve_mlqe_parameter_names
  )

  # Consistency correction:
  # theta_hat = tau_q(theta_star_hat).
  estimate <- .lmve_mlqe_tau(
    estimate_star,
    r = q
  )

  mapped_back <- .lmve_mlqe_tau(
    estimate,
    r = 1 / q
  )

  tau_composition_residual <- max(
    abs(mapped_back - estimate_star) /
    pmax(1, abs(estimate_star))
  )

  fitted_components_star <-
    .lmve_mlqe_model_components(
      theta = estimate_star,
      data = data,
      lambda_scale = q,
      strict = TRUE
    )

  fitted_components <-
    .lmve_mlqe_model_components(
      theta = estimate,
      data = data,
      lambda_scale = 1,
      strict = TRUE
    )

  score_at_estimate_star <- tryCatch(
    score_star(
      estimate_star,
      strict = TRUE
    ),
    error = function(e) {
      warning(
        "The score at the original MLqE estimate could not be computed: ",
        conditionMessage(e)
      )

      setNames(
        rep(NA_real_, p),
        .lmve_mlqe_parameter_names
      )
    }
  )

  score_at_estimate <- as.numeric(
    t(tau_jacobian) %*%
    score_at_estimate_star
  )

  score_at_estimate <- setNames(
    score_at_estimate,
    .lmve_mlqe_parameter_names
  )

  vcov_star <- NULL
  vcov_matrix <- NULL
  vcov_direct <- NULL
  vcov_direct_difference <- NULL
  standard_error_star <- NULL
  standard_error <- NULL
  confidence_interval_star <- NULL
  confidence_interval <- NULL
  J_star <- NULL
  K_star <- NULL
  mean_jacobian_star <- NULL
  covariance_jacobian_star <- NULL
  F_star <- NULL
  vcov_warning <- NULL

  if (compute_vcov) {
    jk <- compute_JK(estimate_star)

    J_star <- jk$J
    K_star <- jk$K
    mean_jacobian_star <-
      jk$mean.jacobian
    covariance_jacobian_star <-
      jk$covariance.jacobian
    F_star <- jk$F

    covariance_result <- tryCatch(
      {
        J_inverse <- chol2inv(
          chol(J_star)
        )

        covariance_star <-
          J_inverse %*%
          K_star %*%
          t(J_inverse)

        covariance_star <- (
          covariance_star +
          t(covariance_star)
        ) / 2

        list(
          value = covariance_star,
          warning = NULL
        )
      },
      error = function(cholesky_error) {
        tryCatch(
          {
            J_inverse <- solve(J_star)

            covariance_star <-
              J_inverse %*%
              K_star %*%
              t(J_inverse)

            covariance_star <- (
              covariance_star +
              t(covariance_star)
            ) / 2

            list(
              value = covariance_star,
              warning = paste(
                "The Cholesky factorization of J_star failed; ",
                "solve(J_star) was used instead:",
                conditionMessage(cholesky_error)
              )
            )
          },
          error = function(inverse_error) {
            list(
              value = matrix(
                NA_real_,
                nrow = p,
                ncol = p,
                dimnames = list(
                  .lmve_mlqe_parameter_names,
                  .lmve_mlqe_parameter_names
                )
              ),
              warning = paste(
                "The original covariance matrix could not be computed:",
                conditionMessage(inverse_error)
              )
            )
          }
        )
      }
    )

    vcov_star <- covariance_result$value

    dimnames(vcov_star) <- list(
      .lmve_mlqe_parameter_names,
      .lmve_mlqe_parameter_names
    )

    vcov_warning <- covariance_result$warning

    if (all(is.finite(vcov_star))) {
      vcov_matrix <-
        tau_jacobian_inverse %*%
        vcov_star %*%
        t(tau_jacobian_inverse)
    } else {
      vcov_matrix <- matrix(
        NA_real_,
        nrow = p,
        ncol = p
      )
    }

    vcov_matrix <- (
      vcov_matrix +
      t(vcov_matrix)
    ) / 2

    dimnames(vcov_matrix) <- list(
      .lmve_mlqe_parameter_names,
      .lmve_mlqe_parameter_names
    )

    vcov_direct <-
      tau_q_jacobian %*%
      vcov_star %*%
      t(tau_q_jacobian)

    vcov_direct <- (
      vcov_direct +
      t(vcov_direct)
    ) / 2

    dimnames(vcov_direct) <- list(
      .lmve_mlqe_parameter_names,
      .lmve_mlqe_parameter_names
    )

    vcov_direct_difference <- max(
      abs(
        vcov_matrix -
        vcov_direct
      )
    )

    diagonal_star <- diag(vcov_star)
    diagonal_corrected <- diag(vcov_matrix)

    standard_error_star <- setNames(
      ifelse(
        is.na(diagonal_star) |
        diagonal_star < 0,
        NA_real_,
        sqrt(diagonal_star)
      ),
      .lmve_mlqe_parameter_names
    )

    standard_error <- setNames(
      ifelse(
        is.na(diagonal_corrected) |
        diagonal_corrected < 0,
        NA_real_,
        sqrt(diagonal_corrected)
      ),
      .lmve_mlqe_parameter_names
    )

    z_value <- qnorm(
      1 - (1 - level) / 2
    )

    confidence_interval_star <- cbind(
      lower =
        estimate_star -
        z_value * standard_error_star,
      upper =
        estimate_star +
        z_value * standard_error_star
    )

    confidence_interval <- cbind(
      lower =
        estimate -
        z_value * standard_error,
      upper =
        estimate +
        z_value * standard_error
    )

    rownames(confidence_interval_star) <-
      .lmve_mlqe_parameter_names

    rownames(confidence_interval) <-
      .lmve_mlqe_parameter_names
  }

  result <- list(
    call = call,
    starting.values = start_corrected,
    starting.values.star = start_star,
    automatic.start = automatic_start,
    coefficients = estimate,
    coefficients.star = estimate_star,
    standard.error = standard_error,
    standard.error.star = standard_error_star,
    vcov = vcov_matrix,
    vcov.star = vcov_star,
    vcov.direct = vcov_direct,
    vcov.direct.difference =
      vcov_direct_difference,
    conf.int = confidence_interval,
    conf.int.star =
      confidence_interval_star,
    J = J_star,
    K = K_star,
    J.star = J_star,
    K.star = K_star,
    mean.jacobian.star =
      mean_jacobian_star,
    covariance.jacobian.star =
      covariance_jacobian_star,
    F.star = F_star,
    tau.jacobian = tau_jacobian,
    tau.jacobian.inverse =
      tau_jacobian_inverse,
    tau.q.jacobian =
      tau_q_jacobian,
    tau.composition.residual =
      tau_composition_residual,
    fitted.values =
      fitted_components$fitted.values,
    fitted.y =
      fitted_components$fitted.y,
    fitted.x =
      fitted_components$fitted.x,
    covariance.array =
      fitted_components$covariance.array,
    variance.y =
      fitted_components$variance.y,
    variance.x =
      fitted_components$variance.x,
    covariance.yx =
      fitted_components$covariance.yx,
    sigma.y =
      fitted_components$sigma.y,
    sigma.x =
      fitted_components$sigma.x,
    residuals =
      fitted_components$residuals,
    residual.y =
      fitted_components$residual.y,
    residual.x =
      fitted_components$residual.x,
    log.density =
      fitted_components$log_density,
    density =
      fitted_components$density,
    fitted.values.star =
      fitted_components_star$fitted.values,
    fitted.y.star =
      fitted_components_star$fitted.y,
    fitted.x.star =
      fitted_components_star$fitted.x,
    covariance.array.star =
      fitted_components_star$covariance.array,
    variance.y.star =
      fitted_components_star$variance.y,
    variance.x.star =
      fitted_components_star$variance.x,
    covariance.yx.star =
      fitted_components_star$covariance.yx,
    sigma.y.star =
      fitted_components_star$sigma.y,
    sigma.x.star =
      fitted_components_star$sigma.x,
    residuals.star =
      fitted_components_star$residuals,
    residual.y.star =
      fitted_components_star$residual.y,
    residual.x.star =
      fitted_components_star$residual.x,
    log.density.star =
      fitted_components_star$log_density,
    density.star =
      fitted_components_star$density,
    objective =
      unname(optimization$value),
    score = score_at_estimate_star,
    score.star =
      score_at_estimate_star,
    score.corrected =
      score_at_estimate,
    convergence =
      optimization$convergence,
    message =
      optimization$message,
    counts =
      optimization$counts,
    q = q,
    level = level,
    nobs = data$n,
    method = method,
    Y = data$Y,
    X = data$X,
    tau_y = data$tau_y,
    tau_x = data$tau_x,
    tau_y.star = q * data$tau_y,
    tau_x.star = q * data$tau_x,
    tau = .lmve_mlqe_tau,
    vcov.warning =
      vcov_warning,
    optimization.control =
      control,
    optim = optimization
  )

  class(result) <- "lmve_mlqe"

  result
}


# ============================================================
# Public fitting interface
# ============================================================

fit_lmve_mlqe <- function(Y,
                          X,
                          tau_y,
                          tau_x,
                          start = NULL,
                          q = 1,
                          level = 0.95,
                          method = c(
                            "BFGS",
                            "L-BFGS-B",
                            "Nelder-Mead"
                          ),
                          lower = NULL,
                          upper = NULL,
                          control = list(),
                          use_score = FALSE,
                          compute_vcov = TRUE,
                          q_control = list()) {
  call <- match.call()

  # Support q as the fifth positional argument when start is omitted.
  if (!missing(start) &&
      missing(q) &&
      (
        (
          is.character(start) &&
          length(start) == 1L &&
          identical(toupper(start), "SQV")
        ) ||
        (
          is.numeric(start) &&
          length(start) == 1L &&
          is.finite(start) &&
          start > 0 &&
          start <= 1
        )
      )) {
    q <- start
    start <- NULL
  }

  sqv_requested <-
    is.character(q) &&
    length(q) == 1L &&
    identical(toupper(q), "SQV")

  if (is.character(q) &&
      !sqv_requested) {
    stop(
      "'q' must be numeric or equal to 'SQV'."
    )
  }

  if (!sqv_requested) {
    fit <- .fit_lmve_mlqe_fixed(
      Y = Y,
      X = X,
      tau_y = tau_y,
      tau_x = tau_x,
      start = start,
      q = q,
      level = level,
      method = method,
      lower = lower,
      upper = upper,
      control = control,
      use_score = use_score,
      compute_vcov = compute_vcov
    )

    fit$call <- call

    return(fit)
  }

  if (!isTRUE(compute_vcov)) {
    stop(
      "SQV selection requires 'compute_vcov = TRUE' because ",
      "corrected standard errors are used to construct the ",
      "standardized estimates."
    )
  }

  fit <- .select_q_lmve_mlqe_sqv(
    Y = Y,
    X = X,
    tau_y = tau_y,
    tau_x = tau_x,
    start = start,
    level = level,
    method = method,
    lower = lower,
    upper = upper,
    control = control,
    use_score = use_score,
    q_control = q_control
  )

  fit$call <- call

  fit
}


# ============================================================
# SQV selection algorithm
# ============================================================

.select_q_lmve_mlqe_sqv <- function(Y,
                                     X,
                                     tau_y,
                                     tau_x,
                                     start = NULL,
                                     level,
                                     method,
                                     lower,
                                     upper,
                                     control,
                                     use_score,
                                     q_control) {
  if (!is.list(q_control)) {
    stop("'q_control' must be a list.")
  }

  method <- match.arg(
    method,
    c(
      "BFGS",
      "L-BFGS-B",
      "Nelder-Mead"
    )
  )

  automatic_start <- is.null(start)

  prepared_data <- .lmve_mlqe_prepare_data(
    Y = Y,
    X = X,
    tau_y = tau_y,
    tau_x = tau_x
  )

  start_corrected <- if (automatic_start) {
    .lmve_mlqe_default_start(
      prepared_data
    )
  } else {
    .lmve_mlqe_validate_start(
      start
    )
  }

  default_q_control <- list(
    m0 = 21L,
    m = 3L,
    q_min = 0.5,
    L = 0.01,
    verbose = FALSE,
    keep_fits = FALSE
  )

  q_control <- modifyList(
    default_q_control,
    q_control
  )

  # Retain only the controls used by the original SQV algorithm
  # and the two output controls.
  q_control <- q_control[
    names(default_q_control)
  ]

  scalar_numeric <- function(x) {
    is.numeric(x) &&
      length(x) == 1L &&
      !is.na(x) &&
      is.finite(x)
  }

  scalar_logical <- function(x) {
    is.logical(x) &&
      length(x) == 1L &&
      !is.na(x)
  }

  integer_control <- function(x,
                              minimum,
                              name) {
    if (!scalar_numeric(x) ||
        x < minimum ||
        abs(x - round(x)) >
          sqrt(.Machine$double.eps)) {
      stop(
        "'", name,
        "' must be an integer greater than or equal to ",
        minimum,
        "."
      )
    }

    as.integer(round(x))
  }

  m0 <- integer_control(
    q_control$m0,
    2L,
    "q_control$m0"
  )

  m <- integer_control(
    q_control$m,
    2L,
    "q_control$m"
  )

  if (m > m0) {
    stop(
      "'q_control$m' cannot exceed 'q_control$m0' because the ",
      "restart endpoint is determined by the initial grid."
    )
  }

  if (!scalar_numeric(q_control$q_min) ||
      q_control$q_min <= 0 ||
      q_control$q_min >= 0.8) {
    stop(
      "'q_control$q_min' must be positive and smaller than 0.8."
    )
  }

  if (!scalar_numeric(q_control$L) ||
      q_control$L <= 0) {
    stop(
      "'q_control$L' must be positive."
    )
  }

  if (!scalar_logical(q_control$verbose)) {
    stop(
      "'q_control$verbose' must be TRUE or FALSE."
    )
  }

  if (!scalar_logical(q_control$keep_fits)) {
    stop(
      "'q_control$keep_fits' must be TRUE or FALSE."
    )
  }

  n <- prepared_data$n
  p <- length(start_corrected)
  q_min <- as.numeric(q_control$q_min)
  L <- as.numeric(q_control$L)

  # Step 1: initial grid from 1 to 0.8.
  initial_grid <- seq(
    from = 1,
    to = 0.8,
    length.out = m0
  )

  initial_spacing <-
    initial_grid[1L] -
    initial_grid[2L]

  tolerance <- max(
    100 * .Machine$double.eps,
    initial_spacing * 1e-10
  )

  q_key <- function(q_value) {
    sprintf(
      "%.15f",
      as.numeric(q_value)
    )
  }

  fit_cache <- new.env(
    parent = emptyenv()
  )

  evaluation_order <- character(0L)
  grid_history <- list()
  sqv_history <- list()

  fit_one_q <- function(q_value) {
    key <- q_key(q_value)

    if (exists(
      key,
      envir = fit_cache,
      inherits = FALSE
    )) {
      return(
        get(
          key,
          envir = fit_cache,
          inherits = FALSE
        )
      )
    }

    fit <- tryCatch(
      .fit_lmve_mlqe_fixed(
        Y = prepared_data$Y,
        X = prepared_data$X,
        tau_y = prepared_data$tau_y,
        tau_x = prepared_data$tau_x,
        start = start_corrected,
        q = q_value,
        level = level,
        method = method,
        lower = lower,
        upper = upper,
        control = control,
        use_score = use_score,
        compute_vcov = TRUE
      ),
      error = function(e) {
        stop(
          "The SQV fit failed at q = ",
          format(q_value, digits = 15),
          ": ",
          conditionMessage(e)
        )
      }
    )

    if (!identical(
      as.integer(fit$convergence),
      0L
    )) {
      stop(
        "The SQV fit did not converge at q = ",
        format(q_value, digits = 15),
        "."
      )
    }

    standard_error <- as.numeric(
      fit$standard.error
    )

    valid_standard_error <-
      length(standard_error) == p &&
      !anyNA(standard_error) &&
      all(is.finite(standard_error)) &&
      all(standard_error > 0)

    if (!valid_standard_error) {
      stop(
        "The SQV standardized estimate cannot be computed at q = ",
        format(q_value, digits = 15),
        " because its standard errors are not finite and strictly positive."
      )
    }

    # Step 2: standardized estimates.
    fit$sqv.standardized <- setNames(
      as.numeric(fit$coefficients) /
        (sqrt(n) * standard_error),
      names(fit$coefficients)
    )

    assign(
      key,
      fit,
      envir = fit_cache
    )

    evaluation_order <<- c(
      evaluation_order,
      key
    )

    fit
  }

  normalize_grid <- function(q_grid) {
    q_grid <- as.numeric(q_grid)
    q_grid[q_grid < q_min] <- q_min

    q_grid <- q_grid[
      c(
        TRUE,
        diff(q_grid) < -tolerance
      )
    ]

    q_grid
  }

  evaluate_grid <- function(q_grid,
                            pass,
                            stage) {
    q_grid <- normalize_grid(q_grid)

    if (length(q_grid) < 2L) {
      return(
        list(
          reached_minimum = TRUE,
          stable = FALSE,
          q_grid = q_grid,
          fits = list(),
          z = NULL,
          sqv = numeric(0L),
          violation = NA_integer_
        )
      )
    }

    if (any(diff(q_grid) >= -tolerance)) {
      stop(
        "Every SQV grid must be strictly decreasing."
      )
    }

    fits <- lapply(
      q_grid,
      fit_one_q
    )

    z_matrix <- do.call(
      rbind,
      lapply(
        fits,
        function(fit) {
          as.numeric(
            fit$sqv.standardized
          )
        }
      )
    )

    parameter_names <- names(start_corrected)

    if (is.null(parameter_names) ||
        any(parameter_names == "")) {
      parameter_names <- paste0(
        "theta",
        seq_len(p)
      )
    }

    colnames(z_matrix) <- make.unique(
      parameter_names
    )

    rownames(z_matrix) <- format(
      q_grid,
      digits = 12
    )

    # Step 3: standardized quadratic variation.
    sqv <- vapply(
      seq_len(length(q_grid) - 1L),
      function(index) {
        sqrt(
          sum(
            (
              z_matrix[index, ] -
                z_matrix[index + 1L, ]
            )^2
          )
        ) / p
      },
      numeric(1L)
    )

    stable_pairs <- sqv < L
    stable_grid <- all(stable_pairs)

    # Step 4: identify the smallest q_k for which SQV_qk >= L.
    # Since the grid is descending, this is the last violating pair.
    violation <- if (stable_grid) {
      NA_integer_
    } else {
      max(which(!stable_pairs))
    }

    grid_history[[length(grid_history) + 1L]] <<- list(
      pass = pass,
      stage = stage,
      q = q_grid,
      standardized = z_matrix,
      sqv = sqv,
      stable = stable_grid,
      violation = violation
    )

    sqv_history[[length(sqv_history) + 1L]] <<- data.frame(
      pass = rep(pass, length(sqv)),
      stage = rep(stage, length(sqv)),
      q_current = q_grid[-length(q_grid)],
      q_next = q_grid[-1L],
      SQV = sqv,
      L = rep(L, length(sqv)),
      stable = stable_pairs,
      stringsAsFactors = FALSE
    )

    if (isTRUE(q_control$verbose)) {
      print(
        sqv_history[[length(sqv_history)]]
      )
    }

    list(
      reached_minimum =
        min(q_grid) <= q_min + tolerance,
      stable = stable_grid,
      q_grid = q_grid,
      fits = fits,
      z = z_matrix,
      sqv = sqv,
      violation = violation
    )
  }

  make_next_grid <- function(grid_result) {
    violation <- grid_result$violation

    if (is.na(violation)) {
      stop(
        "The SQV violation index could not be determined."
      )
    }

    spacing <-
      grid_result$q_grid[1L] -
      grid_result$q_grid[2L]

    # Step 5: the next grid begins at q_(k*+1) and preserves
    # the spacing of the current grid.
    q_start <-
      grid_result$q_grid[violation + 1L]

    q_start -
      spacing *
      seq.int(0L, m - 1L)
  }

  run_pass <- function(pass,
                       pass_initial_grid) {
    q_grid <- pass_initial_grid
    stage <- 0L

    repeat {
      grid_result <- evaluate_grid(
        q_grid = q_grid,
        pass = pass,
        stage = stage
      )

      if (grid_result$stable) {
        return(
          list(
            success = TRUE,
            selected_q = grid_result$q_grid[1L],
            pass = pass,
            stage = stage,
            reason = "stability"
          )
        )
      }

      if (grid_result$reached_minimum) {
        return(
          list(
            success = FALSE,
            selected_q = NA_real_,
            pass = pass,
            stage = stage,
            reason = "q_min_reached_without_stability"
          )
        )
      }

      q_grid <- make_next_grid(
        grid_result
      )

      stage <- stage + 1L
    }
  }

  # Steps 1--6: first pass.
  first_pass <- run_pass(
    pass = 1L,
    pass_initial_grid = initial_grid
  )

  if (first_pass$success) {
    selection_result <- first_pass
  } else {
    # Step 7: restart with m0 <- m and an initial grid shortened
    # from 1 to q_(m-1) of the original initial grid.
    restart_lower <- initial_grid[m]

    restart_grid <- seq(
      from = 1,
      to = restart_lower,
      length.out = m
    )

    second_pass <- run_pass(
      pass = 2L,
      pass_initial_grid = restart_grid
    )

    if (second_pass$success) {
      selection_result <- second_pass
    } else {
      # Step 8: final assignment q = 1.
      selection_result <- list(
        success = FALSE,
        selected_q = 1,
        pass = 2L,
        stage = second_pass$stage,
        reason = "q_min_reached_after_restart_fallback_to_q_1"
      )
    }
  }

  selected_q <- as.numeric(
    selection_result$selected_q
  )

  selected_key <- q_key(selected_q)

  if (!exists(
    selected_key,
    envir = fit_cache,
    inherits = FALSE
  )) {
    fit_one_q(selected_q)
  }

  selected_fit <- get(
    selected_key,
    envir = fit_cache,
    inherits = FALSE
  )

  selected_fit$sqv.returned.fit.source <-
    "fit_used_in_SQV"

  sqv_table <- if (length(sqv_history) == 0L) {
    data.frame(
      pass = integer(0L),
      stage = integer(0L),
      q_current = numeric(0L),
      q_next = numeric(0L),
      SQV = numeric(0L),
      L = numeric(0L),
      stable = logical(0L)
    )
  } else {
    do.call(
      rbind,
      sqv_history
    )
  }

  evaluation_rows <- lapply(
    evaluation_order,
    function(key) {
      fit <- get(
        key,
        envir = fit_cache,
        inherits = FALSE
      )

      corrected_estimates <- setNames(
        as.list(
          as.numeric(
            fit$coefficients
          )
        ),
        paste0(
          "estimate.",
          names(fit$coefficients)
        )
      )

      original_estimates <- setNames(
        as.list(
          as.numeric(
            fit$coefficients.star
          )
        ),
        paste0(
          "estimate.star.",
          names(fit$coefficients.star)
        )
      )

      standard_errors <- setNames(
        as.list(
          as.numeric(
            fit$standard.error
          )
        ),
        paste0(
          "se.",
          names(fit$coefficients)
        )
      )

      standardized <- setNames(
        as.list(
          as.numeric(
            fit$sqv.standardized
          )
        ),
        paste0(
          "z.",
          names(fit$coefficients)
        )
      )

      as.data.frame(
        c(
          list(
            q = fit$q,
            convergence = fit$convergence,
            objective = fit$objective
          ),
          corrected_estimates,
          original_estimates,
          standard_errors,
          standardized
        ),
        check.names = FALSE
      )
    }
  )

  evaluations <- do.call(
    rbind,
    evaluation_rows
  )

  rownames(evaluations) <- NULL

  evaluations <- evaluations[
    order(
      evaluations$q,
      decreasing = TRUE
    ),
    ,
    drop = FALSE
  ]

  cached_fits <- NULL

  if (isTRUE(q_control$keep_fits)) {
    cached_fits <- setNames(
      lapply(
        evaluation_order,
        function(key) {
          get(
            key,
            envir = fit_cache,
            inherits = FALSE
          )
        }
      ),
      vapply(
        evaluation_order,
        function(key) {
          format(
            get(
              key,
              envir = fit_cache,
              inherits = FALSE
            )$q,
            digits = 15
          )
        },
        character(1L)
      )
    )
  }

  selected_fit$starting.values <- start_corrected
  selected_fit$sqv.initial.starting.values <- start_corrected
  selected_fit$automatic.start <- automatic_start
  selected_fit$implementation.version <-
    "LMVE-MLqE-SQV-original-algorithm-2026-07-30"

  selected_fit$selected.q <- selected_q

  selected_fit$q.selection <- list(
    method = "SQV",
    selected.q = selected_q,
    reason = selection_result$reason,
    pass = selection_result$pass,
    stage = selection_result$stage,
    configuration = list(
      initial_lower = 0.8,
      m0 = m0,
      m = m,
      q_min = q_min,
      L = L,
      initial_spacing = initial_spacing,
      restart_lower = initial_grid[m]
    ),
    returned.fit.source =
      selected_fit$sqv.returned.fit.source,
    history = sqv_table,
    grids = grid_history,
    evaluations = evaluations,
    fits = cached_fits
  )

  class(selected_fit) <- c(
    "lmve_mlqe_sqv",
    class(selected_fit)
  )

  selected_fit
}


# ============================================================
# S3 methods
# ============================================================

print.lmve_mlqe <- function(x,
                            digits = max(
                              3L,
                              getOption("digits") - 3L
                            ),
                            ...) {
  cat(
    "Linear normal measurement-error model fitted by MLqE\n"
  )

  cat(
    "q =",
    format(x$q, digits = digits),
    "| n =",
    x$nobs,
    "| optimization method =",
    x$method,
    "\n"
  )

  cat(
    "optim convergence code:",
    x$convergence
  )

  if (!is.null(x$message)) {
    cat(" -", x$message)
  }

  cat(
    "\nConsistency correction: tau_q(theta_star)",
    "\n\n"
  )

  coefficient_table <- cbind(
    Estimate = x$coefficients
  )

  if (!is.null(x$standard.error)) {
    coefficient_table <- cbind(
      coefficient_table,
      `Std. Error` = x$standard.error
    )
  }

  printCoefmat(
    coefficient_table,
    digits = digits,
    na.print = "NA"
  )

  if (!is.null(x$vcov.warning)) {
    cat(
      "\nWarning:",
      x$vcov.warning,
      "\n"
    )
  }

  invisible(x)
}


summary.lmve_mlqe <- function(object,
                              ...) {
  coefficient_table <- cbind(
    Estimate = object$coefficients
  )

  if (!is.null(object$standard.error)) {
    z_value <-
      object$coefficients /
      object$standard.error

    p_value <- 2 *
      pnorm(
        abs(z_value),
        lower.tail = FALSE
      )

    coefficient_table <- cbind(
      coefficient_table,
      `Std. Error` =
        object$standard.error,
      `z value` = z_value,
      `Pr(>|z|)` = p_value
    )
  }

  result <- list(
    call = object$call,
    coefficients = coefficient_table,
    coefficients.star =
      object$coefficients.star,
    objective =
      object$objective,
    score =
      object$score.star,
    convergence =
      object$convergence,
    message =
      object$message,
    q = object$q,
    nobs = object$nobs,
    method = object$method,
    vcov.warning =
      object$vcov.warning
  )

  class(result) <-
    "summary.lmve_mlqe"

  result
}


print.summary.lmve_mlqe <- function(x,
                                    digits = max(
                                      3L,
                                      getOption("digits") - 3L
                                    ),
                                    ...) {
  cat(
    "Linear normal measurement-error model fitted by MLqE\n"
  )

  cat(
    "q =",
    format(x$q, digits = digits),
    "| n =",
    x$nobs,
    "| optimization method =",
    x$method,
    "\n"
  )

  cat(
    "optim convergence code:",
    x$convergence
  )

  if (!is.null(x$message)) {
    cat(" -", x$message)
  }

  cat(
    "\nConsistency correction: tau_q(theta_star)",
    "\n\nCorrected coefficients:\n"
  )

  printCoefmat(
    x$coefficients,
    digits = digits,
    na.print = "NA"
  )

  cat(
    "\nOriginal MLqE coefficients:\n"
  )

  print(
    x$coefficients.star,
    digits = digits
  )

  cat(
    "\nObjective function:",
    format(x$objective, digits = digits),
    "\n"
  )

  cat(
    "Maximum absolute original score component:",
    format(
      max(
        abs(x$score),
        na.rm = TRUE
      ),
      digits = digits
    ),
    "\n"
  )

  if (!is.null(x$vcov.warning)) {
    cat(
      "\nWarning:",
      x$vcov.warning,
      "\n"
    )
  }

  invisible(x)
}


coef.lmve_mlqe <- function(object,
                           type = c(
                             "corrected",
                             "star"
                           ),
                           ...) {
  type <- match.arg(type)

  if (type == "corrected") {
    object$coefficients
  } else {
    object$coefficients.star
  }
}


vcov.lmve_mlqe <- function(object,
                           type = c(
                             "corrected",
                             "star"
                           ),
                           ...) {
  type <- match.arg(type)

  if (type == "corrected") {
    object$vcov
  } else {
    object$vcov.star
  }
}


confint.lmve_mlqe <- function(object,
                              parm = seq_along(
                                object$coefficients
                              ),
                              level = object$level,
                              type = c(
                                "corrected",
                                "star"
                              ),
                              ...) {
  type <- match.arg(type)

  if (type == "corrected") {
    estimates <- object$coefficients
    standard_errors <-
      object$standard.error
  } else {
    estimates <-
      object$coefficients.star
    standard_errors <-
      object$standard.error.star
  }

  if (is.null(standard_errors)) {
    stop(
      "Standard errors were not computed."
    )
  }

  parameter_indices <- if (
    is.character(parm)
  ) {
    match(
      parm,
      names(estimates)
    )
  } else {
    parm
  }

  if (anyNA(parameter_indices)) {
    stop(
      "At least one requested parameter name was not found."
    )
  }

  z_value <- qnorm(
    1 - (1 - level) / 2
  )

  selected_estimates <-
    estimates[parameter_indices]

  selected_standard_errors <-
    standard_errors[parameter_indices]

  interval <- cbind(
    lower =
      selected_estimates -
      z_value *
      selected_standard_errors,
    upper =
      selected_estimates +
      z_value *
      selected_standard_errors
  )

  rownames(interval) <-
    names(selected_estimates)

  interval
}


fitted.lmve_mlqe <- function(object,
                             type = c(
                               "corrected",
                               "star"
                             ),
                             ...) {
  type <- match.arg(type)

  if (type == "corrected") {
    object$fitted.values
  } else {
    object$fitted.values.star
  }
}


residuals.lmve_mlqe <- function(object,
                                type = c(
                                  "corrected",
                                  "star"
                                ),
                                ...) {
  type <- match.arg(type)

  if (type == "corrected") {
    object$residuals
  } else {
    object$residuals.star
  }
}


predict.lmve_mlqe <- function(object,
                              newdata = NULL,
                              type = c(
                                "mean",
                                "mean_y",
                                "mean_x",
                                "structural_mean",
                                "covariance",
                                "variance_y",
                                "variance_x",
                                "covariance_yx",
                                "sigma_y",
                                "sigma_x"
                              ),
                              parameterization = c(
                                "corrected",
                                "star"
                              ),
                              ...) {
  type <- match.arg(type)

  parameterization <- match.arg(
    parameterization
  )

  theta <- if (
    parameterization == "corrected"
  ) {
    object$coefficients
  } else {
    object$coefficients.star
  }

  lambda_scale <- if (
    parameterization == "corrected"
  ) {
    1
  } else {
    object$q
  }

  beta0 <- theta["beta0"]
  beta1 <- theta["beta1"]
  mu_x <- theta["mu_x"]
  sigma2_x <- theta["sigma2_x"]
  sigma2 <- theta["sigma2"]

  if (type == "structural_mean") {
    if (is.null(newdata) ||
        (
          is.null(newdata$x) &&
          is.null(newdata$latent_x)
        )) {
      stop(
        "For type = 'structural_mean', 'newdata' must contain ",
        "'x' or 'latent_x'."
      )
    }

    x_value <- if (
      !is.null(newdata$latent_x)
    ) {
      as.numeric(newdata$latent_x)
    } else {
      as.numeric(newdata$x)
    }

    return(
      beta0 +
      beta1 * x_value
    )
  }

  if (is.null(newdata)) {
    tau_y <- object$tau_y
    tau_x <- object$tau_x
  } else {
    if (is.null(newdata$tau_y) ||
        is.null(newdata$tau_x)) {
      stop(
        "'newdata' must contain 'tau_y' and 'tau_x'."
      )
    }

    new_n <- max(
      length(newdata$tau_y),
      length(newdata$tau_x)
    )

    tau_y <- rep_len(
      as.numeric(newdata$tau_y),
      new_n
    )

    tau_x <- rep_len(
      as.numeric(newdata$tau_x),
      new_n
    )

    if (anyNA(tau_y) ||
        anyNA(tau_x) ||
        any(!is.finite(tau_y)) ||
        any(!is.finite(tau_x)) ||
        any(tau_y < 0) ||
        any(tau_x < 0)) {
      stop(
        "The new measurement-error variances must be finite ",
        "and non-negative."
      )
    }
  }

  prediction_n <- max(
    length(tau_y),
    length(tau_x)
  )

  tau_y <- rep_len(
    tau_y,
    prediction_n
  )

  tau_x <- rep_len(
    tau_x,
    prediction_n
  )

  mean_y <- beta0 + beta1 * mu_x
  mean_x <- mu_x

  if (type == "mean") {
    return(
      cbind(
        Y = rep(mean_y, prediction_n),
        X = rep(mean_x, prediction_n)
      )
    )
  }

  if (type == "mean_y") {
    return(
      rep(mean_y, prediction_n)
    )
  }

  if (type == "mean_x") {
    return(
      rep(mean_x, prediction_n)
    )
  }

  variance_y <-
    beta1^2 * sigma2_x +
    sigma2 +
    lambda_scale * tau_y

  variance_x <-
    sigma2_x +
    lambda_scale * tau_x

  covariance_yx <- rep(
    beta1 * sigma2_x,
    prediction_n
  )

  if (type == "variance_y") {
    return(variance_y)
  }

  if (type == "variance_x") {
    return(variance_x)
  }

  if (type == "covariance_yx") {
    return(covariance_yx)
  }

  if (type == "sigma_y") {
    return(sqrt(variance_y))
  }

  if (type == "sigma_x") {
    return(sqrt(variance_x))
  }

  covariance_array <- array(
    NA_real_,
    dim = c(2L, 2L, prediction_n),
    dimnames = list(
      c("Y", "X"),
      c("Y", "X"),
      NULL
    )
  )

  covariance_array[1L, 1L, ] <-
    variance_y

  covariance_array[1L, 2L, ] <-
    covariance_yx

  covariance_array[2L, 1L, ] <-
    covariance_yx

  covariance_array[2L, 2L, ] <-
    variance_x

  covariance_array
}


print.lmve_mlqe_sqv <- function(x,
                                digits = max(
                                  3L,
                                  getOption("digits") - 3L
                                ),
                                ...) {
  cat(
    "Tuning parameter selected by SQV: q =",
    format(
      x$selected.q,
      digits = digits
    ),
    "\n"
  )

  cat(
    "Selection result:",
    x$q.selection$reason,
    "| pass =",
    x$q.selection$pass,
    "| stage =",
    x$q.selection$stage,
    "\n\n"
  )

  NextMethod("print")
}


summary.lmve_mlqe_sqv <- function(object,
                                  ...) {
  result <- NextMethod("summary")

  result$selected.q <-
    object$selected.q

  result$q.selection.reason <-
    object$q.selection$reason

  result$q.selection.pass <-
    object$q.selection$pass

  result$q.selection.stage <-
    object$q.selection$stage

  class(result) <- c(
    "summary.lmve_mlqe_sqv",
    class(result)
  )

  result
}


print.summary.lmve_mlqe_sqv <- function(x,
                                        digits = max(
                                          3L,
                                          getOption("digits") - 3L
                                        ),
                                        ...) {
  cat(
    "Tuning parameter selected by SQV: q =",
    format(
      x$selected.q,
      digits = digits
    ),
    "\n"
  )

  cat(
    "Selection result:",
    x$q.selection.reason,
    "| pass =",
    x$q.selection.pass,
    "| stage =",
    x$q.selection.stage,
    "\n\n"
  )

  NextMethod("print")
}

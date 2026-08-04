# Generic MDPDE fitting for a linear normal measurement-error model
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
# Therefore, (Y_i, X_i)^T follows a bivariate normal distribution
# with mean
#
#   c(beta0 + beta1 * mu_x, mu_x)
#
# and covariance matrix
#
#   [ beta1^2 * sigma2_x + sigma2 + tau_yi, beta1 * sigma2_x ]
#   [ beta1 * sigma2_x,                       sigma2_x + tau_xi ]
#
# The parameter vector must contain:
#
#   beta0, beta1, mu_x, sigma2_x, sigma2.
#
# Both sigma2_x and sigma2 are estimated on their original variance
# scales and must be strictly positive.


# ============================================================
# Internal utility functions
# ============================================================

.lmve_mdpde_parameter_names <- c(
  "beta0",
  "beta1",
  "mu_x",
  "sigma2_x",
  "sigma2"
)


.lmve_mdpde_validate_start <- function(start) {
  if (!is.numeric(start) || length(start) != 5L) {
    stop(
      "'start' must be a numeric vector of length 5 containing ",
      "beta0, beta1, mu_x, sigma2_x and sigma2."
    )
  }

  if (is.null(names(start))) {
    names(start) <- .lmve_mdpde_parameter_names
  }

  missing_names <- setdiff(
    .lmve_mdpde_parameter_names,
    names(start)
  )

  if (length(missing_names) > 0L) {
    stop(
      "'start' must contain the following named parameters: ",
      paste(
        .lmve_mdpde_parameter_names,
        collapse = ", "
      ),
      "."
    )
  }

  start <- start[
    .lmve_mdpde_parameter_names
  ]

  start <- setNames(
    as.numeric(start),
    .lmve_mdpde_parameter_names
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


.lmve_mdpde_prepare_data <- function(Y,
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


# Compute the default starting values when start = NULL.
#
# These starting values follow the robust moment-based expressions
# used in the original measurement-error simulation scripts.
.lmve_mdpde_default_start <- function(data) {
  n <- data$n
  X <- data$X
  Y <- data$Y
  tau_x <- data$tau_x
  tau_y <- data$tau_y

  Mxy <- sum(
    (
      X -
        median(X)
    ) *
      Y
  ) /
    (
      n -
        1
    )

  Mx <- mad(X)^2
  My <- mad(Y)^2

  Xbar <- median(X)
  Ybar <- median(Y)

  tau_xbar <- median(
    tau_x
  )

  tau_ybar <- median(
    tau_y
  )

  denominator <-
    Mx -
    tau_xbar

  if (!is.finite(denominator) ||
      denominator == 0) {
    stop(
      "The automatic starting values could not be computed because ",
      "mad(X)^2 - median(tau_x) is zero or non-finite. ",
      "Supply 'start' explicitly."
    )
  }

  b1ini <-
    Mxy /
    denominator

  betaini <- c(
    beta0 =
      Ybar -
      Xbar *
      b1ini,
    beta1 =
      b1ini
  )

  mu_xini <- Xbar

  sigma2_xini <-
    Mx -
    tau_xbar

  sigma2ini <-
    My -
    Mxy^2 /
    denominator -
    tau_ybar

  # Protect the automatic variance starting values.
  # Whenever an initial variance is smaller than 0.05,
  # replace it with 0.1 before starting the optimization.
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
    mu_x =
      mu_xini,
    sigma2_x =
      sigma2_xini,
    sigma2 =
      sigma2ini
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
    .lmve_mdpde_parameter_names
  )
}


.lmve_mdpde_model_components <- function(theta,
                                         data,
                                         strict = FALSE) {
  theta <- setNames(
    as.numeric(theta),
    .lmve_mdpde_parameter_names
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

  mean_y <- beta0 +
    beta1 *
    mu_x

  mean_x <- mu_x

  variance_y <-
    beta1^2 *
    sigma2_x +
    sigma2 +
    data$tau_y

  covariance_yx <-
    rep(
      beta1 *
        sigma2_x,
      data$n
    )

  variance_x <-
    sigma2_x +
    data$tau_x

  determinant <-
    variance_y *
    variance_x -
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

  residual_y <-
    data$Y -
    mean_y

  residual_x <-
    data$X -
    mean_x

  # For this bivariate model, the normal log-density is evaluated
  # directly from the analytical determinant and inverse of each
  # 2 by 2 covariance matrix. This is equivalent to lMvn()/dMvn(),
  # but avoids a separate Cholesky decomposition for every observation.
  quadratic_form <-
    (
      variance_x *
        residual_y^2 -
        2 *
        covariance_yx *
        residual_y *
        residual_x +
        variance_y *
        residual_x^2
    ) /
    determinant

  log_density <-
    -log(2 * pi) -
    0.5 *
      log(determinant) -
    0.5 *
      quadratic_form

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
    Y = rep(
      mean_y,
      data$n
    ),
    X = rep(
      mean_x,
      data$n
    )
  )

  residuals <- cbind(
    Y = residual_y,
    X = residual_x
  )

  covariance_array <- array(
    NA_real_,
    dim = c(
      2L,
      2L,
      data$n
    ),
    dimnames = list(
      c(
        "Y",
        "X"
      ),
      c(
        "Y",
        "X"
      ),
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


.lmve_mdpde_derivative_matrix <- function(theta) {
  theta <- setNames(
    as.numeric(theta),
    .lmve_mdpde_parameter_names
  )

  beta1 <- theta["beta1"]
  mu_x <- theta["mu_x"]
  sigma2_x <- theta["sigma2_x"]

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
      c(
        "mean_Y",
        "mean_X"
      ),
      .lmve_mdpde_parameter_names
    )
  )

  covariance_jacobian <- matrix(
    c(
      0,
      2 * sigma2_x * beta1,
      0,
      beta1^2,
      1,
      0,
      sigma2_x,
      0,
      beta1,
      0,
      0,
      sigma2_x,
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
      .lmve_mdpde_parameter_names
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


.lmve_mdpde_block_matrix <- function(top_left,
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
    top_dimension +
      seq_len(bottom_dimension),
    top_dimension +
      seq_len(bottom_dimension)
  ] <- bottom_right

  result
}


# ============================================================
# Fixed-q MDPDE fit
# ============================================================

.fit_lmve_mdpde_fixed <- function(Y,
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

  data <- .lmve_mdpde_prepare_data(
    Y = Y,
    X = X,
    tau_y = tau_y,
    tau_x = tau_x
  )

  automatic_start <- is.null(
    start
  )

  start <- if (automatic_start) {
    .lmve_mdpde_default_start(
      data
    )
  } else {
    .lmve_mdpde_validate_start(
      start
    )
  }

  p <- length(start)

  if (length(q) != 1L ||
      !is.numeric(q) ||
      !is.finite(q) ||
      q <= 0 ||
      q > 1) {
    stop(
      "'q' must be a number in the interval (0, 1]."
    )
  }

  if (length(level) != 1L ||
      !is.numeric(level) ||
      !is.finite(level) ||
      level <= 0 ||
      level >= 1) {
    stop(
      "'level' must be a number in the interval (0, 1)."
    )
  }

  if (!is.list(control)) {
    stop("'control' must be a list.")
  }

  objective <- function(theta) {
    theta <- setNames(
      as.numeric(theta),
      .lmve_mdpde_parameter_names
    )

    sigma2_x <- theta["sigma2_x"]
    sigma2 <- theta["sigma2"]

    # For unconstrained methods such as BFGS, reject parameter
    # values that do not define strictly positive variances.
    if (!is.finite(sigma2_x) ||
        !is.finite(sigma2) ||
        sigma2_x <= 0 ||
        sigma2 <= 0) {
      return(NaN)
    }

    components <- .lmve_mdpde_model_components(
      theta = theta,
      data = data
    )

    if (is.null(components)) {
      return(NaN)
    }

    if (q == 1) {
      value <- sum(
        components$log_density
      )
    } else {
      one_minus_q <- 1 - q

      density_power <- exp(
        one_minus_q *
          components$log_density
      )

      integral_term <- exp(
        -log(2 - q) -
          one_minus_q *
          log(2 * pi) -
          0.5 *
          one_minus_q *
          log(components$determinant)
      )

      value <- sum(
        (2 - q) /
          one_minus_q *
          density_power -
          integral_term
      )
    }

    if (is.finite(value)) {
      value
    } else {
      NaN
    }
  }

  score <- function(theta) {
    components <- .lmve_mdpde_model_components(
      theta = theta,
      data = data,
      strict = TRUE
    )

    derivatives <- .lmve_mdpde_derivative_matrix(
      theta
    )

    F <- derivatives$F
    one_minus_q <- 1 - q

    score_sum <- rep(
      0,
      p
    )

    for (i in seq_len(data$n)) {
      Sigma <- matrix(
        c(
          components$variance.y[i],
          components$covariance.yx[i],
          components$covariance.yx[i],
          components$variance.x[i]
        ),
        nrow = 2L
      )

      Sigma_inverse <- solve(
        Sigma
      )

      H <- .lmve_mdpde_block_matrix(
        top_left = Sigma_inverse,
        bottom_right =
          0.5 *
          kronecker(
            Sigma_inverse,
            Sigma_inverse
          )
      )

      z <- c(
        components$residual.y[i],
        components$residual.x[i]
      )

      centered_second_moment <-
        tcrossprod(z) -
        Sigma

      s <- c(
        z,
        c(centered_second_moment)
      )

      density_power <- if (q == 1) {
        1
      } else {
        exp(
          one_minus_q *
            components$log_density[i]
        )
      }

      correction_scalar <- if (q == 1) {
        0
      } else {
        (q - 1) *
          exp(
            -2 *
              log(2 - q) -
              one_minus_q *
              log(2 * pi) -
              0.5 *
              one_minus_q *
              log(
                components$determinant[i]
              )
          )
      }

      correction_vector <- c(
        0,
        0,
        correction_scalar *
          c(Sigma)
      )

      score_sum <-
        score_sum +
        as.numeric(
          (2 - q) *
            t(F) %*%
            H %*%
            (
              s *
                density_power -
                correction_vector
            )
        )
    }

    setNames(
      score_sum,
      .lmve_mdpde_parameter_names
    )
  }

  compute_JK <- function(theta) {
    components <- .lmve_mdpde_model_components(
      theta = theta,
      data = data,
      strict = TRUE
    )

    derivatives <- .lmve_mdpde_derivative_matrix(
      theta
    )

    F <- derivatives$F
    one_minus_q <- 1 - q

    J <- matrix(
      0,
      nrow = p,
      ncol = p
    )

    K <- matrix(
      0,
      nrow = p,
      ncol = p
    )

    kappa <-
      (2 - 2 * q)^2 -
      one_minus_q^2 *
      (3 - 2 * q)^3 /
      (2 - q)^4

    for (i in seq_len(data$n)) {
      Sigma <- matrix(
        c(
          components$variance.y[i],
          components$covariance.yx[i],
          components$covariance.yx[i],
          components$variance.x[i]
        ),
        nrow = 2L
      )

      Sigma_inverse <- solve(
        Sigma
      )

      Sigma_inverse_vector <- c(
        Sigma_inverse
      )

      covariance_block_J <-
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

      H_J <- .lmve_mdpde_block_matrix(
        top_left = Sigma_inverse,
        bottom_right = covariance_block_J
      )

      c3 <- exp(
        -log(2 - q) -
          one_minus_q *
          log(2 * pi) -
          0.5 *
          one_minus_q *
          log(
            components$determinant[i]
          )
      )

      J <-
        J +
        c3 *
        t(F) %*%
        H_J %*%
        F

      covariance_block_K <-
        (
          0.5 *
            kronecker(
              Sigma_inverse,
              Sigma_inverse
            ) +
            0.25 *
            kappa *
            tcrossprod(
              Sigma_inverse_vector
            )
        ) /
        (3 - 2 * q)

      H_K <- .lmve_mdpde_block_matrix(
        top_left = Sigma_inverse,
        bottom_right = covariance_block_K
      )

      c4 <- exp(
        -log(3 - 2 * q) -
          2 *
          one_minus_q *
          log(2 * pi) -
          one_minus_q *
          log(
            components$determinant[i]
          )
      )

      K <-
        K +
        c4 *
        t(F) %*%
        H_K %*%
        F
    }

    K <-
      K *
      (2 - q)^2 /
      (3 - 2 * q)

    J <- (
      J +
        t(J)
    ) /
      2

    K <- (
      K +
        t(K)
    ) /
      2

    dimnames(J) <- list(
      .lmve_mdpde_parameter_names,
      .lmve_mdpde_parameter_names
    )

    dimnames(K) <- list(
      .lmve_mdpde_parameter_names,
      .lmve_mdpde_parameter_names
    )

    list(
      J = J,
      K = K,
      mean.jacobian =
        derivatives$mean,
      covariance.jacobian =
        derivatives$covariance,
      F = F
    )
  }

  .lmve_mdpde_model_components(
    theta = start,
    data = data,
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
    par = start,
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
      optim_arguments$gr <- score
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
      rep(
        Inf,
        p
      ),
      .lmve_mdpde_parameter_names
    )

    if (is.null(lower)) {
      lower <- default_lower
    }

    if (is.null(upper)) {
      upper <- default_upper
    }

    if (!is.null(names(lower))) {
      lower <- lower[
        .lmve_mdpde_parameter_names
      ]
    }

    if (!is.null(names(upper))) {
      upper <- upper[
        .lmve_mdpde_parameter_names
      ]
    }

    lower <- setNames(
      rep_len(
        as.numeric(lower),
        p
      ),
      .lmve_mdpde_parameter_names
    )

    upper <- setNames(
      rep_len(
        as.numeric(upper),
        p
      ),
      .lmve_mdpde_parameter_names
    )

    if (anyNA(lower) ||
        anyNA(upper)) {
      stop(
        "'lower' and 'upper' must not contain missing values."
      )
    }

    if (any(lower > upper)) {
      stop(
        "Each lower bound must be less than or equal ",
        "to its upper bound."
      )
    }

    if (any(start < lower) ||
        any(start > upper)) {
      stop(
        "All starting values must lie within the ",
        "specified bounds."
      )
    }

    optim_arguments$lower <- lower
    optim_arguments$upper <- upper
  }

  optimization <- do.call(
    optim,
    optim_arguments
  )

  estimate <- setNames(
    as.numeric(
      optimization$par
    ),
    .lmve_mdpde_parameter_names
  )

  fitted_components <-
    .lmve_mdpde_model_components(
      theta = estimate,
      data = data,
      strict = TRUE
    )

  score_at_estimate <- tryCatch(
    score(
      estimate
    ),
    error = function(e) {
      warning(
        "The score at the MDPDE estimate could not be computed: ",
        conditionMessage(e)
      )

      setNames(
        rep(
          NA_real_,
          p
        ),
        .lmve_mdpde_parameter_names
      )
    }
  )

  vcov_matrix <- NULL
  standard_error <- NULL
  confidence_interval <- NULL
  J_matrix <- NULL
  K_matrix <- NULL
  mean_jacobian <- NULL
  covariance_jacobian <- NULL
  F_matrix <- NULL
  vcov_warning <- NULL

  if (compute_vcov) {
    jk <- compute_JK(
      estimate
    )

    J_matrix <- jk$J
    K_matrix <- jk$K
    mean_jacobian <-
      jk$mean.jacobian
    covariance_jacobian <-
      jk$covariance.jacobian
    F_matrix <- jk$F

    covariance_result <- tryCatch(
      {
        J_inverse <- chol2inv(
          chol(
            J_matrix
          )
        )

        covariance <-
          J_inverse %*%
          K_matrix %*%
          t(J_inverse)

        covariance <- (
          covariance +
            t(covariance)
        ) /
          2

        list(
          value = covariance,
          warning = NULL
        )
      },
      error = function(cholesky_error) {
        tryCatch(
          {
            J_inverse <- solve(
              J_matrix
            )

            covariance <-
              J_inverse %*%
              K_matrix %*%
              t(J_inverse)

            covariance <- (
              covariance +
                t(covariance)
            ) /
              2

            list(
              value = covariance,
              warning = paste(
                "The Cholesky factorization of J failed; ",
                "solve(J) was used instead:",
                conditionMessage(
                  cholesky_error
                )
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
                  .lmve_mdpde_parameter_names,
                  .lmve_mdpde_parameter_names
                )
              ),
              warning = paste(
                "The covariance matrix could not be computed:",
                conditionMessage(
                  inverse_error
                )
              )
            )
          }
        )
      }
    )

    vcov_matrix <-
      covariance_result$value

    dimnames(vcov_matrix) <- list(
      .lmve_mdpde_parameter_names,
      .lmve_mdpde_parameter_names
    )

    vcov_warning <-
      covariance_result$warning

    covariance_diagonal <- diag(
      vcov_matrix
    )

    standard_error <- setNames(
      ifelse(
        is.na(
          covariance_diagonal
        ) |
          covariance_diagonal <
          0,
        NA_real_,
        sqrt(
          covariance_diagonal
        )
      ),
      .lmve_mdpde_parameter_names
    )

    z_value <- qnorm(
      1 -
        (
          1 -
            level
        ) /
          2
    )

    confidence_interval <- cbind(
      lower =
        estimate -
        z_value *
          standard_error,
      upper =
        estimate +
        z_value *
          standard_error
    )

    rownames(
      confidence_interval
    ) <- .lmve_mdpde_parameter_names
  }

  result <- list(
    call = call,
    starting.values = start,
    automatic.start = automatic_start,
    coefficients = estimate,
    standard.error = standard_error,
    vcov = vcov_matrix,
    conf.int = confidence_interval,
    J = J_matrix,
    K = K_matrix,
    mean.jacobian = mean_jacobian,
    covariance.jacobian =
      covariance_jacobian,
    F = F_matrix,
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
    objective =
      unname(
        optimization$value
      ),
    score = score_at_estimate,
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
    vcov.warning =
      vcov_warning,
    optimization.control =
      control,
    optim = optimization
  )

  class(result) <- "lmve_mdpde"

  result
}


# ============================================================
# Public fitting interface
# ============================================================

fit_lmve_mdpde <- function(Y,
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

  # Support calls in which q is supplied as the fifth positional
  # argument and no starting values are supplied. For example:
  #
  # fit_lmve_mdpde(Y, X, tau_y, tau_x, "SQV")
  # fit_lmve_mdpde(Y, X, tau_y, tau_x, 0.8)
  #
  # A valid explicit starting vector must have length 5, so a scalar
  # numeric value or the character value "SQV" is unambiguously q.
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
    identical(
      toupper(q),
      "SQV"
    )

  if (is.character(q) &&
      !sqv_requested) {
    stop(
      "'q' must be numeric or equal to 'SQV'."
    )
  }

  if (!sqv_requested) {
    fit <- .fit_lmve_mdpde_fixed(
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

    return(
      fit
    )
  }

  if (!isTRUE(compute_vcov)) {
    stop(
      "SQV selection requires 'compute_vcov = TRUE' because ",
      "standard errors are used to construct the standardized ",
      "estimates."
    )
  }

  fit <- .select_q_lmve_mdpde_sqv(
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

.select_q_lmve_mdpde_sqv <- function(Y,
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

  prepared_data <- .lmve_mdpde_prepare_data(
    Y = Y,
    X = X,
    tau_y = tau_y,
    tau_x = tau_x
  )

  start <- if (automatic_start) {
    .lmve_mdpde_default_start(
      prepared_data
    )
  } else {
    .lmve_mdpde_validate_start(
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
  p <- length(start)
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
      .fit_lmve_mdpde_fixed(
        Y = prepared_data$Y,
        X = prepared_data$X,
        tau_y = prepared_data$tau_y,
        tau_x = prepared_data$tau_x,
        start = start,
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

    parameter_names <- names(start)

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

      estimates <- setNames(
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
          estimates,
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

  selected_fit$starting.values <- start
  selected_fit$sqv.initial.starting.values <- start
  selected_fit$automatic.start <- automatic_start
  selected_fit$implementation.version <-
    "LMVE-MDPDE-SQV-original-algorithm-2026-07-30"

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
    "lmve_mdpde_sqv",
    class(selected_fit)
  )

  selected_fit
}


# ============================================================
# S3 methods
# ============================================================

print.lmve_mdpde <- function(x,
                             digits = max(
                               3L,
                               getOption(
                                 "digits"
                               ) -
                                 3L
                             ),
                             ...) {
  cat(
    "Linear normal measurement-error model fitted by MDPDE\n"
  )

  cat(
    "q =",
    format(
      x$q,
      digits = digits
    ),
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

  if (!is.null(
    x$message
  )) {
    cat(
      " -",
      x$message
    )
  }

  cat(
    "\n\n"
  )

  coefficient_table <- cbind(
    Estimate =
      x$coefficients
  )

  if (!is.null(
    x$standard.error
  )) {
    coefficient_table <- cbind(
      coefficient_table,
      `Std. Error` =
        x$standard.error
    )
  }

  printCoefmat(
    coefficient_table,
    digits = digits,
    na.print = "NA"
  )

  if (!is.null(
    x$vcov.warning
  )) {
    cat(
      "\nWarning:",
      x$vcov.warning,
      "\n"
    )
  }

  invisible(
    x
  )
}


summary.lmve_mdpde <- function(object,
                               ...) {
  coefficient_table <- cbind(
    Estimate =
      object$coefficients
  )

  if (!is.null(
    object$standard.error
  )) {
    z_value <-
      object$coefficients /
      object$standard.error

    p_value <- 2 *
      pnorm(
        abs(
          z_value
        ),
        lower.tail =
          FALSE
      )

    coefficient_table <- cbind(
      coefficient_table,
      `Std. Error` =
        object$standard.error,
      `z value` =
        z_value,
      `Pr(>|z|)` =
        p_value
    )
  }

  result <- list(
    call =
      object$call,
    coefficients =
      coefficient_table,
    objective =
      object$objective,
    score =
      object$score,
    convergence =
      object$convergence,
    message =
      object$message,
    q =
      object$q,
    nobs =
      object$nobs,
    method =
      object$method,
    vcov.warning =
      object$vcov.warning
  )

  class(
    result
  ) <- "summary.lmve_mdpde"

  result
}


print.summary.lmve_mdpde <- function(x,
                                     digits = max(
                                       3L,
                                       getOption(
                                         "digits"
                                       ) -
                                         3L
                                     ),
                                     ...) {
  cat(
    "Linear normal measurement-error model fitted by MDPDE\n"
  )

  cat(
    "q =",
    format(
      x$q,
      digits = digits
    ),
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

  if (!is.null(
    x$message
  )) {
    cat(
      " -",
      x$message
    )
  }

  cat(
    "\n\nCoefficients:\n"
  )

  printCoefmat(
    x$coefficients,
    digits = digits,
    na.print = "NA"
  )

  cat(
    "\nObjective function:",
    format(
      x$objective,
      digits = digits
    ),
    "\n"
  )

  cat(
    "Maximum absolute score component:",
    format(
      max(
        abs(
          x$score
        ),
        na.rm = TRUE
      ),
      digits = digits
    ),
    "\n"
  )

  if (!is.null(
    x$vcov.warning
  )) {
    cat(
      "\nWarning:",
      x$vcov.warning,
      "\n"
    )
  }

  invisible(
    x
  )
}


coef.lmve_mdpde <- function(object,
                            ...) {
  object$coefficients
}


vcov.lmve_mdpde <- function(object,
                            ...) {
  object$vcov
}


confint.lmve_mdpde <- function(object,
                               parm = seq_along(
                                 object$coefficients
                               ),
                               level = object$level,
                               ...) {
  if (is.null(
    object$standard.error
  )) {
    stop(
      "Standard errors were not computed."
    )
  }

  parameter_indices <- if (is.character(
    parm
  )) {
    match(
      parm,
      names(
        object$coefficients
      )
    )
  } else {
    parm
  }

  if (anyNA(
    parameter_indices
  )) {
    stop(
      "At least one requested parameter name was not found."
    )
  }

  z_value <- qnorm(
    1 -
      (
        1 -
          level
      ) /
        2
  )

  estimates <-
    object$coefficients[
      parameter_indices
    ]

  standard_errors <-
    object$standard.error[
      parameter_indices
    ]

  interval <- cbind(
    lower =
      estimates -
      z_value *
        standard_errors,
    upper =
      estimates +
      z_value *
        standard_errors
  )

  rownames(
    interval
  ) <- names(
    estimates
  )

  interval
}


fitted.lmve_mdpde <- function(object,
                              ...) {
  object$fitted.values
}


residuals.lmve_mdpde <- function(object,
                                 ...) {
  object$residuals
}


predict.lmve_mdpde <- function(object,
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
                               ...) {
  type <- match.arg(
    type
  )

  theta <- object$coefficients

  beta0 <- theta["beta0"]
  beta1 <- theta["beta1"]
  mu_x <- theta["mu_x"]
  sigma2_x <- theta["sigma2_x"]
  sigma2 <- theta["sigma2"]

  if (type == "structural_mean") {
    if (is.null(
      newdata
    ) ||
        (
          is.null(
            newdata$x
          ) &&
            is.null(
              newdata$latent_x
            )
        )) {
      stop(
        "For type = 'structural_mean', 'newdata' must contain ",
        "'x' or 'latent_x'."
      )
    }

    x_value <- if (!is.null(
      newdata$latent_x
    )) {
      as.numeric(
        newdata$latent_x
      )
    } else {
      as.numeric(
        newdata$x
      )
    }

    return(
      beta0 +
        beta1 *
        x_value
    )
  }

  if (is.null(
    newdata
  )) {
    tau_y <- object$tau_y
    tau_x <- object$tau_x
  } else {
    if (is.null(
      newdata$tau_y
    ) ||
        is.null(
          newdata$tau_x
        )) {
      stop(
        "'newdata' must contain 'tau_y' and 'tau_x'."
      )
    }

    new_n <- max(
      length(
        newdata$tau_y
      ),
      length(
        newdata$tau_x
      )
    )

    tau_y <- rep_len(
      as.numeric(
        newdata$tau_y
      ),
      new_n
    )

    tau_x <- rep_len(
      as.numeric(
        newdata$tau_x
      ),
      new_n
    )

    if (anyNA(
      tau_y
    ) ||
        anyNA(
          tau_x
        ) ||
        any(!is.finite(
          tau_y
        )) ||
        any(!is.finite(
          tau_x
        )) ||
        any(
          tau_y <
            0
        ) ||
        any(
          tau_x <
            0
        )) {
      stop(
        "The new measurement-error variances must be finite ",
        "and non-negative."
      )
    }
  }

  prediction_n <- length(
    tau_y
  )

  mean_y <- beta0 +
    beta1 *
    mu_x

  mean_x <- mu_x

  if (type == "mean") {
    return(
      cbind(
        Y = rep(
          mean_y,
          prediction_n
        ),
        X = rep(
          mean_x,
          prediction_n
        )
      )
    )
  }

  if (type == "mean_y") {
    return(
      rep(
        mean_y,
        prediction_n
      )
    )
  }

  if (type == "mean_x") {
    return(
      rep(
        mean_x,
        prediction_n
      )
    )
  }

  variance_y <-
    beta1^2 *
    sigma2_x +
    sigma2 +
    tau_y

  variance_x <-
    sigma2_x +
    tau_x

  covariance_yx <- rep(
    beta1 *
      sigma2_x,
    prediction_n
  )

  if (type == "variance_y") {
    return(
      variance_y
    )
  }

  if (type == "variance_x") {
    return(
      variance_x
    )
  }

  if (type == "covariance_yx") {
    return(
      covariance_yx
    )
  }

  if (type == "sigma_y") {
    return(
      sqrt(
        variance_y
      )
    )
  }

  if (type == "sigma_x") {
    return(
      sqrt(
        variance_x
      )
    )
  }

  covariance_array <- array(
    NA_real_,
    dim = c(
      2L,
      2L,
      prediction_n
    ),
    dimnames = list(
      c(
        "Y",
        "X"
      ),
      c(
        "Y",
        "X"
      ),
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


print.lmve_mdpde_sqv <- function(x,
                                 digits = max(
                                   3L,
                                   getOption(
                                     "digits"
                                   ) -
                                     3L
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

  NextMethod(
    "print"
  )
}


summary.lmve_mdpde_sqv <- function(object,
                                   ...) {
  result <- NextMethod(
    "summary"
  )

  result$selected.q <-
    object$selected.q

  result$q.selection.reason <-
    object$q.selection$reason

  result$q.selection.pass <-
    object$q.selection$pass

  result$q.selection.stage <-
    object$q.selection$stage

  class(
    result
  ) <- c(
    "summary.lmve_mdpde_sqv",
    class(
      result
    )
  )

  result
}


print.summary.lmve_mdpde_sqv <- function(x,
                                         digits = max(
                                           3L,
                                           getOption(
                                             "digits"
                                           ) -
                                             3L
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

  NextMethod(
    "print"
  )
}

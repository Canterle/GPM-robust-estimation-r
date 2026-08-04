# ============================================================
# Maximum likelihood estimation for a linear errors-in-variables
# model with a bivariate Student-t distribution
#
# Scalar version of the structural model
#
#   x1_i = beta0 + beta1 * x2_i + q_i,
#   X1_i = x1_i + delta_x1_i,
#   X2_i = x2_i + delta_x2_i.
#
# For each observation i, assume
#
#   (x2_i, q_i, delta_x1_i, delta_x2_i)^T
#     ~ t_4(
#         (mu_x, 0, 0, 0)^T,
#         diag(sigma2_x, sigma2, tau_x1_i, tau_x2_i),
#         nu
#       ),
#
# independently across i. The diagonal block structure refers to the
# Student-t scale matrix. Its components are not marginally independent
# under the Student-t model; they are conditionally independent given
# the common scale-mixture variable for observation i.
#
# By affine closure of the multivariate Student-t family, the observed
# vector
#
#   Y_i = (X1_i, X2_i)^T
#
# has distribution
#
#   Y_i ~ t_2(mu(theta), Sigma_i(theta), nu),
#
# where
#
#   mu(theta) = c(beta0 + beta1 * mu_x, mu_x)
#
# and
#
#   Sigma_i(theta) =
#
#     [ beta1^2 * sigma2_x + sigma2 + tau_x1_i,
#       beta1 * sigma2_x                               ]
#     [ beta1 * sigma2_x,
#       sigma2_x + tau_x2_i                           ].
#
# Here, Sigma_i(theta), sigma2_x, sigma2, tau_x1_i and tau_x2_i
# are Student-t scale quantities. They are not covariance quantities.
# When nu > 2,
#
#   Var(Y_i) = nu / (nu - 2) * Sigma_i(theta).
#
# The parameter vector is
#
#   theta = (beta0, beta1, mu_x, sigma2_x, sigma2),
#
# where sigma2 is the scale parameter of the equation error q_i.
#
# In the function interface, Y and X denote the observed variables X1
# and X2, while tau_y and tau_x denote the known scale quantities
# tau_x1 and tau_x2, respectively.
#
# The degrees of freedom nu are fixed rather than estimated.
# ============================================================


# ============================================================
# Internal utility functions
# ============================================================

.lmve_t_parameter_names <- c(
  "beta0",
  "beta1",
  "mu_x",
  "sigma2_x",
  "sigma2"
)


.lmve_t_prepare_data <- function(Y,
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
      "The known measurement-error scale values 'tau_y' and ",
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


.lmve_t_validate_start <- function(start) {
  if (!is.numeric(start) ||
      length(start) != 5L) {
    stop(
      "'start' must be a numeric vector of length 5 containing ",
      "beta0, beta1, mu_x, sigma2_x and sigma2."
    )
  }

  if (is.null(names(start))) {
    names(start) <- .lmve_t_parameter_names
  }

  missing_names <- setdiff(
    .lmve_t_parameter_names,
    names(start)
  )

  if (length(missing_names) > 0L) {
    stop(
      "'start' must contain the following named parameters: ",
      paste(
        .lmve_t_parameter_names,
        collapse = ", "
      ),
      "."
    )
  }

  start <- setNames(
    as.numeric(
      start[
        .lmve_t_parameter_names
      ]
    ),
    .lmve_t_parameter_names
  )

  if (anyNA(start) ||
      any(!is.finite(start))) {
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



# Robust moment-based starting values. These preserve the rule used in
# the normal errors-in-variables implementation. Because the fitted
# parameters are Student-t scale quantities, empirical second moments
# are converted approximately to the scale parameterization when
# nu > 2. These values are used only to initialize the optimizer.
.lmve_t_default_start <- function(data,
                                  nu) {
  n <- data$n
  X <- data$X
  Y <- data$Y

  X_median <- stats::median(X)
  Y_median <- stats::median(Y)

  Mxy <- sum(
    (X - X_median) *
      (Y - Y_median)
  ) / (n - 1)

  Mx <- stats::mad(X)^2
  My <- stats::mad(Y)^2

  if (nu > 2) {
    moment_to_scale <- (nu - 2) / nu

    Mxy <- moment_to_scale * Mxy
    Mx <- moment_to_scale * Mx
    My <- moment_to_scale * My
  }

  tau_x_center <- stats::median(
    data$tau_x
  )

  tau_y_center <- stats::median(
    data$tau_y
  )

  denominator <- Mx - tau_x_center

  if (!is.finite(denominator) ||
      abs(denominator) <= sqrt(.Machine$double.eps)) {
    stop(
      "The automatic starting values could not be computed because ",
      "the robust scale estimate of X minus median(tau_x) is zero or ",
      "non-finite. Supply 'start' explicitly."
    )
  }

  beta1_start <- Mxy / denominator

  beta0_start <-
    Y_median -
    beta1_start * X_median

  sigma2_x_start <- denominator

  sigma2_start <-
    My -
    Mxy^2 / denominator -
    tau_y_center

  if (!is.finite(sigma2_x_start) ||
      sigma2_x_start < 0.05) {
    sigma2_x_start <- 0.1
  }

  if (!is.finite(sigma2_start) ||
      sigma2_start < 0.05) {
    sigma2_start <- 0.1
  }

  start <- c(
    beta0 = beta0_start,
    beta1 = beta1_start,
    mu_x = X_median,
    sigma2_x = sigma2_x_start,
    sigma2 = sigma2_start
  )

  .lmve_t_validate_start(start)
}


.lmve_t_safe_inverse <- function(M) {
  inverse <- tryCatch(
    chol2inv(
      chol(M)
    ),
    error = function(e) {
      tryCatch(
        solve(M),
        error = function(e2) NULL
      )
    }
  )

  if (is.null(inverse) ||
      anyNA(inverse) ||
      any(!is.finite(inverse))) {
    inverse <- matrix(
      NA_real_,
      nrow = nrow(M),
      ncol = ncol(M),
      dimnames = dimnames(M)
    )
  }

  inverse
}


.lmve_t_derivative_matrices <- function(theta) {
  theta <- setNames(
    as.numeric(theta),
    .lmve_t_parameter_names
  )

  beta1 <- theta["beta1"]
  mu_x <- theta["mu_x"]
  sigma2_x <- theta["sigma2_x"]

  D <- matrix(
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
        "mu_Y",
        "mu_X"
      ),
      .lmve_t_parameter_names
    )
  )

  V <- matrix(
    c(
      0,
      2 * beta1 * sigma2_x,
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
      .lmve_t_parameter_names
    )
  )

  list(
    mean = D,
    scale = V,
    F = rbind(
      D,
      V
    )
  )
}


.lmve_t_model_components <- function(theta,
                                     data,
                                     nu,
                                     strict = FALSE) {
  theta <- setNames(
    as.numeric(theta),
    .lmve_t_parameter_names
  )

  beta0 <- theta["beta0"]
  beta1 <- theta["beta1"]
  mu_x <- theta["mu_x"]
  sigma2_x_parameter <- theta["sigma2_x"]
  sigma2_parameter <- theta["sigma2"]

  invalid_dispersion <-
    !is.finite(sigma2_x_parameter) ||
    !is.finite(sigma2_parameter) ||
    sigma2_x_parameter <= 0 ||
    sigma2_parameter <= 0

  if (invalid_dispersion) {
    if (strict) {
      stop(
        "'sigma2_x' and 'sigma2' must be strictly positive."
      )
    }

    return(NULL)
  }

  sigma2_x_scale <-
    sigma2_x_parameter

  sigma2_scale <-
    sigma2_parameter

  tau_y_scale <-
    data$tau_y

  tau_x_scale <-
    data$tau_x

  mean_y <-
    beta0 +
    beta1 * mu_x

  mean_x <- mu_x

  scale_y <-
    beta1^2 * sigma2_x_scale +
    sigma2_scale +
    tau_y_scale

  scale_yx <- rep(
    beta1 * sigma2_x_scale,
    data$n
  )

  scale_x <-
    sigma2_x_scale +
    tau_x_scale

  determinant <-
    scale_y * scale_x -
    scale_yx^2

  valid <-
    all(is.finite(scale_y)) &&
    all(is.finite(scale_yx)) &&
    all(is.finite(scale_x)) &&
    all(is.finite(determinant)) &&
    all(scale_y > 0) &&
    all(scale_x > 0) &&
    all(determinant > 0)

  if (!valid) {
    if (strict) {
      stop(
        "At least one fitted Student-t scale matrix is not ",
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

  inverse_11 <-
    scale_x /
    determinant

  inverse_12 <-
    -scale_yx /
    determinant

  inverse_22 <-
    scale_y /
    determinant

  quadratic_form <-
    inverse_11 * residual_y^2 +
    2 * inverse_12 * residual_y * residual_x +
    inverse_22 * residual_x^2

  d <- 2

  log_density <-
    lgamma((nu + d) / 2) -
    lgamma(nu / 2) -
    0.5 * d * log(nu * pi) -
    0.5 * log(determinant) -
    0.5 * (nu + d) *
    log1p(quadratic_form / nu)

  if (any(!is.finite(log_density))) {
    if (strict) {
      stop(
        "The model produced non-finite bivariate Student-t ",
        "log-density values."
      )
    }

    return(NULL)
  }

  weight <-
    (nu + d) /
    (nu + quadratic_form)

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

  scale_array <- array(
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

  scale_array[1L, 1L, ] <- scale_y
  scale_array[1L, 2L, ] <- scale_yx
  scale_array[2L, 1L, ] <- scale_yx
  scale_array[2L, 2L, ] <- scale_x

  inverse_array <- array(
    NA_real_,
    dim = c(
      2L,
      2L,
      data$n
    ),
    dimnames = dimnames(scale_array)
  )

  inverse_array[1L, 1L, ] <- inverse_11
  inverse_array[1L, 2L, ] <- inverse_12
  inverse_array[2L, 1L, ] <- inverse_12
  inverse_array[2L, 2L, ] <- inverse_22

  covariance_array <- if (nu > 2) {
    nu /
      (nu - 2) *
      scale_array
  } else {
    NULL
  }

  conditional_location_y_given_x <-
    mean_y +
    scale_yx /
    scale_x *
    residual_x

  conditional_scale2_y_given_x <-
    (
      nu +
      residual_x^2 /
      scale_x
    ) /
    (nu + 1) *
    (
      scale_y -
      scale_yx^2 /
      scale_x
    )

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
    scale.y = scale_y,
    scale.x = scale_x,
    scale.yx = scale_yx,
    scale.array = scale_array,
    scale.inverse.array = inverse_array,
    covariance.array = covariance_array,
    determinant = determinant,
    quadratic.form = quadratic_form,
    weight = weight,
    log.density = log_density,
    density = exp(log_density),
    sigma2.x.scale = unname(sigma2_x_scale),
    sigma2.scale = unname(sigma2_scale),
    tau.y.scale = tau_y_scale,
    tau.x.scale = tau_x_scale,
    conditional.observed.location.y.given.x =
      conditional_location_y_given_x,
    conditional.observed.scale2.y.given.x =
      conditional_scale2_y_given_x,
    conditional.observed.df.y.given.x = nu + 1
  )
}


# ============================================================
# Main fitting function
# ============================================================

fit_lmve_t <- function(Y,
                       X,
                       tau_y,
                       tau_x,
                       start = NULL,
                       nu,
                       level = 0.95,
                       method = c(
                         "BFGS",
                         "L-BFGS-B",
                         "Nelder-Mead",
                         "CG",
                         "BOBYQA"
                       ),
                       lower = NULL,
                       upper = NULL,
                       control = list(),
                       use_score = TRUE,
                       compute_vcov = TRUE,
                       hessian = TRUE,
                       derivative_eps = .Machine$double.eps^(1 / 3)) {
  call <- match.call()

  method <- match.arg(method)

  if (length(nu) != 1L ||
      !is.numeric(nu) ||
      !is.finite(nu) ||
      nu <= 0) {
    stop(
      "'nu' must be a fixed positive number."
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

  if (length(derivative_eps) != 1L ||
      !is.numeric(derivative_eps) ||
      !is.finite(derivative_eps) ||
      derivative_eps <= 0) {
    stop(
      "'derivative_eps' must be positive."
    )
  }

  if (!is.list(control)) {
    stop("'control' must be a list.")
  }

  data <- .lmve_t_prepare_data(
    Y = Y,
    X = X,
    tau_y = tau_y,
    tau_x = tau_x
  )


  automatic_start <- is.null(start)

  start <- if (automatic_start) {
    .lmve_t_default_start(
      data = data,
      nu = nu
    )
  } else {
    .lmve_t_validate_start(start)
  }

  p <- length(start)

  initial_components <- .lmve_t_model_components(
    theta = start,
    data = data,
    nu = nu,
    strict = TRUE
  )

  initial_nll <- -sum(
    initial_components$log.density
  )

  if (!is.finite(initial_nll)) {
    stop(
      "The log-likelihood is not finite at the starting values."
    )
  }

  penalty_base <- max(
    1e8,
    min(
      1e100,
      1e6 *
        max(
          1,
          abs(initial_nll)
        )
    )
  )

  negative_loglikelihood <- function(theta) {
    components <- .lmve_t_model_components(
      theta = theta,
      data = data,
      nu = nu,
      strict = FALSE
    )

    if (is.null(components)) {
      return(
        penalty_base +
          sum(
            (theta - start)^2
          )
      )
    }

    value <- -sum(
      components$log.density
    )

    if (is.finite(value)) {
      value
    } else {
      penalty_base +
        sum(
          (theta - start)^2
        )
    }
  }

  score <- function(theta,
                    strict = FALSE) {
    components <- .lmve_t_model_components(
      theta = theta,
      data = data,
      nu = nu,
      strict = strict
    )

    if (is.null(components)) {
      return(
        setNames(
          rep(NaN, p),
          .lmve_t_parameter_names
        )
      )
    }

    derivatives <- .lmve_t_derivative_matrices(
      theta = theta
    )

    D <- derivatives$mean
    V <- derivatives$scale

    score_sum <- numeric(p)

    for (i in seq_len(data$n)) {
      Sigma_i <- components$scale.array[, , i]
      Sigma_inverse_i <-
        components$scale.inverse.array[, , i]

      z_i <- c(
        components$residual.y[i],
        components$residual.x[i]
      )

      weight_i <- components$weight[i]

      mean_score_i <- drop(
        Sigma_inverse_i %*%
          (weight_i * z_i)
      )

      scale_score_matrix_i <-
        0.5 *
        Sigma_inverse_i %*%
        (
          weight_i *
            tcrossprod(z_i) -
            Sigma_i
        ) %*%
        Sigma_inverse_i

      score_sum <-
        score_sum +
        drop(
          crossprod(
            D,
            mean_score_i
          )
        ) +
        drop(
          crossprod(
            V,
            c(scale_score_matrix_i)
          )
        )
    }

    setNames(
      score_sum,
      .lmve_t_parameter_names
    )
  }

  numerical_gradient <- function(fn,
                                 theta) {
    gradient <- numeric(
      length(theta)
    )

    for (j in seq_along(theta)) {
      step <-
        derivative_eps *
        max(
          1,
          abs(theta[j])
        )

      theta_plus <- theta
      theta_minus <- theta

      theta_plus[j] <-
        theta_plus[j] +
        step

      theta_minus[j] <-
        theta_minus[j] -
        step

      gradient[j] <-
        (
          fn(theta_plus) -
          fn(theta_minus)
        ) /
        (2 * step)
    }

    setNames(
      gradient,
      .lmve_t_parameter_names
    )
  }

  negative_score <- function(theta) {
    components <- .lmve_t_model_components(
      theta = theta,
      data = data,
      nu = nu,
      strict = FALSE
    )

    if (is.null(components)) {
      return(
        numerical_gradient(
          negative_loglikelihood,
          theta
        )
      )
    }

    answer <- tryCatch(
      -score(
        theta,
        strict = TRUE
      ),
      error = function(e) NULL
    )

    if (is.null(answer) ||
        anyNA(answer) ||
        any(!is.finite(answer))) {
      return(
        numerical_gradient(
          negative_loglikelihood,
          theta
        )
      )
    }

    answer
  }

  fisher_information <- function(theta,
                                 strict = TRUE) {
    components <- .lmve_t_model_components(
      theta = theta,
      data = data,
      nu = nu,
      strict = strict
    )

    derivatives <- .lmve_t_derivative_matrices(
      theta = theta
    )

    D <- derivatives$mean
    V <- derivatives$scale

    information <- matrix(
      0,
      nrow = p,
      ncol = p
    )

    d <- 2

    mean_constant <-
      (nu + d) /
      (nu + d + 2)

    for (i in seq_len(data$n)) {
      Sigma_inverse_i <-
        components$scale.inverse.array[, , i]

      scale_information_i <-
        (nu + d) /
        (
          2 *
          (nu + d + 2)
        ) *
        kronecker(
          Sigma_inverse_i,
          Sigma_inverse_i
        ) -
        1 /
        (
          2 *
          (nu + d + 2)
        ) *
        tcrossprod(
          c(Sigma_inverse_i)
        )

      information <-
        information +
        mean_constant *
        crossprod(
          D,
          Sigma_inverse_i %*% D
        ) +
        crossprod(
          V,
          scale_information_i %*% V
        )
    }

    information <-
      (
        information +
        t(information)
      ) / 2

    dimnames(information) <- list(
      .lmve_t_parameter_names,
      .lmve_t_parameter_names
    )

    list(
      information = information,
      mean.jacobian = D,
      scale.jacobian = V,
      F = derivatives$F
    )
  }

  # ----------------------------------------------------------
  # Optimization controls and bounds
  # ----------------------------------------------------------

  if (method == "L-BFGS-B") {
    default_control <- list(
      maxit = 2000,
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
  } else if (method == "BOBYQA") {
    if (!requireNamespace(
      "nloptr",
      quietly = TRUE
    )) {
      stop(
        "Package 'nloptr' is required when method = 'BOBYQA'. ",
        "Install it with install.packages('nloptr')."
      )
    }

    if (!is.null(control$maxit)) {
      if (is.null(control$maxeval)) {
        control$maxeval <- control$maxit
      }

      control$maxit <- NULL
    }

    if (!is.null(control$reltol)) {
      if (is.null(control$ftol_rel)) {
        control$ftol_rel <- control$reltol
      }

      control$reltol <- NULL
    }

    if (!is.null(control$abstol)) {
      if (is.null(control$ftol_abs)) {
        control$ftol_abs <- control$abstol
      }

      control$abstol <- NULL
    }

    control$factr <- NULL
    control$pgtol <- NULL

    default_control <- list(
      maxeval = 100000,
      xtol_rel = 1e-8,
      ftol_rel = 1e-8
    )
  } else {
    default_control <- list(
      maxit = 2000,
      reltol = 1e-9
    )

    control$factr <- NULL
    control$pgtol <- NULL
  }

  control <- utils::modifyList(
    default_control,
    control
  )

  if (method == "L-BFGS-B") {
    control$reltol <- NULL
    control$abstol <- NULL
  } else if (method != "BOBYQA") {
    control$factr <- NULL
    control$pgtol <- NULL
  }

  lower_used <- NULL
  upper_used <- NULL

  if (method %in% c(
    "L-BFGS-B",
    "BOBYQA"
  )) {
    if (is.null(lower)) {
      lower <- c(
        beta0 = -Inf,
        beta1 = -Inf,
        mu_x = -Inf,
        sigma2_x = 1e-10,
        sigma2 = 1e-10
      )
    }

    if (is.null(upper)) {
      upper <- setNames(
        rep(Inf, p),
        .lmve_t_parameter_names
      )
    }

    if (!is.null(names(lower))) {
      missing_lower <- setdiff(
        .lmve_t_parameter_names,
        names(lower)
      )

      if (length(missing_lower) > 0L) {
        stop(
          "'lower' must contain all parameter names."
        )
      }

      lower <- lower[
        .lmve_t_parameter_names
      ]
    }

    if (!is.null(names(upper))) {
      missing_upper <- setdiff(
        .lmve_t_parameter_names,
        names(upper)
      )

      if (length(missing_upper) > 0L) {
        stop(
          "'upper' must contain all parameter names."
        )
      }

      upper <- upper[
        .lmve_t_parameter_names
      ]
    }

    lower <- setNames(
      rep_len(
        as.numeric(lower),
        p
      ),
      .lmve_t_parameter_names
    )

    upper <- setNames(
      rep_len(
        as.numeric(upper),
        p
      ),
      .lmve_t_parameter_names
    )

    if (anyNA(lower) ||
        anyNA(upper)) {
      stop(
        "'lower' and 'upper' must not contain missing values."
      )
    }

    if (any(lower > upper)) {
      stop(
        "Each lower bound must be less than or equal to its upper bound."
      )
    }

    if (any(start < lower) ||
        any(start > upper)) {
      stop(
        "All starting values must lie within the specified bounds."
      )
    }

    lower_used <- lower
    upper_used <- upper
  } else if (!is.null(lower) ||
             !is.null(upper)) {
    stop(
      "Finite bounds require method = 'L-BFGS-B' or method = 'BOBYQA'."
    )
  }

  # ----------------------------------------------------------
  # Optimization
  # ----------------------------------------------------------

  if (method == "BOBYQA") {
    if (isTRUE(use_score)) {
      warning(
        "'use_score' is ignored when method = 'BOBYQA'."
      )
    }

    bobyqa_fit <- nloptr::bobyqa(
      x0 = start,
      fn = negative_loglikelihood,
      lower = lower,
      upper = upper,
      control = control
    )

    raw_convergence <- as.integer(
      bobyqa_fit$convergence
    )

    optimization <- list(
      par = bobyqa_fit$par,
      value = bobyqa_fit$value,
      counts = setNames(
        c(
          as.integer(bobyqa_fit$iter),
          NA_integer_
        ),
        c(
          "function",
          "gradient"
        )
      ),
      convergence = if (
        length(raw_convergence) == 1L &&
        !is.na(raw_convergence) &&
        raw_convergence > 0L
      ) {
        0L
      } else {
        raw_convergence
      },
      message = paste0(
        bobyqa_fit$message,
        " [NLopt status ",
        raw_convergence,
        "]"
      ),
      raw = bobyqa_fit
    )
  } else {
    optim_arguments <- list(
      par = start,
      fn = negative_loglikelihood,
      method = method,
      control = control,
      hessian = FALSE
    )

    gradient_methods <- c(
      "BFGS",
      "CG",
      "L-BFGS-B"
    )

    if (isTRUE(use_score) &&
        method %in% gradient_methods) {
      optim_arguments$gr <- negative_score
    } else if (isTRUE(use_score) &&
               method == "Nelder-Mead") {
      warning(
        "'use_score' is ignored when method = 'Nelder-Mead'."
      )
    }

    if (method == "L-BFGS-B") {
      optim_arguments$lower <- lower
      optim_arguments$upper <- upper
    }

    optimization <- do.call(
      stats::optim,
      optim_arguments
    )
  }

  estimate <- setNames(
    as.numeric(optimization$par),
    .lmve_t_parameter_names
  )

  fitted_components <- .lmve_t_model_components(
    theta = estimate,
    data = data,
    nu = nu,
    strict = TRUE
  )

  information_components <- fisher_information(
    theta = estimate,
    strict = TRUE
  )

  expected_information <-
    information_components$information

  score_at_estimate <- tryCatch(
    score(
      estimate,
      strict = TRUE
    ),
    error = function(e) {
      warning(
        "The score at the estimate could not be computed: ",
        conditionMessage(e)
      )

      setNames(
        rep(NA_real_, p),
        .lmve_t_parameter_names
      )
    }
  )

  # ----------------------------------------------------------
  # Expected and observed information matrices
  # ----------------------------------------------------------

  vcov_fisher <- NULL
  standard_error <- NULL
  confidence_interval <- NULL
  vcov_warning <- NULL

  if (isTRUE(compute_vcov)) {
    vcov_fisher <- .lmve_t_safe_inverse(
      expected_information
    )

    covariance_diagonal <- diag(
      vcov_fisher
    )

    standard_error <- setNames(
      ifelse(
        is.finite(covariance_diagonal) &
        covariance_diagonal >= 0,
        sqrt(covariance_diagonal),
        NA_real_
      ),
      .lmve_t_parameter_names
    )

    if (anyNA(vcov_fisher)) {
      vcov_warning <-
        "The expected Fisher information could not be inverted."
    }

    critical_value <- stats::qnorm(
      (1 + level) / 2
    )

    confidence_interval <- cbind(
      lower =
        estimate -
        critical_value *
        standard_error,
      upper =
        estimate +
        critical_value *
        standard_error
    )

    rownames(confidence_interval) <-
      .lmve_t_parameter_names
  }

  observed_information <- NULL
  vcov_observed <- NULL
  observed_information_warning <- NULL

  if (isTRUE(hessian)) {
    observed_information <- tryCatch(
      {
        Hessian <- stats::optimHess(
          par = estimate,
          fn = negative_loglikelihood,
          gr = if (isTRUE(use_score)) {
            negative_score
          } else {
            NULL
          }
        )

        Hessian <-
          (
            Hessian +
            t(Hessian)
          ) / 2

        dimnames(Hessian) <- list(
          .lmve_t_parameter_names,
          .lmve_t_parameter_names
        )

        Hessian
      },
      error = function(e) {
        observed_information_warning <<-
          paste(
            "The observed information could not be computed:",
            conditionMessage(e)
          )

        NULL
      }
    )

    if (!is.null(observed_information)) {
      vcov_observed <- .lmve_t_safe_inverse(
        observed_information
      )

      if (anyNA(vcov_observed)) {
        observed_information_warning <- paste(
          c(
            observed_information_warning,
            "The observed information could not be inverted."
          ),
          collapse = " "
        )
      }
    }
  }

  maximized_loglikelihood <-
    -unname(optimization$value)

  AIC_value <-
    -2 * maximized_loglikelihood +
    2 * p

  BIC_value <-
    -2 * maximized_loglikelihood +
    log(data$n) * p

  coefficient_table <- cbind(
    Estimate = estimate
  )

  if (!is.null(standard_error)) {
    z_value <-
      estimate /
      standard_error

    p_value <-
      2 *
      stats::pnorm(
        abs(z_value),
        lower.tail = FALSE
      )

    coefficient_table <- cbind(
      coefficient_table,
      `Std. Error` = standard_error,
      `z value` = z_value,
      `Pr(>|z|)` = p_value
    )
  }

  max_abs_score <- if (
    length(score_at_estimate) > 0L &&
    any(is.finite(score_at_estimate))
  ) {
    max(
      abs(
        score_at_estimate[
          is.finite(score_at_estimate)
        ]
      )
    )
  } else {
    NA_real_
  }

  latent_x_variance <- if (nu > 2) {
    nu /
      (nu - 2) *
      fitted_components$sigma2.x.scale
  } else {
    NA_real_
  }

  equation_error_variance <- if (nu > 2) {
    nu /
      (nu - 2) *
      fitted_components$sigma2.scale
  } else {
    NA_real_
  }

  result <- list(
    call = call,
    starting.values = start,
    automatic.start = automatic_start,
    coefficients = estimate,
    beta = estimate[c("beta0", "beta1")],
    mu.x = unname(estimate["mu_x"]),
    sigma2.x = unname(estimate["sigma2_x"]),
    sigma2 = unname(estimate["sigma2"]),
    sigma2.x.scale = fitted_components$sigma2.x.scale,
    sigma2.scale = fitted_components$sigma2.scale,
    latent.x.variance = latent_x_variance,
    equation.error.variance = equation_error_variance,
    standard.error = standard_error,
    conf.int = confidence_interval,
    coefficient.table = coefficient_table,
    vcov = vcov_fisher,
    vcov.fisher = vcov_fisher,
    vcov.observed = vcov_observed,
    fisher.information = expected_information,
    observed.information = observed_information,
    mean.jacobian = information_components$mean.jacobian,
    scale.jacobian = information_components$scale.jacobian,
    F = information_components$F,
    fitted.values = fitted_components$fitted.values,
    fitted.y = fitted_components$fitted.y,
    fitted.x = fitted_components$fitted.x,
    residuals = fitted_components$residuals,
    residual.y = fitted_components$residual.y,
    residual.x = fitted_components$residual.x,
    scale.y = fitted_components$scale.y,
    scale.x = fitted_components$scale.x,
    scale.yx = fitted_components$scale.yx,
    scale.array = fitted_components$scale.array,
    covariance.array = fitted_components$covariance.array,
    scale.inverse.array = fitted_components$scale.inverse.array,
    quadratic.forms = fitted_components$quadratic.form,
    weights = fitted_components$weight,
    log.density = fitted_components$log.density,
    density = fitted_components$density,
    conditional.observed.location.y.given.x =
      fitted_components$conditional.observed.location.y.given.x,
    conditional.observed.scale2.y.given.x =
      fitted_components$conditional.observed.scale2.y.given.x,
    conditional.observed.df.y.given.x =
      fitted_components$conditional.observed.df.y.given.x,
    score = score_at_estimate,
    max.abs.score = max_abs_score,
    logLik = maximized_loglikelihood,
    objective = unname(optimization$value),
    AIC = AIC_value,
    BIC = BIC_value,
    nu = nu,
    level = level,
    nobs = data$n,
    p = p,
    convergence = optimization$convergence,
    message = optimization$message,
    counts = optimization$counts,
    method = method,
    use.score = isTRUE(use_score),
    Y = data$Y,
    X = data$X,
    tau.y = data$tau_y,
    tau.x = data$tau_x,
    tau.y.scale = fitted_components$tau.y.scale,
    tau.x.scale = fitted_components$tau.x.scale,
    vcov.warning = vcov_warning,
    observed.information.warning = observed_information_warning,
    lower = lower_used,
    upper = upper_used,
    optimization.control = control,
    optim = optimization
  )

  class(result) <- "lmve_t_fit"

  if (!identical(
    as.integer(optimization$convergence),
    0L
  )) {
    warning(
      "The optimizer returned convergence code ",
      optimization$convergence,
      "."
    )
  }

  result
}


# ============================================================
# S3 methods
# ============================================================

print.lmve_t_fit <- function(x,
                             digits = max(
                               3L,
                               getOption("digits") - 3L
                             ),
                             ...) {
  cat(
    "Linear errors-in-variables model with a bivariate Student-t distribution\n"
  )

  cat(
    "fixed degrees of freedom =",
    format(x$nu, digits = digits),
    "| observations =",
    x$nobs,
    "| optimization method =",
    x$method,
    "\n"
  )

  cat(
    "optimization convergence code:",
    x$convergence
  )

  if (!is.null(x$message)) {
    cat(
      " -",
      x$message
    )
  }

  cat("\n\n")

  printCoefmat(
    x$coefficient.table,
    digits = digits,
    P.values = ncol(x$coefficient.table) >= 4L,
    has.Pvalue = ncol(x$coefficient.table) >= 4L,
    na.print = "NA"
  )

  cat(
    "\nLog-likelihood:",
    format(x$logLik, digits = digits),
    "\nAIC:",
    format(x$AIC, digits = digits),
    "\nBIC:",
    format(x$BIC, digits = digits),
    "\nMaximum absolute score:",
    format(x$max.abs.score, digits = digits),
    "\n"
  )

  invisible(x)
}


summary.lmve_t_fit <- function(object,
                               ...) {
  result <- list(
    call = object$call,
    coefficients = object$coefficient.table,
    logLik = object$logLik,
    AIC = object$AIC,
    BIC = object$BIC,
    score = object$score,
    max.abs.score = object$max.abs.score,
    convergence = object$convergence,
    message = object$message,
    nu = object$nu,
    nobs = object$nobs,
    method = object$method
  )

  class(result) <-
    "summary.lmve_t_fit"

  result
}


print.summary.lmve_t_fit <- function(x,
                                     digits = max(
                                       3L,
                                       getOption("digits") - 3L
                                     ),
                                     ...) {
  cat(
    "Linear errors-in-variables model fitted by Student-t maximum likelihood\n"
  )

  cat(
    "fixed degrees of freedom =",
    format(x$nu, digits = digits),
    "| n =",
    x$nobs,
    "| optimization method =",
    x$method,
    "\n\n"
  )

  printCoefmat(
    x$coefficients,
    digits = digits,
    P.values = ncol(x$coefficients) >= 4L,
    has.Pvalue = ncol(x$coefficients) >= 4L,
    na.print = "NA"
  )

  cat(
    "\nLog-likelihood:",
    format(x$logLik, digits = digits),
    "\nAIC:",
    format(x$AIC, digits = digits),
    "\nBIC:",
    format(x$BIC, digits = digits),
    "\nMaximum absolute score:",
    format(x$max.abs.score, digits = digits),
    "\nConvergence code:",
    x$convergence,
    "\n"
  )

  if (!is.null(x$message)) {
    cat(
      "Optimizer message:",
      x$message,
      "\n"
    )
  }

  invisible(x)
}


coef.lmve_t_fit <- function(object,
                            ...) {
  object$coefficients
}


vcov.lmve_t_fit <- function(object,
                            type = c(
                              "fisher",
                              "observed"
                            ),
                            ...) {
  type <- match.arg(type)

  if (type == "fisher") {
    if (is.null(object$vcov.fisher)) {
      stop(
        "The expected-information covariance matrix is not available."
      )
    }

    return(object$vcov.fisher)
  }

  if (is.null(object$vcov.observed)) {
    stop(
      "The observed-information covariance matrix is not available."
    )
  }

  object$vcov.observed
}


confint.lmve_t_fit <- function(object,
                               parm = seq_along(
                                 object$coefficients
                               ),
                               level = 0.95,
                               type = c(
                                 "fisher",
                                 "observed"
                               ),
                               ...) {
  type <- match.arg(type)

  if (length(level) != 1L ||
      !is.numeric(level) ||
      !is.finite(level) ||
      level <= 0 ||
      level >= 1) {
    stop(
      "'level' must be between zero and one."
    )
  }

  parameter_indices <- if (is.character(parm)) {
    match(
      parm,
      names(object$coefficients)
    )
  } else {
    parm
  }

  if (anyNA(parameter_indices)) {
    stop(
      "At least one requested parameter name was not found."
    )
  }

  covariance <- vcov(
    object,
    type = type
  )

  estimates <- object$coefficients[
    parameter_indices
  ]

  variances <- diag(covariance)[
    parameter_indices
  ]

  standard_errors <- ifelse(
    is.finite(variances) &
    variances >= 0,
    sqrt(variances),
    NA_real_
  )

  critical_value <- stats::qnorm(
    (1 + level) / 2
  )

  interval <- cbind(
    lower =
      estimates -
      critical_value *
      standard_errors,
    upper =
      estimates +
      critical_value *
      standard_errors
  )

  rownames(interval) <- names(estimates)

  interval
}


fitted.lmve_t_fit <- function(object,
                              component = c(
                                "both",
                                "Y",
                                "X"
                              ),
                              ...) {
  component <- match.arg(component)

  switch(
    component,
    both = object$fitted.values,
    Y = object$fitted.y,
    X = object$fitted.x
  )
}


residuals.lmve_t_fit <- function(object,
                                 component = c(
                                   "both",
                                   "Y",
                                   "X"
                                 ),
                                 ...) {
  component <- match.arg(component)

  switch(
    component,
    both = object$residuals,
    Y = object$residual.y,
    X = object$residual.x
  )
}


predict.lmve_t_fit <- function(object,
                               newx = NULL,
                               type = c(
                                 "structural",
                                 "marginal",
                                 "conditional_observed"
                               ),
                               tau_x = NULL,
                               ...) {
  type <- match.arg(type)

  beta0 <- object$coefficients["beta0"]
  beta1 <- object$coefficients["beta1"]
  mu_x <- object$coefficients["mu_x"]

  if (type == "marginal") {
    n_prediction <- if (is.null(newx)) {
      object$nobs
    } else {
      length(newx)
    }

    return(
      rep(
        beta0 + beta1 * mu_x,
        n_prediction
      )
    )
  }

  if (is.null(newx)) {
    newx <- object$X
  }

  newx <- as.numeric(newx)

  if (anyNA(newx) ||
      any(!is.finite(newx))) {
    stop(
      "'newx' must contain only finite, non-missing values."
    )
  }

  if (type == "structural") {
    return(
      beta0 +
      beta1 * newx
    )
  }

  if (is.null(tau_x)) {
    if (length(newx) == object$nobs) {
      tau_x <- object$tau.x
    } else {
      stop(
        "'tau_x' must be supplied for conditional predictions at ",
        "new observed X values."
      )
    }
  }

  tau_x <- rep_len(
    as.numeric(tau_x),
    length(newx)
  )

  if (anyNA(tau_x) ||
      any(!is.finite(tau_x)) ||
      any(tau_x < 0)) {
    stop(
      "'tau_x' must contain finite non-negative values."
    )
  }

  sigma2_x_scale <-
    object$coefficients["sigma2_x"]

  tau_x_scale <-
    tau_x

  scale_x <-
    sigma2_x_scale +
    tau_x_scale

  scale_yx <-
    beta1 *
    sigma2_x_scale

  beta0 +
    beta1 * mu_x +
    scale_yx /
    scale_x *
    (newx - mu_x)
}


logLik.lmve_t_fit <- function(object,
                              ...) {
  structure(
    object$logLik,
    df = object$p,
    nobs = object$nobs,
    class = "logLik"
  )
}


nobs.lmve_t_fit <- function(object,
                            ...) {
  object$nobs
}

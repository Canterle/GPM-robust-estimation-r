# Generic MLqE fitting for nonlinear normal models
#
# The model function must have the form
#
#   model(theta, data)
#
# and return a list containing:
#
#   mu     : conditional means;
#   sigma2 : conditional variances.
#
# Alternatively, the model may return sigma instead of sigma2.
# A scalar variance is automatically replicated, so homoscedastic
# models are included as a particular case.
#
# The model may also return the optional analytic derivative matrices:
#
#   D : Jacobian of mu with respect to the complete parameter vector;
#   V : Jacobian of sigma2 with respect to the complete parameter vector.
#
# Both matrices must have dimension n by p. When D or V is absent,
# only the missing matrix is computed numerically.
#
# The user supplies a single transformation family:
#
#   tau(theta, r)
#
# which returns tau_r(theta).
#
# If theta_star is the original MLqE estimate, the corrected estimate is
#
#   theta_hat = tau(theta_star, q).
#
# The forward transformation is obtained by replacing q with 1/q:
#
#   theta_star = tau(theta_hat, 1 / q).
#
# The corrected covariance matrix is computed as
#
#   V = A^{-1} V_star A^{-T},
#
# where A is the Jacobian of tau(theta, 1/q), evaluated at theta_hat.
# The equivalent direct form based on the Jacobian of tau(theta_star, q)
# is also computed for verification.


# ============================================================
# Internal utility functions
# ============================================================

.nlm_mlqe_validate_tau_output <- function(value,
                                          parameter_names,
                                          label) {
  value <- as.numeric(value)

  if (length(value) != length(parameter_names) ||
      anyNA(value) ||
      any(!is.finite(value))) {
    stop(
      "'", label, "' must return a finite numeric vector of length ",
      length(parameter_names),
      "."
    )
  }

  setNames(value, parameter_names)
}


.nlm_mlqe_call_tau <- function(tau_function,
                               theta,
                               q,
                               parameter_names,
                               label) {
  theta <- setNames(as.numeric(theta), parameter_names)

  value <- tryCatch(
    tau_function(theta, q),
    error = function(e) {
      stop(
        "'", label, "' failed: ",
        conditionMessage(e)
      )
    }
  )

  .nlm_mlqe_validate_tau_output(
    value = value,
    parameter_names = parameter_names,
    label = label
  )
}


.nlm_mlqe_compute_jacobian <- function(func,
                                       x,
                                       method,
                                       method_args,
                                       label) {
  result <- tryCatch(
    do.call(
      numDeriv::jacobian,
      c(
        list(
          func = func,
          x = x,
          method = method
        ),
        list(method.args = method_args)
      )
    ),
    error = function(e) {
      stop(
        label,
        " could not be computed: ",
        conditionMessage(e)
      )
    }
  )

  result <- as.matrix(result)

  if (anyNA(result) || any(!is.finite(result))) {
    stop(label, " contains non-finite values.")
  }

  result
}


.nlm_mlqe_invert_matrix <- function(matrix_value,
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

.fit_nlm_mlqe_fixed <- function(y,
                                model,
                                start,
                                data = NULL,
                                q = 1,
                                tau = NULL,
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
                                jacobian_method = c(
                                  "Richardson",
                                  "simple"
                                ),
                                jacobian_method_args = list()) {
  force(y)
  force(model)
  force(start)
  force(data)
  force(q)
  force(tau)

  call <- match.call()
  method <- match.arg(method)
  jacobian_method <- match.arg(jacobian_method)

  if (!requireNamespace("numDeriv", quietly = TRUE)) {
    stop(
      "Package 'numDeriv' is required. Install it with ",
      "install.packages('numDeriv')."
    )
  }

  if (!is.function(model)) {
    stop("'model' must be a function.")
  }

  if (length(q) != 1L ||
      !is.numeric(q) ||
      !is.finite(q) ||
      q <= 0 ||
      q > 1) {
    stop("'q' must be in the interval (0, 1].")
  }

  if (q < 1 && !is.function(tau)) {
    stop(
      "For q < 1, 'tau' must be supplied as a function."
    )
  }

  if (!is.null(tau) && !is.function(tau)) {
    stop("'tau' must be NULL or a function.")
  }

  original_parameter_names <- names(start)
  y <- as.numeric(y)
  start <- as.numeric(start)

  n <- length(y)
  p <- length(start)

  if (n == 0L) {
    stop("'y' must not be empty.")
  }
  if (p == 0L) {
    stop("'start' must contain at least one parameter.")
  }
  if (anyNA(y) || any(!is.finite(y))) {
    stop(
      "'y' must contain only finite, non-missing values."
    )
  }
  if (anyNA(start) || any(!is.finite(start))) {
    stop(
      "'start' must contain only finite, non-missing values."
    )
  }
  if (length(level) != 1L ||
      !is.finite(level) ||
      level <= 0 ||
      level >= 1) {
    stop("'level' must be in the interval (0, 1).")
  }
  if (!is.list(control)) {
    stop("'control' must be a list.")
  }
  if (!is.list(jacobian_method_args)) {
    stop("'jacobian_method_args' must be a list.")
  }

  parameter_names <- original_parameter_names

  if (is.null(parameter_names) ||
      any(parameter_names == "")) {
    parameter_names <- paste0(
      "theta",
      seq_len(p)
    )
  }

  parameter_names <- make.unique(parameter_names)
  names(start) <- parameter_names

  call_model <- function(theta) {
    theta <- setNames(
      as.numeric(theta),
      parameter_names
    )

    model(
      theta = theta,
      data = data
    )
  }

  evaluate_model <- function(theta,
                             strict = FALSE) {
    output <- tryCatch(
      call_model(theta),
      error = function(e) e
    )

    if (inherits(output, "error")) {
      if (strict) {
        stop(
          "The model function failed: ",
          conditionMessage(output)
        )
      }

      return(NULL)
    }

    if (!is.list(output) ||
        is.null(output$mu)) {
      if (strict) {
        stop(
          "The model function must return a list ",
          "containing 'mu'."
        )
      }

      return(NULL)
    }

    has_sigma2 <- !is.null(output$sigma2)
    has_sigma <- !is.null(output$sigma)

    if (!has_sigma2 && !has_sigma) {
      if (strict) {
        stop(
          "The model function must return either ",
          "'sigma2' or 'sigma'."
        )
      }

      return(NULL)
    }

    mu <- as.numeric(output$mu)

    if (has_sigma2) {
      sigma2 <- as.numeric(output$sigma2)
    } else {
      sigma2 <- as.numeric(output$sigma)^2
    }

    if (length(mu) == 1L) {
      mu <- rep(mu, n)
    }
    if (length(sigma2) == 1L) {
      sigma2 <- rep(sigma2, n)
    }

    valid <- length(mu) == n &&
      length(sigma2) == n &&
      !anyNA(mu) &&
      !anyNA(sigma2) &&
      all(is.finite(mu)) &&
      all(is.finite(sigma2)) &&
      all(sigma2 > 0)

    if (!valid) {
      if (strict) {
        stop(
          "The model must produce finite vectors 'mu' and ",
          "'sigma2' of length ",
          n,
          ", with all variances strictly positive."
        )
      }

      return(NULL)
    }

    sigma <- sqrt(sigma2)
    residuals <- y - mu
    log_density <- dnorm(
      y,
      mean = mu,
      sd = sigma,
      log = TRUE
    )

    if (any(!is.finite(log_density))) {
      if (strict) {
        stop(
          "The model produced non-finite normal ",
          "log-density values."
        )
      }

      return(NULL)
    }

    list(
      mu = mu,
      sigma2 = sigma2,
      sigma = sigma,
      residuals = residuals,
      log_density = log_density,
      raw = output
    )
  }

  objective <- function(theta_star) {
    components <- evaluate_model(theta_star)

    if (is.null(components)) {
      return(-1e100)
    }

    if (q == 1) {
      value <- sum(components$log_density)
    } else {
      one_minus_q <- 1 - q

      value <- sum(
        expm1(
          one_minus_q *
            components$log_density
        ) /
          one_minus_q
      )
    }

    if (is.finite(value)) {
      value
    } else {
      -1e100
    }
  }

  model_derivatives <- function(theta_star) {
    theta_star <- setNames(
      as.numeric(theta_star),
      parameter_names
    )

    components <- evaluate_model(
      theta_star,
      strict = TRUE
    )

    validate_derivative_matrix <- function(value,
                                           label) {
      value <- as.matrix(value)

      if (!all(dim(value) == c(n, p))) {
        stop(
          label,
          " must be an ",
          n,
          " by ",
          p,
          " matrix."
        )
      }

      if (anyNA(value) ||
          any(!is.finite(value))) {
        stop(
          label,
          " must contain only finite, non-missing values."
        )
      }

      value
    }

    mean_function <- function(par) {
      evaluated <- evaluate_model(
        par,
        strict = TRUE
      )

      evaluated$mu
    }

    variance_function <- function(par) {
      evaluated <- evaluate_model(
        par,
        strict = TRUE
      )

      evaluated$sigma2
    }

    D <- if (!is.null(components$raw$D)) {
      validate_derivative_matrix(
        value = components$raw$D,
        label = "'D' returned by the model"
      )
    } else {
      validate_derivative_matrix(
        value = .nlm_mlqe_compute_jacobian(
          func = mean_function,
          x = theta_star,
          method = jacobian_method,
          method_args = jacobian_method_args,
          label = "The numerical Jacobian of mu"
        ),
        label = "The numerical Jacobian of mu"
      )
    }

    V <- if (!is.null(components$raw$V)) {
      validate_derivative_matrix(
        value = components$raw$V,
        label = "'V' returned by the model"
      )
    } else {
      validate_derivative_matrix(
        value = .nlm_mlqe_compute_jacobian(
          func = variance_function,
          x = theta_star,
          method = jacobian_method,
          method_args = jacobian_method_args,
          label = "The numerical Jacobian of sigma2"
        ),
        label = "The numerical Jacobian of sigma2"
      )
    }

    colnames(D) <- parameter_names
    colnames(V) <- parameter_names

    list(
      D = D,
      V = V
    )
  }

  score <- function(theta_star) {
    components <- evaluate_model(
      theta_star,
      strict = TRUE
    )

    derivatives <- model_derivatives(
      theta_star
    )

    density_power <- if (q == 1) {
      rep(1, n)
    } else {
      exp(
        (1 - q) *
          components$log_density
      )
    }

    mean_weight <-
      components$residuals *
      density_power /
      components$sigma2

    variance_weight <-
      (
        components$residuals^2 -
          components$sigma2
      ) *
      density_power /
      (
        2 *
          components$sigma2^2
      )

    score_value <-
      colSums(
        sweep(
          derivatives$D,
          1L,
          mean_weight,
          FUN = "*"
        )
      ) +
      colSums(
        sweep(
          derivatives$V,
          1L,
          variance_weight,
          FUN = "*"
        )
      )

    setNames(
      score_value,
      parameter_names
    )
  }

  compute_JK <- function(theta_star) {
    components <- evaluate_model(
      theta_star,
      strict = TRUE
    )

    derivatives <- model_derivatives(
      theta_star
    )

    sigma2 <- components$sigma2
    inverse_sigma2 <- 1 / sigma2
    one_minus_q <- 1 - q

    c1 <- exp(
      0.5 * log(q) -
        0.5 *
          one_minus_q *
          log(2 * pi * sigma2)
    )

    mean_weight_J <-
      q *
      c1 *
      inverse_sigma2

    variance_weight_J <-
      q *
      c1 *
      inverse_sigma2^2 /
      2

    J_star <-
      crossprod(
        derivatives$D,
        sweep(
          derivatives$D,
          1L,
          mean_weight_J,
          FUN = "*"
        )
      ) +
      crossprod(
        derivatives$V,
        sweep(
          derivatives$V,
          1L,
          variance_weight_J,
          FUN = "*"
        )
      )

    c_q <- 4 *
      (2 - q) /
      (
        2 +
          one_minus_q^2
      )

    c2 <- c1^2 /
      sqrt(
        q *
          (2 - q)
      )

    mean_weight_K <-
      c2 *
      inverse_sigma2 /
      (2 - q)

    variance_weight_K <-
      c2 *
      inverse_sigma2^2 /
      (
        c_q *
          (2 - q)
      )

    K_star <-
      crossprod(
        derivatives$D,
        sweep(
          derivatives$D,
          1L,
          mean_weight_K,
          FUN = "*"
        )
      ) +
      crossprod(
        derivatives$V,
        sweep(
          derivatives$V,
          1L,
          variance_weight_K,
          FUN = "*"
        )
      )

    dimnames(J_star) <- list(
      parameter_names,
      parameter_names
    )

    dimnames(K_star) <- list(
      parameter_names,
      parameter_names
    )

    list(
      J = J_star,
      K = K_star,
      D = derivatives$D,
      V = derivatives$V
    )
  }

  evaluate_model(
    start,
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
    if (is.null(lower)) {
      lower <- rep(-Inf, p)
    }
    if (is.null(upper)) {
      upper <- rep(Inf, p)
    }

    lower <- rep_len(
      as.numeric(lower),
      p
    )

    upper <- rep_len(
      as.numeric(upper),
      p
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

  estimate_star <- setNames(
    as.numeric(optimization$par),
    parameter_names
  )

  tau_family <- if (is.null(tau)) {
    function(theta, r) theta
  } else {
    tau
  }

  # Consistency correction: theta_hat = tau_q(theta_star_hat).
  estimate <- .nlm_mlqe_call_tau(
    tau_function = tau_family,
    theta = estimate_star,
    q = q,
    parameter_names = parameter_names,
    label = "tau"
  )

  # Check tau_{1/q}(tau_q(theta_star)) = theta_star.
  mapped_back <- .nlm_mlqe_call_tau(
    tau_function = tau_family,
    theta = estimate,
    q = 1 / q,
    parameter_names = parameter_names,
    label = "tau"
  )

  tau_composition_residual <- max(
    abs(mapped_back - estimate_star) /
      pmax(1, abs(estimate_star))
  )

  if (!is.finite(tau_composition_residual) ||
      tau_composition_residual > 1e-6) {
    warning(
      "The supplied transformation does not accurately satisfy ",
      "tau(tau(theta, q), 1/q) = theta. Maximum scaled residual: ",
      format(tau_composition_residual, digits = 8),
      "."
    )
  }

  # A = d tau_{1/q}(theta) / d theta^T,
  # evaluated at the corrected estimate.
  tau_jacobian <- .nlm_mlqe_compute_jacobian(
    func = function(theta) {
      .nlm_mlqe_call_tau(
        tau_function = tau_family,
        theta = theta,
        q = 1 / q,
        parameter_names = parameter_names,
        label = "tau"
      )
    },
    x = estimate,
    method = jacobian_method,
    method_args = jacobian_method_args,
    label = "The Jacobian of tau(theta, 1/q)"
  )

  if (!all(dim(tau_jacobian) == c(p, p))) {
    stop(
      "The Jacobian of tau(theta, 1/q) must be a ",
      p,
      " by ",
      p,
      " matrix."
    )
  }

  dimnames(tau_jacobian) <- list(
    parameter_names,
    parameter_names
  )

  tau_jacobian_inverse <- .nlm_mlqe_invert_matrix(
    matrix_value = tau_jacobian,
    label = "The Jacobian of tau(theta, 1/q)"
  )

  dimnames(tau_jacobian_inverse) <- list(
    parameter_names,
    parameter_names
  )

  # B = d tau_q(theta_star) / d theta_star^T,
  # evaluated at the original MLqE estimate.
  tau_q_jacobian <- .nlm_mlqe_compute_jacobian(
    func = function(theta_star) {
      .nlm_mlqe_call_tau(
        tau_function = tau_family,
        theta = theta_star,
        q = q,
        parameter_names = parameter_names,
        label = "tau"
      )
    },
    x = estimate_star,
    method = jacobian_method,
    method_args = jacobian_method_args,
    label = "The Jacobian of tau(theta, q)"
  )

  if (!all(dim(tau_q_jacobian) == c(p, p))) {
    stop(
      "The Jacobian of tau(theta, q) must be a ",
      p,
      " by ",
      p,
      " matrix."
    )
  }

  dimnames(tau_q_jacobian) <- list(
    parameter_names,
    parameter_names
  )

  fitted_components_star <- evaluate_model(
    estimate_star,
    strict = TRUE
  )

  fitted_components <- evaluate_model(
    estimate,
    strict = TRUE
  )

  score_at_estimate_star <- tryCatch(
    score(estimate_star),
    error = function(e) {
      warning(
        "The score at the original MLqE estimate could not be computed: ",
        conditionMessage(e)
      )

      setNames(
        rep(NA_real_, p),
        parameter_names
      )
    }
  )

  vcov_star <- NULL
  vcov_matrix <- NULL
  standard_error_star <- NULL
  standard_error <- NULL
  confidence_interval_star <- NULL
  confidence_interval <- NULL
  J_star <- NULL
  K_star <- NULL
  D_star <- NULL
  V_star_derivative <- NULL
  vcov_warning <- NULL
  vcov_direct <- NULL
  vcov_direct_difference <- NULL

  if (compute_vcov) {
    jk <- compute_JK(estimate_star)

    J_star <- jk$J
    K_star <- jk$K
    D_star <- jk$D
    V_star_derivative <- jk$V

    covariance_result <- tryCatch(
      {
        J_inverse <- chol2inv(
          chol(J_star)
        )

        covariance_star <-
          J_inverse %*%
          K_star %*%
          t(J_inverse)

        covariance_star <-
          (
            covariance_star +
              t(covariance_star)
          ) /
          2

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

            covariance_star <-
              (
                covariance_star +
                  t(covariance_star)
              ) /
              2

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
                  parameter_names,
                  parameter_names
                )
              ),
              warning = paste(
                "The covariance matrix could not be computed:",
                conditionMessage(inverse_error)
              )
            )
          }
        )
      }
    )

    vcov_star <- covariance_result$value
    dimnames(vcov_star) <- list(
      parameter_names,
      parameter_names
    )

    vcov_warning <- covariance_result$warning

    if (all(is.finite(vcov_star))) {
      vcov_matrix <-
        tau_jacobian_inverse %*%
        vcov_star %*%
        t(tau_jacobian_inverse)

      vcov_matrix <-
        (
          vcov_matrix +
            t(vcov_matrix)
        ) /
        2
    } else {
      vcov_matrix <- matrix(
        NA_real_,
        nrow = p,
        ncol = p
      )
    }

    dimnames(vcov_matrix) <- list(
      parameter_names,
      parameter_names
    )

    vcov_direct <-
      tau_q_jacobian %*%
      vcov_star %*%
      t(tau_q_jacobian)

    vcov_direct <-
      (
        vcov_direct +
          t(vcov_direct)
      ) /
      2

    dimnames(vcov_direct) <- list(
      parameter_names,
      parameter_names
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
      rep(
        NA_real_,
        p
      ),
      parameter_names
    )

    valid_diagonal_star <-
      !is.na(diagonal_star) &
      is.finite(diagonal_star) &
      diagonal_star >= 0

    standard_error_star[
      valid_diagonal_star
    ] <- sqrt(
      diagonal_star[
        valid_diagonal_star
      ]
    )

    standard_error <- setNames(
      rep(
        NA_real_,
        p
      ),
      parameter_names
    )

    valid_diagonal_corrected <-
      !is.na(diagonal_corrected) &
      is.finite(diagonal_corrected) &
      diagonal_corrected >= 0

    standard_error[
      valid_diagonal_corrected
    ] <- sqrt(
      diagonal_corrected[
        valid_diagonal_corrected
      ]
    )

    z_value <- qnorm(
      1 -
        (1 - level) /
          2
    )

    confidence_interval_star <- cbind(
      lower =
        estimate_star -
        z_value *
          standard_error_star,
      upper =
        estimate_star +
        z_value *
          standard_error_star
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

    rownames(confidence_interval_star) <-
      parameter_names

    rownames(confidence_interval) <-
      parameter_names
  }

  result <- list(
    call = call,
    coefficients = estimate,
    coefficients.star = estimate_star,
    standard.error = standard_error,
    standard.error.star = standard_error_star,
    vcov = vcov_matrix,
    vcov.star = vcov_star,
    vcov.direct = vcov_direct,
    vcov.direct.difference = vcov_direct_difference,
    conf.int = confidence_interval,
    conf.int.star = confidence_interval_star,
    J = J_star,
    K = K_star,
    J.star = J_star,
    K.star = K_star,
    mean.jacobian.star = D_star,
    variance.jacobian.star = V_star_derivative,
    tau.jacobian = tau_jacobian,
    tau.jacobian.inverse = tau_jacobian_inverse,
    tau.q.jacobian = tau_q_jacobian,
    tau.composition.residual = tau_composition_residual,
    fitted.values = fitted_components$mu,
    sigma = fitted_components$sigma,
    sigma2 = fitted_components$sigma2,
    residuals = fitted_components$residuals,
    fitted.values.star = fitted_components_star$mu,
    sigma.star = fitted_components_star$sigma,
    sigma2.star = fitted_components_star$sigma2,
    residuals.star = fitted_components_star$residuals,
    objective = unname(optimization$value),
    score = score_at_estimate_star,
    score.star = score_at_estimate_star,
    convergence = optimization$convergence,
    message = optimization$message,
    counts = optimization$counts,
    q = q,
    level = level,
    nobs = n,
    method = method,
    jacobian.method = jacobian_method,
    model = model,
    data = data,
    tau = tau,
    vcov.warning = vcov_warning,
    optimization.control = control,
    optim = optimization,
    implementation.version =
      "MLqE-SQV-analytic-DV-final-2026-07-27"
  )

  class(result) <- "nlm_mlqe"
  result
}


# ============================================================
# Public fitting interface
# ============================================================

fit_nlm_mlqe <- function(y,
                         model,
                         start,
                         data = NULL,
                         q = 1,
                         tau = NULL,
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
                         jacobian_method = c(
                           "Richardson",
                           "simple"
                         ),
                         jacobian_method_args = list(),
                         q_control = list()) {
  force(y)
  force(model)
  force(start)
  force(data)
  force(q)
  force(tau)

  call <- match.call()

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
    fit <- .fit_nlm_mlqe_fixed(
      y = y,
      model = model,
      start = start,
      data = data,
      q = q,
      tau = tau,
      level = level,
      method = method,
      lower = lower,
      upper = upper,
      control = control,
      use_score = use_score,
      compute_vcov = compute_vcov,
      jacobian_method = jacobian_method,
      jacobian_method_args = jacobian_method_args
    )

    fit$call <- call
    return(fit)
  }

  if (!is.function(tau)) {
    stop(
      "SQV selection requires a user-supplied 'tau' function."
    )
  }

  if (!isTRUE(compute_vcov)) {
    stop(
      "SQV selection requires 'compute_vcov = TRUE' because ",
      "standard errors are used in the standardized estimates."
    )
  }

  fit <- .select_q_nlm_mlqe_sqv(
    y = y,
    model = model,
    start = start,
    data = data,
    tau = tau,
    level = level,
    method = method,
    lower = lower,
    upper = upper,
    control = control,
    use_score = use_score,
    jacobian_method = jacobian_method,
    jacobian_method_args = jacobian_method_args,
    q_control = q_control
  )

  fit$call <- call
  fit
}


# ============================================================
# SQV selection algorithm
# ============================================================

.select_q_nlm_mlqe_sqv <- function(y,
                                   model,
                                   start,
                                   data,
                                   tau,
                                   level,
                                   method,
                                   lower,
                                   upper,
                                   control,
                                   use_score,
                                   jacobian_method,
                                   jacobian_method_args,
                                   q_control) {
  force(y)
  force(model)
  force(start)
  force(data)
  force(tau)
  force(level)
  force(method)
  force(lower)
  force(upper)
  force(control)
  force(use_score)
  force(jacobian_method)
  force(jacobian_method_args)
  force(q_control)

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

  jacobian_method <- match.arg(
    jacobian_method,
    c(
      "Richardson",
      "simple"
    )
  )

  original_parameter_names <- names(start)
  start <- as.numeric(start)
  p <- length(start)
  n <- length(y)

  if (p == 0L) {
    stop("'start' must contain at least one parameter.")
  }

  parameter_names <- original_parameter_names

  if (is.null(parameter_names) ||
      any(parameter_names == "")) {
    parameter_names <- paste0(
      "theta",
      seq_len(p)
    )
  }

  parameter_names <- make.unique(
    parameter_names
  )

  names(start) <- parameter_names

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

  q_min <- as.numeric(
    q_control$q_min
  )

  L <- as.numeric(
    q_control$L
  )

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
      .fit_nlm_mlqe_fixed(
        y = y,
        model = model,
        start = start,
        data = data,
        q = q_value,
        tau = tau,
        level = level,
        method = method,
        lower = lower,
        upper = upper,
        control = control,
        use_score = use_score,
        compute_vcov = TRUE,
        jacobian_method = jacobian_method,
        jacobian_method_args =
          jacobian_method_args
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
    # For MLqE, coefficients and standard.error are the corrected
    # estimates and their corrected standard errors.
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

    colnames(z_matrix) <- parameter_names

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

  selected_fit$implementation.version <-
    "NLM-MLqE-SQV-original-algorithm-2026-07-30"

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
    "nlm_mlqe_sqv",
    class(selected_fit)
  )

  selected_fit
}


# ============================================================
# S3 methods
# ============================================================

print.nlm_mlqe <- function(x,
                           digits = max(
                             3L,
                             getOption("digits") - 3L
                           ),
                           ...) {
  cat(
    "Nonlinear normal model fitted by MLqE\n"
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
    "| Jacobian method =",
    x$jacobian.method,
    "\n"
  )

  cat(
    "optim convergence code:",
    x$convergence
  )

  if (!is.null(x$message)) {
    cat(
      " -",
      x$message
    )
  }

  cat(
    "\nConsistency correction: tau(theta_star, q)",
    "\n\n"
  )

  coefficient_table <- cbind(
    Estimate = x$coefficients
  )

  if (!is.null(x$standard.error)) {
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

  if (!is.null(x$vcov.warning)) {
    cat(
      "\nWarning:",
      x$vcov.warning,
      "\n"
    )
  }

  invisible(x)
}


summary.nlm_mlqe <- function(object,
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
      `z value` =
        z_value,
      `Pr(>|z|)` =
        p_value
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
    q =
      object$q,
    nobs =
      object$nobs,
    method =
      object$method,
    jacobian.method =
      object$jacobian.method,
    vcov.warning =
      object$vcov.warning
  )

  class(result) <-
    "summary.nlm_mlqe"

  result
}


print.summary.nlm_mlqe <- function(x,
                                   digits = max(
                                     3L,
                                     getOption("digits") - 3L
                                   ),
                                   ...) {
  cat(
    "Nonlinear normal model fitted by MLqE\n"
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
    "| Jacobian method =",
    x$jacobian.method,
    "\n"
  )

  cat(
    "optim convergence code:",
    x$convergence
  )

  if (!is.null(x$message)) {
    cat(
      " -",
      x$message
    )
  }

  cat(
    "\nConsistency correction: tau(theta_star, q)",
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


coef.nlm_mlqe <- function(object,
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


vcov.nlm_mlqe <- function(object,
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


confint.nlm_mlqe <- function(object,
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
    1 -
      (1 - level) /
      2
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


fitted.nlm_mlqe <- function(object,
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


residuals.nlm_mlqe <- function(object,
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


predict.nlm_mlqe <- function(object,
                             newdata = NULL,
                             type = c(
                               "mean",
                               "sigma",
                               "variance"
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

  prediction_data <- if (
    is.null(newdata)
  ) {
    object$data
  } else {
    newdata
  }

  theta <- if (
    parameterization == "corrected"
  ) {
    object$coefficients
  } else {
    object$coefficients.star
  }

  model_output <- object$model(
    theta,
    prediction_data
  )

  if (!is.list(model_output) ||
      is.null(model_output$mu)) {
    stop(
      "The model function must return a list containing 'mu'."
    )
  }

  if (type == "mean") {
    return(
      as.numeric(model_output$mu)
    )
  }

  if (!is.null(model_output$sigma2)) {
    sigma2 <- as.numeric(
      model_output$sigma2
    )
  } else if (!is.null(model_output$sigma)) {
    sigma2 <- as.numeric(
      model_output$sigma
    )^2
  } else {
    stop(
      "The model function must return either 'sigma2' or 'sigma'."
    )
  }

  if (type == "sigma") {
    sqrt(sigma2)
  } else {
    sigma2
  }
}


print.nlm_mlqe_sqv <- function(x,
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


summary.nlm_mlqe_sqv <- function(object,
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
    "summary.nlm_mlqe_sqv",
    class(result)
  )

  result
}


print.summary.nlm_mlqe_sqv <- function(x,
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

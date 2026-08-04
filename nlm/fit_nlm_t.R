# ============================================================
# Maximum likelihood estimation for Student-t regression
#
# Model:
#   Y_i ~ t_nu(mu_i, sigma2_i)
#
# where sigma2_i is the squared scale parameter of the
# Student-t distribution, not Var(Y_i).
#
# The user-supplied model must have the form:
#
#   model <- function(theta, data) {
#     ...
#     list(
#       mu     = ...,
#       sigma2 = ...
#     )
#   }
#
# Optionally, the model may also return:
#
#   D = d mu / d theta^T
#   V = d sigma2 / d theta^T
#
# Both matrices must have dimension n x p, where p is the
# length of theta. If D or V is not supplied, only the missing
# matrix is computed numerically.
# ============================================================

fit_nlm_t <- function(
    y,
    model,
    start_theta,
    data = NULL,
    nu,
    method = c(
      "BFGS",
      "L-BFGS-B",
      "Nelder-Mead",
      "CG"
    ),
    lower = -Inf,
    upper = Inf,
    control = list(),
    use_gradient = TRUE,
    hessian = TRUE,
    derivative_eps = .Machine$double.eps^(1 / 3),
    conf_level = 0.95
) {
  call <- match.call()
  method <- match.arg(method)

  # ----------------------------------------------------------
  # Basic checks
  # ----------------------------------------------------------

  if (!is.numeric(y) ||
      length(y) == 0L ||
      anyNA(y) ||
      any(!is.finite(y))) {
    stop(
      "'y' must be a finite, nonempty numeric vector."
    )
  }

  y <- as.numeric(y)
  n <- length(y)

  if (!is.function(model)) {
    stop("'model' must be a function.")
  }

  if (length(nu) != 1L ||
      !is.numeric(nu) ||
      !is.finite(nu) ||
      nu <= 0) {
    stop("'nu' must be a fixed positive number.")
  }

  if (!is.numeric(start_theta) ||
      length(start_theta) == 0L ||
      anyNA(start_theta) ||
      any(!is.finite(start_theta))) {
    stop(
      "'start_theta' must be a finite, nonempty numeric vector."
    )
  }

  if (length(conf_level) != 1L ||
      !is.numeric(conf_level) ||
      !is.finite(conf_level) ||
      conf_level <= 0 ||
      conf_level >= 1) {
    stop(
      "'conf_level' must be a number between zero and one."
    )
  }

  if (length(derivative_eps) != 1L ||
      !is.numeric(derivative_eps) ||
      !is.finite(derivative_eps) ||
      derivative_eps <= 0) {
    stop("'derivative_eps' must be positive.")
  }

  if (!is.list(control)) {
    stop("'control' must be a list.")
  }

  original_parameter_names <- names(start_theta)
  start_theta <- as.numeric(start_theta)
  p <- length(start_theta)

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

  start_theta <- setNames(
    start_theta,
    parameter_names
  )

  # ----------------------------------------------------------
  # Bounds
  # ----------------------------------------------------------

  expand_bound <- function(x,
                           p,
                           argument_name) {
    if (!is.numeric(x)) {
      stop(
        "'",
        argument_name,
        "' must be numeric."
      )
    }

    if (length(x) == 1L) {
      x <- rep(
        x,
        p
      )
    }

    if (length(x) != p) {
      stop(
        "'",
        argument_name,
        "' must have length one or length ",
        p,
        "."
      )
    }

    if (anyNA(x)) {
      stop(
        "'",
        argument_name,
        "' must not contain missing values."
      )
    }

    as.numeric(x)
  }

  lower <- expand_bound(
    lower,
    p,
    "lower"
  )

  upper <- expand_bound(
    upper,
    p,
    "upper"
  )

  if (any(lower > upper)) {
    stop(
      "Each lower bound must be less than or equal to ",
      "its upper bound."
    )
  }

  if (any(start_theta < lower) ||
      any(start_theta > upper)) {
    stop(
      "All starting values must lie within the specified bounds."
    )
  }

  has_bounds <-
    any(is.finite(lower)) ||
    any(is.finite(upper))

  if (has_bounds &&
      method != "L-BFGS-B") {
    stop(
      "Finite bounds require method = 'L-BFGS-B'."
    )
  }

  # ----------------------------------------------------------
  # Optimization controls
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
  } else {
    default_control <- list(
      maxit = 2000,
      reltol = 1e-10
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
  } else {
    control$factr <- NULL
    control$pgtol <- NULL
  }

  # ----------------------------------------------------------
  # Standardize derivative matrices
  # ----------------------------------------------------------

  standardize_jacobian <- function(J,
                                   name) {
    J <- as.matrix(J)

    if (!all(dim(J) == c(n, p))) {
      stop(
        "'",
        name,
        "' must have dimension ",
        n,
        " x ",
        p,
        "."
      )
    }

    if (anyNA(J) ||
        any(!is.finite(J))) {
      stop(
        "'",
        name,
        "' must contain only finite, non-missing values."
      )
    }

    colnames(J) <- parameter_names
    J
  }

  # ----------------------------------------------------------
  # Evaluate the user-supplied model
  # ----------------------------------------------------------

  evaluate_model <- function(theta,
                             strict = TRUE) {
    theta <- setNames(
      as.numeric(theta),
      parameter_names
    )

    answer <- tryCatch(
      model(
        theta,
        data
      ),
      error = function(e) e
    )

    if (inherits(answer, "error")) {
      if (strict) {
        stop(
          "The model function failed: ",
          conditionMessage(answer)
        )
      }

      return(NULL)
    }

    if (!is.list(answer)) {
      if (strict) {
        stop(
          "The model function must return a list."
        )
      }

      return(NULL)
    }

    if (is.null(answer$mu)) {
      if (strict) {
        stop(
          "The model function did not return 'mu'."
        )
      }

      return(NULL)
    }

    has_sigma2 <- !is.null(answer$sigma2)
    has_sigma <- !is.null(answer$sigma)

    if (!has_sigma2 &&
        !has_sigma) {
      if (strict) {
        stop(
          "The model function must return either 'sigma2' or 'sigma'."
        )
      }

      return(NULL)
    }

    mu <- as.numeric(
      answer$mu
    )

    sigma2 <- if (has_sigma2) {
      as.numeric(
        answer$sigma2
      )
    } else {
      as.numeric(
        answer$sigma
      )^2
    }

    if (length(mu) == 1L) {
      mu <- rep(
        mu,
        n
      )
    }

    if (length(sigma2) == 1L) {
      sigma2 <- rep(
        sigma2,
        n
      )
    }

    valid <-
      length(mu) == n &&
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
          ", with all scale-squared values strictly positive."
        )
      }

      return(NULL)
    }

    D <- NULL
    V <- NULL

    if (!is.null(answer$D)) {
      D <- standardize_jacobian(
        answer$D,
        "D"
      )
    }

    if (!is.null(answer$V)) {
      V <- standardize_jacobian(
        answer$V,
        "V"
      )
    }

    list(
      mu = mu,
      sigma2 = sigma2,
      D = D,
      V = V,
      raw = answer
    )
  }

  safe_evaluate_model <- function(theta) {
    evaluate_model(
      theta,
      strict = FALSE
    )
  }

  initial_model <- evaluate_model(
    start_theta,
    strict = TRUE
  )

  # ----------------------------------------------------------
  # Negative log-likelihood
  # ----------------------------------------------------------

  nll_from_model <- function(model_value) {
    standardized_y <-
      (
        y -
          model_value$mu
      ) /
      sqrt(
        model_value$sigma2
      )

    log_density <-
      stats::dt(
        standardized_y,
        df = nu,
        log = TRUE
      ) -
      0.5 *
      log(
        model_value$sigma2
      )

    value <- -sum(
      log_density
    )

    if (is.finite(value)) {
      value
    } else {
      Inf
    }
  }

  initial_nll <- nll_from_model(
    initial_model
  )

  if (!is.finite(initial_nll)) {
    stop(
      "The log-likelihood is not finite at the initial values."
    )
  }

  penalty_base <- max(
    1e8,
    min(
      1e100,
      1e6 *
        abs(initial_nll)
    )
  )

  negative_loglikelihood <- function(theta) {
    model_value <- safe_evaluate_model(
      theta
    )

    if (is.null(model_value)) {
      return(
        penalty_base +
          sum(
            (
              theta -
                start_theta
            )^2
          )
      )
    }

    value <- nll_from_model(
      model_value
    )

    if (!is.finite(value)) {
      return(
        penalty_base +
          sum(
            (
              theta -
                start_theta
            )^2
          )
      )
    }

    value
  }

  # ----------------------------------------------------------
  # Analytic or numerical derivatives of mu and sigma2
  # ----------------------------------------------------------

  calculate_jacobians <- function(theta,
                                  model_value = NULL) {
    theta <- setNames(
      as.numeric(theta),
      parameter_names
    )

    if (is.null(model_value)) {
      model_value <- evaluate_model(
        theta,
        strict = TRUE
      )
    }

    D <- model_value$D
    V <- model_value$V

    calculate_numerical_jacobian <- function(component,
                                             label) {
      J <- matrix(
        NA_real_,
        nrow = n,
        ncol = p,
        dimnames = list(
          NULL,
          parameter_names
        )
      )

      for (j in seq_len(p)) {
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

        model_plus <- safe_evaluate_model(
          theta_plus
        )

        model_minus <- safe_evaluate_model(
          theta_minus
        )

        value_current <- model_value[[
          component
        ]]

        if (!is.null(model_plus) &&
            !is.null(model_minus)) {
          J[, j] <-
            (
              model_plus[[component]] -
                model_minus[[component]]
            ) /
            (
              2 *
                step
            )
        } else if (!is.null(model_plus)) {
          J[, j] <-
            (
              model_plus[[component]] -
                value_current
            ) /
            step
        } else if (!is.null(model_minus)) {
          J[, j] <-
            (
              value_current -
                model_minus[[component]]
            ) /
            step
        } else {
          stop(
            "It was not possible to calculate the numerical ",
            label,
            " derivative with respect to '",
            parameter_names[j],
            "'."
          )
        }
      }

      standardize_jacobian(
        J,
        label
      )
    }

    if (is.null(D)) {
      D <- calculate_numerical_jacobian(
        component = "mu",
        label = "D"
      )
    }

    if (is.null(V)) {
      V <- calculate_numerical_jacobian(
        component = "sigma2",
        label = "V"
      )
    }

    list(
      D = D,
      V = V
    )
  }

  # ----------------------------------------------------------
  # Score vector
  # ----------------------------------------------------------

  calculate_score <- function(theta,
                              model_value = NULL) {
    if (is.null(model_value)) {
      model_value <- evaluate_model(
        theta,
        strict = TRUE
      )
    }

    jacobians <- calculate_jacobians(
      theta = theta,
      model_value = model_value
    )

    D <- jacobians$D
    V <- jacobians$V

    residual <- y -
      model_value$mu

    sigma2 <- model_value$sigma2

    denominator <-
      nu *
      sigma2 +
      residual^2

    score_mu <-
      (
        nu +
          1
      ) *
      residual /
      denominator

    score_sigma2 <-
      nu *
      (
        residual^2 -
          sigma2
      ) /
      (
        2 *
        sigma2 *
        denominator
      )

    score <-
      drop(
        crossprod(
          D,
          score_mu
        )
      ) +
      drop(
        crossprod(
          V,
          score_sigma2
        )
      )

    setNames(
      score,
      parameter_names
    )
  }

  # ----------------------------------------------------------
  # Numerical gradient used at invalid trial points
  # ----------------------------------------------------------

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
        (
          2 *
            step
        )
    }

    setNames(
      gradient,
      parameter_names
    )
  }

  negative_score <- function(theta) {
    model_value <- safe_evaluate_model(
      theta
    )

    if (is.null(model_value)) {
      return(
        numerical_gradient(
          negative_loglikelihood,
          theta
        )
      )
    }

    answer <- tryCatch(
      -calculate_score(
        theta,
        model_value
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

  # ----------------------------------------------------------
  # Optimization
  # ----------------------------------------------------------

  optim_arguments <- list(
    par = start_theta,
    fn = negative_loglikelihood,
    method = method,
    control = control,
    hessian = hessian
  )

  gradient_methods <- c(
    "BFGS",
    "CG",
    "L-BFGS-B"
  )

  if (isTRUE(use_gradient) &&
      method %in% gradient_methods) {
    optim_arguments$gr <-
      negative_score
  }

  if (method == "L-BFGS-B") {
    optim_arguments$lower <-
      lower

    optim_arguments$upper <-
      upper
  }

  optimization <- do.call(
    stats::optim,
    optim_arguments
  )

  theta_hat <- setNames(
    as.numeric(
      optimization$par
    ),
    parameter_names
  )

  fitted_model <- evaluate_model(
    theta_hat,
    strict = TRUE
  )

  fitted_jacobians <- calculate_jacobians(
    theta = theta_hat,
    model_value = fitted_model
  )

  D_hat <- fitted_jacobians$D
  V_hat <- fitted_jacobians$V

  score_hat <- calculate_score(
    theta = theta_hat,
    model_value = fitted_model
  )

  # ----------------------------------------------------------
  # Expected Fisher information
  # ----------------------------------------------------------

  sigma2_hat <- fitted_model$sigma2

  mean_constant <-
    (
      nu +
        1
    ) /
    (
      nu +
        3
    )

  variance_constant <-
    nu /
    (
      2 *
      (
        nu +
          3
      )
    )

  mean_weights <-
    mean_constant /
    sigma2_hat

  variance_weights <-
    variance_constant /
    sigma2_hat^2

  fisher_information <-
    crossprod(
      D_hat,
      sweep(
        D_hat,
        1L,
        mean_weights,
        FUN = "*"
      )
    ) +
    crossprod(
      V_hat,
      sweep(
        V_hat,
        1L,
        variance_weights,
        FUN = "*"
      )
    )

  dimnames(
    fisher_information
  ) <- list(
    parameter_names,
    parameter_names
  )

  # ----------------------------------------------------------
  # Matrix inversion
  # ----------------------------------------------------------

  safe_inverse <- function(M) {
    inverse <- tryCatch(
      {
        chol2inv(
          chol(M)
        )
      },
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

  vcov_fisher <- safe_inverse(
    fisher_information
  )

  observed_information <- NULL
  vcov_observed <- NULL

  if (isTRUE(hessian) &&
      !is.null(optimization$hessian)) {
    observed_information <-
      as.matrix(
        optimization$hessian
      )

    dimnames(
      observed_information
    ) <- list(
      parameter_names,
      parameter_names
    )

    vcov_observed <- safe_inverse(
      observed_information
    )
  }

  # ----------------------------------------------------------
  # Standard errors and confidence intervals
  # ----------------------------------------------------------

  fisher_variances <- diag(
    vcov_fisher
  )

  standard_error <- setNames(
    ifelse(
      is.finite(fisher_variances) &
        fisher_variances >= 0,
      sqrt(fisher_variances),
      NA_real_
    ),
    parameter_names
  )

  z_value <-
    theta_hat /
    standard_error

  p_value <-
    2 *
    stats::pnorm(
      abs(z_value),
      lower.tail = FALSE
    )

  critical_value <- stats::qnorm(
    (
      1 +
        conf_level
    ) /
    2
  )

  conf_int <- cbind(
    lower =
      theta_hat -
      critical_value *
      standard_error,
    upper =
      theta_hat +
      critical_value *
      standard_error
  )

  coefficient_table <- cbind(
    Estimate = theta_hat,
    `Std. Error` =
      standard_error,
    `z value` =
      z_value,
    `Pr(>|z|)` =
      p_value
  )

  # ----------------------------------------------------------
  # Residuals and additional quantities
  # ----------------------------------------------------------

  fitted_values <- fitted_model$mu
  residual_values <- y -
    fitted_values

  u_values <-
    residual_values^2 /
    sigma2_hat

  t_weights <-
    (
      nu +
        1
    ) /
    (
      nu +
        u_values
    )

  scale_residuals <-
    residual_values /
    sqrt(
      sigma2_hat
    )

  if (nu > 2) {
    pearson_residuals <-
      residual_values /
      sqrt(
        (
          nu /
          (
            nu -
              2
          )
        ) *
        sigma2_hat
      )
  } else {
    pearson_residuals <- rep(
      NA_real_,
      n
    )
  }

  maximized_loglikelihood <-
    -optimization$value

  AIC_value <-
    -2 *
    maximized_loglikelihood +
    2 *
    p

  BIC_value <-
    -2 *
    maximized_loglikelihood +
    log(n) *
    p

  result <- list(
    call = call,
    coefficients = theta_hat,
    standard.error = standard_error,
    conf.int = conf_int,
    coefficient.table = coefficient_table,
    vcov = vcov_fisher,
    vcov.fisher = vcov_fisher,
    vcov.observed = vcov_observed,
    fisher.information = fisher_information,
    observed.information = observed_information,
    score = score_hat,
    max.abs.score = max(
      abs(score_hat)
    ),
    logLik = maximized_loglikelihood,
    objective = optimization$value,
    AIC = AIC_value,
    BIC = BIC_value,
    fitted.values = fitted_values,
    sigma2 = sigma2_hat,
    scale = sqrt(
      sigma2_hat
    ),
    residuals = residual_values,
    scale.residuals = scale_residuals,
    pearson.residuals = pearson_residuals,
    u = u_values,
    weights = t_weights,
    mean.jacobian = D_hat,
    scale2.jacobian = V_hat,
    nu = nu,
    nobs = n,
    n = n,
    p = p,
    convergence = optimization$convergence,
    message = optimization$message,
    counts = optimization$counts,
    method = method,
    model = model,
    data = data,
    y = y,
    optimization.control = control,
    optim = optimization
  )

  class(result) <- "nlm_t_fit"

  if (!identical(
    as.integer(
      optimization$convergence
    ),
    0L
  )) {
    warning(
      "optim() returned convergence code ",
      optimization$convergence,
      "."
    )
  }

  result
}


# ============================================================
# S3 methods
# ============================================================

print.nlm_t_fit <- function(
    x,
    digits = max(
      3L,
      getOption("digits") -
        3L
    ),
    ...
) {
  cat(
    "\nStudent-t regression model\n"
  )

  cat(
    "Fixed degrees of freedom:",
    x$nu,
    "\n"
  )

  cat(
    "Number of observations:",
    x$nobs,
    "\n\n"
  )

  printCoefmat(
    x$coefficient.table,
    digits = digits,
    P.values = TRUE,
    has.Pvalue = TRUE
  )

  cat(
    "\nLog-likelihood:",
    format(
      x$logLik,
      digits = digits
    )
  )

  cat(
    "\nAIC:",
    format(
      x$AIC,
      digits = digits
    )
  )

  cat(
    "\nBIC:",
    format(
      x$BIC,
      digits = digits
    )
  )

  cat(
    "\nMaximum absolute score:",
    format(
      x$max.abs.score,
      digits = digits
    )
  )

  cat(
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


summary.nlm_t_fit <- function(object,
                              ...) {
  result <- list(
    call = object$call,
    coefficients =
      object$coefficient.table,
    logLik = object$logLik,
    AIC = object$AIC,
    BIC = object$BIC,
    score = object$score,
    max.abs.score =
      object$max.abs.score,
    convergence =
      object$convergence,
    message = object$message,
    nu = object$nu,
    nobs = object$nobs,
    method = object$method
  )

  class(result) <-
    "summary.nlm_t_fit"

  result
}


print.summary.nlm_t_fit <- function(
    x,
    digits = max(
      3L,
      getOption("digits") -
        3L
    ),
    ...
) {
  cat(
    "Student-t regression model fitted by maximum likelihood\n"
  )

  cat(
    "Fixed degrees of freedom:",
    x$nu,
    "| n =",
    x$nobs,
    "| optimization method =",
    x$method,
    "\n\n"
  )

  printCoefmat(
    x$coefficients,
    digits = digits,
    P.values = TRUE,
    has.Pvalue = TRUE
  )

  cat(
    "\nLog-likelihood:",
    format(
      x$logLik,
      digits = digits
    ),
    "\nAIC:",
    format(
      x$AIC,
      digits = digits
    ),
    "\nBIC:",
    format(
      x$BIC,
      digits = digits
    ),
    "\nMaximum absolute score:",
    format(
      x$max.abs.score,
      digits = digits
    ),
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


coef.nlm_t_fit <- function(object,
                           ...) {
  object$coefficients
}


vcov.nlm_t_fit <- function(
    object,
    type = c(
      "fisher",
      "observed"
    ),
    ...
) {
  type <- match.arg(type)

  if (type == "fisher") {
    return(
      object$vcov.fisher
    )
  }

  if (is.null(
    object$vcov.observed
  )) {
    stop(
      "The observed covariance matrix is not available."
    )
  }

  object$vcov.observed
}


confint.nlm_t_fit <- function(
    object,
    parm = seq_along(
      object$coefficients
    ),
    level = 0.95,
    ...
) {
  if (length(level) != 1L ||
      !is.numeric(level) ||
      !is.finite(level) ||
      level <= 0 ||
      level >= 1) {
    stop(
      "'level' must be between zero and one."
    )
  }

  parameter_indices <- if (
    is.character(parm)
  ) {
    match(
      parm,
      names(
        object$coefficients
      )
    )
  } else {
    parm
  }

  if (anyNA(parameter_indices)) {
    stop(
      "At least one requested parameter name was not found."
    )
  }

  estimates <-
    object$coefficients[
      parameter_indices
    ]

  standard_errors <-
    object$standard.error[
      parameter_indices
    ]

  critical_value <- stats::qnorm(
    (
      1 +
        level
    ) /
    2
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

  rownames(interval) <- names(
    estimates
  )

  interval
}


fitted.nlm_t_fit <- function(object,
                             ...) {
  object$fitted.values
}


residuals.nlm_t_fit <- function(
    object,
    type = c(
      "response",
      "scale",
      "pearson"
    ),
    ...
) {
  type <- match.arg(type)

  switch(
    type,
    response =
      object$residuals,
    scale =
      object$scale.residuals,
    pearson =
      object$pearson.residuals
  )
}


predict.nlm_t_fit <- function(
    object,
    newdata = NULL,
    type = c(
      "mean",
      "sigma",
      "variance"
    ),
    ...
) {
  type <- match.arg(type)

  prediction_data <- if (
    is.null(newdata)
  ) {
    object$data
  } else {
    newdata
  }

  model_output <- object$model(
    object$coefficients,
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
      as.numeric(
        model_output$mu
      )
    )
  }

  if (!is.null(
    model_output$sigma2
  )) {
    sigma2 <- as.numeric(
      model_output$sigma2
    )
  } else if (!is.null(
    model_output$sigma
  )) {
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


logLik.nlm_t_fit <- function(object,
                             ...) {
  structure(
    object$logLik,
    df = object$p,
    nobs = object$nobs,
    class = "logLik"
  )
}


nobs.nlm_t_fit <- function(object,
                           ...) {
  object$nobs
}

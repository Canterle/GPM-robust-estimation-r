# ============================================================
# Maximum likelihood estimation for linear mixed models with
# a multivariate Student-t marginal distribution
#
# For independent groups i = 1, ..., n_groups, the model is
#
#   Y_i ~ t_{d_i}(X_i beta, Sigma_i, nu),
#
# where
#
#   Sigma_i = Z_i Delta(gamma) Z_i^T + sigma2 I_{d_i}.
#
# Here, Sigma_i is the scale matrix of the multivariate Student-t
# distribution, not its covariance matrix. When nu > 2,
#
#   Var(Y_i) = nu / (nu - 2) * Sigma_i.
#
# A hierarchical interpretation that leads exactly to this marginal
# model is obtained with one common scale-mixture variable per group:
#
#   lambda_i ~ Gamma(nu / 2, rate = nu / 2),
#   b_i | lambda_i ~ Normal(0, Delta(gamma) / lambda_i),
#   epsilon_i | lambda_i ~ Normal(0, sigma2 I_{d_i} / lambda_i),
#
# conditionally independently, with
#
#   Y_i = X_i beta + Z_i b_i + epsilon_i.
#
# This is not the same model as taking b_i Gaussian and epsilon_i
# Student-t independently. That convolution is not multivariate
# Student-t and does not have the likelihood used below.
#
# The random-effects scale structure is supplied through
#
#   Delta(gamma),
#
# which must return an m by m symmetric positive-definite matrix.
# The optional function Delta_jacobian(gamma) must return
#
#   d vec(Delta(gamma)) / d gamma^T,
#
# either as an m^2 by p_gamma matrix or as an m by m by p_gamma
# array. R column-major vectorization is used. If the derivative is
# omitted, it is computed numerically with numDeriv::jacobian().
#
# The parameter vector is ordered as
#
#   theta = (beta, gamma, sigma2),
#
# and nu is fixed rather than estimated.
# ============================================================


# ============================================================
# Internal utility functions
# ============================================================

.lmm_t_as_group_list <- function(value,
                                 name) {
  if (is.list(value)) {
    return(value)
  }

  if (is.matrix(value) ||
      is.data.frame(value) ||
      is.numeric(value)) {
    return(list(value))
  }

  stop(
    "'", name, "' must be a list or an object representing one group."
  )
}


.lmm_t_prepare_data <- function(Y,
                                X,
                                Z) {
  Y <- .lmm_t_as_group_list(
    Y,
    "Y"
  )

  X <- .lmm_t_as_group_list(
    X,
    "X"
  )

  Z <- .lmm_t_as_group_list(
    Z,
    "Z"
  )

  n_groups <- length(Y)

  if (n_groups == 0L) {
    stop("'Y' must contain at least one group.")
  }

  if (length(X) != n_groups ||
      length(Z) != n_groups) {
    stop(
      "'Y', 'X' and 'Z' must contain the same number of groups."
    )
  }

  Y <- lapply(
    Y,
    as.numeric
  )

  X <- lapply(
    X,
    as.matrix
  )

  Z <- lapply(
    Z,
    as.matrix
  )

  group_sizes <- integer(n_groups)

  for (i in seq_len(n_groups)) {
    group_sizes[i] <- length(Y[[i]])

    if (group_sizes[i] == 0L) {
      stop(
        "Every group in 'Y' must contain at least one observation."
      )
    }

    if (nrow(X[[i]]) != group_sizes[i]) {
      stop(
        "The number of rows of X[[", i,
        "]] must equal length(Y[[", i, "]])."
      )
    }

    if (nrow(Z[[i]]) != group_sizes[i]) {
      stop(
        "The number of rows of Z[[", i,
        "]] must equal length(Y[[", i, "]])."
      )
    }

    if (anyNA(Y[[i]]) ||
        any(!is.finite(Y[[i]])) ||
        anyNA(X[[i]]) ||
        any(!is.finite(X[[i]])) ||
        anyNA(Z[[i]]) ||
        any(!is.finite(Z[[i]]))) {
      stop(
        "'Y', 'X' and 'Z' must contain only finite, non-missing values."
      )
    }
  }

  p_beta <- ncol(X[[1]])
  m <- ncol(Z[[1]])

  if (p_beta == 0L) {
    stop(
      "The fixed-effects design matrices must have at least one column."
    )
  }

  if (m == 0L) {
    stop(
      "The random-effects design matrices must have at least one column."
    )
  }

  if (any(
    vapply(
      X,
      ncol,
      integer(1L)
    ) != p_beta
  )) {
    stop(
      "All fixed-effects design matrices must have the same number of columns."
    )
  }

  if (any(
    vapply(
      Z,
      ncol,
      integer(1L)
    ) != m
  )) {
    stop(
      "All random-effects design matrices must have the same number of columns."
    )
  }

  beta_names <- colnames(X[[1]])

  if (is.null(beta_names) ||
      any(is.na(beta_names)) ||
      any(beta_names == "")) {
    beta_names <- paste0(
      "beta",
      seq_len(p_beta)
    )
  }

  beta_names <- make.unique(beta_names)

  random_effect_names <- colnames(Z[[1]])

  if (is.null(random_effect_names) ||
      any(is.na(random_effect_names)) ||
      any(random_effect_names == "")) {
    random_effect_names <- paste0(
      "b",
      seq_len(m)
    )
  }

  random_effect_names <- make.unique(
    random_effect_names
  )

  list(
    Y = Y,
    X = X,
    Z = Z,
    n.groups = n_groups,
    nobs = sum(group_sizes),
    group.sizes = group_sizes,
    p.beta = p_beta,
    m = m,
    beta.names = beta_names,
    random.effect.names = random_effect_names
  )
}


# Validate the matrix identifying where each gamma parameter appears
# in Delta(gamma). This is needed only for automatic starting values.
# A gamma parameter appearing in a diagonal position receives starting
# value 1; a gamma parameter appearing only off the diagonal receives 0.
.lmm_t_prepare_gamma_structure <- function(gamma_structure,
                                           m) {
  if (is.null(gamma_structure)) {
    stop(
      "'gamma_structure' must be supplied when 'start = NULL'."
    )
  }

  gamma_structure <- as.matrix(
    gamma_structure
  )

  if (nrow(gamma_structure) != m ||
      ncol(gamma_structure) != m) {
    stop(
      "'gamma_structure' must be a square ",
      m, " by ", m, " matrix."
    )
  }

  gamma_structure <- matrix(
    as.character(gamma_structure),
    nrow = m,
    ncol = m,
    dimnames = dimnames(gamma_structure)
  )

  empty_entries <-
    is.na(gamma_structure) |
    gamma_structure == "" |
    gamma_structure == "0" |
    toupper(gamma_structure) == "NA"

  gamma_structure[empty_entries] <- NA_character_

  nonempty_entries <- gamma_structure[
    !is.na(gamma_structure)
  ]

  if (length(nonempty_entries) == 0L) {
    stop(
      "'gamma_structure' must identify at least one gamma parameter."
    )
  }

  if (any(
    !grepl(
      "^[.A-Za-z][.A-Za-z0-9_]*$",
      nonempty_entries
    )
  )) {
    stop(
      "Every non-missing entry of 'gamma_structure' must be ",
      "a valid gamma parameter name."
    )
  }

  for (row_index in seq_len(m)) {
    for (column_index in seq_len(m)) {
      left <- gamma_structure[
        row_index,
        column_index
      ]

      right <- gamma_structure[
        column_index,
        row_index
      ]

      same_missing <-
        is.na(left) &&
        is.na(right)

      same_name <-
        !is.na(left) &&
        !is.na(right) &&
        identical(left, right)

      if (!same_missing &&
          !same_name) {
        stop(
          "'gamma_structure' must be symmetric. The entries at [",
          row_index, ", ", column_index, "] and [",
          column_index, ", ", row_index, "] do not match."
        )
      }
    }
  }

  gamma_names <- unique(
    as.vector(gamma_structure)
  )

  gamma_names <- gamma_names[
    !is.na(gamma_names)
  ]

  diagonal_names <- unique(
    diag(gamma_structure)
  )

  diagonal_names <- diagonal_names[
    !is.na(diagonal_names)
  ]

  gamma_start <- setNames(
    as.numeric(
      gamma_names %in%
      diagonal_names
    ),
    gamma_names
  )

  list(
    structure = gamma_structure,
    names = gamma_names,
    diagonal.names = diagonal_names,
    start = gamma_start
  )
}


# Compute robust default starting values when start = NULL.
# This rule is suitable when gamma directly parameterizes entries of
# Delta(gamma). For log-scale or Cholesky parameterizations, supply
# start explicitly.
.lmm_t_default_start <- function(data,
                                 gamma_structure) {
  if (!requireNamespace(
    "MASS",
    quietly = TRUE
  )) {
    stop(
      "Package 'MASS' is required when 'start = NULL'."
    )
  }

  gamma_information <-
    .lmm_t_prepare_gamma_structure(
      gamma_structure = gamma_structure,
      m = data$m
    )

  Y_stacked <- unlist(
    data$Y,
    use.names = FALSE
  )

  X_stacked <- do.call(
    rbind,
    data$X
  )

  robust_fit <- tryCatch(
    MASS::rlm(
      x = X_stacked,
      y = Y_stacked,
      maxit = 100
    ),
    error = function(e) {
      stop(
        "The robust starting-value fit using MASS::rlm() failed: ",
        conditionMessage(e)
      )
    }
  )

  beta_start <- as.numeric(
    stats::coef(robust_fit)
  )

  if (length(beta_start) != data$p.beta ||
      anyNA(beta_start) ||
      any(!is.finite(beta_start))) {
    stop(
      "MASS::rlm() did not produce a finite fixed-effect ",
      "starting vector of the expected length."
    )
  }

  names(beta_start) <- data$beta.names

  sigma2_start <- robust_fit$s^2

  if (length(sigma2_start) != 1L ||
      !is.finite(sigma2_start) ||
      sigma2_start < 0.05) {
    sigma2_start <- 0.1
  }

  start <- c(
    beta_start,
    gamma_information$start,
    sigma2 = sigma2_start
  )

  list(
    start = start,
    gamma.structure = gamma_information$structure,
    gamma.diagonal.names = gamma_information$diagonal.names,
    robust.fit = robust_fit
  )
}


.lmm_t_prepare_start <- function(start,
                                 p_beta,
                                 beta_names = NULL) {
  if (!is.numeric(start) ||
      length(start) < p_beta + 2L) {
    stop(
      "'start' must be a numeric vector containing the fixed effects, ",
      "at least one parameter of Delta(gamma), and sigma2."
    )
  }

  original_names <- names(start)

  start <- as.numeric(start)
  p <- length(start)
  p_gamma <- p - p_beta - 1L

  if (is.null(beta_names) ||
      length(beta_names) != p_beta) {
    beta_names <- paste0(
      "beta",
      seq_len(p_beta)
    )
  }

  if (is.null(original_names)) {
    parameter_names <- c(
      beta_names,
      paste0(
        "gamma",
        seq_len(p_gamma)
      ),
      "sigma2"
    )
  } else {
    parameter_names <- original_names

    empty_names <- is.na(parameter_names) |
      parameter_names == ""

    if (any(empty_names)) {
      default_names <- c(
        beta_names,
        paste0(
          "gamma",
          seq_len(p_gamma)
        ),
        "sigma2"
      )

      parameter_names[empty_names] <-
        default_names[empty_names]
    }

    parameter_names <- make.unique(
      parameter_names
    )
  }

  names(start) <- parameter_names

  if (anyNA(start) ||
      any(!is.finite(start))) {
    stop(
      "'start' must contain only finite, non-missing values."
    )
  }

  beta_indices <- seq_len(p_beta)

  gamma_indices <- p_beta +
    seq_len(p_gamma)

  sigma2_index <- p

  beta_parameter_names <- parameter_names[
    beta_indices
  ]

  gamma_names <- parameter_names[
    gamma_indices
  ]

  sigma2_name <- parameter_names[
    sigma2_index
  ]

  if (start[sigma2_index] <= 0) {
    stop(
      "The starting value of sigma2 must be strictly positive."
    )
  }

  list(
    start = start,
    p = p,
    p.gamma = p_gamma,
    beta.indices = beta_indices,
    gamma.indices = gamma_indices,
    covariance.indices = c(
      gamma_indices,
      sigma2_index
    ),
    sigma2.index = sigma2_index,
    parameter.names = parameter_names,
    beta.names = beta_parameter_names,
    gamma.names = gamma_names,
    sigma2.name = sigma2_name
  )
}


.lmm_t_call_Delta <- function(Delta,
                              gamma,
                              m,
                              strict = FALSE) {
  output <- tryCatch(
    Delta(gamma),
    error = function(e) e
  )

  if (inherits(output, "error")) {
    if (strict) {
      stop(
        "'Delta' failed: ",
        conditionMessage(output)
      )
    }

    return(NULL)
  }

  output <- as.matrix(output)

  valid_dimensions <-
    nrow(output) == m &&
    ncol(output) == m

  valid_values <-
    !anyNA(output) &&
    all(is.finite(output))

  if (!valid_dimensions ||
      !valid_values) {
    if (strict) {
      stop(
        "'Delta(gamma)' must return a finite ",
        m, " by ", m, " matrix."
      )
    }

    return(NULL)
  }

  symmetry_tolerance <-
    100 *
    sqrt(.Machine$double.eps) *
    max(
      1,
      max(abs(output))
    )

  if (max(abs(output - t(output))) >
      symmetry_tolerance) {
    if (strict) {
      stop(
        "'Delta(gamma)' must be symmetric."
      )
    }

    return(NULL)
  }

  output <- (
    output +
    t(output)
  ) / 2

  cholesky <- tryCatch(
    chol(output),
    error = function(e) NULL
  )

  if (is.null(cholesky)) {
    if (strict) {
      stop(
        "'Delta(gamma)' must be positive definite."
      )
    }

    return(NULL)
  }

  list(
    value = output,
    chol = cholesky
  )
}


.lmm_t_compute_Delta_jacobian <- function(gamma,
                                          Delta,
                                          Delta_jacobian,
                                          m,
                                          method,
                                          method_args) {
  p_gamma <- length(gamma)

  if (!is.null(Delta_jacobian)) {
    if (!is.function(Delta_jacobian)) {
      stop(
        "'Delta_jacobian' must be NULL or a function."
      )
    }

    derivative <- tryCatch(
      Delta_jacobian(gamma),
      error = function(e) {
        stop(
          "'Delta_jacobian' failed: ",
          conditionMessage(e)
        )
      }
    )

    if (length(dim(derivative)) == 3L) {
      expected_dimensions <- c(
        m,
        m,
        p_gamma
      )

      if (!all(
        dim(derivative) ==
        expected_dimensions
      )) {
        stop(
          "An array returned by 'Delta_jacobian' must have dimensions ",
          "m by m by length(gamma)."
        )
      }

      derivative <- do.call(
        cbind,
        lapply(
          seq_len(p_gamma),
          function(index) {
            c(derivative[, , index])
          }
        )
      )
    } else {
      derivative <- as.matrix(derivative)
    }
  } else {
    if (!requireNamespace(
      "numDeriv",
      quietly = TRUE
    )) {
      stop(
        "Package 'numDeriv' is required when 'Delta_jacobian' is omitted. ",
        "Install it with install.packages('numDeriv')."
      )
    }

    gamma_names <- names(gamma)

    delta_vector_function <- function(value) {
      value <- setNames(
        as.numeric(value),
        gamma_names
      )

      delta_result <- .lmm_t_call_Delta(
        Delta = Delta,
        gamma = value,
        m = m,
        strict = TRUE
      )

      c(delta_result$value)
    }

    derivative <- tryCatch(
      do.call(
        numDeriv::jacobian,
        c(
          list(
            func = delta_vector_function,
            x = gamma,
            method = method
          ),
          list(
            method.args = method_args
          )
        )
      ),
      error = function(e) {
        stop(
          "The numerical Jacobian of vec(Delta(gamma)) could not be computed: ",
          conditionMessage(e)
        )
      }
    )

    derivative <- as.matrix(derivative)
  }

  expected_dimensions <- c(
    m^2,
    p_gamma
  )

  if (!all(
    dim(derivative) ==
    expected_dimensions
  )) {
    stop(
      "The Jacobian of vec(Delta(gamma)) must have dimensions ",
      m^2, " by ", p_gamma, "."
    )
  }

  if (anyNA(derivative) ||
      any(!is.finite(derivative))) {
    stop(
      "The Jacobian of vec(Delta(gamma)) contains non-finite values."
    )
  }

  # Because Delta(gamma) is symmetric, every derivative matrix must
  # also be symmetric. Symmetrizing here removes only numerical noise
  # and ensures compatibility with the Fisher-information formula.
  for (j in seq_len(p_gamma)) {
    derivative_matrix <- matrix(
      derivative[, j],
      nrow = m,
      ncol = m
    )

    asymmetry_tolerance <-
      100 *
      sqrt(.Machine$double.eps) *
      max(
        1,
        max(abs(derivative_matrix))
      )

    if (max(
      abs(
        derivative_matrix -
        t(derivative_matrix)
      )
    ) > asymmetry_tolerance) {
      stop(
        "Each column of the Jacobian of vec(Delta(gamma)) must ",
        "represent a symmetric derivative matrix."
      )
    }

    derivative[, j] <- c(
      (
        derivative_matrix +
        t(derivative_matrix)
      ) / 2
    )
  }

  colnames(derivative) <- names(gamma)

  rownames(derivative) <- paste0(
    "Delta",
    seq_len(m^2)
  )

  derivative
}


.lmm_t_derivative_matrices <- function(theta,
                                       data,
                                       structure,
                                       Delta,
                                       Delta_jacobian,
                                       jacobian_method,
                                       jacobian_method_args) {
  theta <- setNames(
    as.numeric(theta),
    structure$parameter.names
  )

  gamma <- theta[
    structure$gamma.indices
  ]

  delta_derivative <-
    .lmm_t_compute_Delta_jacobian(
      gamma = gamma,
      Delta = Delta,
      Delta_jacobian = Delta_jacobian,
      m = data$m,
      method = jacobian_method,
      method_args = jacobian_method_args
    )

  mean_jacobians <- vector(
    "list",
    data$n.groups
  )

  covariance_jacobians <- vector(
    "list",
    data$n.groups
  )

  F_matrices <- vector(
    "list",
    data$n.groups
  )

  for (i in seq_len(data$n.groups)) {
    d_i <- data$group.sizes[i]

    D_i <- data$X[[i]]

    V_gamma_i <-
      kronecker(
        data$Z[[i]],
        data$Z[[i]]
      ) %*%
      delta_derivative

    V_i <- cbind(
      V_gamma_i,
      c(diag(d_i))
    )

    F_i <- rbind(
      cbind(
        D_i,
        matrix(
          0,
          nrow = d_i,
          ncol = structure$p.gamma + 1L
        )
      ),
      cbind(
        matrix(
          0,
          nrow = d_i^2,
          ncol = data$p.beta
        ),
        V_i
      )
    )

    colnames(D_i) <- structure$beta.names

    colnames(V_i) <- c(
      structure$gamma.names,
      structure$sigma2.name
    )

    colnames(F_i) <- structure$parameter.names

    mean_jacobians[[i]] <- D_i
    covariance_jacobians[[i]] <- V_i
    F_matrices[[i]] <- F_i
  }

  list(
    Delta.jacobian = delta_derivative,
    mean.jacobians = mean_jacobians,
    covariance.jacobians = covariance_jacobians,
    F = F_matrices
  )
}


.lmm_t_evaluate_model <- function(theta,
                                  data,
                                  structure,
                                  Delta,
                                  nu,
                                  strict = FALSE) {
  theta <- setNames(
    as.numeric(theta),
    structure$parameter.names
  )

  beta <- theta[
    structure$beta.indices
  ]

  gamma <- theta[
    structure$gamma.indices
  ]

  sigma2 <- theta[
    structure$sigma2.index
  ]

  if (!is.finite(sigma2) ||
      sigma2 <= 0) {
    if (strict) {
      stop(
        "'sigma2' must be strictly positive."
      )
    }

    return(NULL)
  }

  delta_result <- .lmm_t_call_Delta(
    Delta = Delta,
    gamma = gamma,
    m = data$m,
    strict = strict
  )

  if (is.null(delta_result)) {
    return(NULL)
  }

  delta_matrix <- delta_result$value

  means <- vector(
    "list",
    data$n.groups
  )

  scale_matrices <- vector(
    "list",
    data$n.groups
  )

  scale_inverses <- vector(
    "list",
    data$n.groups
  )

  variance_matrices <- if (nu > 2) {
    vector(
      "list",
      data$n.groups
    )
  } else {
    NULL
  }

  residuals <- vector(
    "list",
    data$n.groups
  )

  standardized_residuals <- vector(
    "list",
    data$n.groups
  )

  log_determinants <- numeric(
    data$n.groups
  )

  quadratic_forms <- numeric(
    data$n.groups
  )

  t_weights <- numeric(
    data$n.groups
  )

  log_densities <- numeric(
    data$n.groups
  )

  for (i in seq_len(data$n.groups)) {
    d_i <- data$group.sizes[i]

    mean_i <- drop(
      data$X[[i]] %*%
      beta
    )

    Sigma_i <-
      data$Z[[i]] %*%
      delta_matrix %*%
      t(data$Z[[i]]) +
      sigma2 *
      diag(d_i)

    chol_Sigma <- tryCatch(
      chol(Sigma_i),
      error = function(e) NULL
    )

    if (is.null(chol_Sigma)) {
      if (strict) {
        stop(
          "The marginal scale matrix is not positive definite ",
          "for group ", i, "."
        )
      }

      return(NULL)
    }

    root_inverse <- backsolve(
      chol_Sigma,
      diag(d_i)
    )

    Sigma_inverse <- tcrossprod(
      root_inverse
    )

    residual_i <-
      data$Y[[i]] -
      mean_i

    standardized_residual_i <- drop(
      crossprod(
        root_inverse,
        residual_i
      )
    )

    quadratic_form <- sum(
      standardized_residual_i^2
    )

    log_determinant <-
      2 *
      sum(
        log(diag(chol_Sigma))
      )

    log_density <-
      lgamma((nu + d_i) / 2) -
      lgamma(nu / 2) -
      0.5 * d_i * log(nu * pi) -
      0.5 * log_determinant -
      0.5 * (nu + d_i) *
      log1p(quadratic_form / nu)

    if (!is.finite(log_density)) {
      if (strict) {
        stop(
          "The model produced a non-finite log-density for group ",
          i, "."
        )
      }

      return(NULL)
    }

    means[[i]] <- mean_i
    scale_matrices[[i]] <- Sigma_i
    scale_inverses[[i]] <- Sigma_inverse

    if (nu > 2) {
      variance_matrices[[i]] <-
        nu /
        (nu - 2) *
        Sigma_i
    }

    residuals[[i]] <- residual_i
    standardized_residuals[[i]] <-
      standardized_residual_i

    log_determinants[i] <-
      log_determinant

    quadratic_forms[i] <-
      quadratic_form

    t_weights[i] <-
      (nu + d_i) /
      (nu + quadratic_form)

    log_densities[i] <-
      log_density
  }

  list(
    beta = beta,
    gamma = gamma,
    sigma2 = unname(sigma2),
    Delta = delta_matrix,
    Delta.variance = if (nu > 2) {
      nu /
      (nu - 2) *
      delta_matrix
    } else {
      NULL
    },
    residual.variance = if (nu > 2) {
      nu /
      (nu - 2) *
      unname(sigma2)
    } else {
      NA_real_
    },
    means = means,
    scale.matrices = scale_matrices,
    scale.inverses = scale_inverses,
    variance.matrices = variance_matrices,
    residuals = residuals,
    standardized.residuals = standardized_residuals,
    log.determinants = log_determinants,
    quadratic.forms = quadratic_forms,
    weights = t_weights,
    log.densities = log_densities,
    densities = exp(log_densities)
  )
}


.lmm_t_compute_random_effects <- function(components,
                                          data) {
  random_effects <- vector(
    "list",
    data$n.groups
  )

  conditional_means <- vector(
    "list",
    data$n.groups
  )

  conditional_residuals <- vector(
    "list",
    data$n.groups
  )

  for (i in seq_len(data$n.groups)) {
    b_hat_i <-
      components$Delta %*%
      t(data$Z[[i]]) %*%
      components$scale.inverses[[i]] %*%
      components$residuals[[i]]

    b_hat_i <- drop(b_hat_i)

    names(b_hat_i) <-
      data$random.effect.names

    conditional_mean_i <-
      components$means[[i]] +
      drop(
        data$Z[[i]] %*%
        b_hat_i
      )

    random_effects[[i]] <- b_hat_i
    conditional_means[[i]] <-
      conditional_mean_i
    conditional_residuals[[i]] <-
      data$Y[[i]] -
      conditional_mean_i
  }

  list(
    random.effects = random_effects,
    conditional.means = conditional_means,
    conditional.residuals = conditional_residuals
  )
}


.lmm_t_safe_inverse <- function(M) {
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


# ============================================================
# Main fitting function
# ============================================================

fit_lmm_t <- function(Y,
                      X,
                      Z,
                      Delta,
                      start = NULL,
                      gamma_structure = NULL,
                      Delta_jacobian = NULL,
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
                      jacobian_method = c(
                        "Richardson",
                        "simple"
                      ),
                      jacobian_method_args = list(),
                      derivative_eps = .Machine$double.eps^(1 / 3)) {
  call <- match.call()

  method <- match.arg(method)

  jacobian_method <- match.arg(
    jacobian_method
  )

  if (!is.function(Delta)) {
    stop(
      "'Delta' must be a function of gamma that returns ",
      "the random-effects scale matrix."
    )
  }

  if (!is.null(Delta_jacobian) &&
      !is.function(Delta_jacobian)) {
    stop(
      "'Delta_jacobian' must be NULL or a function."
    )
  }

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

  if (!is.list(jacobian_method_args)) {
    stop(
      "'jacobian_method_args' must be a list."
    )
  }

  data <- .lmm_t_prepare_data(
    Y = Y,
    X = X,
    Z = Z
  )

  automatic_start <- is.null(start)
  automatic_start_information <- NULL

  if (automatic_start) {
    automatic_start_information <-
      .lmm_t_default_start(
        data = data,
        gamma_structure = gamma_structure
      )

    start <- automatic_start_information$start
  }

  structure <- .lmm_t_prepare_start(
    start = start,
    p_beta = data$p.beta,
    beta_names = data$beta.names
  )

  start <- structure$start
  p <- structure$p

  initial_components <-
    .lmm_t_evaluate_model(
      theta = start,
      data = data,
      structure = structure,
      Delta = Delta,
      nu = nu,
      strict = TRUE
    )

  # Validate the covariance derivative at the starting values.
  .lmm_t_derivative_matrices(
    theta = start,
    data = data,
    structure = structure,
    Delta = Delta,
    Delta_jacobian = Delta_jacobian,
    jacobian_method = jacobian_method,
    jacobian_method_args = jacobian_method_args
  )

  initial_nll <- -sum(
    initial_components$log.densities
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
    components <- .lmm_t_evaluate_model(
      theta = theta,
      data = data,
      structure = structure,
      Delta = Delta,
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
      components$log.densities
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
    components <- .lmm_t_evaluate_model(
      theta = theta,
      data = data,
      structure = structure,
      Delta = Delta,
      nu = nu,
      strict = strict
    )

    if (is.null(components)) {
      return(
        setNames(
          rep(NaN, p),
          structure$parameter.names
        )
      )
    }

    derivatives <-
      .lmm_t_derivative_matrices(
        theta = theta,
        data = data,
        structure = structure,
        Delta = Delta,
        Delta_jacobian = Delta_jacobian,
        jacobian_method = jacobian_method,
        jacobian_method_args = jacobian_method_args
      )

    score_sum <- numeric(p)

    for (i in seq_len(data$n.groups)) {
      residual_i <-
        components$residuals[[i]]

      Sigma_i <-
        components$scale.matrices[[i]]

      Sigma_inverse_i <-
        components$scale.inverses[[i]]

      weight_i <-
        components$weights[i]

      mean_score_i <-
        drop(
          Sigma_inverse_i %*%
          (weight_i * residual_i)
        )

      covariance_score_matrix_i <-
        0.5 *
        Sigma_inverse_i %*%
        (
          weight_i *
          tcrossprod(residual_i) -
          Sigma_i
        ) %*%
        Sigma_inverse_i

      score_sum[
        structure$beta.indices
      ] <-
        score_sum[
          structure$beta.indices
        ] +
        drop(
          crossprod(
            derivatives$mean.jacobians[[i]],
            mean_score_i
          )
        )

      score_sum[
        structure$covariance.indices
      ] <-
        score_sum[
          structure$covariance.indices
        ] +
        drop(
          crossprod(
            derivatives$covariance.jacobians[[i]],
            c(covariance_score_matrix_i)
          )
        )
    }

    setNames(
      score_sum,
      structure$parameter.names
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
      structure$parameter.names
    )
  }

  negative_score <- function(theta) {
    components <- .lmm_t_evaluate_model(
      theta = theta,
      data = data,
      structure = structure,
      Delta = Delta,
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
    components <- .lmm_t_evaluate_model(
      theta = theta,
      data = data,
      structure = structure,
      Delta = Delta,
      nu = nu,
      strict = strict
    )

    derivatives <-
      .lmm_t_derivative_matrices(
        theta = theta,
        data = data,
        structure = structure,
        Delta = Delta,
        Delta_jacobian = Delta_jacobian,
        jacobian_method = jacobian_method,
        jacobian_method_args = jacobian_method_args
      )

    information <- matrix(
      0,
      nrow = p,
      ncol = p
    )

    for (i in seq_len(data$n.groups)) {
      d_i <- data$group.sizes[i]

      Sigma_inverse_i <-
        components$scale.inverses[[i]]

      D_i <-
        derivatives$mean.jacobians[[i]]

      V_i <-
        derivatives$covariance.jacobians[[i]]

      mean_constant_i <-
        (nu + d_i) /
        (nu + d_i + 2)

      covariance_information_i <-
        (nu + d_i) /
        (
          2 *
          (nu + d_i + 2)
        ) *
        kronecker(
          Sigma_inverse_i,
          Sigma_inverse_i
        ) -
        1 /
        (
          2 *
          (nu + d_i + 2)
        ) *
        tcrossprod(
          c(Sigma_inverse_i)
        )

      information[
        structure$beta.indices,
        structure$beta.indices
      ] <-
        information[
          structure$beta.indices,
          structure$beta.indices
        ] +
        mean_constant_i *
        crossprod(
          D_i,
          Sigma_inverse_i %*%
          D_i
        )

      information[
        structure$covariance.indices,
        structure$covariance.indices
      ] <-
        information[
          structure$covariance.indices,
          structure$covariance.indices
        ] +
        crossprod(
          V_i,
          covariance_information_i %*%
          V_i
        )
    }

    information <-
      (
        information +
        t(information)
      ) / 2

    dimnames(information) <- list(
      structure$parameter.names,
      structure$parameter.names
    )

    list(
      information = information,
      Delta.jacobian = derivatives$Delta.jacobian,
      mean.jacobians = derivatives$mean.jacobians,
      covariance.jacobians = derivatives$covariance.jacobians,
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
      lower <- setNames(
        rep(-Inf, p),
        structure$parameter.names
      )

      lower[
        structure$sigma2.index
      ] <- 1e-10
    }

    if (is.null(upper)) {
      upper <- setNames(
        rep(Inf, p),
        structure$parameter.names
      )
    }

    if (!is.null(names(lower))) {
      missing_lower <- setdiff(
        structure$parameter.names,
        names(lower)
      )

      if (length(missing_lower) > 0L) {
        stop(
          "'lower' must contain all parameter names."
        )
      }

      lower <- lower[
        structure$parameter.names
      ]
    }

    if (!is.null(names(upper))) {
      missing_upper <- setdiff(
        structure$parameter.names,
        names(upper)
      )

      if (length(missing_upper) > 0L) {
        stop(
          "'upper' must contain all parameter names."
        )
      }

      upper <- upper[
        structure$parameter.names
      ]
    }

    lower <- setNames(
      rep_len(
        as.numeric(lower),
        p
      ),
      structure$parameter.names
    )

    upper <- setNames(
      rep_len(
        as.numeric(upper),
        p
      ),
      structure$parameter.names
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
      optim_arguments$gr <-
        negative_score
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
    structure$parameter.names
  )

  fitted_components <-
    .lmm_t_evaluate_model(
      theta = estimate,
      data = data,
      structure = structure,
      Delta = Delta,
      nu = nu,
      strict = TRUE
    )

  information_components <-
    fisher_information(
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
        structure$parameter.names
      )
    }
  )

  conditional_components <-
    .lmm_t_compute_random_effects(
      components = fitted_components,
      data = data
    )

  # ----------------------------------------------------------
  # Expected and observed information matrices
  # ----------------------------------------------------------

  vcov_fisher <- NULL
  standard_error <- NULL
  confidence_interval <- NULL
  vcov_warning <- NULL

  if (isTRUE(compute_vcov)) {
    vcov_fisher <- .lmm_t_safe_inverse(
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
      structure$parameter.names
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
      structure$parameter.names
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

        Hessian <- (
          Hessian +
          t(Hessian)
        ) / 2

        dimnames(Hessian) <- list(
          structure$parameter.names,
          structure$parameter.names
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
      vcov_observed <- .lmm_t_safe_inverse(
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
    -2 *
    maximized_loglikelihood +
    2 *
    p

  BIC_groups_value <-
    -2 *
    maximized_loglikelihood +
    log(data$n.groups) *
    p

  BIC_observations_value <-
    -2 *
    maximized_loglikelihood +
    log(data$nobs) *
    p

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

  result <- list(
    call = call,
    starting.values = start,
    automatic.start = automatic_start,
    gamma.structure = if (automatic_start) {
      automatic_start_information$gamma.structure
    } else {
      gamma_structure
    },
    gamma.diagonal.names = if (automatic_start) {
      automatic_start_information$gamma.diagonal.names
    } else {
      NULL
    },
    robust.start.fit = if (automatic_start) {
      automatic_start_information$robust.fit
    } else {
      NULL
    },
    coefficients = estimate,
    beta = fitted_components$beta,
    gamma = fitted_components$gamma,
    sigma2 = fitted_components$sigma2,
    Delta = fitted_components$Delta,
    Delta.variance = fitted_components$Delta.variance,
    residual.variance = fitted_components$residual.variance,
    standard.error = standard_error,
    conf.int = confidence_interval,
    coefficient.table = coefficient_table,
    vcov = vcov_fisher,
    vcov.fisher = vcov_fisher,
    vcov.observed = vcov_observed,
    fisher.information = expected_information,
    observed.information = observed_information,
    Delta.jacobian = information_components$Delta.jacobian,
    mean.jacobians = information_components$mean.jacobians,
    covariance.jacobians = information_components$covariance.jacobians,
    F = information_components$F,
    fitted.values = fitted_components$means,
    marginal.fitted.values = fitted_components$means,
    conditional.fitted.values = conditional_components$conditional.means,
    random.effects = conditional_components$random.effects,
    residuals = fitted_components$residuals,
    marginal.residuals = fitted_components$residuals,
    conditional.residuals = conditional_components$conditional.residuals,
    standardized.residuals = fitted_components$standardized.residuals,
    scale.matrices = fitted_components$scale.matrices,
    covariance.matrices = fitted_components$variance.matrices,
    scale.inverses = fitted_components$scale.inverses,
    quadratic.forms = fitted_components$quadratic.forms,
    weights = fitted_components$weights,
    log.density = fitted_components$log.densities,
    density = fitted_components$densities,
    score = score_at_estimate,
    max.abs.score = max_abs_score,
    logLik = maximized_loglikelihood,
    objective = unname(optimization$value),
    AIC = AIC_value,
    BIC = BIC_groups_value,
    BIC.groups = BIC_groups_value,
    BIC.observations = BIC_observations_value,
    nu = nu,
    level = level,
    n.groups = data$n.groups,
    nobs = data$nobs,
    group.sizes = data$group.sizes,
    p = p,
    p.beta = data$p.beta,
    p.gamma = structure$p.gamma,
    m = data$m,
    beta.names = structure$beta.names,
    gamma.names = structure$gamma.names,
    sigma2.name = structure$sigma2.name,
    random.effect.names = data$random.effect.names,
    convergence = optimization$convergence,
    message = optimization$message,
    counts = optimization$counts,
    method = method,
    use.score = isTRUE(use_score),
    jacobian.method = jacobian_method,
    Y = data$Y,
    X = data$X,
    Z = data$Z,
    Delta.function = Delta,
    Delta.jacobian.function = Delta_jacobian,
    vcov.warning = vcov_warning,
    observed.information.warning = observed_information_warning,
    lower = lower_used,
    upper = upper_used,
    optimization.control = control,
    optim = optimization
  )

  class(result) <- "lmm_t_fit"

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

print.lmm_t_fit <- function(x,
                            digits = max(
                              3L,
                              getOption("digits") - 3L
                            ),
                            ...) {
  cat(
    "Linear mixed model with multivariate Student-t marginal distribution\n"
  )

  cat(
    "fixed degrees of freedom =",
    format(x$nu, digits = digits),
    "| groups =",
    x$n.groups,
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
    "\nEstimated Delta(gamma) scale matrix:\n"
  )

  print(
    x$Delta,
    digits = digits
  )

  cat(
    "\nLog-likelihood:",
    format(x$logLik, digits = digits),
    "\nAIC:",
    format(x$AIC, digits = digits),
    "\nBIC based on independent groups:",
    format(x$BIC.groups, digits = digits),
    "\nMaximum absolute score component:",
    format(x$max.abs.score, digits = digits),
    "\n"
  )

  if (!is.null(x$vcov.warning)) {
    cat(
      "\nWarning:",
      x$vcov.warning,
      "\n"
    )
  }

  if (!is.null(x$observed.information.warning)) {
    cat(
      "\nObserved-information warning:",
      x$observed.information.warning,
      "\n"
    )
  }

  invisible(x)
}


summary.lmm_t_fit <- function(object,
                              ...) {
  result <- list(
    call = object$call,
    coefficients = object$coefficient.table,
    Delta = object$Delta,
    Delta.variance = object$Delta.variance,
    sigma2 = object$sigma2,
    residual.variance = object$residual.variance,
    logLik = object$logLik,
    AIC = object$AIC,
    BIC.groups = object$BIC.groups,
    BIC.observations = object$BIC.observations,
    score = object$score,
    max.abs.score = object$max.abs.score,
    convergence = object$convergence,
    message = object$message,
    nu = object$nu,
    n.groups = object$n.groups,
    nobs = object$nobs,
    method = object$method,
    jacobian.method = object$jacobian.method,
    vcov.warning = object$vcov.warning,
    observed.information.warning = object$observed.information.warning
  )

  class(result) <-
    "summary.lmm_t_fit"

  result
}


print.summary.lmm_t_fit <- function(x,
                                    digits = max(
                                      3L,
                                      getOption("digits") - 3L
                                    ),
                                    ...) {
  cat(
    "Linear mixed model fitted by Student-t maximum likelihood\n"
  )

  cat(
    "fixed degrees of freedom =",
    format(x$nu, digits = digits),
    "| groups =",
    x$n.groups,
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

  cat(
    "\n\nCoefficients:\n"
  )

  printCoefmat(
    x$coefficients,
    digits = digits,
    P.values = ncol(x$coefficients) >= 4L,
    has.Pvalue = ncol(x$coefficients) >= 4L,
    na.print = "NA"
  )

  cat(
    "\nEstimated Delta(gamma) scale matrix:\n"
  )

  print(
    x$Delta,
    digits = digits
  )

  if (!is.null(x$Delta.variance)) {
    cat(
      "\nImplied random-effects covariance matrix:\n"
    )

    print(
      x$Delta.variance,
      digits = digits
    )
  }

  cat(
    "\nResidual squared scale sigma2:",
    format(x$sigma2, digits = digits),
    "\n"
  )

  if (is.finite(x$residual.variance)) {
    cat(
      "Implied residual variance:",
      format(x$residual.variance, digits = digits),
      "\n"
    )
  }

  cat(
    "\nLog-likelihood:",
    format(x$logLik, digits = digits),
    "\nAIC:",
    format(x$AIC, digits = digits),
    "\nBIC based on independent groups:",
    format(x$BIC.groups, digits = digits),
    "\nBIC based on total observations:",
    format(x$BIC.observations, digits = digits),
    "\nMaximum absolute score component:",
    format(x$max.abs.score, digits = digits),
    "\n"
  )

  if (!is.null(x$vcov.warning)) {
    cat(
      "\nWarning:",
      x$vcov.warning,
      "\n"
    )
  }

  if (!is.null(x$observed.information.warning)) {
    cat(
      "\nObserved-information warning:",
      x$observed.information.warning,
      "\n"
    )
  }

  invisible(x)
}


coef.lmm_t_fit <- function(object,
                           ...) {
  object$coefficients
}


vcov.lmm_t_fit <- function(object,
                           type = c(
                             "fisher",
                             "observed"
                           ),
                           ...) {
  type <- match.arg(type)

  if (type == "fisher") {
    if (is.null(object$vcov.fisher)) {
      stop(
        "The Fisher covariance matrix was not computed."
      )
    }

    return(object$vcov.fisher)
  }

  if (is.null(object$vcov.observed)) {
    stop(
      "The observed covariance matrix is not available."
    )
  }

  object$vcov.observed
}


confint.lmm_t_fit <- function(object,
                              parm = seq_along(
                                object$coefficients
                              ),
                              level = object$level,
                              ...) {
  if (is.null(object$standard.error)) {
    stop(
      "Standard errors were not computed."
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

  estimates <- object$coefficients[
    parameter_indices
  ]

  standard_errors <- object$standard.error[
    parameter_indices
  ]

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


fitted.lmm_t_fit <- function(object,
                             type = c(
                               "marginal",
                               "conditional"
                             ),
                             ...) {
  type <- match.arg(type)

  if (type == "marginal") {
    object$marginal.fitted.values
  } else {
    object$conditional.fitted.values
  }
}


residuals.lmm_t_fit <- function(object,
                                type = c(
                                  "marginal",
                                  "conditional",
                                  "standardized"
                                ),
                                ...) {
  type <- match.arg(type)

  switch(
    type,
    marginal = object$marginal.residuals,
    conditional = object$conditional.residuals,
    standardized = object$standardized.residuals
  )
}


predict.lmm_t_fit <- function(object,
                              newdata = NULL,
                              type = c(
                                "mean",
                                "scale",
                                "variance",
                                "random_effects",
                                "conditional_mean"
                              ),
                              ...) {
  type <- match.arg(type)

  if (is.null(newdata)) {
    if (type == "mean") {
      return(
        object$marginal.fitted.values
      )
    }

    if (type == "scale") {
      return(
        object$scale.matrices
      )
    }

    if (type == "variance") {
      if (object$nu <= 2) {
        stop(
          "The multivariate Student-t variance does not exist when nu <= 2."
        )
      }

      return(
        object$covariance.matrices
      )
    }

    if (type == "random_effects") {
      return(
        object$random.effects
      )
    }

    return(
      object$conditional.fitted.values
    )
  }

  if (!is.list(newdata) ||
      is.null(newdata$X) ||
      is.null(newdata$Z)) {
    stop(
      "'newdata' must be a list containing 'X' and 'Z'."
    )
  }

  X_new <- .lmm_t_as_group_list(
    newdata$X,
    "newdata$X"
  )

  Z_new <- .lmm_t_as_group_list(
    newdata$Z,
    "newdata$Z"
  )

  if (length(X_new) != length(Z_new)) {
    stop(
      "'newdata$X' and 'newdata$Z' must contain the same number of groups."
    )
  }

  X_new <- lapply(
    X_new,
    as.matrix
  )

  Z_new <- lapply(
    Z_new,
    as.matrix
  )

  if (any(
    vapply(
      X_new,
      ncol,
      integer(1L)
    ) != object$p.beta
  )) {
    stop(
      "The new fixed-effects design matrices have an incorrect number of columns."
    )
  }

  if (any(
    vapply(
      Z_new,
      ncol,
      integer(1L)
    ) != object$m
  )) {
    stop(
      "The new random-effects design matrices have an incorrect number of columns."
    )
  }

  if (any(
    vapply(
      seq_along(X_new),
      function(index) {
        nrow(X_new[[index]]) !=
        nrow(Z_new[[index]])
      },
      logical(1L)
    )
  )) {
    stop(
      "Each new X matrix must have the same number of rows as its corresponding Z matrix."
    )
  }

  marginal_means <- lapply(
    X_new,
    function(X_i) {
      drop(
        X_i %*%
        object$beta
      )
    }
  )

  scale_matrices <- lapply(
    Z_new,
    function(Z_i) {
      Z_i %*%
      object$Delta %*%
      t(Z_i) +
      object$sigma2 *
      diag(nrow(Z_i))
    }
  )

  if (type == "mean") {
    return(marginal_means)
  }

  if (type == "scale") {
    return(scale_matrices)
  }

  if (type == "variance") {
    if (object$nu <= 2) {
      stop(
        "The multivariate Student-t variance does not exist when nu <= 2."
      )
    }

    return(
      lapply(
        scale_matrices,
        function(Sigma_i) {
          object$nu /
          (object$nu - 2) *
          Sigma_i
        }
      )
    )
  }

  if (is.null(newdata$Y)) {
    stop(
      "'newdata$Y' is required for random-effects or conditional-mean predictions."
    )
  }

  Y_new <- .lmm_t_as_group_list(
    newdata$Y,
    "newdata$Y"
  )

  Y_new <- lapply(
    Y_new,
    as.numeric
  )

  if (length(Y_new) != length(X_new)) {
    stop(
      "'newdata$Y', 'newdata$X' and 'newdata$Z' must contain the same number of groups."
    )
  }

  random_effects <- vector(
    "list",
    length(Y_new)
  )

  conditional_means <- vector(
    "list",
    length(Y_new)
  )

  for (i in seq_along(Y_new)) {
    if (length(Y_new[[i]]) !=
        length(marginal_means[[i]])) {
      stop(
        "The response and design matrices have incompatible dimensions in new group ",
        i, "."
      )
    }

    scale_inverse_i <- chol2inv(
      chol(scale_matrices[[i]])
    )

    random_effects[[i]] <- drop(
      object$Delta %*%
      t(Z_new[[i]]) %*%
      scale_inverse_i %*%
      (
        Y_new[[i]] -
        marginal_means[[i]]
      )
    )

    names(random_effects[[i]]) <-
      object$random.effect.names

    conditional_means[[i]] <-
      marginal_means[[i]] +
      drop(
        Z_new[[i]] %*%
        random_effects[[i]]
      )
  }

  if (type == "random_effects") {
    random_effects
  } else {
    conditional_means
  }
}


logLik.lmm_t_fit <- function(object,
                             ...) {
  structure(
    object$logLik,
    df = object$p,
    nobs = object$n.groups,
    class = "logLik"
  )
}


nobs.lmm_t_fit <- function(object,
                           ...) {
  object$nobs
}

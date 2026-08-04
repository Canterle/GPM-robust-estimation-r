# Generic MDPDE fitting for Gaussian linear mixed models
#
# For group i, the model is
#
#   Y_i = X_i beta + Z_i b_i + epsilon_i,
#
# where
#
#   b_i       ~ Normal(0, Delta(gamma)),
#   epsilon_i ~ Normal(0, sigma2 I_{d_i}).
#
# Consequently,
#
#   Y_i ~ Normal(
#     X_i beta,
#     Z_i Delta(gamma) Z_i^T + sigma2 I_{d_i}
#   ).
#
# The covariance structure of the random effects is supplied through
# the argument
#
#   Delta(gamma),
#
# which must be a function returning an m by m symmetric positive
# definite matrix. This permits, for example:
#
#   - a diagonal covariance matrix;
#   - an unstructured covariance matrix;
#   - a compound-symmetry covariance matrix;
#   - a covariance matrix defined through a Cholesky factor;
#   - any other differentiable parameterization.
#
# The optional function Delta_jacobian(gamma) may return
#
#   d vec(Delta(gamma)) / d gamma^T,
#
# as either:
#
#   - an m^2 by p_gamma matrix; or
#   - an m by m by p_gamma array.
#
# If Delta_jacobian is omitted, the derivative is computed numerically
# with numDeriv::jacobian().
#
# When start = NULL, robust initial values are computed as follows:
#
#   - beta is obtained from MASS::rlm() applied to the stacked data;
#   - sigma2 is the squared robust residual scale returned by rlm();
#   - gamma parameters appearing on the diagonal of gamma_structure
#     receive starting value 1;
#   - gamma parameters appearing only off the diagonal receive 0.
#
# Optimization may be performed with BFGS, L-BFGS-B, Nelder-Mead,
# or BOBYQA. The BOBYQA option uses nloptr::bobyqa(), corresponding
# to the NLOPT_LN_BOBYQA algorithm, and requires the nloptr package.
#
# The parameter vector supplied in start is ordered as
#
#   theta = (beta, gamma, sigma2),
#
# where:
#
#   - the first ncol(X[[1]]) values are fixed-effect parameters;
#   - the final value is the residual variance sigma2;
#   - the intermediate values parameterize Delta(gamma).


# ============================================================
# Internal utility functions
# ============================================================

.lmm_mdpde_as_group_list <- function(value,
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


.lmm_mdpde_prepare_data <- function(Y,
                                    X,
                                    Z) {
  Y <- .lmm_mdpde_as_group_list(
    Y,
    "Y"
  )

  X <- .lmm_mdpde_as_group_list(
    X,
    "X"
  )

  Z <- .lmm_mdpde_as_group_list(
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
    function(value) {
      as.numeric(value)
    }
  )

  X <- lapply(
    X,
    function(value) {
      as.matrix(value)
    }
  )

  Z <- lapply(
    Z,
    function(value) {
      as.matrix(value)
    }
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

  list(
    Y = Y,
    X = X,
    Z = Z,
    n.groups = n_groups,
    nobs = sum(group_sizes),
    group.sizes = group_sizes,
    p.beta = p_beta,
    m = m
  )
}


# Validate the matrix identifying where each gamma parameter appears
# in Delta(gamma).
#
# Entries must be gamma names or NA. A gamma appearing in at least
# one diagonal position receives starting value 1. A gamma appearing
# only outside the main diagonal receives starting value 0.
.lmm_mdpde_prepare_gamma_structure <- function(gamma_structure,
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

  nonempty_entries <-
    gamma_structure[
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
    as.vector(
      gamma_structure
    )
  )

  gamma_names <- gamma_names[
    !is.na(gamma_names)
  ]

  diagonal_names <- unique(
    diag(
      gamma_structure
    )
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
#
# The fixed effects and residual scale are obtained from MASS::rlm()
# fitted to the stacked marginal model, initially ignoring the
# within-group covariance structure.
#
# Gamma starting values are defined by gamma_structure:
#
#   - 1 for gamma parameters appearing on the main diagonal;
#   - 0 for gamma parameters appearing only off the diagonal.
.lmm_mdpde_default_start <- function(data,
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
    .lmm_mdpde_prepare_gamma_structure(
      gamma_structure =
        gamma_structure,
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
    stats::coef(
      robust_fit
    )
  )

  if (length(beta_start) !=
      data$p.beta ||
      anyNA(beta_start) ||
      any(!is.finite(beta_start))) {
    stop(
      "MASS::rlm() did not produce a finite fixed-effect ",
      "starting vector of the expected length."
    )
  }

  beta_names <- colnames(
    X_stacked
  )

  if (is.null(beta_names) ||
      any(beta_names == "")) {
    beta_names <- paste0(
      "beta",
      seq_len(
        data$p.beta
      )
    )
  }

  beta_names <- make.unique(
    beta_names
  )

  names(beta_start) <-
    beta_names

  sigma2_start <- robust_fit$s^2

  if (length(sigma2_start) != 1L ||
      !is.finite(sigma2_start) ||
      sigma2_start < 0.05) {
    sigma2_start <- 0.1
  }

  start <- c(
    beta_start,
    gamma_information$start,
    sigma2 =
      sigma2_start
  )

  list(
    start = start,
    gamma.structure =
      gamma_information$structure,
    gamma.diagonal.names =
      gamma_information$diagonal.names,
    robust.fit =
      robust_fit
  )
}


.lmm_mdpde_prepare_start <- function(start,
                                     p_beta) {
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

  if (is.null(original_names)) {
    parameter_names <- c(
      paste0(
        "beta",
        seq_len(p_beta)
      ),
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
        paste0(
          "beta",
          seq_len(p_beta)
        ),
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

  beta_names <- parameter_names[
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
      "The starting value of the residual variance must be strictly positive."
    )
  }

  list(
    start = start,
    p = p,
    p.gamma = p_gamma,
    beta.indices = beta_indices,
    gamma.indices = gamma_indices,
    sigma2.index = sigma2_index,
    parameter.names = parameter_names,
    beta.names = beta_names,
    gamma.names = gamma_names,
    sigma2.name = sigma2_name
  )
}


.lmm_mdpde_call_Delta <- function(Delta,
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

  dimnames(output) <- list(
    colnames(output),
    colnames(output)
  )

  list(
    value = output,
    chol = cholesky
  )
}


.lmm_mdpde_compute_Delta_jacobian <- function(gamma,
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
            c(
              derivative[, , index]
            )
          }
        )
      )
    } else {
      derivative <- as.matrix(
        derivative
      )
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

      delta_result <- .lmm_mdpde_call_Delta(
        Delta = Delta,
        gamma = value,
        m = m,
        strict = TRUE
      )

      c(
        delta_result$value
      )
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
            method.args =
              method_args
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

    derivative <- as.matrix(
      derivative
    )
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

  colnames(derivative) <-
    names(gamma)

  rownames(derivative) <-
    paste0(
      "Delta",
      seq_len(m^2)
    )

  derivative
}


.lmm_mdpde_block_matrix <- function(top_left,
                                    bottom_right) {
  top_dimension <- nrow(top_left)
  bottom_dimension <- nrow(bottom_right)

  result <- matrix(
    0,
    nrow =
      top_dimension +
      bottom_dimension,
    ncol =
      top_dimension +
      bottom_dimension
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


.lmm_mdpde_evaluate_model <- function(theta,
                                      data,
                                      structure,
                                      Delta,
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

  delta_result <- .lmm_mdpde_call_Delta(
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

  covariance_matrices <- vector(
    "list",
    data$n.groups
  )

  covariance_inverses <- vector(
    "list",
    data$n.groups
  )

  residuals <- vector(
    "list",
    data$n.groups
  )

  log_determinants <- numeric(
    data$n.groups
  )

  log_densities <- numeric(
    data$n.groups
  )

  for (i in seq_len(
    data$n.groups
  )) {
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
          "The marginal covariance matrix is not positive definite ",
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

    quadratic_form <- sum(
      crossprod(
        root_inverse,
        residual_i
      )^2
    )

    log_determinant <- 2 *
      sum(
        log(
          diag(chol_Sigma)
        )
      )

    log_density <-
      -0.5 *
      d_i *
      log(2 * pi) -
      0.5 *
      log_determinant -
      0.5 *
      quadratic_form

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
    covariance_matrices[[i]] <-
      Sigma_i

    covariance_inverses[[i]] <-
      Sigma_inverse

    residuals[[i]] <- residual_i
    log_determinants[i] <-
      log_determinant

    log_densities[i] <-
      log_density
  }

  list(
    beta = beta,
    gamma = gamma,
    sigma2 = unname(sigma2),
    Delta = delta_matrix,
    means = means,
    covariance.matrices =
      covariance_matrices,
    covariance.inverses =
      covariance_inverses,
    residuals = residuals,
    log.determinants =
      log_determinants,
    log.densities =
      log_densities,
    densities =
      exp(log_densities)
  )
}


.lmm_mdpde_derivative_matrices <- function(theta,
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
    .lmm_mdpde_compute_Delta_jacobian(
      gamma = gamma,
      Delta = Delta,
      Delta_jacobian =
        Delta_jacobian,
      m = data$m,
      method =
        jacobian_method,
      method_args =
        jacobian_method_args
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

  for (i in seq_len(
    data$n.groups
  )) {
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
      c(
        diag(d_i)
      )
    )

    F_i <- rbind(
      cbind(
        D_i,
        matrix(
          0,
          nrow = d_i,
          ncol =
            structure$p.gamma +
            1L
        )
      ),
      cbind(
        matrix(
          0,
          nrow = d_i^2,
          ncol =
            data$p.beta
        ),
        V_i
      )
    )

    colnames(D_i) <-
      structure$beta.names

    colnames(V_i) <- c(
      structure$gamma.names,
      structure$sigma2.name
    )

    colnames(F_i) <-
      structure$parameter.names

    mean_jacobians[[i]] <-
      D_i

    covariance_jacobians[[i]] <-
      V_i

    F_matrices[[i]] <-
      F_i
  }

  list(
    Delta.jacobian =
      delta_derivative,
    mean.jacobians =
      mean_jacobians,
    covariance.jacobians =
      covariance_jacobians,
    F = F_matrices
  )
}


.lmm_mdpde_compute_random_effects <- function(components,
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

  for (i in seq_len(
    data$n.groups
  )) {
    b_hat_i <-
      components$Delta %*%
      t(data$Z[[i]]) %*%
      components$covariance.inverses[[i]] %*%
      components$residuals[[i]]

    b_hat_i <- drop(
      b_hat_i
    )

    names(b_hat_i) <-
      colnames(
        data$Z[[i]]
      )

    conditional_mean_i <-
      components$means[[i]] +
      drop(
        data$Z[[i]] %*%
        b_hat_i
      )

    random_effects[[i]] <-
      b_hat_i

    conditional_means[[i]] <-
      conditional_mean_i

    conditional_residuals[[i]] <-
      data$Y[[i]] -
      conditional_mean_i
  }

  list(
    random.effects =
      random_effects,
    conditional.means =
      conditional_means,
    conditional.residuals =
      conditional_residuals
  )
}


# ============================================================
# Fixed-q MDPDE fit
# ============================================================

.fit_lmm_mdpde_fixed <- function(Y,
                                 X,
                                 Z,
                                 Delta,
                                 start = NULL,
                                 gamma_structure = NULL,
                                 Delta_jacobian = NULL,
                                 q = 1,
                                 level = 0.95,
                                 method = c(
                                   "BFGS",
                                   "L-BFGS-B",
                                   "Nelder-Mead",
                                   "BOBYQA"
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
  call <- match.call()
  method <- match.arg(method)
  jacobian_method <-
    match.arg(jacobian_method)

  if (!is.function(Delta)) {
    stop(
      "'Delta' must be a function of gamma that returns ",
      "the random-effects covariance matrix."
    )
  }

  if (!is.null(Delta_jacobian) &&
      !is.function(Delta_jacobian)) {
    stop(
      "'Delta_jacobian' must be NULL or a function."
    )
  }

  if (!is.list(control)) {
    stop("'control' must be a list.")
  }

  if (!is.list(
    jacobian_method_args
  )) {
    stop(
      "'jacobian_method_args' must be a list."
    )
  }

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

  data <- .lmm_mdpde_prepare_data(
    Y = Y,
    X = X,
    Z = Z
  )

  automatic_start <- is.null(
    start
  )

  automatic_start_information <- NULL

  if (automatic_start) {
    automatic_start_information <-
      .lmm_mdpde_default_start(
        data = data,
        gamma_structure =
          gamma_structure
      )

    start <-
      automatic_start_information$start
  }

  structure <- .lmm_mdpde_prepare_start(
    start = start,
    p_beta = data$p.beta
  )

  start <- structure$start
  p <- structure$p

  initial_components <-
    .lmm_mdpde_evaluate_model(
      theta = start,
      data = data,
      structure = structure,
      Delta = Delta,
      strict = TRUE
    )

  .lmm_mdpde_derivative_matrices(
    theta = start,
    data = data,
    structure = structure,
    Delta = Delta,
    Delta_jacobian =
      Delta_jacobian,
    jacobian_method =
      jacobian_method,
    jacobian_method_args =
      jacobian_method_args
  )

  objective <- function(theta) {
    components <-
      .lmm_mdpde_evaluate_model(
        theta = theta,
        data = data,
        structure = structure,
        Delta = Delta,
        strict = FALSE
      )

    if (is.null(components)) {
      return(-Inf)
    }

    if (q == 1) {
      value <- sum(
        components$log.densities
      )
    } else {
      one_minus_q <- 1 - q

      density_power <- exp(
        one_minus_q *
        components$log.densities
      )

      correction <- exp(
        -0.5 *
        data$group.sizes *
        log(2 - q) -
        0.5 *
        one_minus_q *
        data$group.sizes *
        log(2 * pi) -
        0.5 *
        one_minus_q *
        components$log.determinants
      )

      value <- sum(
        (2 - q) /
        one_minus_q *
        density_power -
        correction
      )
    }

    if (is.finite(value)) {
      value
    } else {
      -Inf
    }
  }

  score <- function(theta,
                    strict = FALSE) {
    components <-
      .lmm_mdpde_evaluate_model(
        theta = theta,
        data = data,
        structure = structure,
        Delta = Delta,
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
      .lmm_mdpde_derivative_matrices(
        theta = theta,
        data = data,
        structure = structure,
        Delta = Delta,
        Delta_jacobian =
          Delta_jacobian,
        jacobian_method =
          jacobian_method,
        jacobian_method_args =
          jacobian_method_args
      )

    one_minus_q <- 1 - q

    score_sum <- rep(
      0,
      p
    )

    for (i in seq_len(
      data$n.groups
    )) {
      d_i <- data$group.sizes[i]

      Sigma_i <-
        components$covariance.matrices[[i]]

      Sigma_inverse_i <-
        components$covariance.inverses[[i]]

      H_i <- .lmm_mdpde_block_matrix(
        top_left =
          Sigma_inverse_i,
        bottom_right =
          0.5 *
          kronecker(
            Sigma_inverse_i,
            Sigma_inverse_i
          )
      )

      residual_i <-
        components$residuals[[i]]

      s_i <- c(
        residual_i,
        c(
          tcrossprod(
            residual_i
          ) -
          Sigma_i
        )
      )

      density_power_i <- if (
        q == 1
      ) {
        1
      } else {
        exp(
          one_minus_q *
          components$log.densities[i]
        )
      }

      correction_scalar_i <- if (
        q == 1
      ) {
        0
      } else {
        (q - 1) *
        exp(
          -0.5 *
          (d_i + 2) *
          log(2 - q) -
          0.5 *
          one_minus_q *
          d_i *
          log(2 * pi) -
          0.5 *
          one_minus_q *
          components$log.determinants[i]
        )
      }

      correction_vector_i <- c(
        rep(0, d_i),
        correction_scalar_i *
        c(Sigma_i)
      )

      score_sum <-
        score_sum +
        as.numeric(
          (2 - q) *
          t(
            derivatives$F[[i]]
          ) %*%
          H_i %*%
          (
            s_i *
            density_power_i -
            correction_vector_i
          )
        )
    }

    setNames(
      score_sum,
      structure$parameter.names
    )
  }

  compute_JK <- function(theta) {
    components <-
      .lmm_mdpde_evaluate_model(
        theta = theta,
        data = data,
        structure = structure,
        Delta = Delta,
        strict = TRUE
      )

    derivatives <-
      .lmm_mdpde_derivative_matrices(
        theta = theta,
        data = data,
        structure = structure,
        Delta = Delta,
        Delta_jacobian =
          Delta_jacobian,
        jacobian_method =
          jacobian_method,
        jacobian_method_args =
          jacobian_method_args
      )

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

    for (i in seq_len(
      data$n.groups
    )) {
      d_i <- data$group.sizes[i]

      Sigma_inverse_i <-
        components$covariance.inverses[[i]]

      Sigma_inverse_vector_i <- c(
        Sigma_inverse_i
      )

      covariance_block_J <-
        (
          0.5 *
          kronecker(
            Sigma_inverse_i,
            Sigma_inverse_i
          ) +
          0.25 *
          one_minus_q^2 *
          tcrossprod(
            Sigma_inverse_vector_i
          )
        ) /
        (2 - q)

      H_J_i <- .lmm_mdpde_block_matrix(
        top_left =
          Sigma_inverse_i,
        bottom_right =
          covariance_block_J
      )

      c3_i <- exp(
        -0.5 *
        d_i *
        log(2 - q) -
        0.5 *
        one_minus_q *
        d_i *
        log(2 * pi) -
        0.5 *
        one_minus_q *
        components$log.determinants[i]
      )

      J <-
        J +
        c3_i *
        t(
          derivatives$F[[i]]
        ) %*%
        H_J_i %*%
        derivatives$F[[i]]

      kappa_i <-
        (2 - 2 * q)^2 -
        one_minus_q^2 *
        (3 - 2 * q)^(
          2 +
          d_i / 2
        ) /
        (2 - q)^(
          d_i +
          2
        )

      covariance_block_K <-
        (
          0.5 *
          kronecker(
            Sigma_inverse_i,
            Sigma_inverse_i
          ) +
          0.25 *
          kappa_i *
          tcrossprod(
            Sigma_inverse_vector_i
          )
        ) /
        (3 - 2 * q)

      H_K_i <- .lmm_mdpde_block_matrix(
        top_left =
          Sigma_inverse_i,
        bottom_right =
          covariance_block_K
      )

      c4_i <- exp(
        -0.5 *
        d_i *
        log(3 - 2 * q) -
        one_minus_q *
        d_i *
        log(2 * pi) -
        one_minus_q *
        components$log.determinants[i]
      )

      K <-
        K +
        c4_i *
        t(
          derivatives$F[[i]]
        ) %*%
        H_K_i %*%
        derivatives$F[[i]]
    }

    K <-
      (2 - q)^2 /
      (3 - 2 * q) *
      K

    J <- (
      J +
      t(J)
    ) / 2

    K <- (
      K +
      t(K)
    ) / 2

    dimnames(J) <- list(
      structure$parameter.names,
      structure$parameter.names
    )

    dimnames(K) <- list(
      structure$parameter.names,
      structure$parameter.names
    )

    list(
      J = J,
      K = K,
      Delta.jacobian =
        derivatives$Delta.jacobian,
      mean.jacobians =
        derivatives$mean.jacobians,
      covariance.jacobians =
        derivatives$covariance.jacobians,
      F = derivatives$F
    )
  }

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

    # Translate common optim-style controls to NLopt controls.
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
    control$fnscale <- NULL

    default_control <- list(
      maxeval = 10000,
      xtol_rel = 1e-8,
      ftol_rel = 1e-8
    )
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
  }

  if (method == "BOBYQA") {
    if (use_score) {
      warning(
        "'use_score' is ignored when method = 'BOBYQA'."
      )
    }

    bobyqa_objective <- function(theta) {
      value <- objective(theta)

      if (is.finite(value)) {
        -value
      } else {
        1e100
      }
    }

    bobyqa_fit <- nloptr::bobyqa(
      x0 = start,
      fn = bobyqa_objective,
      lower = lower,
      upper = upper,
      control = control
    )

    raw_convergence <- as.integer(
      bobyqa_fit$convergence
    )

    optimization <- list(
      par = bobyqa_fit$par,
      value = if (is.finite(bobyqa_fit$value)) {
        -bobyqa_fit$value
      } else {
        NA_real_
      },
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
      optim_arguments$lower <- lower
      optim_arguments$upper <- upper
    }

    optimization <- do.call(
      stats::optim,
      optim_arguments
    )
  }

  estimate <- setNames(
    as.numeric(
      optimization$par
    ),
    structure$parameter.names
  )

  fitted_components <-
    .lmm_mdpde_evaluate_model(
      theta = estimate,
      data = data,
      structure = structure,
      Delta = Delta,
      strict = TRUE
    )

  derivative_components <-
    .lmm_mdpde_derivative_matrices(
      theta = estimate,
      data = data,
      structure = structure,
      Delta = Delta,
      Delta_jacobian =
        Delta_jacobian,
      jacobian_method =
        jacobian_method,
      jacobian_method_args =
        jacobian_method_args
    )

  conditional_components <-
    .lmm_mdpde_compute_random_effects(
      components =
        fitted_components,
      data = data
    )

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

  vcov_matrix <- NULL
  standard_error <- NULL
  confidence_interval <- NULL
  J_matrix <- NULL
  K_matrix <- NULL
  vcov_warning <- NULL

  if (compute_vcov) {
    jk <- compute_JK(
      estimate
    )

    J_matrix <- jk$J
    K_matrix <- jk$K

    derivative_components <-
      list(
        Delta.jacobian =
          jk$Delta.jacobian,
        mean.jacobians =
          jk$mean.jacobians,
        covariance.jacobians =
          jk$covariance.jacobians,
        F = jk$F
      )

    covariance_result <- tryCatch(
      {
        J_inverse <- chol2inv(
          chol(J_matrix)
        )

        covariance <-
          J_inverse %*%
          K_matrix %*%
          t(J_inverse)

        covariance <- (
          covariance +
          t(covariance)
        ) / 2

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
            ) / 2

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
                  structure$parameter.names,
                  structure$parameter.names
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
      structure$parameter.names,
      structure$parameter.names
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
        covariance_diagonal < 0,
        NA_real_,
        sqrt(
          covariance_diagonal
        )
      ),
      structure$parameter.names
    )

    z_value <- qnorm(
      1 -
      (1 - level) / 2
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
    ) <- structure$parameter.names
  }

  result <- list(
    call = call,
    starting.values = start,
    automatic.start =
      automatic_start,
    gamma.structure =
      if (automatic_start) {
        automatic_start_information$gamma.structure
      } else {
        gamma_structure
      },
    gamma.diagonal.names =
      if (automatic_start) {
        automatic_start_information$gamma.diagonal.names
      } else {
        NULL
      },
    robust.start.fit =
      if (automatic_start) {
        automatic_start_information$robust.fit
      } else {
        NULL
      },
    coefficients = estimate,
    beta =
      fitted_components$beta,
    gamma =
      fitted_components$gamma,
    sigma2 =
      fitted_components$sigma2,
    Delta =
      fitted_components$Delta,
    standard.error =
      standard_error,
    vcov =
      vcov_matrix,
    conf.int =
      confidence_interval,
    J = J_matrix,
    K = K_matrix,
    Delta.jacobian =
      derivative_components$Delta.jacobian,
    mean.jacobians =
      derivative_components$mean.jacobians,
    covariance.jacobians =
      derivative_components$covariance.jacobians,
    F =
      derivative_components$F,
    fitted.values =
      fitted_components$means,
    marginal.fitted.values =
      fitted_components$means,
    conditional.fitted.values =
      conditional_components$conditional.means,
    random.effects =
      conditional_components$random.effects,
    residuals =
      fitted_components$residuals,
    marginal.residuals =
      fitted_components$residuals,
    conditional.residuals =
      conditional_components$conditional.residuals,
    covariance.matrices =
      fitted_components$covariance.matrices,
    covariance.inverses =
      fitted_components$covariance.inverses,
    log.density =
      fitted_components$log.densities,
    density =
      fitted_components$densities,
    objective =
      unname(
        optimization$value
      ),
    score =
      score_at_estimate,
    convergence =
      optimization$convergence,
    message =
      optimization$message,
    counts =
      optimization$counts,
    q = q,
    level = level,
    n.groups =
      data$n.groups,
    nobs =
      data$nobs,
    group.sizes =
      data$group.sizes,
    p.beta =
      data$p.beta,
    p.gamma =
      structure$p.gamma,
    m =
      data$m,
    beta.names =
      structure$beta.names,
    gamma.names =
      structure$gamma.names,
    sigma2.name =
      structure$sigma2.name,
    method = method,
    jacobian.method =
      jacobian_method,
    Y = data$Y,
    X = data$X,
    Z = data$Z,
    Delta.function =
      Delta,
    Delta.jacobian.function =
      Delta_jacobian,
    vcov.warning =
      vcov_warning,
    lower = lower_used,
    upper = upper_used,
    optimization.control =
      control,
    optim = optimization
  )

  class(result) <- "lmm_mdpde"

  result
}


# ============================================================
# Public fitting interface
# ============================================================

fit_lmm_mdpde <- function(Y,
                          X,
                          Z,
                          Delta,
                          start = NULL,
                          gamma_structure = NULL,
                          Delta_jacobian = NULL,
                          q = 1,
                          level = 0.95,
                          method = c(
                            "BFGS",
                            "L-BFGS-B",
                            "Nelder-Mead",
                            "BOBYQA"
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
    fit <- .fit_lmm_mdpde_fixed(
      Y = Y,
      X = X,
      Z = Z,
      Delta = Delta,
      start = start,
      gamma_structure =
        gamma_structure,
      Delta_jacobian =
        Delta_jacobian,
      q = q,
      level = level,
      method = method,
      lower = lower,
      upper = upper,
      control = control,
      use_score = use_score,
      compute_vcov =
        compute_vcov,
      jacobian_method =
        jacobian_method,
      jacobian_method_args =
        jacobian_method_args
    )

    fit$call <- call

    return(fit)
  }

  if (!isTRUE(
    compute_vcov
  )) {
    stop(
      "SQV selection requires 'compute_vcov = TRUE' because ",
      "standard errors are used to construct the standardized estimates."
    )
  }

  fit <- .select_q_lmm_mdpde_sqv(
    Y = Y,
    X = X,
    Z = Z,
    Delta = Delta,
    start = start,
    gamma_structure =
      gamma_structure,
    Delta_jacobian =
      Delta_jacobian,
    level = level,
    method = method,
    lower = lower,
    upper = upper,
    control = control,
    use_score = use_score,
    jacobian_method =
      jacobian_method,
    jacobian_method_args =
      jacobian_method_args,
    q_control = q_control
  )

  fit$call <- call

  fit
}


# ============================================================
# SQV selection algorithm
# ============================================================

.select_q_lmm_mdpde_sqv <- function(Y,
                                    X,
                                    Z,
                                    Delta,
                                    start = NULL,
                                    gamma_structure = NULL,
                                    Delta_jacobian,
                                    level,
                                    method,
                                    lower,
                                    upper,
                                    control,
                                    use_score,
                                    jacobian_method,
                                    jacobian_method_args,
                                    q_control) {
  if (!is.list(q_control)) {
    stop("'q_control' must be a list.")
  }

  method <- match.arg(
    method,
    c(
      "BFGS",
      "L-BFGS-B",
      "Nelder-Mead",
      "BOBYQA"
    )
  )

  if (method == "BOBYQA") {
    use_score <- FALSE
  }

  automatic_start <- is.null(start)

  prepared_data <- .lmm_mdpde_prepare_data(
    Y = Y,
    X = X,
    Z = Z
  )

  automatic_start_information <- NULL

  if (automatic_start) {
    automatic_start_information <-
      .lmm_mdpde_default_start(
        data = prepared_data,
        gamma_structure = gamma_structure
      )

    start <- automatic_start_information$start
  }

  start_structure <- .lmm_mdpde_prepare_start(
    start = start,
    p_beta = prepared_data$p.beta
  )

  start <- start_structure$start
  p <- start_structure$p

  default_q_control <- list(
    m0 = 21L,
    m = 3L,
    q_min = 0.5,
    L = 0.005,
    verbose = FALSE,
    keep_fits = FALSE
  )

  q_control <- modifyList(
    default_q_control,
    q_control
  )

  # Only the controls appearing in the original SQV algorithm
  # and output controls are retained. Legacy entries have no effect.
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
        "'", name, "' must be an integer greater than or equal to ",
        minimum, "."
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
      "restart endpoint is q_(m - 1) from the initial grid."
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
    stop("'q_control$L' must be positive.")
  }

  if (!scalar_logical(q_control$verbose)) {
    stop("'q_control$verbose' must be TRUE or FALSE.")
  }

  if (!scalar_logical(q_control$keep_fits)) {
    stop("'q_control$keep_fits' must be TRUE or FALSE.")
  }

  n <- prepared_data$n.groups
  q_min <- as.numeric(q_control$q_min)
  L <- as.numeric(q_control$L)

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
      .fit_lmm_mdpde_fixed(
          Y = prepared_data$Y,
          X = prepared_data$X,
          Z = prepared_data$Z,
          Delta = Delta,
          start = start,
          gamma_structure = gamma_structure,
          Delta_jacobian = Delta_jacobian,
          q = q_value,
          level = level,
          method = method,
          lower = lower,
          upper = upper,
          control = control,
          use_score = use_score,
          compute_vcov = TRUE,
          jacobian_method = jacobian_method,
          jacobian_method_args = jacobian_method_args
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
      stop("Every SQV grid must be strictly decreasing.")
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

    # The original algorithm uses the smallest q_k at which
    # SQV_qk >= L. Because the grid is descending, this is the
    # last violating pair in the grid.
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
      stop("The SQV violation index could not be determined.")
    }

    spacing <-
      grid_result$q_grid[1L] -
      grid_result$q_grid[2L]

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

  first_pass <- run_pass(
    pass = 1L,
    pass_initial_grid = initial_grid
  )

  if (first_pass$success) {
    selection_result <- first_pass
  } else {
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

  selected_fit$gamma.structure <-
    if (automatic_start) {
      automatic_start_information$gamma.structure
    } else {
      gamma_structure
    }

  selected_fit$gamma.diagonal.names <-
    if (automatic_start) {
      automatic_start_information$gamma.diagonal.names
    } else {
      NULL
    }

  selected_fit$robust.start.fit <-
    if (automatic_start) {
      automatic_start_information$robust.fit
    } else {
      NULL
    }

  selected_fit$implementation.version <-
    "LMM-MDPDE-SQV-original-algorithm-2026-07-30"

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
    "lmm_mdpde_sqv",
    class(selected_fit)
  )

  selected_fit
}


# ============================================================
# S3 methods
# ============================================================

print.lmm_mdpde <- function(x,
                            digits = max(
                              3L,
                              getOption(
                                "digits"
                              ) -
                              3L
                            ),
                            ...) {
  cat(
    "Gaussian linear mixed model fitted by MDPDE\n"
  )

  cat(
    "q =",
    format(
      x$q,
      digits = digits
    ),
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

  if (!is.null(
    x$message
  )) {
    cat(
      " -",
      x$message
    )
  }

  cat("\n\n")

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

  cat(
    "\nEstimated Delta(gamma):\n"
  )

  print(
    x$Delta,
    digits = digits
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

  invisible(x)
}


summary.lmm_mdpde <- function(object,
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
    call =
      object$call,
    coefficients =
      coefficient_table,
    Delta =
      object$Delta,
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
    n.groups =
      object$n.groups,
    nobs =
      object$nobs,
    method =
      object$method,
    jacobian.method =
      object$jacobian.method,
    vcov.warning =
      object$vcov.warning
  )

  class(
    result
  ) <- "summary.lmm_mdpde"

  result
}


print.summary.lmm_mdpde <- function(x,
                                    digits = max(
                                      3L,
                                      getOption(
                                        "digits"
                                      ) -
                                      3L
                                    ),
                                    ...) {
  cat(
    "Gaussian linear mixed model fitted by MDPDE\n"
  )

  cat(
    "q =",
    format(
      x$q,
      digits = digits
    ),
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
    "\nEstimated Delta(gamma):\n"
  )

  print(
    x$Delta,
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

  invisible(x)
}


coef.lmm_mdpde <- function(object,
                           ...) {
  object$coefficients
}


vcov.lmm_mdpde <- function(object,
                           ...) {
  object$vcov
}


confint.lmm_mdpde <- function(object,
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

  if (anyNA(
    parameter_indices
  )) {
    stop(
      "At least one requested parameter name was not found."
    )
  }

  z_value <- qnorm(
    1 -
    (1 - level) / 2
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

  rownames(interval) <-
    names(estimates)

  interval
}


fitted.lmm_mdpde <- function(object,
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


residuals.lmm_mdpde <- function(object,
                                type = c(
                                  "marginal",
                                  "conditional"
                                ),
                                ...) {
  type <- match.arg(type)

  if (type == "marginal") {
    object$marginal.residuals
  } else {
    object$conditional.residuals
  }
}


predict.lmm_mdpde <- function(object,
                              newdata = NULL,
                              type = c(
                                "mean",
                                "covariance",
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

    if (type == "covariance") {
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

  X_new <- .lmm_mdpde_as_group_list(
    newdata$X,
    "newdata$X"
  )

  Z_new <- .lmm_mdpde_as_group_list(
    newdata$Z,
    "newdata$Z"
  )

  if (length(X_new) !=
      length(Z_new)) {
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

  covariance_matrices <- lapply(
    Z_new,
    function(Z_i) {
      Z_i %*%
      object$Delta %*%
      t(Z_i) +
      object$sigma2 *
      diag(
        nrow(Z_i)
      )
    }
  )

  if (type == "mean") {
    return(
      marginal_means
    )
  }

  if (type == "covariance") {
    return(
      covariance_matrices
    )
  }

  if (is.null(newdata$Y)) {
    stop(
      "'newdata$Y' is required for random-effects or conditional-mean predictions."
    )
  }

  Y_new <- .lmm_mdpde_as_group_list(
    newdata$Y,
    "newdata$Y"
  )

  Y_new <- lapply(
    Y_new,
    as.numeric
  )

  if (length(Y_new) !=
      length(X_new)) {
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

  for (i in seq_along(
    Y_new
  )) {
    if (length(Y_new[[i]]) !=
        length(marginal_means[[i]])) {
      stop(
        "The response and design matrices have incompatible dimensions in new group ",
        i, "."
      )
    }

    Sigma_inverse_i <- chol2inv(
      chol(
        covariance_matrices[[i]]
      )
    )

    random_effects[[i]] <-
      drop(
        object$Delta %*%
        t(Z_new[[i]]) %*%
        Sigma_inverse_i %*%
        (
          Y_new[[i]] -
          marginal_means[[i]]
        )
      )

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


print.lmm_mdpde_sqv <- function(x,
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

  NextMethod("print")
}


summary.lmm_mdpde_sqv <- function(object,
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

  class(result) <- c(
    "summary.lmm_mdpde_sqv",
    class(result)
  )

  result
}


print.summary.lmm_mdpde_sqv <- function(x,
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

  NextMethod("print")
}

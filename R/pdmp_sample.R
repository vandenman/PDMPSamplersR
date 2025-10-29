# Internal function to validate and prepare common PDMP parameters
validate_pdmp_params <- function(d, flow, algorithm, T, t0 = 0.0,
                                flow_mean = NULL, flow_cov = NULL, c0 = 1e-2,
                                x0 = NULL, theta0 = NULL, show_progress = TRUE,
                                discretize_dt = NULL,
                                sticky = FALSE, can_stick = NULL, model_prior = NULL, parameter_prior = NULL,
                                grid_n = 30, grid_t_max = 2.0) {

  # Validate basic parameters
  d <- cast_integer(d, n = 1)
  validate_type(d, type = "integer", n = 1, positive = TRUE)
  validate_type(T, type = "double", n = 1, positive = TRUE)

  flow <- match.arg(flow, c("ZigZag", "BouncyParticle", "Boomerang"))
  algorithm <- match.arg(algorithm, c("ThinningStrategy", "GridThinningStrategy", "RootsPoissonStrategy"))

  validate_type(c0, type = "double", n = 1, positive = TRUE)
  validate_type(show_progress, type = "logical", n = 1)

  # Handle and validate flow_mean and flow_cov
  if (is.null(flow_mean)) {
    flow_mean <- rep(0, d)
  } else {
    validate_type(flow_mean, type = "double", n = d)
  }

  if (is.null(flow_cov)) {
    flow_cov <- diag(1, nrow = d, ncol = d)
  } else {
    validate_type(flow_cov, type = "double", dims = c(d, d))
    if (!isSymmetric(flow_cov)) {
      cli::cli_abort("Argument 'flow_cov' must be a symmetric matrix.")
    }
  }

  # Validate x0
  if (is.null(x0)) {
    x0 <- rnorm(d)
  } else {
    validate_type(x0, type = "double", n = d)
  }

  # Validate theta0
  if (!is.null(theta0)) {
    validate_type(theta0, type = "double", n = d)
    if (flow == "ZigZag") {
      if (!all(theta0 %in% c(-1, 1))) {
        cli::cli_abort("Argument 'theta0' should only contain -1 or 1 when using ZigZag dynamics.")
      }
    }
  }

  # Validate discretize_dt
  if (!is.null(discretize_dt))
    validate_type(discretize_dt, type = "double", n = 1, positive = TRUE)

  # Validate sticky parameters
  validate_type(sticky, type = "logical", n = 1)
  if (sticky) {
    if (is.null(can_stick)) {
      can_stick <- rep(FALSE, d)
    } else {
      validate_type(can_stick, type = "logical", n = d)
    }
    if (is.null(model_prior) || !(is.bernoulli(model_prior) || is.betabernoulli(model_prior))) {
      cli::cli_abort("Argument 'model_prior' must be provided when 'sticky' is TRUE. It should be an object of class 'bernoulli' or 'beta-bernoulli'.")
    } else if (is.bernoulli(model_prior)) {
      if (length(model_prior$prob) == 1) {
        model_prior <- bernoulli(prob = rep(model_prior$prob, d))
      } else if (length(model_prior$prob) != d) {
        cli::cli_abort("For model_prior 'bernoulli', 'prob' must be a single value or a vector of length d.")
      }
    }

    if (is.null(parameter_prior))
      cli::cli_abort("Argument 'parameter_prior' must be provided when 'sticky' is TRUE. It must be a single value or a vector of length d.")

    validate_type(parameter_prior, type = "double", n = d, positive = TRUE)

  }

  grid_n <- cast_integer(grid_n, n = 1)
  validate_type(grid_n, type = "integer", n = 1, positive = TRUE)
  validate_type(grid_t_max, type = "double", n = 1, positive = TRUE)

  return(list(
    d = d, flow = flow, algorithm = algorithm, T = T, t0 = t0,
    flow_mean = flow_mean, flow_cov = flow_cov, c0 = c0,
    x0 = x0, theta0 = theta0, show_progress = show_progress,
    discretize_dt = discretize_dt,
    sticky = sticky, can_stick = can_stick, model_prior = model_prior, parameter_prior = parameter_prior,
    grid_n = grid_n, grid_t_max = grid_t_max
  ))
}

#' PDMP Sampling
#'
#' Performs Piecewise Deterministic Markov Process (PDMP) sampling using the
#' PDMPSamplers.jl Julia package.
#'
#' @param f Function that computes the NEGATIVE gradient. Should take a numeric vector
#'   and return a numeric vector of the same length.
#' @param d Integer, dimension of the problem.
#' @param flow Character string specifying the flow type. One of "ZigZag",
#'   "BouncyParticle", or "Boomerang".
#' @param algorithm Character string specifying the algorithm. One of
#'   "ThinningStrategy", "GridThinningStrategy", or "RootsPoissonStrategy".
#' @param T Numeric, total sampling time (default: 50000).
#' @param t0 Numeric, initial time (default: 0.0).
#' @param flow_mean Numeric vector of length d, mean vector for the flow (default: zero vector).
#' @param flow_cov Numeric matrix of size d x d, covariance matrix for the flow (default: identity matrix).
#' @param c0 Numeric, bound parameter (default: 1e-2).
#' @param x0 Numeric vector of length d, initial position (default: random normal).
#' @param theta0 Numeric vector of length d, initial velocity (default: random based on the flow).
#' @param hessian Function that returns the negative Hessian matrix of size d x d (default: NULL).
#' @param sticky Logical, whether to use sticky sampling (default: FALSE).
#' @param can_stick Logical vector of length d, which coordinates can stick (default: all FALSE).
#' @param kappa Numeric, stickiness parameter (default: NULL).
#' @param show_progress Logical, whether to show progress bar (default: TRUE).
#' @param discretize_dt Numeric, discretization time step. If NULL, uses mean
#'   inter-event time (default: NULL).
#'
#' @return A list containing:
#'   \item{samples}{Matrix of discretized samples}
#'   \item{trace}{The full trace object}
#'   \item{stats}{Sampling statistics}
#'
#' @export
pdmp_sample <- function(f, d,
                        flow = c("ZigZag", "BouncyParticle", "Boomerang"),
                        algorithm = c("ThinningStrategy", "GridThinningStrategy"), #"RootsPoissonStrategy"),
                        T = 50000, t0 = 0.0,
                        flow_mean = NULL, flow_cov = NULL, c0 = 1e-2,
                        x0 = NULL, theta0 = NULL,
                        hessian = NULL,
                        sticky = FALSE, can_stick = NULL, model_prior = NULL, parameter_prior = NULL,
                        show_progress = TRUE, discretize_dt = NULL) {

  check_for_julia_setup()

  # Validate function argument
  if (!rlang::is_function(f)) {
    cli::cli_abort(c(
      "{.fn f} must be a function",
      "x" = "You've supplied an object of {.cls {class(f)}}."
    ))
  }

  # Use common validation function
  params <- validate_pdmp_params(d, flow, algorithm, T, t0, flow_mean, flow_cov,
                                c0, x0, theta0, show_progress, discretize_dt,
                                sticky, can_stick, model_prior, parameter_prior)

  # Test the function with a sample input
  tryCatch({
    test_output <- f(params$x0)
    if (!is.numeric(test_output) || length(test_output) != params$d) {
      cli::cli_abort(paste0("Function 'f' must return a numeric vector of length ", params$d, "."))
    }
  }, error = function(e) {
    cli::cli_abort(paste0("Function 'f' failed on test input: ", e$message))
  })

  # Test the hessian with a sample input
  if (!is.null(hessian)) {
    tryCatch({
      test_output <- hessian(params$x0)
      if (!is.numeric(test_output) || !is.matrix(test_output) || all(dim(test_output) != c(params$d, params$d))) {
        cli::cli_abort(paste0("Function 'hessian' must return a numeric matrix of dimensions ", params$d, "x", params$d, "."))
      }
    }, error = function(e) {
      cli::cli_abort(paste0("Function 'hessian' failed on test input: ", e$message))
    })
  }

    # Pass arguments to Julia
  for (nm in names(params))
    JuliaCall::julia_assign(nm, params[[nm]])

  # Pass arguments to Julia
  JuliaCall::julia_assign("f",             f)
  JuliaCall::julia_command("grad!(out, x) = out .= f(x)")
  JuliaCall::julia_assign("hessian",       hessian)
  JuliaCall::julia_assign("hvp!",          NULL)

  result <- run_r_interface_function()
  return(result)

}

#' PDMP Sampling from Stan Model
#'
#' Performs Piecewise Deterministic Markov Process (PDMP) sampling from a Stan model
#' using the PDMPSamplers.jl Julia package.
#'
#' @param path_to_stanmodel Character, path to the compiled Stan model file.
#' @param path_to_standata Character, path to the Stan data file (JSON format).
#' @param flow Character string specifying the flow type. One of "ZigZag",
#'   "BouncyParticle", or "Boomerang".
#' @param algorithm Character string specifying the algorithm. One of
#'   "ThinningStrategy", "GridThinningStrategy", or "RootsPoissonStrategy".
#' @param T Numeric, total sampling time (default: 50000).
#' @param t0 Numeric, initial time (default: 0.0).
#' @param flow_mean Numeric vector of length d, mean vector for the flow (default: zero vector).
#' @param flow_cov Numeric matrix of size d x d, covariance matrix for the flow (default: identity matrix).
#' @param c0 Numeric, bound parameter (default: 1e-2).
#' @param x0 Numeric vector of length d, initial position (default: random normal).
#' @param theta0 Numeric vector of length d, initial velocity (default: random based on the flow).
#' @param hessian Function that returns the negative Hessian matrix of size d x d (default: NULL).
#' @param sticky Logical, whether to use sticky sampling (default: FALSE).
#' @param can_stick Logical vector of length d, which coordinates can stick (default: all FALSE).
#' @param model_prior Prior distribution object for model selection. Should be of class 'bernoulli' or 'beta-bernoulli' (default: NULL).
#' @param parameter_prior Numeric vector of length d, prior parameters for sticky sampling (default: NULL).
#' @param grid_n Integer, number of grid points for GridThinningStrategy (default: 30).
#' @param grid_t_max Numeric, maximum time for grid in GridThinningStrategy (default: 2.0).
#' @param show_progress Logical, whether to show progress bar (default: TRUE).
#' @param discretize_dt Numeric, discretization time step. If NULL, uses mean
#'   inter-event time (default: NULL).
#'
#' @return A list containing:
#'   \item{samples}{Matrix of discretized samples}
#'   \item{trace}{The full trace object}
#'   \item{stats}{Sampling statistics}
#'
#' @export
pdmp_sample_from_stanmodel <- function(path_to_stanmodel, path_to_standata,
                        flow = c("ZigZag", "BouncyParticle", "Boomerang"),
                        algorithm = c("ThinningStrategy", "GridThinningStrategy", "RootsPoissonStrategy"),
                        T = 50000, t0 = 0.0,
                        flow_mean = NULL, flow_cov = NULL, c0 = 1e-2,
                        x0 = NULL, theta0 = NULL,
                        hessian = NULL,
                        sticky = FALSE, can_stick = NULL, model_prior = NULL, parameter_prior = NULL,
                        grid_n = 30, grid_t_max = 2.0,
                        show_progress = TRUE, discretize_dt = NULL) {

  check_for_julia_setup()

  # Setup Stan model and get dimension
  JuliaCall::julia_assign("path_to_stan_model", path_to_stanmodel)
  JuliaCall::julia_assign("path_to_stan_data",  path_to_standata)
  JuliaCall::julia_command("grad!, hvp!, d, sm = grad_and_hvp_from_stanmodel(path_to_stan_model, path_to_stan_data);")
  JuliaCall::julia_command("hessian = nothing")
  d <- JuliaCall::julia_eval("d")

  # Use common validation function
  params <- validate_pdmp_params(d, flow, algorithm, T, t0, flow_mean, flow_cov,
                                 c0, x0, theta0, show_progress, discretize_dt,
                                 sticky, can_stick, model_prior, parameter_prior, grid_n, grid_t_max)

  # Pass arguments to Julia
  for (nm in names(params)) {
    JuliaCall::julia_assign(nm, params[[nm]])
  }

  result <- run_r_interface_function()

  return(result)
}

run_r_interface_function <- function() {
  JuliaCall::julia_eval("r_interface_function(
    grad!, d, x0, flow, algorithm, flow_mean, flow_cov, nothing,
    c0, grid_n, grid_t_max,
    t0, T, discretize_dt, hessian, hvp!,
    sticky, can_stick, model_prior, parameter_prior,
    show_progress
  )")
}

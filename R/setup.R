get_julia_project_dir <- function() {
  cache_dir <- fs::path_expand(rappdirs::user_cache_dir("PDMPSamplersR"))
  project_dir <- fs::path(cache_dir, "_juliaProject")
  return(project_dir)
}

julia_project_exists <- function() {
  project_dir <- get_julia_project_dir()
  return(fs::dir_exists(project_dir) && fs::file_exists(fs::path(project_dir, "Project.toml")))
}

setup_julia_project <- function() {
  project_dir <- get_julia_project_dir()
  if (!fs::dir_exists(project_dir)) fs::dir_create(project_dir, recurse = TRUE)

  JuliaCall::julia_command("import Pkg")
  JuliaCall::julia_command(sprintf('Pkg.activate("%s")', project_dir))
  pmdpsamplersjl_location <- getOption("PDMPSamplersR.pmdpsamplersjl_location", NULL)
  if (!is.null(pmdpsamplersjl_location)) {
    JuliaCall::julia_command(sprintf('Pkg.develop(path="%s")', pmdpsamplersjl_location))
  } else {
    pmdpsamplersjl_url <- "https://github.com/vandenman/PDMPSamplers.jl"
    JuliaCall::julia_command(sprintf('Pkg.add(url="%s")', pmdpsamplersjl_url))
  }
  # BridgeStan is needed to trigger the PDMPSamplers BridgeStan extension
  JuliaCall::julia_command('Pkg.add("BridgeStan")')
  JuliaCall::julia_command('Pkg.update()')
  JuliaCall::julia_command('Pkg.precompile()')

}

load_julia_project <- function() {
  project_dir <- get_julia_project_dir()
  JuliaCall::julia_command(sprintf('Base.active_project() == joinpath("%1$s", "Project.toml") || Pkg.activate("%1$s");', project_dir))
}

#' Update the Julia project for PDMPSamplers.jl
#'
#' Updates the Julia package dependencies used by PDMPSamplersR. If the
#' Julia project does not exist yet, it is created automatically.
#'
#' @returns Invisible \code{NULL}. Called for its side effects.
#' @seealso \code{\link{pdmpsamplers_setup}} to perform the initial setup.
#' @export
pdmpsamplers_update <- function() {
  project_dir <- get_julia_project_dir()
  if (fs::dir_exists(project_dir) && fs::file_exists(fs::path(project_dir, "Project.toml"))) {
    JuliaCall::julia_command(sprintf('Pkg.activate("%s")', project_dir))
    JuliaCall::julia_command('Pkg.update()')
    JuliaCall::julia_command('Pkg.precompile()')
  } else {
    cli::cli_inform("Julia project not found. Setting up a new project.")
    setup_julia_project()
  }
}

#' Use a local version of PDMPSamplers.jl
#'
#' Switch the Julia project to use a local (development) copy of
#' PDMPSamplers.jl instead of the GitHub version, or switch back to the
#' GitHub version.
#'
#' @param path Character string giving the path to the local PDMPSamplers.jl
#'   directory. If \code{NULL} (default), switches back to the GitHub version.
#'
#' @returns Invisible \code{NULL}. Called for its side effects.
#' @export
use_local_pdmpsamplers <- function(path = NULL) {
  project_dir <- get_julia_project_dir()
  if (!julia_project_exists()) {
    cli::cli_abort("Julia project not found! Run {.fn setup_julia_project} first.")
  }
  JuliaCall::julia_command("import Pkg")
  JuliaCall::julia_command(sprintf('Pkg.activate("%s")', project_dir))
  if (!is.null(path)) {
    path <- normalizePath(path, mustWork = TRUE)
    cli::cli_alert_info("Switching to local PDMPSamplers.jl at {.path {path}}")
    JuliaCall::julia_command(sprintf('Pkg.develop(path="%s")', path))
  } else {
    cli::cli_alert_info("Switching to GitHub version of PDMPSamplers.jl")
    pmdpsamplersjl_url <- "https://github.com/vandenman/PDMPSamplers.jl"
    JuliaCall::julia_command('Pkg.rm("PDMPSamplers")')
    JuliaCall::julia_command(sprintf('Pkg.add(url="%s")', pmdpsamplersjl_url))
  }
  JuliaCall::julia_command('Pkg.precompile()')
  invisible(NULL)
}

load_interface_function <- function() {
  interface_file <- system.file("julia", "main_interface_function.jl", package = "PDMPSamplersR")
  if (fs::file_exists(interface_file)) {
    JuliaCall::julia_source(interface_file)
  } else {
    cli::cli_abort("Interface file not found! Reinstall the package and try again.")
  }
}

#' Check if the Julia Project is Setup Properly
#'
#' @param error_if_not_exists logical, if TRUE, will throw an error if the Julia project is not found. If FALSE, will attempt to set up the project if setup_if_not_exists is TRUE.
#' @param setup_if_not_exists logical, if FALSE and `error_if_not_exists` is FALSE, will not throw an error but printif the Julia project is not found. If TRUE, will attempt to set up the project if `error_if_not_exists` is FALSE.
#'
#' @returns TRUE if the Julia project is successfully loaded or set up, \code{FALSE} otherwise (if `error_if_not_exists` is also \code{FALSE}). If `error_if_not_exists` is \code{TRUE}, will throw an error if the project cannot be loaded or set up.
#' @details This function checks for the existence of the Julia project for PDMPSamplers.jl, and if it exists, loads it and the required packages. If it does not exist, it can either throw an error or attempt to set up the project, depending on the parameters. It is mostly useful for internal use to ensure the Julia environment is ready before calling any functions that depend on it, but can also be called directly by users if they want to check or set up the Julia project.
#' @export
check_for_julia_setup <- function(error_if_not_exists = FALSE, setup_if_not_exists = TRUE) {

  # TODO: some check on Julia installation?
  if (julia_project_exists()) {
    load_julia_project()
    load_required_julia_packages()
    load_interface_function()
    return(TRUE)
  } else if (!error_if_not_exists) {
    if (setup_if_not_exists) {
      cli::cli_alert_info("Setting up Julia project for PDMPSamplers.jl...")
      setup_julia_project()
      check_for_julia_setup(error_if_not_exists = TRUE)
      return(TRUE)
    } else {
      cli::cli_abort("Julia project for PDMPSamplers.jl not found! Run PDMPSamplersR::setup_julia_project() to create it.")
      return(FALSE)
    }
  } else {
    cli::cli_abort("Setting up Julia project for PDMPSamplers.jl failed!")
  }
}

load_required_julia_packages <- function() {
  # TODO: get this from github!
  # JuliaCall::julia_install_package_if_needed("PDMPSamplers")
  JuliaCall::julia_install_package_if_needed("LinearAlgebra")
  JuliaCall::julia_install_package_if_needed("BridgeStan")
}

#' Set up the Julia environment for PDMPSamplers.jl
#'
#' Initializes the Julia runtime via \pkg{JuliaCall} and sets up the
#' Julia project that PDMPSamplersR depends on. This installs all
#' required Julia packages and precompiles them. Call this once before
#' using any sampling functions.
#'
#' @param verbose Logical, if \code{TRUE} (default), progress messages are
#'   printed.
#'
#' @returns Invisible \code{TRUE} if the setup succeeds.
#' @seealso \code{\link{pdmpsamplers_update}} to update Julia dependencies.
#' @export
pdmpsamplers_setup <- function(verbose = TRUE) {
  if (verbose) cli::cli_inform("Setting up Julia")
  JuliaCall::julia_setup(verbose = verbose)#, install=FALSE, useRCall=FALSE)
  if (verbose) cli::cli_inform("Setting up Julia project for PDMPSamplers.jl")
  check_for_julia_setup()
}

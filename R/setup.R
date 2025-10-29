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

  # pmdpsamplersjl_location <- "/home/don/hdd/surfdrive/Postdoc/PDMPSamplers.jl"
    JuliaCall::julia_command(sprintf('Pkg.develop(path="%s")', pmdpsamplersjl_location))
  } else {
    pmdpsamplersjl_url <- "https://github.com/vandenman/PDMPSamplers.jl"
    JuliaCall::julia_command(sprintf('Pkg.add(url="%s")', pmdpsamplersjl_url))
  }
  JuliaCall::julia_command('Pkg.update()')
  JuliaCall::julia_command('Pkg.precompile()')

}

load_julia_project <- function() {
  project_dir <- get_julia_project_dir()
  JuliaCall::julia_command(sprintf('Base.active_project() == joinpath("%1$s", "Project.toml") || Pkg.activate("%1$s")', project_dir))
}

#'@export
update_julia_project <- function() {
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

load_interface_function <- function() {
  interface_file <- system.file("julia", "main_interface_function.jl", package = "PDMPSamplersR")
  if (fs::file_exists(interface_file)) {
    JuliaCall::julia_source(interface_file)
  } else {
    cli::cli_abort("Interface file not found! Reinstall the package and try again.")
  }
}

check_for_julia_setup <- function(error_if_not_exists = FALSE) {

  # TODO: some check on Julia installation?
  if (julia_project_exists()) {
    load_julia_project()
    load_required_julia_packages()
    load_interface_function()
  } else if (!error_if_not_exists) {
    cli::cli_alert_info("Setting up Julia project for PDMPSamplers.jl...")
    setup_julia_project()
    check_for_julia_setup(error_if_not_exists = TRUE)
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

#'@export
setup_cmdstanr <- function() {
  if (!require("cmdstanr"))
    install.packages("cmdstanr", repos = c('https://stan-dev.r-universe.dev', getOption("repos")))
  cmdstanr::check_cmdstan_toolchain()
  cmdstanr::install_cmdstan(cores = 2) # TODO: make cores a parameter with default
}

#'@export
setup <- function(verbose = TRUE) {
  if (verbose) cli::cli_inform("Setting up Julia")
  JuliaCall::julia_setup(verbose = verbose)#, install=FALSE, useRCall=FALSE)
  if (verbose) cli::cli_inform("Setting up Julia project for PDMPSampler.jl")
  check_for_julia_setup()
  if (verbose) cli::cli_inform("Setting up cmdstanr")
  setup_cmdstanr()
}

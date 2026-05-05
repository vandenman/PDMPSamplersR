get_julia_project_dir <- function() {
  cache_dir <- fs::path_expand(rappdirs::user_cache_dir("PDMPSamplersR"))
  project_dir <- fs::path(cache_dir, "_juliaProject")
  return(project_dir)
}

julia_project_exists <- function() {
  project_dir <- get_julia_project_dir()
  return(fs::dir_exists(project_dir) && fs::file_exists(fs::path(project_dir, "Project.toml")))
}

ensure_julia_runtime <- function(verbose = FALSE) {
  if (verbose) {
    JuliaCall::julia_setup(verbose = TRUE)
  } else {
    suppressMessages(JuliaCall::julia_setup(verbose = FALSE))
  }

  invisible(TRUE)
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

  # Verify that PDMPSamplers actually loads; if precompilation silently failed
  # (e.g., OOM-killed worker on CI), `using` will trigger precompilation again
  # and surface the real error.
  tryCatch(
    JuliaCall::julia_command("using PDMPSamplers"),
    error = function(e) {
      cli::cli_warn(c(
        "!" = "PDMPSamplers failed to load after precompilation.",
        "i" = "This is sometimes caused by insufficient memory during precompilation on CI.",
        "i" = "Retrying precompilation with verbose output..."
      ))
      # Show the precompile log if available
      tryCatch({
        JuliaCall::julia_command('
          for (uuid, info) in Pkg.Types.Context().env.manifest.deps
            if info.name == "PDMPSamplers"
              compiled_dir = joinpath(first(DEPOT_PATH), "compiled", "v$(VERSION.major).$(VERSION.minor)", "PDMPSamplers")
              if isdir(compiled_dir)
                for f in readdir(compiled_dir)
                  if endswith(f, ".log")
                    logpath = joinpath(compiled_dir, f)
                    println("=== Precompile log: $logpath ===")
                    println(read(logpath, String))
                  end
                end
              end
            end
          end
        ')
      }, error = function(e2) NULL)

      # Retry once
      tryCatch({
        JuliaCall::julia_command('Pkg.precompile()')
        JuliaCall::julia_command("using PDMPSamplers")
      }, error = function(e2) {
        cli::cli_abort(c(
          "PDMPSamplers.jl failed to precompile.",
          "x" = conditionMessage(e),
          "i" = "Check the Julia precompile logs above for details."
        ))
      })
    }
  )

}

load_julia_project <- function() {
  project_dir <- get_julia_project_dir()
  JuliaCall::julia_command(sprintf('Base.active_project() == joinpath("%1$s", "Project.toml") || Pkg.activate("%1$s");', project_dir))
}

.pdmp_env <- new.env(parent = emptyenv())
.pdmp_env$version_printed <- FALSE

print_pdmpsamplers_version <- function() {
  if (.pdmp_env$version_printed) return(invisible(NULL))
  version_str <- tryCatch(
    JuliaCall::julia_eval('
      let deps = Pkg.dependencies()
        entry = [dep for (_, dep) in deps if dep.name == "PDMPSamplers"]
        if isempty(entry)
          "unknown"
        else
          d = first(entry)
          d.is_tracking_path || d.is_tracking_repo ? "$(d.version) (dev)" : string(d.version)
        end
      end
    '),
    error = function(e) "unknown"
  )
  cli::cli_inform("Using PDMPSamplers.jl v{version_str}")
  .pdmp_env$version_printed <- TRUE
  invisible(NULL)
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

#' Show the status of the PDMPSamplers.jl installation
#'
#' Displays the current version of PDMPSamplers.jl and whether it is installed
#' from a local path or from GitHub.
#'
#' @returns Invisible \code{NULL}. Called for its side effects.
#' @seealso \code{\link{pdmpsamplers_setup}} to perform the initial setup,
#'   \code{\link{pdmpsamplers_update}} to update Julia dependencies.
#' @export
pdmpsamplers_status <- function() {
  if (!julia_project_exists()) {
    cli::cli_alert_warning("Julia project not found. Run {.fn pdmpsamplers_setup} first.")
    return(invisible(NULL))
  }
  load_julia_project()
  info <- tryCatch(
    JuliaCall::julia_eval('
      let deps = Pkg.dependencies()
        entries = [dep for (_, dep) in deps if dep.name == "PDMPSamplers"]
        if isempty(entries)
          "not installed|||none|||"
        else
          d = first(entries)
          ver = d.version === nothing ? "unknown" : string(d.version)
          git_src = isnothing(d.git_source) ? "" : d.git_source
          source_type = d.is_tracking_path ? "local" : (d.is_tracking_repo ? "GitHub" : "registry")
          join([ver, source_type, git_src], "|||")
        end
      end
    '),
    error = function(e) "unknown|||unknown|||"
  )
  parts <- strsplit(info, "|||", fixed = TRUE)[[1]]
  version     <- parts[1]
  source_type <- if (length(parts) >= 2) parts[2] else "unknown"
  source_path <- if (length(parts) >= 3) parts[3] else ""

  source_label <- switch(source_type,
    "local"    = paste0("local (", source_path, ")"),
    "GitHub"   = "GitHub",
    "registry" = "registry",
    source_type
  )

  cli::cli_inform(c(
    "PDMPSamplers.jl status:",
    "*" = "Version: {version}",
    "*" = "Source:  {source_label}"
  ))
  invisible(NULL)
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
  if (isTRUE(JuliaCall::julia_eval("isdefined(Main, :PDMPSamplersRBridge)"))) return(invisible(NULL))
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
  ensure_julia_runtime()

  if (julia_project_exists()) {
    load_julia_project()
    load_required_julia_packages()
    load_interface_function()
    print_pdmpsamplers_version()
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
  ensure_julia_runtime(verbose = verbose)
  if (verbose) cli::cli_inform("Setting up Julia project for PDMPSamplers.jl")
  check_for_julia_setup()
}

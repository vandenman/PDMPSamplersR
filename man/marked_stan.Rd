\name{stan_marked_subsampling}
\alias{stan_marked_subsampling}
\alias{stan_residual_envelope}
\alias{omrf_residual_envelope}
\alias{pdmp_subsample_hpp_path}
\alias{compile_pdmp_stan_model}
\alias{stan_parameter_mapping}
\title{Exact marked custom-Stan sampling interfaces}
\usage{
stan_residual_envelope(weights, growth_rates = NULL)

stan_marked_subsampling(n_observations, subsample_size, prior_standata,
  residual_envelope, anchor = NULL, n_anchor_updates = 0L,
  use_anchor_bank = FALSE, bank_capacity = 1L)

omrf_residual_envelope(X, seen, thresholds, interactions,
  backend = c("stan_header", "analytic"))

pdmp_subsample_hpp_path()

compile_pdmp_stan_model(path)

stan_parameter_mapping(path_to_stanmodel, standata, parameters = NULL)
}
\arguments{
\item{weights}{Finite nonnegative component-by-observation envelope weights.}
\item{growth_rates}{Optional nonnegative growth rate per component.}
\item{n_observations}{Total number of marked observations.}
\item{subsample_size}{Number of fresh observations per marked proposal.}
\item{prior_standata}{Prior-only data for the same compiled Stan model.}
\item{residual_envelope}{A certified package residual-envelope object.}
\item{anchor}{Optional unconstrained control-variate anchor.}
\item{n_anchor_updates}{Number of requested anchor updates.}
\item{use_anchor_bank}{Whether to use an analytic provider's anchor bank.}
\item{bank_capacity}{Positive anchor-bank capacity.}
\item{X}{Integer person-by-node OMRF response matrix encoded from zero.}
\item{seen}{Number of observed categories per OMRF node.}
\item{thresholds}{Stan threshold block name or prefix.}
\item{interactions}{Stan interaction block name or prefix.}
\item{backend}{Selected-gradient backend.}
\item{path}{Stan source path to compile with the package subset hook.}
\item{path_to_stanmodel}{Stan source or compiled-library path.}
\item{standata}{Named Stan data list or JSON path.}
\item{parameters}{Optional names, block prefixes, or 1-based indices.}
}
\value{A marked/envelope specification, installed header or library path, or a parameter mapping.}
\description{
Public interfaces for persistent BridgeStan selected-gradient contexts and
certified exact marked thinning. Julia owns minibatch scaling and draws a fresh
subset at every marked proposal.
}

# Persistent selected-gradient backend for public custom-Stan subsampling sampling.

mutable struct StanSubsamplingCallCounts
    full_gradient::Int
    prior_gradient::Int
    selected_gradient::Int
    persons_evaluated::Int
    anchor_cache_gradient::Int
    anchor_cache_persons::Int
    analytic_residual::Int
    analytic_persons::Int
    analytic_factors::Int
    analytic_node_conditionals::Int
    model_constructions::Int
    data_constructions::Int
    native_full_data_scans::Int
    event_bound_observation_visits::Int
    native_full_gradient_seconds::Float64
    analytic_residual_seconds::Float64
    envelope_point_calls::Int
    envelope_point_seconds::Float64
    envelope_cell_calls::Int
    envelope_cell_seconds::Float64
    proposal_deterministic_bound_sum::Float64
    proposal_aggregate_residual_bound_sum::Float64
    proposal_subset_bound_sum::Float64
    proposal_deterministic_actual_sum::Float64
    proposal_exact_residual_rate_sum::Float64
    proposal_actual_rate_sum::Float64
    accepted_proposal_bound_sum::Float64
    accepted_actual_rate_sum::Float64
    proposal_count::Int
    proposal_acceptance_probability_sum::Float64
    proposal_certificate_ratio_sum::Float64
    proposal_intrinsic_ratio_sum::Float64
    lazy_branch_evaluations::Int
    lazy_branch_rejections::Int
    lazy_leaf_evaluations::Int
    lazy_sampling_seconds::Float64
    cell_cached_point_screen_seconds::Float64
end
StanSubsamplingCallCounts() = StanSubsamplingCallCounts(
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0.0, 0.0, 0, 0.0, 0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0, 0.0, 0.0, 0.0,
    0, 0, 0, 0.0, 0.0)

mutable struct StanSubsamplingContext
    full::PDMPModel
    prior::PDMPModel
    selected::PDMPModel
    set_fn::Ptr{Nothing}
    clear_fn::Ptr{Nothing}
    subset0::Vector{Int32}
    m::Int
    counts::StanSubsamplingCallCounts
end

function _resolve_subset_hooks(sm::BridgeStan.StanModel)
    set_fn = Libdl.dlsym(sm.lib, :pdmp_set_subsample_indices; throw_error=false)
    clear_fn = Libdl.dlsym(sm.lib, :pdmp_clear_subsample_indices; throw_error=false)
    (set_fn == C_NULL || clear_fn == C_NULL) && throw(ArgumentError(
        "the compiled Stan model does not export the PDMPSamplersR subset hook; compile it with compile_pdmp_stan_model()"))
    return set_fn, clear_fn
end

@inline function _clear_stan_subset!(ctx::StanSubsamplingContext)
    @ccall $(ctx.clear_fn)()::Cvoid
    return nothing
end

@inline function _install_stan_subset!(ctx::StanSubsamplingContext, subset;
        require_configured_size::Bool=true)
    require_configured_size && length(subset) != ctx.m &&
        throw(DimensionMismatch("selected subset has the wrong size"))
    length(subset) <= length(ctx.subset0) || throw(DimensionMismatch(
        "selected subset exceeds the configured buffer size"))
    @inbounds for j in eachindex(subset)
        ctx.subset0[j] = Int32(subset[j] - 1)
    end
    subset0 = ctx.subset0
    set_fn = ctx.set_fn
    m = Cint(length(subset))
    GC.@preserve subset0 begin
        @ccall $set_fn(pointer(subset0)::Ptr{Int32}, m::Cint)::Cvoid
    end
    return nothing
end

function _stan_subset_gradient!(ctx::StanSubsamplingContext, out, x, subset;
        require_configured_size::Bool=true)
    lock(_BRIDGESTAN_CALL_LOCK) do
        _install_stan_subset!(ctx, subset; require_configured_size)
        try
            compute_gradient!(ctx.selected.grad, x, out)
        finally
            _clear_stan_subset!(ctx)
        end
    end
    return out
end

function _clear_gradient!(ctx::StanSubsamplingContext, which::Symbol, out, x)
    which in (:full, :prior) || throw(ArgumentError(
        "gradient selector must be :full or :prior"))
    lock(_BRIDGESTAN_CALL_LOCK) do
        _clear_stan_subset!(ctx)
        model = which === :full ? ctx.full : ctx.prior
        compute_gradient!(model.grad, x, out)
    end
    if which === :full
        ctx.counts.full_gradient += 1
    else
        ctx.counts.prior_gradient += 1
    end
    return out
end

function _selected_gradient!(ctx::StanSubsamplingContext, out, x, subset)
    _stan_subset_gradient!(ctx, out, x, subset)
    ctx.counts.selected_gradient += 1
    ctx.counts.persons_evaluated += ctx.m
    return out
end


function _anchor_observation_gradient!(ctx::StanSubsamplingContext, out, x, i::Int)
    _stan_subset_gradient!(ctx, out, x, (i,); require_configured_size=false)
    ctx.counts.anchor_cache_gradient += 1
    ctx.counts.anchor_cache_persons += 1
    return out
end

mutable struct StanSubsamplingOracle
    ctx::StanSubsamplingContext
    full_anchor::Vector{Float64}
    prior_anchor::Vector{Float64}
    prior_x::Vector{Float64}
    prior_x_position::Vector{Float64}
    prior_x_valid::Bool
    selected_x::Vector{Float64}
    selected_anchor::Vector{Float64}
    anchor_likelihood_cache::Matrix{Float64}
end

function StanSubsamplingOracle(ctx::StanSubsamplingContext,
        anchor::Vector{Float64}, N::Int)
    d = length(anchor)
    oracle = StanSubsamplingOracle(ctx, zeros(d), zeros(d), zeros(d),
        zeros(d), false, zeros(d), zeros(d), Matrix{Float64}(undef, d, N))
    _clear_gradient!(ctx, :full, oracle.full_anchor, anchor)
    _clear_gradient!(ctx, :prior, oracle.prior_anchor, anchor)
    for i in 1:N
        _anchor_observation_gradient!(ctx, oracle.selected_anchor, anchor, i)
        @views @. oracle.anchor_likelihood_cache[:, i] =
            oracle.selected_anchor - oracle.prior_anchor
    end
    return oracle
end

function _prior_gradient_at_x!(oracle::StanSubsamplingOracle, x)
    if !oracle.prior_x_valid || oracle.prior_x_position != x
        _clear_gradient!(oracle.ctx, :prior, oracle.prior_x, x)
        copyto!(oracle.prior_x_position, x)
        oracle.prior_x_valid = true
    end
    return oracle.prior_x
end

function (oracle::StanSubsamplingOracle)(out, x)
    _prior_gradient_at_x!(oracle, x)
    @. out = oracle.full_anchor + oracle.prior_x - oracle.prior_anchor
    return out
end

function (oracle::StanSubsamplingOracle)(out, x, subset, anchor)
    _selected_gradient!(oracle.ctx, oracle.selected_x, x, subset)
    copyto!(oracle.selected_anchor, oracle.prior_anchor)
    @inbounds for i in subset, j in eachindex(oracle.selected_anchor)
        oracle.selected_anchor[j] += oracle.anchor_likelihood_cache[j, i]
    end
    _prior_gradient_at_x!(oracle, x)
    @. out = (oracle.selected_x - oracle.prior_x) -
        (oracle.selected_anchor - oracle.prior_anchor)
    return out
end

"""
Workspace for the exact OMRF person-likelihood residual.

The threshold and interaction coordinates refer to the full unconstrained
state. Only those coordinates are written by `omrf_residual!`; all nuisance
and prior-only coordinates remain zero. Thresholds are ordered node-major and
interactions use the upper-triangular order `(1,2), (1,3), ..., (P-1,P)`.
"""
mutable struct OMRFResidualContext
    X::Matrix{Int}
    X_float::Matrix{Float64}
    seen::Vector{Int}
    threshold_indices::Vector{Int}
    interaction_indices::Vector{Int}
    global_scale_index::Int
    node_scale_indices::Vector{Int}
    threshold_starts::Vector{Int}
    incident_edges::Vector{Vector{Int}}
    incident_neighbours::Vector{Vector{Int}}
    anchor_probabilities::Matrix{Float64}
    anchor_expected_scores::Matrix{Float64}
    expected_scores_workspace::Matrix{Float64}
    interaction_matrix::Matrix{Float64}
    linear_fields::Matrix{Float64}
    score_crossproducts::Matrix{Float64}
    observed_threshold_counts::Vector{Float64}
    observed_interaction_sums::Vector{Float64}
    probability_buffer::Vector{Float64}
    hcv_buffer::Vector{Float64}
    likelihood_hessian::Union{Nothing,Matrix{Float64}}
    use_hcv::Bool
    hcv_active::Bool
    hcv_damping::Float64
    displacement::Vector{Float64}
    hcv_product::Vector{Float64}
    factorization::Symbol
    touched_positions::Vector{Int}
    touched_mask::BitVector
    counts::StanSubsamplingCallCounts
end

struct OMRFIndependentAnalyticPrior
    threshold_indices::Vector{Int}
    interaction_indices::Vector{Int}
    threshold_alpha::Float64
    threshold_beta::Float64
    interaction_type::Symbol
    interaction_scale::Float64
end

"""Analytic prior for the node-shared OMRF slab.

The production Stan model declares positive global and node scales.  BridgeStan
therefore represents them on log coordinates, and the `-log(scale)` Jacobian
terms are part of the unconstrained negative-gradient target.
"""
mutable struct OMRFNodeScaleAnalyticPrior
    threshold_indices::Vector{Int}
    interaction_indices::Vector{Int}
    global_scale_index::Int
    node_scale_indices::Vector{Int}
    threshold_alpha::Float64
    threshold_beta::Float64
    interaction_scale::Float64
    global_scale_sd::Float64
    node_scale_sd::Float64
    node_sums::Vector{Float64}
    node_scales::Vector{Float64}
end

const OMRFAnalyticPrior = Union{OMRFIndependentAnalyticPrior,
    OMRFNodeScaleAnalyticPrior}

function _omrf_analytic_prior(envelope_spec, context::OMRFResidualContext)
    _haskey(envelope_spec, :analytic_prior) || return nothing
    spec = _rget(envelope_spec, :analytic_prior)
    α = Float64(_rget(spec, :threshold_alpha))
    β = Float64(_rget(spec, :threshold_beta))
    scale = Float64(_rget(spec, :interaction_scale))
    interaction_type = Symbol(String(_rget(spec, :type)))
    all(isfinite, (α, β, scale)) && α > 0 && β > 0 && scale > 0 ||
        throw(ArgumentError("analytic OMRF prior parameters must be finite and positive"))
    if interaction_type in (:cauchy, :gaussian)
        return OMRFIndependentAnalyticPrior(
            context.threshold_indices, context.interaction_indices,
            α, β, interaction_type, scale)
    elseif interaction_type === :node_scale_gaussian
        global_scale_sd = Float64(_rget(spec, :global_scale_sd))
        node_scale_sd = Float64(_rget(spec, :node_scale_sd))
        all(isfinite, (global_scale_sd, node_scale_sd)) &&
            global_scale_sd > 0 && node_scale_sd > 0 ||
            throw(ArgumentError(
                "node-scale OMRF prior parameters must be finite and positive"))
        return OMRFNodeScaleAnalyticPrior(
            context.threshold_indices, context.interaction_indices,
            context.global_scale_index, context.node_scale_indices,
            α, β, scale, global_scale_sd, node_scale_sd,
            zeros(length(context.node_scale_indices)),
            zeros(length(context.node_scale_indices)))
    end
    throw(ArgumentError("unsupported analytic OMRF interaction prior: $interaction_type"))
end

@inline _omrf_logistic(x::Float64) = LogExpFunctions.logistic(x)

function _omrf_prior_gradient!(out, prior::OMRFIndependentAnalyticPrior, x)
    fill!(out, 0.0)
    α, β = prior.threshold_alpha, prior.threshold_beta
    @inbounds for j in prior.threshold_indices
        out[j] = (α + β) * _omrf_logistic(x[j]) - α
    end
    scale2 = prior.interaction_scale^2
    if prior.interaction_type === :cauchy
        @inbounds for j in prior.interaction_indices
            out[j] = 2x[j] / (scale2 + x[j]^2)
        end
    else
        @inbounds for j in prior.interaction_indices
            out[j] = x[j] / scale2
        end
    end
    return out
end

function _omrf_prior_gradient!(out, prior::OMRFNodeScaleAnalyticPrior, x)
    fill!(out, 0.0)
    α, β = prior.threshold_alpha, prior.threshold_beta
    @inbounds for j in prior.threshold_indices
        out[j] = (α + β) * _omrf_logistic(x[j]) - α
    end

    global_index = prior.global_scale_index
    global_scale = exp(x[global_index])
    global_sum = global_scale^2 / prior.global_scale_sd^2 - 1.0
    node_sums = prior.node_sums
    node_scales = prior.node_scales
    @inbounds for node in eachindex(node_scales)
        node_scales[node] = exp(x[prior.node_scale_indices[node]])
        node_sums[node] = 0.0
    end
    edge = 1
    @inbounds for left in 1:(length(node_scales) - 1), right in (left + 1):length(node_scales)
        interaction = x[prior.interaction_indices[edge]]
        scale2 = prior.interaction_scale^2 * global_scale^2 *
            node_scales[left] * node_scales[right]
        ratio = interaction^2 / scale2
        out[prior.interaction_indices[edge]] = interaction / scale2
        global_sum += 1.0 - ratio
        node_sums[left] += 0.5 * (1.0 - ratio)
        node_sums[right] += 0.5 * (1.0 - ratio)
        edge += 1
    end
    out[global_index] = global_sum
    @inbounds for node in eachindex(node_sums)
        scale = node_scales[node]
        out[prior.node_scale_indices[node]] = node_sums[node] +
            scale^2 / prior.node_scale_sd^2 - 1.0
    end
    return out
end

"""
Allocation-free negative gradient for a complete independent-prior OMRF.

Unlike `OMRFAnalyticSubsamplingOracle`, this provider evaluates every
person-node conditional and therefore contains no factorized event dynamics.
The threshold and interaction blocks must cover the complete unconstrained
Stan state.  Its sign follows `PDMPModel`: it returns the gradient of the
negative log posterior.
"""
mutable struct OMRFFullGradientProvider <: Function
    context::OMRFResidualContext
    prior::OMRFAnalyticPrior
end

function Base.copy(provider::OMRFFullGradientProvider)
    return deepcopy(provider)
end

@inline function _omrf_fill_interaction_matrix!(context::OMRFResidualContext, x)
    matrix = context.interaction_matrix
    fill!(matrix, 0.0)
    edge = 1
    P = size(matrix, 1)
    @inbounds for left in 1:(P - 1), right in (left + 1):P
        value = x[context.interaction_indices[edge]]
        matrix[left, right] = value
        matrix[right, left] = value
        edge += 1
    end
    return matrix
end

"""
Add the complete OMRF pseudo-likelihood gradient to `out`.

The observed-category numerator is a fixed sufficient statistic prepared in
`_omrf_residual_context`. Retained evaluations therefore perform two dense
matrix products (`X * beta` and `X' * E[X | X_-j]`) and one allocation-free
categorical-probability pass; they never revisit the numerator person by
person. When `cache_anchor` is true the same pass also refreshes the anchor
probability tables used by the residual oracle.
"""
function _omrf_add_likelihood_gradient!(out, context::OMRFResidualContext, x;
        cache_anchor::Bool=false)
    X = context.X_float
    N, P = size(X)
    _omrf_fill_interaction_matrix!(context, x)
    mul!(context.linear_fields, X, context.interaction_matrix)
    expected_scores = context.expected_scores_workspace

    @inbounds for node in 1:P
        q = context.seen[node] - 1
        start = context.threshold_starts[node]
        for n in 1:N
            field = context.linear_fields[n, node]
            max_eta = 0.0
            for u in 1:q
                eta = x[context.threshold_indices[start + u - 1]] + u * field
                context.probability_buffer[u] = eta
                max_eta = max(max_eta, eta)
            end
            denominator = exp(-max_eta)
            for u in 1:q
                probability = exp(context.probability_buffer[u] - max_eta)
                context.probability_buffer[u] = probability
                denominator += probability
            end
            expected = 0.0
            for u in 1:q
                probability = context.probability_buffer[u] / denominator
                out[context.threshold_indices[start + u - 1]] += probability
                expected += u * probability
                if cache_anchor
                    context.anchor_probabilities[start + u - 1, n] = probability
                end
            end
            expected_scores[n, node] = expected
            cache_anchor && (context.anchor_expected_scores[node, n] = expected)
        end
        for u in 1:q
            out[context.threshold_indices[start + u - 1]] -=
                context.observed_threshold_counts[start + u - 1]
        end
    end

    mul!(context.score_crossproducts, transpose(X), expected_scores)
    edge = 1
    @inbounds for left in 1:(P - 1), right in (left + 1):P
        out[context.interaction_indices[edge]] +=
            context.score_crossproducts[right, left] +
            context.score_crossproducts[left, right] -
            context.observed_interaction_sums[edge]
        edge += 1
    end
    return out
end

function (provider::OMRFFullGradientProvider)(out, x)
    started = time_ns()
    context = provider.context
    _omrf_prior_gradient!(out, provider.prior, x)
    _omrf_add_likelihood_gradient!(out, context, x)
    context.counts.native_full_data_scans += 1
    context.counts.native_full_gradient_seconds +=
        (time_ns() - started) * 1e-9
    return out
end

"""Refresh an OMRF anchor cache and its exact full gradient in one pass."""
function _cache_omrf_anchor_and_full_gradient!(context::OMRFResidualContext,
        full, prior::OMRFAnalyticPrior, anchor)
    _omrf_prior_gradient!(full, prior, anchor)
    _omrf_add_likelihood_gradient!(full, context, anchor; cache_anchor=true)
    context.counts.native_full_data_scans += 1
    return full
end

"""Exact Hessian-vector product paired with `OMRFFullGradientProvider`."""
mutable struct OMRFFullHVPProvider <: Function
    context::OMRFResidualContext
    prior::OMRFAnalyticPrior
end

function Base.copy(provider::OMRFFullHVPProvider)
    return deepcopy(provider)
end

function _omrf_prior_hvp!(out, prior::OMRFIndependentAnalyticPrior, x, direction)
    fill!(out, 0.0)
    αβ = prior.threshold_alpha + prior.threshold_beta
    @inbounds for j in prior.threshold_indices
        p = _omrf_logistic(x[j])
        out[j] = αβ * p * (1 - p) * direction[j]
    end
    scale2 = prior.interaction_scale^2
    if prior.interaction_type === :cauchy
        @inbounds for j in prior.interaction_indices
            denominator = scale2 + x[j]^2
            out[j] = 2 * (scale2 - x[j]^2) / denominator^2 * direction[j]
        end
    else
        @inbounds for j in prior.interaction_indices
            out[j] = direction[j] / scale2
        end
    end

    return out
end

function _omrf_prior_hvp!(out, prior::OMRFNodeScaleAnalyticPrior, x, direction)
    fill!(out, 0.0)
    αβ = prior.threshold_alpha + prior.threshold_beta
    @inbounds for j in prior.threshold_indices
        p = _omrf_logistic(x[j])
        out[j] = αβ * p * (1 - p) * direction[j]
    end

    global_index = prior.global_scale_index
    global_scale = exp(x[global_index])
    global_direction = direction[global_index]
    out[global_index] = 2 * global_scale^2 /
        prior.global_scale_sd^2 * global_direction

    node_scales = prior.node_scales
    @inbounds for node in eachindex(node_scales)
        node_scales[node] = exp(x[prior.node_scale_indices[node]])
        node_index = prior.node_scale_indices[node]
        out[node_index] = 2 * node_scales[node]^2 /
            prior.node_scale_sd^2 * direction[node_index]
    end

    edge = 1
    @inbounds for left in 1:(length(node_scales) - 1), right in (left + 1):length(node_scales)
        interaction_index = prior.interaction_indices[edge]
        interaction = x[interaction_index]
        scale2 = prior.interaction_scale^2 * global_scale^2 *
            node_scales[left] * node_scales[right]
        inv_scale2 = inv(scale2)
        ratio = interaction^2 * inv_scale2
        scale_direction = global_direction +
            direction[prior.node_scale_indices[left]] +
            direction[prior.node_scale_indices[right]]
        out[interaction_index] = inv_scale2 * direction[interaction_index] -
            2 * interaction * inv_scale2 * scale_direction
        dc = -2 * interaction * inv_scale2 * direction[interaction_index] +
            2 * ratio * scale_direction
        out[global_index] += dc
        out[prior.node_scale_indices[left]] += 0.5 * dc
        out[prior.node_scale_indices[right]] += 0.5 * dc
        edge += 1
    end
    return out
end

function (provider::OMRFFullHVPProvider)(out, x, direction)
    context = provider.context
    prior = provider.prior
    _omrf_prior_hvp!(out, prior, x, direction)

    N, P = size(context.X)
    @inbounds for n in 1:N, node in 1:P
        _omrf_node_probabilities!(context, x, n, node)
        q = context.seen[node] - 1
        start = context.threshold_starts[node]
        directional_field = 0.0
        for k in eachindex(context.incident_edges[node])
            edge = context.incident_edges[node][k]
            neighbour = context.incident_neighbours[node][k]
            directional_field += context.X[n, neighbour] *
                direction[context.interaction_indices[edge]]
        end
        mean_direction = 0.0
        for u in 1:q
            local_index = start + u - 1
            local_direction = direction[context.threshold_indices[local_index]] +
                u * directional_field
            context.hcv_buffer[u] = local_direction
            mean_direction += context.probability_buffer[u] * local_direction
        end
        expected_derivative = 0.0
        for u in 1:q
            local_index = start + u - 1
            probability_derivative = context.probability_buffer[u] *
                (context.hcv_buffer[u] - mean_direction)
            out[context.threshold_indices[local_index]] += probability_derivative
            expected_derivative += u * probability_derivative
        end
        for k in eachindex(context.incident_edges[node])
            edge = context.incident_edges[node][k]
            neighbour = context.incident_neighbours[node][k]
            out[context.interaction_indices[edge]] +=
                context.X[n, neighbour] * expected_derivative
        end
    end
    context.counts.native_full_data_scans += 1
    return out
end

"""
    build_omrf_full_model(spec, unconstrained_names, d)

Build a non-factorized `PDMPModel` whose complete OMRF negative gradient and
Hessian-vector product are evaluated in native Julia.  The analytic prior and
all of its declared parameter blocks must cover all `d` unconstrained coordinates;
this intentionally rejects models with unimplemented nuisance parameters.
"""
function build_omrf_full_model(spec, unconstrained_names, d::Integer)
    d = Int(d)
    anchor = zeros(d)
    gradient_context = _omrf_residual_context(
        spec, unconstrained_names, anchor, StanSubsamplingCallCounts())
    prior_indices = gradient_context.global_scale_index > 0 ?
        [gradient_context.global_scale_index; gradient_context.node_scale_indices] : Int[]
    covered = sort!(vcat(copy(gradient_context.threshold_indices),
        copy(gradient_context.interaction_indices), prior_indices))
    covered == collect(1:d) || throw(ArgumentError(
        "native OMRF full gradients require all analytic parameter blocks " *
        "to cover every unconstrained coordinate; covered $covered of 1:$d"))
    gradient_prior = _omrf_analytic_prior(spec, gradient_context)
    gradient_prior === nothing && throw(ArgumentError(
        "native OMRF full gradients require a supported analytic OMRF prior"))
    gradient = OMRFFullGradientProvider(gradient_context, gradient_prior)

    hvp_context = _omrf_residual_context(
        spec, unconstrained_names, anchor, StanSubsamplingCallCounts())
    hvp_prior = something(_omrf_analytic_prior(spec, hvp_context))
    hvp = OMRFFullHVPProvider(hvp_context, hvp_prior)
    return PDMPModel(d, FullGradient(gradient), hvp)
end

function _omrf_deterministic_gradient!(out,
        oracle, prior::OMRFIndependentAnalyticPrior, x)
    @. out = oracle.full_anchor - oracle.prior_anchor
    α, β = prior.threshold_alpha, prior.threshold_beta
    @inbounds for j in prior.threshold_indices
        out[j] += (α + β) * _omrf_logistic(x[j]) - α
    end
    scale2 = prior.interaction_scale^2
    if prior.interaction_type === :cauchy
        @inbounds for j in prior.interaction_indices
            out[j] += 2x[j] / (scale2 + x[j]^2)
        end
    else
        @inbounds for j in prior.interaction_indices
            out[j] += x[j] / scale2
        end
    end
    return out
end

function _omrf_deterministic_gradient!(out,
        oracle, prior::OMRFNodeScaleAnalyticPrior, x)
    @. out = oracle.full_anchor - oracle.prior_anchor
    _omrf_prior_gradient!(oracle.prior_x, prior, x)
    @. out += oracle.prior_x
    return out
end

function _omrf_residual_context(envelope_spec, unc_names, anchor,
        counts::StanSubsamplingCallCounts)
    X = Matrix{Int}(_rget(envelope_spec, :X))
    seen = _as_int_vector(_rget(envelope_spec, :seen))
    N, P = size(X)
    length(seen) == P || throw(DimensionMismatch(
        "OMRF seen must contain one category count per node"))
    all(>=(1), seen) || throw(ArgumentError(
        "OMRF category counts must be positive"))
    @inbounds for j in 1:P, n in 1:N
        0 <= X[n, j] < seen[j] || throw(ArgumentError(
            "OMRF column $j must be encoded in 0:$(seen[j] - 1)"))
    end
    threshold_indices = _resolve_unc_spec(
        _rget(envelope_spec, :thresholds), unc_names, "OMRF thresholds")
    interaction_indices = _resolve_unc_spec(
        _rget(envelope_spec, :interactions), unc_names, "OMRF interactions")
    global_scale_index = 0
    node_scale_indices = Int[]
    if _haskey(envelope_spec, :analytic_prior)
        prior_spec = _rget(envelope_spec, :analytic_prior)
        if String(_rget(prior_spec, :type)) == "node_scale_gaussian"
            global_scale_index = only(_resolve_unc_spec(
                _rget(prior_spec, :global_scale), unc_names,
                "OMRF global scale"))
            node_scale_indices = _resolve_unc_spec(
                _rget(prior_spec, :node_scale), unc_names,
                "OMRF node scales")
            length(node_scale_indices) == P || throw(DimensionMismatch(
                "resolved OMRF node-scale block has $(length(node_scale_indices)) " *
                "coordinates; expected $P"))
            isempty(intersect(node_scale_indices,
                [global_scale_index; threshold_indices; interaction_indices])) ||
                throw(ArgumentError("resolved OMRF node-scale coordinates overlap another block"))
        end
    end
    no_thresholds = sum(seen .- 1)
    no_edges = P * (P - 1) ÷ 2
    length(threshold_indices) == no_thresholds || throw(DimensionMismatch(
        "resolved OMRF threshold block has $(length(threshold_indices)) coordinates; expected $no_thresholds"))
    length(interaction_indices) == no_edges || throw(DimensionMismatch(
        "resolved OMRF interaction block has $(length(interaction_indices)) coordinates; expected $no_edges"))
    isempty(intersect(threshold_indices, interaction_indices)) ||
        throw(ArgumentError("resolved OMRF threshold and interaction blocks overlap"))

    threshold_starts = Vector{Int}(undef, P)
    next_threshold = 1
    @inbounds for j in 1:P
        threshold_starts[j] = next_threshold
        next_threshold += seen[j] - 1
    end
    incident_edges = [Int[] for _ in 1:P]
    incident_neighbours = [Int[] for _ in 1:P]
    edge = 1
    if P > 1
        @inbounds for left in 1:(P - 1), right in (left + 1):P
            push!(incident_edges[left], edge)
            push!(incident_neighbours[left], right)
            push!(incident_edges[right], edge)
            push!(incident_neighbours[right], left)
            edge += 1
        end
    end
    use_hcv = _haskey(envelope_spec, :use_hcv) &&
        Bool(_rget(envelope_spec, :use_hcv))
    hcv_after_warmup = _haskey(envelope_spec, :hcv_after_warmup) &&
        Bool(_rget(envelope_spec, :hcv_after_warmup))
    hcv_damping = use_hcv ? Float64(_rget(envelope_spec, :hcv_damping)) : 1.0
    !isnan(hcv_damping) && hcv_damping > 0 || throw(ArgumentError(
        "OMRF HCV damping must be positive"))
    d = length(anchor)
    factorization = _haskey(envelope_spec, :factorization) ?
        Symbol(String(_rget(envelope_spec, :factorization))) : :person
    factorization in (:person, :person_node) || throw(ArgumentError(
        "unknown OMRF factorization: $factorization"))
    X_float = Matrix{Float64}(X)
    observed_threshold_counts = zeros(no_thresholds)
    @inbounds for node in 1:P, n in 1:N
        observed = X[n, node]
        if observed > 0
            observed_threshold_counts[
                threshold_starts[node] + observed - 1] += 1.0
        end
    end
    observed_interaction_sums = zeros(no_edges)
    edge = 1
    @inbounds for left in 1:(P - 1), right in (left + 1):P
        observed_interaction_sums[edge] =
            2.0 * dot(view(X_float, :, left), view(X_float, :, right))
        edge += 1
    end
    context = OMRFResidualContext(X, X_float, seen, threshold_indices,
        interaction_indices, global_scale_index, node_scale_indices,
        threshold_starts, incident_edges,
        incident_neighbours,
        Matrix{Float64}(undef, no_thresholds, N),
        Matrix{Float64}(undef, P, N), zeros(N, P), zeros(P, P),
        zeros(N, P), zeros(P, P), observed_threshold_counts,
        observed_interaction_sums, zeros(maximum(seen) - 1),
        zeros(maximum(seen) - 1), nothing, use_hcv,
        use_hcv && !hcv_after_warmup, hcv_damping,
        zeros(d), zeros(d), factorization, Int[], falses(d), counts)
    _cache_omrf_anchor!(context, anchor)
    use_hcv && (context.likelihood_hessian =
        _omrf_likelihood_hessian(context, anchor))
    return context
end

function _omrf_likelihood_hessian(context::OMRFResidualContext, anchor)
    d = length(anchor)
    hessian = zeros(d, d)
    N, P = size(context.X)
    max_support = maximum(context.seen .- 1) + max(P - 1, 0)
    design = zeros(maximum(context.seen) - 1, max_support)
    support = zeros(Int, max_support)
    mean_design = zeros(max_support)
    @inbounds for n in 1:N, node in 1:P
        q = context.seen[node] - 1
        iszero(q) && continue
        fill!(design, 0.0)
        fill!(mean_design, 0.0)
        start = context.threshold_starts[node]
        s = 0
        for u in 1:q
            s += 1
            support[s] = context.threshold_indices[start + u - 1]
            design[u, s] = 1.0
        end
        edges = context.incident_edges[node]
        neighbours = context.incident_neighbours[node]
        for k in eachindex(edges)
            s += 1
            support[s] = context.interaction_indices[edges[k]]
            neighbour_value = context.X[n, neighbours[k]]
            for u in 1:q
                design[u, s] = u * neighbour_value
            end
        end
        for u in 1:q
            p = context.anchor_probabilities[start + u - 1, n]
            for col in 1:s
                mean_design[col] += p * design[u, col]
            end
        end
        for row in 1:s, col in 1:s
            value = -mean_design[row] * mean_design[col]
            for u in 1:q
                value += context.anchor_probabilities[start + u - 1, n] *
                    design[u, row] * design[u, col]
            end
            hessian[support[row], support[col]] += value
        end
    end
    context.counts.native_full_data_scans += 1
    return hessian
end

@inline function _omrf_node_probabilities!(context::OMRFResidualContext,
        x, n::Int, node::Int)
    q = context.seen[node] - 1
    iszero(q) && return 0.0
    field = 0.0
    edges = context.incident_edges[node]
    neighbours = context.incident_neighbours[node]
    @inbounds for k in eachindex(edges)
        field += context.X[n, neighbours[k]] *
            x[context.interaction_indices[edges[k]]]
    end
    start = context.threshold_starts[node]
    max_eta = 0.0
    @inbounds for u in 1:q
        eta = x[context.threshold_indices[start + u - 1]] + u * field
        context.probability_buffer[u] = eta
        max_eta = max(max_eta, eta)
    end
    denominator = exp(-max_eta)
    @inbounds for u in 1:q
        probability = exp(context.probability_buffer[u] - max_eta)
        context.probability_buffer[u] = probability
        denominator += probability
    end
    expected = 0.0
    @inbounds for u in 1:q
        probability = context.probability_buffer[u] / denominator
        context.probability_buffer[u] = probability
        expected += u * probability
    end
    return expected
end

function _cache_omrf_anchor!(context::OMRFResidualContext, anchor)
    N, P = size(context.X)
    @inbounds for n in 1:N, node in 1:P
        expected = _omrf_node_probabilities!(context, anchor, n, node)
        context.anchor_expected_scores[node, n] = expected
        q = context.seen[node] - 1
        start = context.threshold_starts[node]
        for u in 1:q
            context.anchor_probabilities[start + u - 1, n] =
                context.probability_buffer[u]
        end
    end
    context.counts.native_full_data_scans += 1
    return context
end

function _omrf_node_residual!(context::OMRFResidualContext, out, x,
        n::Int, node::Int, α::Float64)
    expected = _omrf_node_probabilities!(context, x, n, node)
    expected_residual = expected - context.anchor_expected_scores[node, n]
    q = context.seen[node] - 1
    start = context.threshold_starts[node]
    correction_mean = 0.0
    if context.hcv_active
        field_delta = 0.0
        edges = context.incident_edges[node]
        neighbours = context.incident_neighbours[node]
        @inbounds for k in eachindex(edges)
            field_delta += context.X[n, neighbours[k]] *
                context.displacement[context.interaction_indices[edges[k]]]
        end
        @inbounds for u in 1:q
            delta_eta = context.displacement[
                context.threshold_indices[start + u - 1]] + u * field_delta
            context.hcv_buffer[u] = delta_eta
            correction_mean += context.anchor_probabilities[
                start + u - 1, n] * delta_eta
        end
    end
    correction_expected = 0.0
    @inbounds for u in 1:q
        local_threshold = start + u - 1
        correction = context.hcv_active ?
            context.anchor_probabilities[local_threshold, n] *
            (context.hcv_buffer[u] - correction_mean) : 0.0
        out[context.threshold_indices[local_threshold]] +=
            context.probability_buffer[u] -
            context.anchor_probabilities[local_threshold, n] - α * correction
        correction_expected += u * correction
    end
    edges = context.incident_edges[node]
    neighbours = context.incident_neighbours[node]
    @inbounds for k in eachindex(edges)
        out[context.interaction_indices[edges[k]]] +=
            context.X[n, neighbours[k]] *
            (expected_residual - α * correction_expected)
    end
    return out
end

function omrf_residual!(context::OMRFResidualContext, out, x, subset)
    started = time_ns()
    fill!(out, 0.0)
    @inbounds for index in context.touched_positions
        context.touched_mask[index] = false
    end
    empty!(context.touched_positions)
    N, P = size(context.X)
    α = if context.hcv_active
        displacement2 = dot(context.displacement, context.displacement)
        isinf(context.hcv_damping) ? 1.0 :
            context.hcv_damping / (context.hcv_damping + displacement2)
    else
        0.0
    end
    if context.factorization === :person
        @inbounds for n in subset
            1 <= n <= N || throw(BoundsError(context.X, (n, :)))
            for node in 1:P
                _record_omrf_touched_support!(context, node)
                _omrf_node_residual!(context, out, x, n, node, α)
            end
        end
        context.counts.analytic_persons += length(subset)
        context.counts.analytic_node_conditionals += P * length(subset)
    else
        @inbounds for factor in subset
            1 <= factor <= N * P || throw(BoundsError(1:(N * P), factor))
            node0, n = divrem(factor - 1, N)
            _record_omrf_touched_support!(context, node0 + 1)
            _omrf_node_residual!(context, out, x, n + 1, node0 + 1, α)
        end
        context.counts.analytic_factors += length(subset)
        context.counts.analytic_node_conditionals += length(subset)
    end
    context.counts.analytic_residual += 1
    context.counts.analytic_residual_seconds += (time_ns() - started) * 1e-9
    return out
end

"""
Frozen person-node control variate with independently chosen local references.

`anchor_logits[:, node, person]` stores the reference categorical logits for
one observed node conditional.  The references do not have to be projections
of a single parameter vector.  This is nevertheless an ordinary exact
control variate: `full_anchor_likelihood` is the complete sum of the factor
scores at those references and `omrf_residual!` returns the corresponding
factor-score differences.  Preparation is O(NP) work, while later residual
queries touch only the requested persons or person-node factors.

This small type deliberately excludes HCV.  A Hessian control variate shared
across incompatible local references would require a separate factor-local
Hessian representation.
"""
struct OMRFFactorLocalFrozenCV
    context::OMRFResidualContext
    anchor_logits::Array{Float64,3}
    full_anchor_likelihood::Vector{Float64}
    source_indices::Matrix{Int}
end

@inline function _omrf_node_logits!(out, context::OMRFResidualContext,
        x, n::Int, node::Int)
    q = context.seen[node] - 1
    field = 0.0
    @inbounds for k in eachindex(context.incident_edges[node])
        edge = context.incident_edges[node][k]
        neighbour = context.incident_neighbours[node][k]
        field += context.X[n, neighbour] *
            x[context.interaction_indices[edge]]
    end
    start = context.threshold_starts[node]
    @inbounds for u in 1:q
        out[u] = x[context.threshold_indices[start + u - 1]] + u * field
    end
    return q
end

function _prepare_omrf_factor_local_cv(envelope_spec, unc_names,
        candidates::AbstractMatrix,
        counts::StanSubsamplingCallCounts=StanSubsamplingCallCounts())
    d, K = size(candidates)
    K > 0 || throw(ArgumentError(
        "factor-local OMRF control variates need at least one candidate"))
    context = _omrf_residual_context(
        envelope_spec, unc_names, view(candidates, :, 1), counts)
    context.use_hcv && throw(ArgumentError(
        "factor-local OMRF control variates do not yet support HCV"))
    N, P = size(context.X)
    qmax = maximum(context.seen) - 1
    candidate_logits = zeros(qmax, K)
    anchor_logits = zeros(qmax, P, N)
    source_indices = Matrix{Int}(undef, P, N)
    full_anchor_likelihood = zeros(d)

    # Each factor selects the candidate closest to the centre of its own
    # projected warmup path.  This is a medoid, so the stored probabilities
    # are always probabilities evaluated at an observed warmup state.
    @inbounds for n in 1:N, node in 1:P
        q = context.seen[node] - 1
        iszero(q) && continue
        for candidate in 1:K
            _omrf_node_logits!(view(candidate_logits, :, candidate),
                context, view(candidates, :, candidate), n, node)
        end
        best = 1
        best_loss = Inf
        for candidate in 1:K
            loss = 0.0
            for other in 1:K, u in 1:q
                delta = candidate_logits[u, candidate] -
                    candidate_logits[u, other]
                loss += delta * delta
            end
            if loss < best_loss
                best = candidate
                best_loss = loss
            end
        end
        source_indices[node, n] = best
        copyto!(view(anchor_logits, 1:q, node, n),
            view(candidate_logits, 1:q, best))

        expected = _omrf_node_probabilities!(
            context, view(candidates, :, best), n, node)
        context.anchor_expected_scores[node, n] = expected
        observed = context.X[n, node]
        start = context.threshold_starts[node]
        for u in 1:q
            probability = context.probability_buffer[u]
            context.anchor_probabilities[start + u - 1, n] = probability
            full_anchor_likelihood[
                context.threshold_indices[start + u - 1]] +=
                probability - (observed == u)
        end
        for k in eachindex(context.incident_edges[node])
            edge = context.incident_edges[node][k]
            neighbour = context.incident_neighbours[node][k]
            full_anchor_likelihood[context.interaction_indices[edge]] +=
                context.X[n, neighbour] * (expected - observed)
        end
    end
    counts.native_full_data_scans += 1
    return OMRFFactorLocalFrozenCV(context, anchor_logits,
        full_anchor_likelihood, source_indices)
end

function (cv::OMRFFactorLocalFrozenCV)(out, x, subset)
    return omrf_residual!(cv.context, out, x, subset)
end

"""Pointwise scalar-rate certificate for one locally anchored node factor."""
function _omrf_factor_local_bound(cv::OMRFFactorLocalFrozenCV,
        x, velocity, n::Int, node::Int)
    context = cv.context
    q = context.seen[node] - 1
    iszero(q) && return 0.0
    field = 0.0
    field_velocity = 0.0
    @inbounds for k in eachindex(context.incident_edges[node])
        edge = context.incident_edges[node][k]
        neighbour = context.incident_neighbours[node][k]
        covariate = context.X[n, neighbour]
        field += covariate * x[context.interaction_indices[edge]]
        field_velocity += covariate *
            velocity[context.interaction_indices[edge]]
    end
    minimum_displacement = 0.0
    maximum_displacement = 0.0
    minimum_velocity = 0.0
    maximum_velocity = 0.0
    start = context.threshold_starts[node]
    @inbounds for u in 1:q
        displacement = x[context.threshold_indices[start + u - 1]] +
            u * field - cv.anchor_logits[u, node, n]
        logit_velocity = velocity[
            context.threshold_indices[start + u - 1]] + u * field_velocity
        minimum_displacement = min(minimum_displacement, displacement)
        maximum_displacement = max(maximum_displacement, displacement)
        minimum_velocity = min(minimum_velocity, logit_velocity)
        maximum_velocity = max(maximum_velocity, logit_velocity)
    end
    return 0.25 * (maximum_displacement - minimum_displacement) *
        (maximum_velocity - minimum_velocity)
end

"""
Prepared node-local bank whose references may be recombined without data work.

For node `j`, source `k` contains the exact cached sum
`sum_n grad U_nj(anchor_k)`. At a PDMP event boundary a selector may choose a
different `k_j` for every node. The deterministic likelihood term is then
`sum_j baseline[:,j,k_j]` and a sampled person's residual uses the matching
cached node probabilities. Selection and recombination cost O(PK + Pd), not
O(NP), and all `KNP` factor evaluations occur during preparation.
"""
struct OMRFNodeLocalFrozenBank
    contexts::Vector{OMRFResidualContext}
    node_baselines::Array{Float64,3}
    candidates::Matrix{Float64}
end

"""Warmup-prepared node-specific gradient/Hessian control-variate bank."""
mutable struct OMRFNodeLocalHCVBank
    contexts::Vector{OMRFResidualContext}
    candidates::Matrix{Float64}
    node_gradients::Array{Float64,3}
    node_hessians::Array{Float64,4}
    node_supports::Vector{Vector{Int}}
    sources::Vector{Int}
    requested_sources::Vector{Int}
    active_sources::BitVector
    proxy_persons::Vector{Int}
    node_pattern_persons::Vector{Vector{Int}}
    node_pattern_counts::Vector{Vector{Int}}
    intercept::Vector{Float64}
    hessian::Matrix{Float64}
    product::Vector{Float64}
    selector_positions::Matrix{Float64}
    selector_velocities::Matrix{Float64}
    selector_sources::Matrix{Int}
end

function _prepare_omrf_node_local_hcv_bank(entries, counts)
    K = length(entries)
    K > 0 || throw(ArgumentError("node-local HCV bank cannot be empty"))
    contexts = [entry.state.residual_context for entry in entries]
    all(context -> context.use_hcv && isinf(context.hcv_damping), contexts) ||
        throw(ArgumentError("node-local HCV requires exact undamped contexts"))
    d = length(entries[1].state.anchor)
    N, P = size(contexts[1].X)
    candidates = reduce(hcat, (entry.state.anchor for entry in entries))
    gradients = zeros(d, P, K)
    hessians = zeros(d, d, P, K)
    max_support = maximum(contexts[1].seen .- 1) + max(P - 1, 0)
    design = zeros(maximum(contexts[1].seen) - 1, max_support)
    support = zeros(Int, max_support)
    mean_design = zeros(max_support)
    @inbounds for candidate in 1:K
        context = contexts[candidate]
        for n in 1:N, node in 1:P
            expected = context.anchor_expected_scores[node, n]
            observed = context.X[n, node]
            q = context.seen[node] - 1
            start = context.threshold_starts[node]
            fill!(design, 0.0)
            fill!(mean_design, 0.0)
            s = 0
            for u in 1:q
                s += 1
                index = context.threshold_indices[start + u - 1]
                support[s] = index
                design[u, s] = 1.0
                gradients[index, node, candidate] +=
                    context.anchor_probabilities[start + u - 1, n] -
                    (observed == u)
            end
            for k in eachindex(context.incident_edges[node])
                s += 1
                edge = context.incident_edges[node][k]
                neighbour = context.incident_neighbours[node][k]
                index = context.interaction_indices[edge]
                support[s] = index
                value = context.X[n, neighbour]
                gradients[index, node, candidate] += value *
                    (expected - observed)
                for u in 1:q
                    design[u, s] = u * value
                end
            end
            for u in 1:q
                probability = context.anchor_probabilities[start + u - 1, n]
                for col in 1:s
                    mean_design[col] += probability * design[u, col]
                end
            end
            for row in 1:s, col in 1:s
                value = -mean_design[row] * mean_design[col]
                for u in 1:q
                    value += context.anchor_probabilities[start + u - 1, n] *
                        design[u, row] * design[u, col]
                end
                hessians[support[row], support[col], node, candidate] += value
            end
        end
    end
    counts.native_full_data_scans += K
    proxy_count = min(16, N)
    proxy_persons = unique(round.(Int, range(1, N; length=proxy_count)))
    node_pattern_persons = [Int[] for _ in 1:P]
    node_pattern_counts = [Int[] for _ in 1:P]
    @inbounds for node in 1:P
        neighbours = contexts[1].incident_neighbours[node]
        pattern_indices = Dict{Tuple,Int}()
        for person in 1:N
            pattern = Tuple(contexts[1].X[person, neighbour]
                for neighbour in neighbours)
            pattern_index = get(pattern_indices, pattern, 0)
            if iszero(pattern_index)
                push!(node_pattern_persons[node], person)
                push!(node_pattern_counts[node], 1)
                pattern_indices[pattern] = length(node_pattern_persons[node])
            else
                node_pattern_counts[node][pattern_index] += 1
            end
        end
    end
    node_supports = _omrf_node_supports(contexts[1])
    bank = OMRFNodeLocalHCVBank(contexts, candidates, gradients, hessians,
        node_supports, ones(Int, P), ones(Int, P), falses(K), proxy_persons,
        node_pattern_persons, node_pattern_counts,
        zeros(d), zeros(d, d), zeros(d), zeros(d, 0), zeros(d, 0),
        zeros(Int, P, 0))
    return _install_omrf_node_local_sources!(bank, bank.sources)
end

function _install_omrf_node_local_sources!(bank::OMRFNodeLocalHCVBank,
        requested)
    length(requested) == length(bank.sources) || throw(DimensionMismatch(
        "node-local HCV needs one source per node"))
    fill!(bank.intercept, 0.0)
    fill!(bank.hessian, 0.0)
    _, P, K = size(bank.node_gradients)
    @inbounds for node in 1:P
        source = Int(requested[node])
        1 <= source <= K || throw(BoundsError(1:K, source))
        bank.sources[node] = source
        support = bank.node_supports[node]
        for row in support
            value = bank.node_gradients[row, node, source]
            for col in support
                h = bank.node_hessians[row, col, node, source]
                bank.hessian[row, col] += h
                value -= h * bank.candidates[col, source]
            end
            bank.intercept[row] += value
        end
    end
    fill!(bank.active_sources, false)
    @inbounds for source in bank.sources
        bank.active_sources[source] = true
    end
    return bank
end

function _install_omrf_node_local_sources!(manager, cv, requested)
    bank = something(manager.oracle.node_local_hcv)
    changed = requested != bank.sources
    changed && _install_omrf_node_local_sources!(bank, requested)
    scales = cv.envelope.component_scales!
    if scales isa OMRFBatchedFactorHCVScales
        _install_omrf_batched_node_sources!(scales, bank)
    elseif scales isa Union{OMRFGroupedPersonHCVBoxScales,
            OMRFNodeGroupedPersonHCVBoxScales}
        _install_omrf_grouped_box_node_sources!(scales, bank)
    elseif scales isa OMRFCellCachedNodeGroupedPersonHCVBoxScales
        _install_omrf_grouped_box_node_sources!(scales.base, bank)
    elseif scales isa OMRFLazyNodeGroupedPersonHCVBoxScales
        _install_omrf_grouped_box_node_sources!(scales, bank)
    elseif scales isa OMRFPatternLocalHCVScales
        _install_omrf_pattern_node_sources!(scales, bank)
    else
        throw(ArgumentError(
            "node-local HCV sources require a compatible HCV envelope"))
    end
    return changed
end

@inline function _omrf_node_hcv_bound_at(context::OMRFResidualContext,
        anchor, state, flow, person::Int, node::Int)
    q = context.seen[node] - 1
    iszero(q) && return 0.0
    field_displacement = 0.0
    field_velocity = 0.0
    @inbounds for k in eachindex(context.incident_edges[node])
        covariate = context.X[person, context.incident_neighbours[node][k]]
        iszero(covariate) && continue
        index = context.interaction_indices[context.incident_edges[node][k]]
        displacement, velocity = _omrf_coordinate_at(
            state, flow, anchor, index, 0.0)
        field_displacement += covariate * displacement
        field_velocity += covariate * velocity
    end
    min_displacement = 0.0
    max_displacement = 0.0
    min_velocity = 0.0
    max_velocity = 0.0
    start = context.threshold_starts[node]
    @inbounds for u in 1:q
        index = context.threshold_indices[start + u - 1]
        displacement, velocity = _omrf_coordinate_at(
            state, flow, anchor, index, 0.0)
        logit_displacement = displacement + u * field_displacement
        logit_velocity = velocity + u * field_velocity
        min_displacement = min(min_displacement, logit_displacement)
        max_displacement = max(max_displacement, logit_displacement)
        min_velocity = min(min_velocity, logit_velocity)
        max_velocity = max(max_velocity, logit_velocity)
    end
    return (max_velocity - min_velocity) *
        (max_displacement - min_displacement)^2 / 8
end

function _select_omrf_node_local_sources!(manager, cv, state, flow)
    selection_started = time_ns()
    bank = something(manager.oracle.node_local_hcv)
    context = bank.contexts[1]
    P = size(context.X, 2)
    K = size(bank.candidates, 2)
    requested = bank.requested_sources
    scales = cv.envelope.component_scales!
    @inbounds for node in 1:P
        best_source = 1
        best_score = Inf
        for source in 1:K
            score = 0.0
            anchor = view(bank.candidates, :, source)
            if scales isa OMRFGroupedPersonHCVBoxScales
                # The scalable clock owns a setup-time summary of each group.
                # Score candidate sources with those same certified boxes so
                # retained anchor selection does not revisit even a proxy
                # panel of person rows.
                for group in axes(scales.minima, 1)
                    score += cv.envelope.totals[group] *
                        _omrf_grouped_box_node_bound_at(
                            scales, state, flow, group, node, 0.0, anchor)
                end
            elseif scales isa OMRFNodeGroupedPersonHCVBoxScales
                source_anchor = view(bank.candidates, :, source)
                _cache_omrf_node_grouped_point_node!(
                    scales, state, flow, node, source_anchor)
                for component in scales.node_component_ranges[node]
                    score += cv.envelope.totals[component] *
                        _omrf_node_grouped_point_bound(
                            scales, component, node)
                end
            elseif scales isa OMRFCellCachedNodeGroupedPersonHCVBoxScales
                base = scales.base
                source_anchor = view(bank.candidates, :, source)
                _cache_omrf_node_grouped_point_node!(
                    base, state, flow, node, source_anchor)
                for component in base.node_component_ranges[node]
                    score += cv.envelope.totals[component] *
                        _omrf_node_grouped_point_bound(
                            base, component, node)
                end
            elseif scales isa OMRFLazyNodeGroupedPersonHCVBoxScales
                base = scales.base
                source_anchor = view(bank.candidates, :, source)
                _cache_omrf_node_grouped_point_node!(
                    base, state, flow, node, source_anchor)
                for root in scales.roots
                    base.component_nodes[root] == node || continue
                    score += length(scales.members[root]) *
                        _omrf_node_grouped_point_bound(base, root, node)
                end
            elseif scales isa OMRFPatternLocalHCVScales
                _cache_omrf_pattern_point_coordinates!(
                    scales, state, flow, 0.0, bank.candidates, source)
                for component in scales.node_component_ranges[node]
                    score += cv.envelope.totals[component] *
                        _omrf_pattern_hcv_scale_cached(scales, component)
                end
            else
                # Deprecated batched-person comparison path.
                for person in bank.proxy_persons
                    score += _omrf_node_hcv_bound_at(
                        context, anchor, state, flow, person, node)
                end
            end
            if score < best_score
                best_source = source
                best_score = score
            end
        end
        requested[node] = best_source
    end
    changed = _install_omrf_node_local_sources!(manager, cv, requested)
    if changed
        manager.main_activations += 1
        manager.node_source_changes += 1
    end
    manager.main_selections += 1
    manager.selection_seconds += (time_ns() - selection_started) * 1e-9
    if scales isa Union{OMRFNodeGroupedPersonHCVBoxScales,
            OMRFCellCachedNodeGroupedPersonHCVBoxScales,
            OMRFLazyNodeGroupedPersonHCVBoxScales,
            OMRFPatternLocalHCVScales}
        component_count = if scales isa OMRFLazyNodeGroupedPersonHCVBoxScales
            length(scales.roots)
        elseif scales isa OMRFCellCachedNodeGroupedPersonHCVBoxScales
            length(scales.base.component_nodes)
        else
            length(scales.component_nodes)
        end
        manager.selection_component_evaluations +=
            size(bank.candidates, 2) * component_count
    end
    return changed
end

@inline function _omrf_node_source_moved_enough(
        manager, x)
    threshold = manager.node_source_check_distance
    iszero(threshold) && return true
    distance2 = 0.0
    @inbounds for j in eachindex(x, manager.node_source_last_position)
        delta = (x[j] - manager.node_source_last_position[j]) /
            manager.distance_scales[j]
        distance2 += delta * delta
    end
    return distance2 >= threshold^2
end

function _omrf_node_local_hcv_residual!(bank::OMRFNodeLocalHCVBank,
        out, x, subset, counts; factorization::Symbol=:person)
    fill!(out, 0.0)
    N, P = size(bank.contexts[1].X)
    @inbounds for source in eachindex(bank.active_sources)
        bank.active_sources[source] || continue
        context = bank.contexts[source]
        for j in eachindex(context.displacement, x)
            context.displacement[j] = x[j] - bank.candidates[j, source]
        end
    end
    if factorization === :person
        @inbounds for n in subset, node in 1:P
            _omrf_node_residual!(bank.contexts[bank.sources[node]],
                out, x, n, node, 1.0)
        end
        counts.analytic_persons += length(subset)
        counts.analytic_node_conditionals += P * length(subset)
    else
        @inbounds for factor in subset
            node0, n0 = divrem(factor - 1, N)
            node, n = node0 + 1, n0 + 1
            _omrf_node_residual!(bank.contexts[bank.sources[node]],
                out, x, n, node, 1.0)
        end
        counts.analytic_factors += length(subset)
        counts.analytic_node_conditionals += length(subset)
    end
    counts.analytic_residual += 1
    return out
end

function _prepare_omrf_node_local_bank(envelope_spec, unc_names,
        candidates::AbstractMatrix,
        counts::StanSubsamplingCallCounts=StanSubsamplingCallCounts())
    d, K = size(candidates)
    K > 0 || throw(ArgumentError(
        "node-local OMRF banks need at least one candidate"))
    contexts = [_omrf_residual_context(envelope_spec, unc_names,
        view(candidates, :, candidate), counts) for candidate in 1:K]
    any(context -> context.use_hcv, contexts) && throw(ArgumentError(
        "node-local OMRF banks do not yet support HCV"))
    N, P = size(contexts[1].X)
    baselines = zeros(d, P, K)
    @inbounds for candidate in 1:K
        context = contexts[candidate]
        for n in 1:N, node in 1:P
            expected = context.anchor_expected_scores[node, n]
            observed = context.X[n, node]
            q = context.seen[node] - 1
            start = context.threshold_starts[node]
            for u in 1:q
                probability = context.anchor_probabilities[start + u - 1, n]
                baselines[context.threshold_indices[start + u - 1],
                    node, candidate] += probability - (observed == u)
            end
            for k in eachindex(context.incident_edges[node])
                edge = context.incident_edges[node][k]
                neighbour = context.incident_neighbours[node][k]
                baselines[context.interaction_indices[edge], node, candidate] +=
                    context.X[n, neighbour] * (expected - observed)
            end
        end
    end
    counts.native_full_data_scans += K
    return OMRFNodeLocalFrozenBank(
        contexts, baselines, Matrix{Float64}(candidates))
end

function _omrf_node_local_baseline!(out, bank::OMRFNodeLocalFrozenBank,
        sources::AbstractVector{<:Integer})
    P = size(bank.node_baselines, 2)
    K = size(bank.node_baselines, 3)
    length(sources) == P || throw(DimensionMismatch(
        "node-local source vector must have one entry per node"))
    fill!(out, 0.0)
    @inbounds for node in 1:P
        source = Int(sources[node])
        1 <= source <= K || throw(BoundsError(1:K, source))
        for j in eachindex(out)
            out[j] += bank.node_baselines[j, node, source]
        end
    end
    return out
end

function _omrf_node_local_residual!(bank::OMRFNodeLocalFrozenBank, out, x,
        subset, sources::AbstractVector{<:Integer}; factorization::Symbol=:person)
    fill!(out, 0.0)
    context = bank.contexts[1]
    N, P = size(context.X)
    if factorization === :person
        @inbounds for n in subset, node in 1:P
            _omrf_node_residual!(bank.contexts[sources[node]],
                out, x, n, node, 0.0)
        end
    elseif factorization === :person_node
        @inbounds for factor in subset
            node0, n0 = divrem(factor - 1, N)
            node, n = node0 + 1, n0 + 1
            _omrf_node_residual!(bank.contexts[sources[node]],
                out, x, n, node, 0.0)
        end
    else
        throw(ArgumentError("unknown OMRF factorization: $factorization"))
    end
    return out
end

function _record_omrf_touched_support!(context::OMRFResidualContext, node::Int)
    q = context.seen[node] - 1
    start = context.threshold_starts[node]
    @inbounds for u in 1:q
        _record_omrf_touched!(context,
            context.threshold_indices[start + u - 1])
    end
    @inbounds for edge in context.incident_edges[node]
        _record_omrf_touched!(context, context.interaction_indices[edge])
    end
    return nothing
end

@inline function _record_omrf_touched!(context::OMRFResidualContext, index::Int)
    context.touched_mask[index] && return nothing
    context.touched_mask[index] = true
    push!(context.touched_positions, index)
    return nothing
end

mutable struct OMRFAnalyticSubsamplingOracle
    ctx::StanSubsamplingContext
    anchor::Vector{Float64}
    full_anchor::Vector{Float64}
    prior_anchor::Vector{Float64}
    prior_x::Vector{Float64}
    prior_x_position::Vector{Float64}
    prior_x_valid::Bool
    residual_context::OMRFResidualContext
    analytic_prior::Union{Nothing,OMRFAnalyticPrior}
    node_local_hcv::Union{Nothing,OMRFNodeLocalHCVBank}
    lazy_residual_subset_bound::Float64
end

function OMRFAnalyticSubsamplingOracle(ctx::StanSubsamplingContext,
        envelope_spec, unc_names, anchor::Vector{Float64})
    d = length(anchor)
    residual_context = _omrf_residual_context(
        envelope_spec, unc_names, anchor, ctx.counts)
    analytic_prior = _omrf_analytic_prior(envelope_spec, residual_context)
    oracle = OMRFAnalyticSubsamplingOracle(ctx, copy(anchor), zeros(d), zeros(d),
        zeros(d), zeros(d), false, residual_context, analytic_prior, nothing,
        NaN)
    _clear_gradient!(ctx, :full, oracle.full_anchor, anchor)
    if analytic_prior === nothing
        _clear_gradient!(ctx, :prior, oracle.prior_anchor, anchor)
    else
        _omrf_prior_gradient!(oracle.prior_anchor, analytic_prior, anchor)
    end
    return oracle
end

function _prior_gradient_at_x!(oracle::OMRFAnalyticSubsamplingOracle, x)
    if !oracle.prior_x_valid || oracle.prior_x_position != x
        if oracle.analytic_prior === nothing
            _clear_gradient!(oracle.ctx, :prior, oracle.prior_x, x)
        else
            _omrf_prior_gradient!(oracle.prior_x, oracle.analytic_prior, x)
        end
        copyto!(oracle.prior_x_position, x)
        oracle.prior_x_valid = true
    end
    return oracle.prior_x
end

function (oracle::OMRFAnalyticSubsamplingOracle)(out, x)
    if oracle.node_local_hcv !== nothing
        prior = something(oracle.analytic_prior)
        _omrf_prior_gradient!(out, prior, x)
        bank = oracle.node_local_hcv
        mul!(bank.product, bank.hessian, x)
        @. out += bank.intercept + bank.product
        return out
    end
    if oracle.analytic_prior === nothing
        _prior_gradient_at_x!(oracle, x)
        @. out = oracle.full_anchor + oracle.prior_x - oracle.prior_anchor
    else
        _omrf_deterministic_gradient!(out, oracle, oracle.analytic_prior, x)
    end
    context = oracle.residual_context
    if context.hcv_active
        @. context.displacement = x - oracle.anchor
        displacement2 = dot(context.displacement, context.displacement)
        α = isinf(context.hcv_damping) ? 1.0 :
            context.hcv_damping / (context.hcv_damping + displacement2)
        mul!(context.hcv_product, something(context.likelihood_hessian),
            context.displacement)
        @. out += α * context.hcv_product
    end
    return out
end

function (oracle::OMRFAnalyticSubsamplingOracle)(out, x, subset, anchor)
    if oracle.node_local_hcv !== nothing
        started = time_ns()
        _omrf_node_local_hcv_residual!(oracle.node_local_hcv, out, x, subset,
            oracle.residual_context.counts;
            factorization=oracle.residual_context.factorization)
        oracle.residual_context.counts.analytic_residual_seconds +=
            (time_ns() - started) * 1e-9
        return out
    end
    context = oracle.residual_context
    context.hcv_active && (@. context.displacement = x - oracle.anchor)
    return omrf_residual!(oracle.residual_context, out, x, subset)
end

function PDMPSamplers.subsampling_candidate_rate!(
        oracle::OMRFAnalyticSubsamplingOracle, state, gradient, residual,
        scale, flow::ZigZag, deterministic_rate, subset)
    actual = deterministic_rate
    @inbounds for index in oracle.residual_context.touched_positions
        velocity = state.ξ.θ[index]
        old_gradient = gradient[index]
        new_gradient = old_gradient + scale * residual[index]
        gradient[index] = new_gradient
        actual += max(0.0, velocity * new_gradient) -
            max(0.0, velocity * old_gradient)
    end
    return max(0.0, actual)
end

function PDMPSamplers.record_subsampling_proposal!(
        oracle::OMRFAnalyticSubsamplingOracle,
        D, B, subset_bound, deterministic_actual, residual_rate,
        actual, accepted)
    counts = oracle.residual_context.counts
    counts.proposal_deterministic_bound_sum += D
    counts.proposal_aggregate_residual_bound_sum += B
    counts.proposal_subset_bound_sum += subset_bound
    counts.proposal_deterministic_actual_sum += deterministic_actual
    counts.proposal_exact_residual_rate_sum += residual_rate
    counts.proposal_actual_rate_sum += actual
    certified_exact = deterministic_actual + residual_rate
    counts.proposal_count += 1
    counts.proposal_acceptance_probability_sum +=
        actual / max(subset_bound, eps(Float64))
    counts.proposal_certificate_ratio_sum +=
        certified_exact / max(subset_bound, eps(Float64))
    counts.proposal_intrinsic_ratio_sum +=
        actual / max(certified_exact, eps(Float64))
    if accepted
        counts.accepted_proposal_bound_sum += subset_bound
        counts.accepted_actual_rate_sum += actual
    end
    return nothing
end

_anchor_cache_bytes(oracle::StanSubsamplingOracle) =
    Base.summarysize(oracle.anchor_likelihood_cache)
_anchor_cache_bytes(oracle::OMRFAnalyticSubsamplingOracle) =
    Base.summarysize(oracle.residual_context.anchor_probabilities) +
    Base.summarysize(oracle.residual_context.anchor_expected_scores) +
    Base.summarysize(oracle.residual_context.likelihood_hessian)

"""
Certified node-local OMRF trajectory geometry.

Each component uses only its node thresholds and incident interactions. The
pointwise scale is `norm(v[S]) * norm(x[S] - anchor[S])`; multiplying it by
`0.5 * opnorm(B_ni)^2` dominates the node-conditional residual event rate for
both coordinatewise Zig-Zag rates and scalar BPS/Boomerang rates. Cell methods
bound every coordinate over the complete interval before taking both norms.
"""
struct OMRFNodeLocalScales
    anchor::Vector{Float64}
    supports::Vector{Vector{Int}}
end

"""
Certified pattern-local OMRF residual-rate scales.

For a node with non-reference logit displacement `z` and `q` non-reference
categories, write `R = maximum((0, z...)) - minimum((0, z...))`.  Along the
line from the anchor to the current state, categorical covariance gives
`abs(Δp_u) <= R/4` and Hoeffding's covariance inequality gives
`abs(Δμ) <= qR/4`.  Consequently the complete node-factor Zig-Zag rate, and
also the absolute scalar rate used by BPS/Boomerang, is bounded by

`R/4 * (sum(abs(v_threshold)) + q * sum(X_neighbour * abs(v_edge)))`.

Components group factors with identical node-local covariates.  Cell methods
bound the logit range and velocity factor over the whole closed interval.
"""
mutable struct OMRFPatternLocalScales
    anchor::Vector{Float64}
    component_nodes::Vector{Int}
    component_covariates::Matrix{Int}
    context::OMRFResidualContext
    offsets::Vector{Float64}
    cosines::Vector{Float64}
    sines::Vector{Float64}
    coordinate_offsets::Vector{Float64}
    coordinate_cosines::Vector{Float64}
    coordinate_sines::Vector{Float64}
end

"""
Certified pattern-local scalar certificate for the undamped categorical HCV.

For one node let `a = B*(x-anchor)` and `b = B*v`, including the reference
category with logit zero. The directional second-order Taylor remainder of
`logsumexp` is bounded by `range(b) * range(a)^2 / 8`. Components group equal
observed node patterns exactly as `OMRFPatternLocalScales` does.
"""
mutable struct OMRFPatternLocalHCVScales
    anchor::Vector{Float64}
    node_anchors::Union{Nothing,Matrix{Float64}}
    node_sources::Vector{Int}
    component_nodes::Vector{Int}
    node_component_ranges::Vector{UnitRange{Int}}
    component_covariates::Matrix{Int}
    context::OMRFResidualContext
    offsets::Vector{Float64}
    cosines::Vector{Float64}
    sines::Vector{Float64}
    coordinate_offsets::Vector{Float64}
    coordinate_cosines::Vector{Float64}
    coordinate_sines::Vector{Float64}
end

"""
Legacy allocation-free person-node HCV certificates retained for comparison.

The three `field_*` matrices hold `X * beta_*`, where the interaction
coefficient matrices encode displacement/velocity harmonic coefficients.
Consequently each callback uses three small BLAS matrix products followed by
one contiguous person-node/category pass. `node_anchors` permits every node to
own an independently selected warmup-bank source; the resulting output entry
is the complete mass of that person-node/source residual factor.

This exact callback materializes one state-dependent certificate per
person-node at every point and cell query. It therefore performs `O(NP)`
retained work and is deprecated for new scalable runs. Use
`geometry = "grouped_box"`, backed by `OMRFGroupedPersonHCVBoxScales`.
"""
mutable struct OMRFBatchedFactorHCVScales
    anchor::Vector{Float64}
    node_anchors::Union{Nothing,Matrix{Float64}}
    node_sources::Vector{Int}
    context::OMRFResidualContext
    interaction_offsets::Matrix{Float64}
    interaction_cosines::Matrix{Float64}
    interaction_sines::Matrix{Float64}
    field_offsets::Matrix{Float64}
    field_cosines::Matrix{Float64}
    field_sines::Matrix{Float64}
end

"""
Exact whole-person HCV clock over fixed setup-time covariate boxes.

Each output is a complete bound for one group of whole-person residuals. The
group extrema and membership are immutable after setup; retained callbacks
inspect only `G` boxes and model coordinates, never the `N` person rows.
"""
mutable struct OMRFGroupedPersonHCVBoxScales
    anchor::Vector{Float64}
    node_anchors::Union{Nothing,Matrix{Float64}}
    node_sources::Vector{Int}
    minima::Matrix{Float64}
    maxima::Matrix{Float64}
    context::OMRFResidualContext
end

"""
Exact node-aware whole-person HCV clock over node-local support boxes.

Components are node/box pairs. Each person belongs to one box in every node
partition, so summing its `P` component bounds certifies its complete residual.
Trajectory coefficients are cached once per node/category pair and reused by
all boxes of that node.
"""
mutable struct OMRFNodeGroupedPersonHCVBoxScales
    anchor::Vector{Float64}
    node_anchors::Union{Nothing,Matrix{Float64}}
    node_sources::Vector{Int}
    component_nodes::Vector{Int}
    node_component_ranges::Vector{UnitRange{Int}}
    minima::Matrix{Float64}
    maxima::Matrix{Float64}
    context::OMRFResidualContext
    pair_offsets::Vector{Int}
    edge_offsets_by_node::Vector{Int}
    point_displacement_centers::Vector{Float64}
    point_velocity_centers::Vector{Float64}
    edge_displacements::Vector{Float64}
    edge_velocities::Vector{Float64}
    cell_offsets::Vector{Float64}
    cell_cosines::Vector{Float64}
    cell_sines::Vector{Float64}
    edge_offsets::Vector{Float64}
    edge_cosines::Vector{Float64}
    edge_sines::Vector{Float64}
    edge_displacement_bounds::Vector{Float64}
    edge_velocity_bounds::Vector{Float64}
end

"""
Exact cell-cached node-aware whole-person HCV clock.

Grid construction has already certified every node/support component over the
cell from which a proposal is drawn. This wrapper retains those masses for
marked-person selection, then evaluates only the selected person's pointwise
certificate. Thinning from its cell mass to its point mass recovers the same
tight pointwise marked clock as the eager implementation without traversing
all boxes at every proposal.
"""
mutable struct OMRFCellCachedNodeGroupedPersonHCVBoxScales
    base::OMRFNodeGroupedPersonHCVBoxScales
    cell_lefts::Vector{Float64}
    cell_rights::Vector{Float64}
    cell_scales::Matrix{Float64}
    n_cached_cells::Int
    cached_state::Any
    cached_flow::Any
    cached_time::Float64
    selected_person_point_bound::Float64
end

"""
Exact lazy support hierarchy for the node-grouped whole-person HCV clock.

The embedded `base` provider owns every internal and leaf support box. Only
the root frontier is exposed to the grid clock. Pointwise marked-person draws
then descend one certified branch at a time by thinning a parent box into its
children. This avoids materializing all leaf masses while preserving their
exact size-biased mixture.
"""
mutable struct OMRFLazyNodeGroupedPersonHCVBoxScales
    base::OMRFNodeGroupedPersonHCVBoxScales
    parents::Vector{Int}
    children::Vector{NTuple{2,Int}}
    roots::Vector{Int}
    leaves::BitVector
    members::Vector{Vector{Int}}
    root_totals::Vector{Float64}
    root_scales::Vector{Float64}
    root_cumulative_masses::Vector{Float64}
    leaf_scales::Vector{Float64}
    leaf_groups::Matrix{Int}
    cached_state::Any
    cached_flow::Any
    cached_time::Float64
    cache_valid::Bool
    person_leaf_bounds::Vector{Float64}
end

function _install_omrf_batched_node_sources!(
        scales::OMRFBatchedFactorHCVScales,
        bank::OMRFNodeLocalHCVBank)
    scales.node_anchors = bank.candidates
    resize!(scales.node_sources, length(bank.sources))
    copyto!(scales.node_sources, bank.sources)
    return scales
end

function _install_omrf_grouped_box_node_sources!(
        scales::OMRFGroupedPersonHCVBoxScales,
        bank::OMRFNodeLocalHCVBank)
    scales.node_anchors = bank.candidates
    resize!(scales.node_sources, length(bank.sources))
    copyto!(scales.node_sources, bank.sources)
    return scales
end

function _install_omrf_grouped_box_node_sources!(
        scales::OMRFNodeGroupedPersonHCVBoxScales,
        bank::OMRFNodeLocalHCVBank)
    scales.node_anchors = bank.candidates
    resize!(scales.node_sources, length(bank.sources))
    copyto!(scales.node_sources, bank.sources)
    return scales
end

function _install_omrf_grouped_box_node_sources!(
        scales::OMRFLazyNodeGroupedPersonHCVBoxScales,
        bank::OMRFNodeLocalHCVBank)
    _install_omrf_grouped_box_node_sources!(scales.base, bank)
    scales.cache_valid = false
    return scales
end

function _install_omrf_pattern_node_sources!(
        scales::OMRFPatternLocalHCVScales,
        bank::OMRFNodeLocalHCVBank)
    scales.node_anchors = bank.candidates
    resize!(scales.node_sources, length(bank.sources))
    copyto!(scales.node_sources, bank.sources)
    return scales
end

"""Low-rank node-local certificate using neighbour-vector Euclidean norms."""
struct OMRFNormLocalScales
    anchor::Vector{Float64}
    context::OMRFResidualContext
end

"""
Certified empirical-covariate expansion of the scalar OMRF residual rate.

For node `j`, let `A,C` bound the ranges of the threshold displacement and
velocity, and let `d_k,v_k` bound an incident edge's displacement and
velocity.  For person `n`, the node contribution is bounded by

`(A + q*sum_k X[n,k]d_k) * (C + q*sum_k X[n,k]v_k) / 4`.

The expansion uses static person covariates as envelope weights and shared
trajectory terms as scales.  It therefore retains substantially more design
information than the norm-local certificate without evaluating every
observed node pattern during each grid build.
"""
struct OMRFCovariateLocalScales
    anchor::Vector{Float64}
    context::OMRFResidualContext
    edge_displacements::Vector{Float64}
    edge_velocities::Vector{Float64}
    coordinate_offsets::Vector{Float64}
    coordinate_cosines::Vector{Float64}
    coordinate_sines::Vector{Float64}
end

"""Certified node-local scales for the damped categorical HCV remainder."""
struct OMRFDampedHCVNodeLocalScales
    anchor::Vector{Float64}
    supports::Vector{Vector{Int}}
    context::OMRFResidualContext
end

@inline _omrf_coordinate_free(state::StickyPDMPState, j) = state.free[j]
@inline _omrf_coordinate_free(state::PDMPSamplers.AbstractPDMPState, j) = true

@inline function _omrf_coordinate_at(state, ::Union{ZigZag,BouncyParticle},
        anchor, j, t)
    displacement = state.ξ.x[j] + t * state.ξ.θ[j] - anchor[j]
    velocity = _omrf_coordinate_free(state, j) ? state.ξ.θ[j] : 0.0
    return displacement, velocity
end

@inline function _omrf_coordinate_at_anchor_value(
        state, ::Union{ZigZag,BouncyParticle}, anchor_value, j, t)
    displacement = state.ξ.x[j] + t * state.ξ.θ[j] - anchor_value
    velocity = _omrf_coordinate_free(state, j) ? state.ξ.θ[j] : 0.0
    return displacement, velocity
end

@inline function _omrf_coordinate_at(state, flow::PDMPSamplers.AnyBoomerang,
        anchor, j, t)
    return _omrf_coordinate_at_anchor_value(
        state, flow, anchor[j], j, t)
end

@inline function _omrf_coordinate_at_anchor_value(
        state, flow::PDMPSamplers.AnyBoomerang, anchor_value, j, t)
    if !_omrf_coordinate_free(state, j)
        return state.ξ.x[j] - anchor_value, 0.0
    end
    s, c = sincos(t)
    centered = state.ξ.x[j] - flow.μ[j]
    displacement = flow.μ[j] - anchor_value +
        centered * c + state.ξ.θ[j] * s
    velocity = -centered * s + state.ξ.θ[j] * c
    return displacement, velocity
end

@inline _omrf_coordinate_at(state, flow::PDMPSamplers.PreconditionedDynamics,
    anchor, j, t) = _omrf_coordinate_at(state, flow.dynamics, anchor, j, t)

@inline _omrf_coordinate_at_anchor_value(state,
    flow::PDMPSamplers.PreconditionedDynamics, anchor_value, j, t) =
    _omrf_coordinate_at_anchor_value(
        state, flow.dynamics, anchor_value, j, t)

function _omrf_norm_node_terms_at(scales::OMRFNormLocalScales,
        state, flow, node::Int, t)
    context = scales.context
    q = context.seen[node] - 1
    iszero(q) && return 0.0, 0.0, 0.0, 0.0
    threshold_start = context.threshold_starts[node]
    minimum_threshold = 0.0
    maximum_threshold = 0.0
    threshold_velocity = 0.0
    @inbounds for u in 1:q
        index = context.threshold_indices[threshold_start + u - 1]
        displacement, velocity = _omrf_coordinate_at(
            state, flow, scales.anchor, index, t)
        minimum_threshold = min(minimum_threshold, displacement)
        maximum_threshold = max(maximum_threshold, displacement)
        threshold_velocity += abs(velocity)
    end
    edge_displacement2 = 0.0
    edge_velocity2 = 0.0
    @inbounds for edge in context.incident_edges[node]
        index = context.interaction_indices[edge]
        displacement, velocity = _omrf_coordinate_at(
            state, flow, scales.anchor, index, t)
        edge_displacement2 += displacement^2
        edge_velocity2 += velocity^2
    end
    A = maximum_threshold - minimum_threshold
    B = q * sqrt(edge_displacement2)
    C = threshold_velocity
    D = q * sqrt(edge_velocity2)
    return A, B, C, D
end

@inline function _fill_omrf_norm_scales!(out, node, P, A, B, C, D)
    out[node] = 0.25 * A * C
    out[P + node] = 0.25 * (A * D + B * C)
    out[2P + node] = 0.25 * B * D
    return nothing
end

function (scales::OMRFNormLocalScales)(out, state, flow, t)
    P = length(scales.context.seen)
    @inbounds for node in 1:P
        A, B, C, D = _omrf_norm_node_terms_at(
            scales, state, flow, node, t)
        _fill_omrf_norm_scales!(out, node, P, A, B, C, D)
    end
    return out
end

function (scales::OMRFNormLocalScales)(out, state,
        flow::Union{ZigZag,BouncyParticle}, left, right)
    P = length(scales.context.seen)
    @inbounds for node in 1:P
        Al, Bl, C, D = _omrf_norm_node_terms_at(
            scales, state, flow, node, left)
        Ar, Br, _, _ = _omrf_norm_node_terms_at(
            scales, state, flow, node, right)
        _fill_omrf_norm_scales!(out, node, P,
            max(Al, Ar), max(Bl, Br), C, D)
    end
    return out
end

function (scales::OMRFNormLocalScales)(out, state,
        flow::PDMPSamplers.AnyBoomerang, left, right)
    context = scales.context
    P = length(context.seen)
    @inbounds for node in 1:P
        q = context.seen[node] - 1
        if iszero(q)
            _fill_omrf_norm_scales!(out, node, P, 0.0, 0.0, 0.0, 0.0)
            continue
        end
        threshold_start = context.threshold_starts[node]
        threshold_coefficients = Vector{NTuple{3,Float64}}(undef, q + 1)
        threshold_coefficients[1] = (0.0, 0.0, 0.0)
        C = 0.0
        for u in 1:q
            index = context.threshold_indices[threshold_start + u - 1]
            coefficients = _omrf_harmonic_coefficients(
                state, flow, scales.anchor, index)
            threshold_coefficients[u + 1] = coefficients
            C += PDMPSamplers._harmonic_abs_max(
                0.0, coefficients[3], -coefficients[2], left, right)
        end
        A = 0.0
        for u in 2:(q + 1), v in 1:(u - 1)
            cu, cv = threshold_coefficients[u], threshold_coefficients[v]
            A = max(A, PDMPSamplers._harmonic_abs_max(
                cu[1] - cv[1], cu[2] - cv[2], cu[3] - cv[3], left, right))
        end
        displacement2 = 0.0
        velocity2 = 0.0
        for edge in context.incident_edges[node]
            index = context.interaction_indices[edge]
            offset, cosine, sine = _omrf_harmonic_coefficients(
                state, flow, scales.anchor, index)
            displacement = PDMPSamplers._harmonic_abs_max(
                offset, cosine, sine, left, right)
            velocity = PDMPSamplers._harmonic_abs_max(
                0.0, sine, -cosine, left, right)
            displacement2 += displacement^2
            velocity2 += velocity^2
        end
        _fill_omrf_norm_scales!(out, node, P, A,
            q * sqrt(displacement2), C, q * sqrt(velocity2))
    end
    return out
end

@inline function _fill_omrf_covariate_node_scales!(out, first, q, A, C,
        edge_displacements, edge_velocities, degree)
    destination = first
    out[destination] = 0.25 * A * C
    destination += 1
    @inbounds for k in 1:degree
        out[destination] = 0.25 * q * C * edge_displacements[k]
        destination += 1
    end
    @inbounds for k in 1:degree
        out[destination] = 0.25 * q * A * edge_velocities[k]
        destination += 1
    end
    @inbounds for k in 1:degree, l in 1:degree
        out[destination] = 0.25 * q^2 *
            edge_displacements[k] * edge_velocities[l]
        destination += 1
    end
    return destination
end

function (scales::OMRFCovariateLocalScales)(out, state, flow, t)
    context = scales.context
    _cache_omrf_pattern_point_coordinates!(scales, state, flow, t)
    destination = 1
    @inbounds for node in eachindex(context.seen)
        q = context.seen[node] - 1
        minimum_displacement = 0.0
        maximum_displacement = 0.0
        minimum_velocity = 0.0
        maximum_velocity = 0.0
        threshold_start = context.threshold_starts[node]
        if q > 0
            for u in 1:q
                index = context.threshold_indices[threshold_start + u - 1]
                displacement = scales.coordinate_offsets[index]
                velocity = scales.coordinate_cosines[index]
                minimum_displacement = min(minimum_displacement, displacement)
                maximum_displacement = max(maximum_displacement, displacement)
                minimum_velocity = min(minimum_velocity, velocity)
                maximum_velocity = max(maximum_velocity, velocity)
            end
        end
        edges = context.incident_edges[node]
        degree = length(edges)
        edge_displacements = scales.edge_displacements
        edge_velocities = scales.edge_velocities
        for k in 1:degree
            index = context.interaction_indices[edges[k]]
            displacement = scales.coordinate_offsets[index]
            velocity = scales.coordinate_cosines[index]
            edge_displacements[k] = abs(displacement)
            edge_velocities[k] = abs(velocity)
        end
        destination = _fill_omrf_covariate_node_scales!(out, destination, q,
            maximum_displacement - minimum_displacement,
            maximum_velocity - minimum_velocity,
            edge_displacements, edge_velocities, degree)
    end
    return out
end

function (scales::OMRFCovariateLocalScales)(out, state,
        flow::Union{ZigZag,BouncyParticle}, left, right)
    context = scales.context
    destination = 1
    @inbounds for node in eachindex(context.seen)
        q = context.seen[node] - 1
        threshold_start = context.threshold_starts[node]
        minimum_velocity = 0.0
        maximum_velocity = 0.0
        A = 0.0
        if q > 0
            for u in 1:q, v in 0:(u - 1)
                index = context.threshold_indices[threshold_start + u - 1]
                displacement_left, velocity = _omrf_coordinate_at(
                    state, flow, scales.anchor, index, left)
                displacement_right, _ = _omrf_coordinate_at(
                    state, flow, scales.anchor, index, right)
                minimum_velocity = min(minimum_velocity, velocity)
                maximum_velocity = max(maximum_velocity, velocity)
                if v > 0
                    index = context.threshold_indices[threshold_start + v - 1]
                    other_left, _ = _omrf_coordinate_at(
                        state, flow, scales.anchor, index, left)
                    other_right, _ = _omrf_coordinate_at(
                        state, flow, scales.anchor, index, right)
                    displacement_left -= other_left
                    displacement_right -= other_right
                end
                A = max(A, abs(displacement_left), abs(displacement_right))
            end
        end
        edges = context.incident_edges[node]
        degree = length(edges)
        edge_displacements = scales.edge_displacements
        edge_velocities = scales.edge_velocities
        for k in 1:degree
            index = context.interaction_indices[edges[k]]
            left_displacement, velocity = _omrf_coordinate_at(
                state, flow, scales.anchor, index, left)
            right_displacement, _ = _omrf_coordinate_at(
                state, flow, scales.anchor, index, right)
            edge_displacements[k] = max(
                abs(left_displacement), abs(right_displacement))
            edge_velocities[k] = abs(velocity)
        end
        destination = _fill_omrf_covariate_node_scales!(out, destination, q,
            A, maximum_velocity - minimum_velocity,
            edge_displacements, edge_velocities, degree)
    end
    return out
end

function PDMPSamplers.screening_residual_bound(
        envelope::PDMPSamplers.BlockSeparableResidualEnvelope{F}, state,
        flow::Union{ZigZag,BouncyParticle}, t) where
        {F<:OMRFCovariateLocalScales}
    scales = envelope.component_scales!
    context = scales.context
    _cache_omrf_pattern_point_coordinates!(scales, state, flow, t)
    totals = envelope.totals
    destination = 1
    result = 0.0
    @inbounds for node in eachindex(context.seen)
        q = context.seen[node] - 1
        minimum_displacement = 0.0
        maximum_displacement = 0.0
        minimum_velocity = 0.0
        maximum_velocity = 0.0
        threshold_start = context.threshold_starts[node]
        for u in 1:q
            index = context.threshold_indices[threshold_start + u - 1]
            displacement = scales.coordinate_offsets[index]
            velocity = scales.coordinate_cosines[index]
            minimum_displacement = min(minimum_displacement, displacement)
            maximum_displacement = max(maximum_displacement, displacement)
            minimum_velocity = min(minimum_velocity, velocity)
            maximum_velocity = max(maximum_velocity, velocity)
        end
        A = maximum_displacement - minimum_displacement
        C = maximum_velocity - minimum_velocity
        result += 0.25 * A * C * totals[destination]
        destination += 1
        edges = context.incident_edges[node]
        degree = length(edges)
        edge_displacements = scales.edge_displacements
        edge_velocities = scales.edge_velocities
        for k in 1:degree
            index = context.interaction_indices[edges[k]]
            edge_displacements[k] = abs(scales.coordinate_offsets[index])
            edge_velocities[k] = abs(scales.coordinate_cosines[index])
            result += 0.25 * q * C * edge_displacements[k] *
                totals[destination]
            destination += 1
        end
        for k in 1:degree
            result += 0.25 * q * A * edge_velocities[k] *
                totals[destination]
            destination += 1
        end
        for k in 1:degree, l in 1:degree
            result += 0.25 * q^2 * edge_displacements[k] *
                edge_velocities[l] * totals[destination]
            destination += 1
        end
    end
    return nextfloat(result)
end

function PDMPSamplers.prepare_residual_sampling!(
        envelope::PDMPSamplers.BlockSeparableResidualEnvelope{F}, state,
        flow::Union{ZigZag,BouncyParticle}, t) where
        {F<:OMRFCovariateLocalScales}
    PDMPSamplers.total_residual_bound(envelope, state, flow, t)
    return nothing
end

function (scales::OMRFCovariateLocalScales)(out, state,
        flow::PDMPSamplers.AnyBoomerang, left, right)
    context = scales.context
    _cache_omrf_pattern_harmonic_coefficients!(scales, state, flow)
    midpoint = (left + right) / 2
    midpoint_sine, midpoint_cosine = sincos(midpoint)
    half_width = (right - left) / 2
    destination = 1
    @inbounds for node in eachindex(context.seen)
        q = context.seen[node] - 1
        threshold_start = context.threshold_starts[node]
        A = 0.0
        C = 0.0
        if q > 0
            for u in 1:q, v in 0:(u - 1)
                index = context.threshold_indices[threshold_start + u - 1]
                offset = scales.coordinate_offsets[index]
                cosine = scales.coordinate_cosines[index]
                sine = scales.coordinate_sines[index]
                if v > 0
                    index = context.threshold_indices[threshold_start + v - 1]
                    offset -= scales.coordinate_offsets[index]
                    cosine -= scales.coordinate_cosines[index]
                    sine -= scales.coordinate_sines[index]
                end
                A = max(A, _omrf_harmonic_abs_cell_bound(
                    offset, cosine, sine, midpoint_sine,
                    midpoint_cosine, half_width))
                C = max(C, _omrf_harmonic_abs_cell_bound(
                    0.0, sine, -cosine, midpoint_sine,
                    midpoint_cosine, half_width))
            end
        end
        edges = context.incident_edges[node]
        degree = length(edges)
        edge_displacements = scales.edge_displacements
        edge_velocities = scales.edge_velocities
        for k in 1:degree
            index = context.interaction_indices[edges[k]]
            offset = scales.coordinate_offsets[index]
            cosine = scales.coordinate_cosines[index]
            sine = scales.coordinate_sines[index]
            edge_displacements[k] = _omrf_harmonic_abs_cell_bound(
                offset, cosine, sine, midpoint_sine,
                midpoint_cosine, half_width)
            edge_velocities[k] = _omrf_harmonic_abs_cell_bound(
                0.0, sine, -cosine, midpoint_sine,
                midpoint_cosine, half_width)
        end
        destination = _fill_omrf_covariate_node_scales!(out, destination, q,
            A, C, edge_displacements, edge_velocities, degree)
    end
    return out
end

(scales::OMRFCovariateLocalScales)(out, state,
    flow::PDMPSamplers.PreconditionedDynamics, args...) =
    scales(out, state, flow.dynamics, args...)

(scales::OMRFNormLocalScales)(out, state,
    flow::PDMPSamplers.PreconditionedDynamics, args...) =
    scales(out, state, flow.dynamics, args...)

function _omrf_pattern_scale_at(scales::OMRFPatternLocalScales,
        state, flow, component::Int, t)
    context = scales.context
    node = scales.component_nodes[component]
    q = context.seen[node] - 1
    iszero(q) && return 0.0
    threshold_start = context.threshold_starts[node]
    field_displacement = 0.0
    field_velocity = 0.0
    weighted_edge_velocity = 0.0
    @inbounds for k in eachindex(context.incident_edges[node])
        covariate = scales.component_covariates[
            component, context.incident_neighbours[node][k]]
        iszero(covariate) && continue
        edge_index = context.interaction_indices[context.incident_edges[node][k]]
        displacement, velocity = _omrf_coordinate_at(
            state, flow, scales.anchor, edge_index, t)
        field_displacement += covariate * displacement
        field_velocity += covariate * velocity
        weighted_edge_velocity += covariate * abs(velocity)
    end
    minimum_logit = 0.0
    maximum_logit = 0.0
    minimum_logit_velocity = 0.0
    maximum_logit_velocity = 0.0
    threshold_velocity = 0.0
    @inbounds for u in 1:q
        index = context.threshold_indices[threshold_start + u - 1]
        displacement, velocity = _omrf_coordinate_at(
            state, flow, scales.anchor, index, t)
        logit = displacement + u * field_displacement
        minimum_logit = min(minimum_logit, logit)
        maximum_logit = max(maximum_logit, logit)
        logit_velocity = velocity + u * field_velocity
        minimum_logit_velocity = min(minimum_logit_velocity, logit_velocity)
        maximum_logit_velocity = max(maximum_logit_velocity, logit_velocity)
        threshold_velocity += abs(velocity)
    end
    velocity_range = if flow isa ZigZag
        # Coordinatewise Zig-Zag rates require an L1 bound; signed
        # directional cancellation is available only to scalar-rate flows.
        threshold_velocity + q * weighted_edge_velocity
    else
        maximum_logit_velocity - minimum_logit_velocity
    end
    return 0.25 * (maximum_logit - minimum_logit) * velocity_range
end

function _cache_omrf_pattern_point_coordinates!(
        scales::Union{OMRFPatternLocalScales,OMRFPatternLocalHCVScales,
            OMRFCovariateLocalScales},
        state, flow, t)
    context = scales.context
    displacements = scales.coordinate_offsets
    velocities = scales.coordinate_cosines
    @inbounds for index in context.threshold_indices
        displacements[index], velocities[index] = _omrf_coordinate_at(
            state, flow, scales.anchor, index, t)
    end
    @inbounds for index in context.interaction_indices
        displacements[index], velocities[index] = _omrf_coordinate_at(
            state, flow, scales.anchor, index, t)
    end
    return nothing
end

function _cache_omrf_pattern_point_coordinates!(
        scales::OMRFPatternLocalHCVScales, state,
        flow::PDMPSamplers.AnyBoomerang, t,
        anchors::Matrix{Float64}, source::Int)
    context = scales.context
    displacements = scales.coordinate_offsets
    velocities = scales.coordinate_cosines
    sine, cosine = sincos(t)
    @inbounds for index in Iterators.flatten(
            (context.threshold_indices, context.interaction_indices))
        if _omrf_coordinate_free(state, index)
            centered = state.ξ.x[index] - flow.μ[index]
            displacements[index] = flow.μ[index] - anchors[index, source] +
                centered * cosine + state.ξ.θ[index] * sine
            velocities[index] = -centered * sine + state.ξ.θ[index] * cosine
        else
            displacements[index] = state.ξ.x[index] - anchors[index, source]
            velocities[index] = 0.0
        end
    end
    return nothing
end

function _cache_omrf_pattern_point_coordinates!(
        scales::Union{OMRFPatternLocalScales,OMRFPatternLocalHCVScales,
            OMRFCovariateLocalScales},
        state, flow::PDMPSamplers.AnyBoomerang, t)
    context = scales.context
    displacements = scales.coordinate_offsets
    velocities = scales.coordinate_cosines
    sine, cosine = sincos(t)
    @inbounds for index in context.threshold_indices
        if _omrf_coordinate_free(state, index)
            centered = state.ξ.x[index] - flow.μ[index]
            displacements[index] = flow.μ[index] - scales.anchor[index] +
                centered * cosine + state.ξ.θ[index] * sine
            velocities[index] = -centered * sine + state.ξ.θ[index] * cosine
        else
            displacements[index] = state.ξ.x[index] - scales.anchor[index]
            velocities[index] = 0.0
        end
    end
    @inbounds for index in context.interaction_indices
        if _omrf_coordinate_free(state, index)
            centered = state.ξ.x[index] - flow.μ[index]
            displacements[index] = flow.μ[index] - scales.anchor[index] +
                centered * cosine + state.ξ.θ[index] * sine
            velocities[index] = -centered * sine + state.ξ.θ[index] * cosine
        else
            displacements[index] = state.ξ.x[index] - scales.anchor[index]
            velocities[index] = 0.0
        end
    end
    return nothing
end

function _omrf_pattern_scale_cached(scales::OMRFPatternLocalScales,
        flow, component::Int)
    context = scales.context
    node = scales.component_nodes[component]
    q = context.seen[node] - 1
    iszero(q) && return 0.0
    displacements = scales.coordinate_offsets
    velocities = scales.coordinate_cosines
    field_displacement = 0.0
    field_velocity = 0.0
    weighted_edge_velocity = 0.0
    @inbounds for k in eachindex(context.incident_edges[node])
        covariate = scales.component_covariates[
            component, context.incident_neighbours[node][k]]
        iszero(covariate) && continue
        index = context.interaction_indices[context.incident_edges[node][k]]
        field_displacement += covariate * displacements[index]
        field_velocity += covariate * velocities[index]
        weighted_edge_velocity += covariate * abs(velocities[index])
    end
    minimum_logit = 0.0
    maximum_logit = 0.0
    minimum_logit_velocity = 0.0
    maximum_logit_velocity = 0.0
    threshold_velocity = 0.0
    threshold_start = context.threshold_starts[node]
    @inbounds for u in 1:q
        index = context.threshold_indices[threshold_start + u - 1]
        logit = displacements[index] + u * field_displacement
        minimum_logit = min(minimum_logit, logit)
        maximum_logit = max(maximum_logit, logit)
        logit_velocity = velocities[index] + u * field_velocity
        minimum_logit_velocity = min(minimum_logit_velocity, logit_velocity)
        maximum_logit_velocity = max(maximum_logit_velocity, logit_velocity)
        threshold_velocity += abs(velocities[index])
    end
    velocity_range = if flow isa ZigZag
        threshold_velocity + q * weighted_edge_velocity
    else
        maximum_logit_velocity - minimum_logit_velocity
    end
    return 0.25 * (maximum_logit - minimum_logit) * velocity_range
end

@inline function _omrf_person_pattern_bound_at(context::OMRFResidualContext,
        anchor, state, flow, person::Int, t)
    bound = 0.0
    @inbounds for node in eachindex(context.seen)
        q = context.seen[node] - 1
        iszero(q) && continue
        field_displacement = 0.0
        field_velocity = 0.0
        for k in eachindex(context.incident_edges[node])
            covariate = context.X[
                person, context.incident_neighbours[node][k]]
            iszero(covariate) && continue
            edge_index = context.interaction_indices[
                context.incident_edges[node][k]]
            displacement, velocity = _omrf_coordinate_at(
                state, flow, anchor, edge_index, t)
            field_displacement += covariate * displacement
            field_velocity += covariate * velocity
        end
        minimum_logit = 0.0
        maximum_logit = 0.0
        minimum_velocity = 0.0
        maximum_velocity = 0.0
        threshold_start = context.threshold_starts[node]
        for u in 1:q
            index = context.threshold_indices[threshold_start + u - 1]
            displacement, velocity = _omrf_coordinate_at(
                state, flow, anchor, index, t)
            logit = displacement + u * field_displacement
            logit_velocity = velocity + u * field_velocity
            minimum_logit = min(minimum_logit, logit)
            maximum_logit = max(maximum_logit, logit)
            minimum_velocity = min(minimum_velocity, logit_velocity)
            maximum_velocity = max(maximum_velocity, logit_velocity)
        end
        bound += 0.25 * (maximum_logit - minimum_logit) *
            (maximum_velocity - minimum_velocity)
    end
    return bound
end

@inline function _omrf_person_node_pattern_bound_at(
        context::OMRFResidualContext, anchor, state, flow,
        person::Int, node::Int, t)
    q = context.seen[node] - 1
    iszero(q) && return 0.0
    field_displacement = 0.0
    field_velocity = 0.0
    weighted_edge_velocity = 0.0
    @inbounds for k in eachindex(context.incident_edges[node])
        covariate = context.X[
            person, context.incident_neighbours[node][k]]
        iszero(covariate) && continue
        edge_index = context.interaction_indices[
            context.incident_edges[node][k]]
        displacement, velocity = _omrf_coordinate_at(
            state, flow, anchor, edge_index, t)
        field_displacement += covariate * displacement
        field_velocity += covariate * velocity
        weighted_edge_velocity += covariate * abs(velocity)
    end
    minimum_logit = 0.0
    maximum_logit = 0.0
    minimum_logit_velocity = 0.0
    maximum_logit_velocity = 0.0
    threshold_velocity = 0.0
    threshold_start = context.threshold_starts[node]
    @inbounds for u in 1:q
        index = context.threshold_indices[threshold_start + u - 1]
        displacement, velocity = _omrf_coordinate_at(
            state, flow, anchor, index, t)
        logit = displacement + u * field_displacement
        minimum_logit = min(minimum_logit, logit)
        maximum_logit = max(maximum_logit, logit)
        logit_velocity = velocity + u * field_velocity
        minimum_logit_velocity = min(
            minimum_logit_velocity, logit_velocity)
        maximum_logit_velocity = max(
            maximum_logit_velocity, logit_velocity)
        threshold_velocity += abs(velocity)
    end
    zigzag_dynamics = flow isa ZigZag ||
        (flow isa PDMPSamplers.PreconditionedDynamics &&
         flow.dynamics isa ZigZag)
    velocity_range = zigzag_dynamics ?
        threshold_velocity + q * weighted_edge_velocity :
        maximum_logit_velocity - minimum_logit_velocity
    return 0.25 * (maximum_logit - minimum_logit) *
        velocity_range
end

function PDMPSamplers.subsampling_residual_subset_bound(
        oracle::OMRFAnalyticSubsamplingOracle, state, flow,
        D, M, subset, scale)
    if isfinite(oracle.lazy_residual_subset_bound)
        bound = D + oracle.lazy_residual_subset_bound
        oracle.lazy_residual_subset_bound = NaN
        return min(M, bound)
    end
    context = oracle.residual_context
    factor_pattern_dynamics = flow isa Union{ZigZag,BouncyParticle} ||
        (flow isa PDMPSamplers.PreconditionedDynamics &&
         flow.dynamics isa Union{ZigZag,BouncyParticle})
    if factor_pattern_dynamics && context.factorization === :person_node &&
            !context.use_hcv
        N = size(context.X, 1)
        pattern_bound = 0.0
        @inbounds for factor in subset
            node0, person0 = divrem(factor - 1, N)
            pattern_bound += _omrf_person_node_pattern_bound_at(
                context, oracle.anchor, state, flow,
                person0 + 1, node0 + 1, 0.0)
        end
        residual_bound = min(max(0.0, M - D), scale * pattern_bound)
        return D + residual_bound
    end
    scalar_dynamics = flow isa Union{BouncyParticle,PDMPSamplers.AnyBoomerang} ||
        (flow isa PDMPSamplers.PreconditionedDynamics &&
         flow.dynamics isa Union{BouncyParticle,PDMPSamplers.AnyBoomerang})
    if !scalar_dynamics || context.factorization !== :person || context.use_hcv
        return M
    end
    pattern_bound = 0.0
    @inbounds for person in subset
        pattern_bound += _omrf_person_pattern_bound_at(
            context, oracle.anchor, state, flow, person, 0.0)
    end
    residual_bound = min(max(0.0, M - D), scale * pattern_bound)
    return D + residual_bound
end

function (scales::OMRFPatternLocalScales)(out, state, flow, t)
    _cache_omrf_pattern_point_coordinates!(scales, state, flow, t)
    @inbounds for component in eachindex(out, scales.component_nodes)
        out[component] = _omrf_pattern_scale_cached(
            scales, flow, component)
    end
    return out
end

function (scales::OMRFPatternLocalScales)(out, state,
        flow::Union{ZigZag,BouncyParticle}, left, right)
    @inbounds for component in eachindex(out, scales.component_nodes)
        # The logit range is convex for a linear trajectory; its interval
        # maximum is therefore attained at an endpoint.  Velocity is constant.
        out[component] = max(
            _omrf_pattern_scale_at(scales, state, flow, component, left),
            _omrf_pattern_scale_at(scales, state, flow, component, right))
    end
    return out
end

@inline function _omrf_harmonic_coefficients(state, flow, anchor, j)
    return _omrf_harmonic_coefficients_anchor_value(
        state, flow, anchor[j], j)
end

@inline function _omrf_harmonic_coefficients_anchor_value(
        state, flow, anchor_value, j)
    if !_omrf_coordinate_free(state, j)
        return state.ξ.x[j] - anchor_value, 0.0, 0.0
    end
    return flow.μ[j] - anchor_value,
        state.ξ.x[j] - flow.μ[j], state.ξ.θ[j]
end

@inline function _omrf_harmonic_abs_cell_bound(c, a, b, left, right)
    midpoint = (left + right) / 2
    sine, cosine = sincos(midpoint)
    return _omrf_harmonic_abs_cell_bound(
        c, a, b, sine, cosine, (right - left) / 2)
end

@inline function _omrf_harmonic_abs_cell_bound(
        c, a, b, midpoint_sine, midpoint_cosine, half_width)
    amplitude = hypot(a, b)
    # For f(t) = c + a*cos(t) + b*sin(t), |f'| <= hypot(a, b).
    # The midpoint Lipschitz roof covers the complete cell.  The global
    # amplitude roof is also certified and prevents wide cells from inflating
    # the result unnecessarily.
    midpoint_roof = abs(c + a * midpoint_cosine + b * midpoint_sine) +
        half_width * amplitude
    return nextfloat(min(abs(c) + amplitude, midpoint_roof))
end

function _cache_omrf_pattern_harmonic_coefficients!(
        scales::Union{OMRFPatternLocalScales,OMRFPatternLocalHCVScales,
            OMRFCovariateLocalScales},
        state, flow)
    context = scales.context
    offsets = scales.coordinate_offsets
    cosines = scales.coordinate_cosines
    sines = scales.coordinate_sines
    @inbounds for index in context.threshold_indices
        offsets[index], cosines[index], sines[index] =
            _omrf_harmonic_coefficients(state, flow, scales.anchor, index)
    end
    @inbounds for index in context.interaction_indices
        offsets[index], cosines[index], sines[index] =
            _omrf_harmonic_coefficients(state, flow, scales.anchor, index)
    end
    return nothing
end

function _cache_omrf_pattern_harmonic_coefficients!(
        scales::OMRFPatternLocalHCVScales, state, flow,
        anchors::Matrix{Float64}, source::Int)
    context = scales.context
    @inbounds for index in Iterators.flatten(
            (context.threshold_indices, context.interaction_indices))
        if !_omrf_coordinate_free(state, index)
            scales.coordinate_offsets[index] =
                state.ξ.x[index] - anchors[index, source]
            scales.coordinate_cosines[index] = 0.0
            scales.coordinate_sines[index] = 0.0
        else
            scales.coordinate_offsets[index] =
                flow.μ[index] - anchors[index, source]
            scales.coordinate_cosines[index] = state.ξ.x[index] - flow.μ[index]
            scales.coordinate_sines[index] = state.ξ.θ[index]
        end
    end
    return nothing
end

function (scales::OMRFPatternLocalScales)(out, state,
        flow::PDMPSamplers.AnyBoomerang, left, right)
    context = scales.context
    _cache_omrf_pattern_harmonic_coefficients!(scales, state, flow)
    midpoint_sine, midpoint_cosine = sincos((left + right) / 2)
    half_width = (right - left) / 2
    coordinate_offsets = scales.coordinate_offsets
    coordinate_cosines = scales.coordinate_cosines
    coordinate_sines = scales.coordinate_sines
    @inbounds for component in eachindex(out, scales.component_nodes)
        node = scales.component_nodes[component]
        q = context.seen[node] - 1
        if iszero(q)
            out[component] = 0.0
            continue
        end
        offsets = scales.offsets
        cosines = scales.cosines
        sines = scales.sines
        offsets[1] = 0.0
        cosines[1] = 0.0
        sines[1] = 0.0
        threshold_start = context.threshold_starts[node]
        for u in 1:q
            index = context.threshold_indices[threshold_start + u - 1]
            offsets[u + 1] = coordinate_offsets[index]
            cosines[u + 1] = coordinate_cosines[index]
            sines[u + 1] = coordinate_sines[index]
        end
        for k in eachindex(context.incident_edges[node])
            covariate = scales.component_covariates[
                component, context.incident_neighbours[node][k]]
            iszero(covariate) && continue
            edge_index = context.interaction_indices[context.incident_edges[node][k]]
            offset = coordinate_offsets[edge_index]
            cosine = coordinate_cosines[edge_index]
            sine = coordinate_sines[edge_index]
            for u in 1:q
                offsets[u + 1] += u * covariate * offset
                cosines[u + 1] += u * covariate * cosine
                sines[u + 1] += u * covariate * sine
            end
        end
        range_bound = 0.0
        velocity_range_bound = 0.0
        for u in 1:(q + 1), v in 1:(u - 1)
            cosine_difference = cosines[u] - cosines[v]
            sine_difference = sines[u] - sines[v]
            range_bound = max(range_bound, _omrf_harmonic_abs_cell_bound(
                offsets[u] - offsets[v], cosine_difference,
                sine_difference, midpoint_sine, midpoint_cosine, half_width))
            velocity_range_bound = max(velocity_range_bound,
                _omrf_harmonic_abs_cell_bound(
                    0.0, sine_difference, -cosine_difference,
                    midpoint_sine, midpoint_cosine, half_width))
        end
        out[component] = 0.25 * range_bound * velocity_range_bound
    end
    return out
end

function _omrf_pattern_hcv_scale_cached(
        scales::OMRFPatternLocalHCVScales, component::Int)
    context = scales.context
    node = scales.component_nodes[component]
    q = context.seen[node] - 1
    iszero(q) && return 0.0
    displacements = scales.coordinate_offsets
    velocities = scales.coordinate_cosines
    field_displacement = 0.0
    field_velocity = 0.0
    @inbounds for k in eachindex(context.incident_edges[node])
        covariate = scales.component_covariates[
            component, context.incident_neighbours[node][k]]
        iszero(covariate) && continue
        index = context.interaction_indices[context.incident_edges[node][k]]
        field_displacement += covariate * displacements[index]
        field_velocity += covariate * velocities[index]
    end
    minimum_logit = 0.0
    maximum_logit = 0.0
    minimum_velocity = 0.0
    maximum_velocity = 0.0
    threshold_start = context.threshold_starts[node]
    @inbounds for u in 1:q
        index = context.threshold_indices[threshold_start + u - 1]
        logit = displacements[index] + u * field_displacement
        logit_velocity = velocities[index] + u * field_velocity
        minimum_logit = min(minimum_logit, logit)
        maximum_logit = max(maximum_logit, logit)
        minimum_velocity = min(minimum_velocity, logit_velocity)
        maximum_velocity = max(maximum_velocity, logit_velocity)
    end
    displacement_range = maximum_logit - minimum_logit
    velocity_range = maximum_velocity - minimum_velocity
    return scales.context.hcv_active ?
        velocity_range * displacement_range^2 / 8 :
        velocity_range * displacement_range / 4
end

function _omrf_pattern_hcv_scale_at(scales::OMRFPatternLocalHCVScales,
        state, flow, component::Int, anchor, t)
    context = scales.context
    node = scales.component_nodes[component]
    q = context.seen[node] - 1
    iszero(q) && return 0.0
    field_displacement = 0.0
    field_velocity = 0.0
    @inbounds for k in eachindex(context.incident_edges[node])
        covariate = scales.component_covariates[
            component, context.incident_neighbours[node][k]]
        iszero(covariate) && continue
        index = context.interaction_indices[context.incident_edges[node][k]]
        displacement, velocity = _omrf_coordinate_at(
            state, flow, anchor, index, t)
        field_displacement += covariate * displacement
        field_velocity += covariate * velocity
    end
    minimum_logit = 0.0
    maximum_logit = 0.0
    minimum_velocity = 0.0
    maximum_velocity = 0.0
    start = context.threshold_starts[node]
    @inbounds for u in 1:q
        index = context.threshold_indices[start + u - 1]
        displacement, velocity = _omrf_coordinate_at(
            state, flow, anchor, index, t)
        logit = displacement + u * field_displacement
        logit_velocity = velocity + u * field_velocity
        minimum_logit = min(minimum_logit, logit)
        maximum_logit = max(maximum_logit, logit)
        minimum_velocity = min(minimum_velocity, logit_velocity)
        maximum_velocity = max(maximum_velocity, logit_velocity)
    end
    displacement_range = maximum_logit - minimum_logit
    velocity_range = maximum_velocity - minimum_velocity
    return context.hcv_active ?
        velocity_range * displacement_range^2 / 8 :
        velocity_range * displacement_range / 4
end

function _omrf_pattern_hcv_cell_scale_at(scales::OMRFPatternLocalHCVScales,
        state, flow, component::Int, anchor, left, right)
    context = scales.context
    node = scales.component_nodes[component]
    q = context.seen[node] - 1
    iszero(q) && return 0.0
    midpoint_sine, midpoint_cosine = sincos((left + right) / 2)
    half_width = (right - left) / 2
    offsets = scales.offsets
    cosines = scales.cosines
    sines = scales.sines
    offsets[1] = 0.0
    cosines[1] = 0.0
    sines[1] = 0.0
    start = context.threshold_starts[node]
    @inbounds for u in 1:q
        index = context.threshold_indices[start + u - 1]
        offsets[u + 1], cosines[u + 1], sines[u + 1] =
            _omrf_harmonic_coefficients(state, flow, anchor, index)
    end
    @inbounds for k in eachindex(context.incident_edges[node])
        covariate = scales.component_covariates[
            component, context.incident_neighbours[node][k]]
        iszero(covariate) && continue
        index = context.interaction_indices[context.incident_edges[node][k]]
        offset, cosine, sine =
            _omrf_harmonic_coefficients(state, flow, anchor, index)
        for u in 1:q
            offsets[u + 1] += u * covariate * offset
            cosines[u + 1] += u * covariate * cosine
            sines[u + 1] += u * covariate * sine
        end
    end
    displacement_range = 0.0
    velocity_range = 0.0
    @inbounds for u in 1:(q + 1), v in 1:(u - 1)
        cosine_difference = cosines[u] - cosines[v]
        sine_difference = sines[u] - sines[v]
        displacement_range = max(displacement_range,
            _omrf_harmonic_abs_cell_bound(
                offsets[u] - offsets[v], cosine_difference, sine_difference,
                midpoint_sine, midpoint_cosine, half_width))
        velocity_range = max(velocity_range,
            _omrf_harmonic_abs_cell_bound(
                0.0, sine_difference, -cosine_difference,
                midpoint_sine, midpoint_cosine, half_width))
    end
    return context.hcv_active ?
        velocity_range * displacement_range^2 / 8 :
        velocity_range * displacement_range / 4
end

function _omrf_pattern_hcv_cell_scale_cached(
        scales::OMRFPatternLocalHCVScales, component::Int,
        midpoint_sine, midpoint_cosine, half_width)
    context = scales.context
    node = scales.component_nodes[component]
    q = context.seen[node] - 1
    iszero(q) && return 0.0
    offsets = scales.offsets
    cosines = scales.cosines
    sines = scales.sines
    offsets[1] = 0.0
    cosines[1] = 0.0
    sines[1] = 0.0
    start = context.threshold_starts[node]
    @inbounds for u in 1:q
        index = context.threshold_indices[start + u - 1]
        offsets[u + 1] = scales.coordinate_offsets[index]
        cosines[u + 1] = scales.coordinate_cosines[index]
        sines[u + 1] = scales.coordinate_sines[index]
    end
    @inbounds for k in eachindex(context.incident_edges[node])
        covariate = scales.component_covariates[
            component, context.incident_neighbours[node][k]]
        iszero(covariate) && continue
        index = context.interaction_indices[context.incident_edges[node][k]]
        for u in 1:q
            offsets[u + 1] += u * covariate * scales.coordinate_offsets[index]
            cosines[u + 1] += u * covariate * scales.coordinate_cosines[index]
            sines[u + 1] += u * covariate * scales.coordinate_sines[index]
        end
    end
    displacement_range = 0.0
    velocity_range = 0.0
    @inbounds for u in 1:(q + 1), v in 1:(u - 1)
        cosine_difference = cosines[u] - cosines[v]
        sine_difference = sines[u] - sines[v]
        displacement_range = max(displacement_range,
            _omrf_harmonic_abs_cell_bound(
                offsets[u] - offsets[v], cosine_difference, sine_difference,
                midpoint_sine, midpoint_cosine, half_width))
        velocity_range = max(velocity_range,
            _omrf_harmonic_abs_cell_bound(
                0.0, sine_difference, -cosine_difference,
                midpoint_sine, midpoint_cosine, half_width))
    end
    return context.hcv_active ?
        velocity_range * displacement_range^2 / 8 :
        velocity_range * displacement_range / 4
end

function (scales::OMRFPatternLocalHCVScales)(out, state,
        flow::PDMPSamplers.AnyBoomerang, t)
    started = time_ns()
    if scales.node_anchors !== nothing
        P = length(scales.context.seen)
        @inbounds for node in 1:P
            source = scales.node_sources[node]
            _cache_omrf_pattern_point_coordinates!(
                scales, state, flow, t, scales.node_anchors, source)
            for component in scales.node_component_ranges[node]
                out[component] = _omrf_pattern_hcv_scale_cached(
                    scales, component)
            end
        end
        scales.context.counts.envelope_point_calls += 1
        scales.context.counts.envelope_point_seconds +=
            (time_ns() - started) * 1e-9
        return out
    end
    _cache_omrf_pattern_point_coordinates!(scales, state, flow, t)
    @inbounds for component in eachindex(out, scales.component_nodes)
        out[component] = _omrf_pattern_hcv_scale_cached(scales, component)
    end
    scales.context.counts.envelope_point_calls += 1
    scales.context.counts.envelope_point_seconds +=
        (time_ns() - started) * 1e-9
    return out
end

function (scales::OMRFPatternLocalHCVScales)(out, state,
        flow::PDMPSamplers.AnyBoomerang, left, right)
    started = time_ns()
    context = scales.context
    if scales.node_anchors !== nothing
        midpoint_sine, midpoint_cosine = sincos((left + right) / 2)
        half_width = (right - left) / 2
        @inbounds for node in eachindex(context.seen)
            source = scales.node_sources[node]
            _cache_omrf_pattern_harmonic_coefficients!(
                scales, state, flow, scales.node_anchors, source)
            for component in scales.node_component_ranges[node]
                out[component] = _omrf_pattern_hcv_cell_scale_cached(
                    scales, component, midpoint_sine, midpoint_cosine,
                    half_width)
            end
        end
        context.counts.envelope_cell_calls += 1
        context.counts.envelope_cell_seconds +=
            (time_ns() - started) * 1e-9
        return out
    end
    _cache_omrf_pattern_harmonic_coefficients!(scales, state, flow)
    midpoint_sine, midpoint_cosine = sincos((left + right) / 2)
    half_width = (right - left) / 2
    coordinate_offsets = scales.coordinate_offsets
    coordinate_cosines = scales.coordinate_cosines
    coordinate_sines = scales.coordinate_sines
    @inbounds for component in eachindex(out, scales.component_nodes)
        node = scales.component_nodes[component]
        q = context.seen[node] - 1
        if iszero(q)
            out[component] = 0.0
            continue
        end
        offsets = scales.offsets
        cosines = scales.cosines
        sines = scales.sines
        offsets[1] = 0.0
        cosines[1] = 0.0
        sines[1] = 0.0
        threshold_start = context.threshold_starts[node]
        for u in 1:q
            index = context.threshold_indices[threshold_start + u - 1]
            offsets[u + 1] = coordinate_offsets[index]
            cosines[u + 1] = coordinate_cosines[index]
            sines[u + 1] = coordinate_sines[index]
        end
        for k in eachindex(context.incident_edges[node])
            covariate = scales.component_covariates[
                component, context.incident_neighbours[node][k]]
            iszero(covariate) && continue
            edge_index = context.interaction_indices[
                context.incident_edges[node][k]]
            offset = coordinate_offsets[edge_index]
            cosine = coordinate_cosines[edge_index]
            sine = coordinate_sines[edge_index]
            for u in 1:q
                offsets[u + 1] += u * covariate * offset
                cosines[u + 1] += u * covariate * cosine
                sines[u + 1] += u * covariate * sine
            end
        end
        displacement_range = 0.0
        velocity_range = 0.0
        for u in 1:(q + 1), v in 1:(u - 1)
            cosine_difference = cosines[u] - cosines[v]
            sine_difference = sines[u] - sines[v]
            displacement_range = max(displacement_range,
                _omrf_harmonic_abs_cell_bound(
                    offsets[u] - offsets[v], cosine_difference,
                    sine_difference, midpoint_sine, midpoint_cosine,
                    half_width))
            velocity_range = max(velocity_range,
                _omrf_harmonic_abs_cell_bound(
                    0.0, sine_difference, -cosine_difference,
                    midpoint_sine, midpoint_cosine, half_width))
        end
        out[component] = scales.context.hcv_active ?
            velocity_range * displacement_range^2 / 8 :
            velocity_range * displacement_range / 4
    end
    scales.context.counts.envelope_cell_calls += 1
    scales.context.counts.envelope_cell_seconds +=
        (time_ns() - started) * 1e-9
    return out
end

(scales::OMRFPatternLocalHCVScales)(out, state,
    flow::PDMPSamplers.PreconditionedDynamics, args...) =
    scales(out, state, flow.dynamics, args...)

@inline function _omrf_factor_anchor(scales::OMRFBatchedFactorHCVScales,
        index::Int, node::Int)
    anchors = scales.node_anchors
    return anchors === nothing ? scales.anchor[index] :
        anchors[index, scales.node_sources[node]]
end

@inline function _omrf_factor_anchor(scales::OMRFGroupedPersonHCVBoxScales,
        index::Int, node::Int)
    anchors = scales.node_anchors
    return anchors === nothing ? scales.anchor[index] :
        anchors[index, scales.node_sources[node]]
end

@inline function _omrf_grouped_anchor_value(
        scales::OMRFGroupedPersonHCVBoxScales, index::Int, node::Int)
    anchors = scales.node_anchors
    return anchors === nothing ? scales.anchor[index] :
        anchors[index, scales.node_sources[node]]
end

@inline function _omrf_grouped_anchor_value(
        scales::OMRFNodeGroupedPersonHCVBoxScales, index::Int, node::Int)
    anchors = scales.node_anchors
    return anchors === nothing ? scales.anchor[index] :
        anchors[index, scales.node_sources[node]]
end

function _omrf_grouped_box_node_bound_at(
        scales::Union{OMRFGroupedPersonHCVBoxScales,
            OMRFNodeGroupedPersonHCVBoxScales},
        state, flow, group::Int, node::Int, t::Real)
    return _omrf_grouped_box_node_bound_at_impl(
        scales, state, flow, group, node, t, nothing)
end

function _omrf_grouped_box_node_bound_at(
        scales::Union{OMRFGroupedPersonHCVBoxScales,
            OMRFNodeGroupedPersonHCVBoxScales},
        state, flow, group::Int, node::Int, t::Real, anchor)
    return _omrf_grouped_box_node_bound_at_impl(
        scales, state, flow, group, node, t, anchor)
end

function _omrf_grouped_box_node_bound_at_impl(
        scales::Union{OMRFGroupedPersonHCVBoxScales,
            OMRFNodeGroupedPersonHCVBoxScales},
        state, flow, group::Int, node::Int, t::Real, anchor_override)
    context = scales.context
    q = context.seen[node] - 1
    iszero(q) && return 0.0
    edges = context.incident_edges[node]
    neighbours = context.incident_neighbours[node]
    displacement_range = 0.0
    velocity_range = 0.0
    start = context.threshold_starts[node]
    @inbounds for u in 1:q, v in 0:(u - 1)
        u_index = context.threshold_indices[start + u - 1]
        u_anchor = anchor_override === nothing ?
            _omrf_grouped_anchor_value(scales, u_index, node) :
            anchor_override[u_index]
        u_displacement, u_velocity = _omrf_coordinate_at_anchor_value(
            state, flow, u_anchor, u_index, t)
        if iszero(v)
            v_displacement = 0.0
            v_velocity = 0.0
        else
            v_index = context.threshold_indices[start + v - 1]
            v_anchor = anchor_override === nothing ?
                _omrf_grouped_anchor_value(scales, v_index, node) :
                anchor_override[v_index]
            v_displacement, v_velocity = _omrf_coordinate_at_anchor_value(
                state, flow, v_anchor, v_index, t)
        end
        displacement_center = u_displacement - v_displacement
        velocity_center = u_velocity - v_velocity
        displacement_uncertainty = 0.0
        velocity_uncertainty = 0.0
        multiplier = u - v
        for k in eachindex(edges)
            index = context.interaction_indices[edges[k]]
            edge_anchor = anchor_override === nothing ?
                _omrf_grouped_anchor_value(scales, index, node) :
                anchor_override[index]
            displacement, velocity = _omrf_coordinate_at_anchor_value(
                state, flow, edge_anchor, index, t)
            lower = scales.minima[group, neighbours[k]]
            upper = scales.maxima[group, neighbours[k]]
            midpoint = (lower + upper) / 2
            radius = (upper - lower) / 2
            displacement_center += multiplier * midpoint * displacement
            velocity_center += multiplier * midpoint * velocity
            displacement_uncertainty += multiplier * radius * abs(displacement)
            velocity_uncertainty += multiplier * radius * abs(velocity)
        end
        displacement_range = max(displacement_range,
            abs(displacement_center) + displacement_uncertainty)
        velocity_range = max(velocity_range,
            abs(velocity_center) + velocity_uncertainty)
    end
    return context.hcv_active ?
        velocity_range * displacement_range^2 / 8 :
        velocity_range * displacement_range / 4
end

function _omrf_grouped_box_node_cell_bound(
        scales::OMRFGroupedPersonHCVBoxScales,
        state, flow::PDMPSamplers.AnyBoomerang,
        group::Int, node::Int, left::Real, right::Real)
    context = scales.context
    q = context.seen[node] - 1
    iszero(q) && return 0.0
    edges = context.incident_edges[node]
    neighbours = context.incident_neighbours[node]
    displacement_range = 0.0
    velocity_range = 0.0
    start = context.threshold_starts[node]
    @inbounds for u in 1:q, v in 0:(u - 1)
        u_index = context.threshold_indices[start + u - 1]
        u_offset, u_cosine, u_sine =
            _omrf_harmonic_coefficients_anchor_value(
                state, flow,
                _omrf_grouped_anchor_value(scales, u_index, node), u_index)
        if iszero(v)
            offset = u_offset
            cosine = u_cosine
            sine = u_sine
        else
            v_index = context.threshold_indices[start + v - 1]
            v_offset, v_cosine, v_sine =
                _omrf_harmonic_coefficients_anchor_value(
                    state, flow,
                    _omrf_grouped_anchor_value(scales, v_index, node), v_index)
            offset = u_offset - v_offset
            cosine = u_cosine - v_cosine
            sine = u_sine - v_sine
        end
        displacement_uncertainty = 0.0
        velocity_uncertainty = 0.0
        multiplier = u - v
        for k in eachindex(edges)
            index = context.interaction_indices[edges[k]]
            edge_offset, edge_cosine, edge_sine =
                _omrf_harmonic_coefficients_anchor_value(
                    state, flow,
                    _omrf_grouped_anchor_value(scales, index, node), index)
            lower = scales.minima[group, neighbours[k]]
            upper = scales.maxima[group, neighbours[k]]
            midpoint = (lower + upper) / 2
            radius = (upper - lower) / 2
            offset += multiplier * midpoint * edge_offset
            cosine += multiplier * midpoint * edge_cosine
            sine += multiplier * midpoint * edge_sine
            displacement_uncertainty += multiplier * radius *
                _omrf_harmonic_abs_cell_bound(
                    edge_offset, edge_cosine, edge_sine, left, right)
            velocity_uncertainty += multiplier * radius *
                _omrf_harmonic_abs_cell_bound(
                    0.0, edge_sine, -edge_cosine, left, right)
        end
        displacement_range = max(displacement_range,
            _omrf_harmonic_abs_cell_bound(
                offset, cosine, sine, left, right) +
            displacement_uncertainty)
        velocity_range = max(velocity_range,
            _omrf_harmonic_abs_cell_bound(
                0.0, sine, -cosine, left, right) + velocity_uncertainty)
    end
    return context.hcv_active ?
        velocity_range * displacement_range^2 / 8 :
        velocity_range * displacement_range / 4
end

function (scales::OMRFGroupedPersonHCVBoxScales)(out, state,
        flow::PDMPSamplers.AnyBoomerang, t)
    started = time_ns()
    @inbounds for group in eachindex(out)
        bound = 0.0
        for node in eachindex(scales.context.seen)
            bound += _omrf_grouped_box_node_bound_at(
                scales, state, flow, group, node, t)
        end
        out[group] = bound
    end
    scales.context.counts.envelope_point_calls += 1
    scales.context.counts.envelope_point_seconds +=
        (time_ns() - started) * 1e-9
    return out
end

function (scales::OMRFGroupedPersonHCVBoxScales)(out, state,
        flow::PDMPSamplers.AnyBoomerang, left, right)
    started = time_ns()
    @inbounds for group in eachindex(out)
        bound = 0.0
        for node in eachindex(scales.context.seen)
            bound += _omrf_grouped_box_node_cell_bound(
                scales, state, flow, group, node, left, right)
        end
        out[group] = bound
    end
    scales.context.counts.envelope_cell_calls += 1
    scales.context.counts.envelope_cell_seconds +=
        (time_ns() - started) * 1e-9
    return out
end

(scales::OMRFGroupedPersonHCVBoxScales)(out, state,
    flow::PDMPSamplers.PreconditionedDynamics, args...) =
    scales(out, state, flow.dynamics, args...)

function _cache_omrf_node_grouped_point!(
        scales::OMRFNodeGroupedPersonHCVBoxScales,
        state, flow::PDMPSamplers.AnyBoomerang, t::Real)
    context = scales.context
    edge_destination = 1
    @inbounds for node in eachindex(context.seen)
        q = context.seen[node] - 1
        pair_destination = scales.pair_offsets[node]
        start = context.threshold_starts[node]
        for u in 1:q, v in 0:(u - 1)
            u_index = context.threshold_indices[start + u - 1]
            u_displacement, u_velocity = _omrf_coordinate_at_anchor_value(
                state, flow,
                _omrf_grouped_anchor_value(scales, u_index, node), u_index, t)
            if iszero(v)
                v_displacement = 0.0
                v_velocity = 0.0
            else
                v_index = context.threshold_indices[start + v - 1]
                v_displacement, v_velocity = _omrf_coordinate_at_anchor_value(
                    state, flow,
                    _omrf_grouped_anchor_value(scales, v_index, node),
                    v_index, t)
            end
            scales.point_displacement_centers[pair_destination] =
                u_displacement - v_displacement
            scales.point_velocity_centers[pair_destination] =
                u_velocity - v_velocity
            pair_destination += 1
        end
        for edge in context.incident_edges[node]
            index = context.interaction_indices[edge]
            displacement, velocity = _omrf_coordinate_at_anchor_value(
                state, flow,
                _omrf_grouped_anchor_value(scales, index, node), index, t)
            scales.edge_displacements[edge_destination] = displacement
            scales.edge_velocities[edge_destination] = velocity
            edge_destination += 1
        end
    end
    return scales
end

function _cache_omrf_node_grouped_point_node!(
        scales::OMRFNodeGroupedPersonHCVBoxScales,
        state, flow::PDMPSamplers.AnyBoomerang, node::Int, anchor)
    context = scales.context
    q = context.seen[node] - 1
    pair_destination = scales.pair_offsets[node]
    start = context.threshold_starts[node]
    @inbounds for u in 1:q, v in 0:(u - 1)
        u_index = context.threshold_indices[start + u - 1]
        u_displacement, u_velocity = _omrf_coordinate_at_anchor_value(
            state, flow, anchor[u_index], u_index, 0.0)
        if iszero(v)
            v_displacement = 0.0
            v_velocity = 0.0
        else
            v_index = context.threshold_indices[start + v - 1]
            v_displacement, v_velocity = _omrf_coordinate_at_anchor_value(
                state, flow, anchor[v_index], v_index, 0.0)
        end
        scales.point_displacement_centers[pair_destination] =
            u_displacement - v_displacement
        scales.point_velocity_centers[pair_destination] =
            u_velocity - v_velocity
        pair_destination += 1
    end
    edge_destination = scales.edge_offsets_by_node[node] + 1
    @inbounds for edge in context.incident_edges[node]
        index = context.interaction_indices[edge]
        displacement, velocity = _omrf_coordinate_at_anchor_value(
            state, flow, anchor[index], index, 0.0)
        scales.edge_displacements[edge_destination] = displacement
        scales.edge_velocities[edge_destination] = velocity
        edge_destination += 1
    end
    return scales
end

@inline function _omrf_node_grouped_point_bound(
        scales::OMRFNodeGroupedPersonHCVBoxScales,
        component::Int, node::Int)
    context = scales.context
    q = context.seen[node] - 1
    iszero(q) && return 0.0
    edges = context.incident_edges[node]
    neighbours = context.incident_neighbours[node]
    pair_destination = scales.pair_offsets[node]
    edge_start = scales.edge_offsets_by_node[node]
    displacement_range = 0.0
    velocity_range = 0.0
    @inbounds for u in 1:q, v in 0:(u - 1)
        displacement_center = scales.point_displacement_centers[pair_destination]
        velocity_center = scales.point_velocity_centers[pair_destination]
        displacement_uncertainty = 0.0
        velocity_uncertainty = 0.0
        multiplier = u - v
        for k in eachindex(edges)
            lower = scales.minima[component, neighbours[k]]
            upper = scales.maxima[component, neighbours[k]]
            midpoint = (lower + upper) / 2
            radius = (upper - lower) / 2
            displacement = scales.edge_displacements[edge_start + k]
            velocity = scales.edge_velocities[edge_start + k]
            displacement_center += multiplier * midpoint * displacement
            velocity_center += multiplier * midpoint * velocity
            displacement_uncertainty += multiplier * radius * abs(displacement)
            velocity_uncertainty += multiplier * radius * abs(velocity)
        end
        displacement_range = max(displacement_range,
            abs(displacement_center) + displacement_uncertainty)
        velocity_range = max(velocity_range,
            abs(velocity_center) + velocity_uncertainty)
        pair_destination += 1
    end
    return velocity_range * displacement_range^2 / 8
end

"""
Return the pointwise HCV certificate for one observed person-node factor.

Unlike `_omrf_node_grouped_point_bound`, this uses the selected person's
actual covariates rather than the extrema of its setup-time support box.  It
is therefore both no larger than the corresponding box certificate and
requires no data traversal: a marked proposal reads only the selected row.
The shared trajectory geometry must already have been populated by
`_cache_omrf_node_grouped_point!`.
"""
@inline function _omrf_node_grouped_person_point_bound(
        scales::OMRFNodeGroupedPersonHCVBoxScales,
        person::Int, node::Int)
    context = scales.context
    q = context.seen[node] - 1
    iszero(q) && return 0.0
    edges = context.incident_edges[node]
    neighbours = context.incident_neighbours[node]
    pair_destination = scales.pair_offsets[node]
    edge_start = scales.edge_offsets_by_node[node]
    displacement_range = 0.0
    velocity_range = 0.0
    @inbounds for u in 1:q, v in 0:(u - 1)
        displacement = scales.point_displacement_centers[pair_destination]
        velocity = scales.point_velocity_centers[pair_destination]
        multiplier = u - v
        for k in eachindex(edges)
            covariate = context.X[person, neighbours[k]]
            iszero(covariate) && continue
            displacement += multiplier * covariate *
                scales.edge_displacements[edge_start + k]
            velocity += multiplier * covariate *
                scales.edge_velocities[edge_start + k]
        end
        displacement_range = max(displacement_range, abs(displacement))
        velocity_range = max(velocity_range, abs(velocity))
        pair_destination += 1
    end
    return velocity_range * displacement_range^2 / 8
end

function (scales::OMRFNodeGroupedPersonHCVBoxScales)(out, state,
        flow::PDMPSamplers.AnyBoomerang, t)
    started = time_ns()
    _cache_omrf_node_grouped_point!(scales, state, flow, t)
    @inbounds for component in eachindex(out)
        node = scales.component_nodes[component]
        out[component] = _omrf_node_grouped_point_bound(
            scales, component, node)
    end
    scales.context.counts.envelope_point_calls += 1
    scales.context.counts.envelope_point_seconds +=
        (time_ns() - started) * 1e-9
    return out
end

function _cache_omrf_node_grouped_cell!(
        scales::OMRFNodeGroupedPersonHCVBoxScales,
        state, flow::PDMPSamplers.AnyBoomerang,
        midpoint_sine::Real, midpoint_cosine::Real, half_width::Real)
    context = scales.context
    edge_destination = 1
    @inbounds for node in eachindex(context.seen)
        q = context.seen[node] - 1
        pair_destination = scales.pair_offsets[node]
        start = context.threshold_starts[node]
        for u in 1:q, v in 0:(u - 1)
            u_index = context.threshold_indices[start + u - 1]
            u_offset, u_cosine, u_sine =
                _omrf_harmonic_coefficients_anchor_value(
                    state, flow,
                    _omrf_grouped_anchor_value(scales, u_index, node), u_index)
            if iszero(v)
                offset, cosine, sine = u_offset, u_cosine, u_sine
            else
                v_index = context.threshold_indices[start + v - 1]
                v_offset, v_cosine, v_sine =
                    _omrf_harmonic_coefficients_anchor_value(
                        state, flow,
                        _omrf_grouped_anchor_value(scales, v_index, node),
                        v_index)
                offset = u_offset - v_offset
                cosine = u_cosine - v_cosine
                sine = u_sine - v_sine
            end
            scales.cell_offsets[pair_destination] = offset
            scales.cell_cosines[pair_destination] = cosine
            scales.cell_sines[pair_destination] = sine
            pair_destination += 1
        end
        for edge in context.incident_edges[node]
            index = context.interaction_indices[edge]
            offset, cosine, sine = _omrf_harmonic_coefficients_anchor_value(
                state, flow,
                _omrf_grouped_anchor_value(scales, index, node), index)
            scales.edge_offsets[edge_destination] = offset
            scales.edge_cosines[edge_destination] = cosine
            scales.edge_sines[edge_destination] = sine
            scales.edge_displacement_bounds[edge_destination] =
                _omrf_harmonic_abs_cell_bound(
                    offset, cosine, sine,
                    midpoint_sine, midpoint_cosine, half_width)
            scales.edge_velocity_bounds[edge_destination] =
                _omrf_harmonic_abs_cell_bound(
                    0.0, sine, -cosine,
                    midpoint_sine, midpoint_cosine, half_width)
            edge_destination += 1
        end
    end
    return scales
end

@inline function _omrf_node_grouped_cell_bound(
        scales::OMRFNodeGroupedPersonHCVBoxScales,
        component::Int, node::Int, midpoint_sine::Real,
        midpoint_cosine::Real, half_width::Real)
    context = scales.context
    q = context.seen[node] - 1
    iszero(q) && return 0.0
    edges = context.incident_edges[node]
    neighbours = context.incident_neighbours[node]
    pair_destination = scales.pair_offsets[node]
    edge_start = scales.edge_offsets_by_node[node]
    displacement_range = 0.0
    velocity_range = 0.0
    @inbounds for u in 1:q, v in 0:(u - 1)
        offset = scales.cell_offsets[pair_destination]
        cosine = scales.cell_cosines[pair_destination]
        sine = scales.cell_sines[pair_destination]
        displacement_uncertainty = 0.0
        velocity_uncertainty = 0.0
        multiplier = u - v
        for k in eachindex(edges)
            lower = scales.minima[component, neighbours[k]]
            upper = scales.maxima[component, neighbours[k]]
            midpoint = (lower + upper) / 2
            radius = (upper - lower) / 2
            destination = edge_start + k
            edge_offset = scales.edge_offsets[destination]
            edge_cosine = scales.edge_cosines[destination]
            edge_sine = scales.edge_sines[destination]
            offset += multiplier * midpoint * edge_offset
            cosine += multiplier * midpoint * edge_cosine
            sine += multiplier * midpoint * edge_sine
            displacement_uncertainty += multiplier * radius *
                scales.edge_displacement_bounds[destination]
            velocity_uncertainty += multiplier * radius *
                scales.edge_velocity_bounds[destination]
        end
        displacement_range = max(displacement_range,
            _omrf_harmonic_abs_cell_bound(
                offset, cosine, sine, midpoint_sine, midpoint_cosine,
                half_width) + displacement_uncertainty)
        velocity_range = max(velocity_range,
            _omrf_harmonic_abs_cell_bound(
                0.0, sine, -cosine, midpoint_sine, midpoint_cosine,
                half_width) + velocity_uncertainty)
        pair_destination += 1
    end
    return velocity_range * displacement_range^2 / 8
end

function (scales::OMRFNodeGroupedPersonHCVBoxScales)(out, state,
        flow::PDMPSamplers.AnyBoomerang, left, right)
    started = time_ns()
    midpoint = (left + right) / 2
    midpoint_sine, midpoint_cosine = sincos(midpoint)
    half_width = (right - left) / 2
    _cache_omrf_node_grouped_cell!(scales, state, flow,
        midpoint_sine, midpoint_cosine, half_width)
    @inbounds for component in eachindex(out)
        node = scales.component_nodes[component]
        out[component] = _omrf_node_grouped_cell_bound(
            scales, component, node, midpoint_sine, midpoint_cosine,
            half_width)
    end
    scales.context.counts.envelope_cell_calls += 1
    scales.context.counts.envelope_cell_seconds +=
        (time_ns() - started) * 1e-9
    return out
end

(scales::OMRFNodeGroupedPersonHCVBoxScales)(out, state,
    flow::PDMPSamplers.PreconditionedDynamics, args...) =
    scales(out, state, flow.dynamics, args...)

function _omrf_cell_cache_slot!(
        scales::OMRFCellCachedNodeGroupedPersonHCVBoxScales,
        left::Float64, right::Float64, n_components::Int)
    slot = scales.n_cached_cells + 1
    if slot > length(scales.cell_lefts)
        new_capacity = max(slot, 2 * length(scales.cell_lefts), 8)
        resize!(scales.cell_lefts, new_capacity)
        resize!(scales.cell_rights, new_capacity)
        old = scales.cell_scales
        replacement = zeros(n_components, new_capacity)
        isempty(old) || copyto!(view(replacement, :, 1:size(old, 2)), old)
        scales.cell_scales = replacement
    end
    scales.cell_lefts[slot] = left
    scales.cell_rights[slot] = right
    scales.n_cached_cells = slot
    return slot
end

function (scales::OMRFCellCachedNodeGroupedPersonHCVBoxScales)(out, state,
        flow::PDMPSamplers.AnyBoomerang, left, right)
    iszero(left) && (scales.n_cached_cells = 0)
    scales.base(out, state, flow, left, right)
    slot = _omrf_cell_cache_slot!(
        scales, Float64(left), Float64(right), length(out))
    copyto!(view(scales.cell_scales, :, slot), out)
    return out
end

function (scales::OMRFCellCachedNodeGroupedPersonHCVBoxScales)(out, state,
        flow::PDMPSamplers.AnyBoomerang, t)
    time = Float64(t)
    # The proposal time lies in a cell already appended to the active clock.
    # Search backwards because the newest cells are proposed from first.
    @inbounds for slot in scales.n_cached_cells:-1:1
        if scales.cell_lefts[slot] <= time <= scales.cell_rights[slot]
            copyto!(out, view(scales.cell_scales, :, slot))
            scales.cached_state = state
            scales.cached_flow = flow
            scales.cached_time = time
            return out
        end
    end
    # Defensive direct-API fallback. The event process always has a cached
    # cell, but keeping the callback total makes diagnostics less surprising.
    scales.cached_state = state
    scales.cached_flow = flow
    scales.cached_time = time
    return scales.base(out, state, flow, t)
end

(scales::OMRFCellCachedNodeGroupedPersonHCVBoxScales)(out, state,
    flow::PDMPSamplers.PreconditionedDynamics, args...) =
    scales(out, state, flow.dynamics, args...)

function PDMPSamplers.draw_subset!(rng::Random.AbstractRNG,
        cv::PDMPSamplers.SubsampledControlVariate{<:Any,<:Any,
            <:PDMPSamplers.GroupedResidualEnvelope{
                <:OMRFCellCachedNodeGroupedPersonHCVBoxScales}},
        D::Real, B::Real)
    envelope = cv.envelope
    scales = envelope.component_scales!
    started = time_ns()
    N = PDMPSamplers.n_observations(envelope)
    m = cv.m
    total = D + B
    distinguished = 0
    if !iszero(B) && (iszero(D) || rand(rng) * total > D)
        component = PDMPSamplers._draw_component(rng, envelope, B)
        distinguished = PDMPSamplers._draw_distinguished(
            rng, envelope, component)
        PDMPSamplers._draw_base_subset_conditional!(
            rng, cv, cv.subset_design, distinguished)
    else
        PDMPSamplers._draw_base_subset!(rng, cv, cv.subset_design)
    end

    cell_weight = 0.0
    base_flow = _omrf_lazy_base_flow(scales.cached_flow)
    @inbounds for observation in cv.subset
        for partition in axes(envelope.groups, 1)
            component = envelope.groups[partition, observation]
            cell_weight += envelope.scales[component]
        end
    end
    # Populate the shared node/category geometry only for the selected point;
    # no component loop is performed.
    _cache_omrf_node_grouped_point!(scales.base, scales.cached_state,
        base_flow, scales.cached_time)
    point_weight = 0.0
    @inbounds for observation in cv.subset
        for partition in axes(envelope.groups, 1)
            point_weight += _omrf_node_grouped_person_point_bound(
                scales.base, observation, partition)
        end
    end
    scales.selected_person_point_bound = (N / m) * point_weight
    oracle = cv.residual_oracle
    oracle isa PDMPSamplers.WithResidualStats && (oracle = oracle.f)
    oracle isa OMRFAnalyticSubsamplingOracle &&
        (oracle.lazy_residual_subset_bound = scales.selected_person_point_bound)
    scales.base.context.counts.cell_cached_point_screen_seconds +=
        (time_ns() - started) * 1e-9
    return D + (N / m) * cell_weight
end

@inline function _omrf_lazy_base_flow(flow)
    return flow isa PDMPSamplers.PreconditionedDynamics ? flow.dynamics : flow
end

function _cache_omrf_lazy_point!(
        scales::OMRFLazyNodeGroupedPersonHCVBoxScales,
        state, flow, t::Real)
    base_flow = _omrf_lazy_base_flow(flow)
    _cache_omrf_node_grouped_point!(scales.base, state, base_flow, t)
    scales.cached_state = state
    scales.cached_flow = flow
    scales.cached_time = Float64(t)
    scales.cache_valid = true
    return scales
end

@inline function _omrf_lazy_component_point_bound(
        scales::OMRFLazyNodeGroupedPersonHCVBoxScales, component::Int)
    node = scales.base.component_nodes[component]
    return _omrf_node_grouped_point_bound(scales.base, component, node)
end

function (scales::OMRFLazyNodeGroupedPersonHCVBoxScales)(out, state,
        flow::PDMPSamplers.AnyBoomerang, t)
    started = time_ns()
    _cache_omrf_lazy_point!(scales, state, flow, t)
    @inbounds for root_index in eachindex(scales.roots, out)
        out[root_index] = _omrf_lazy_component_point_bound(
            scales, scales.roots[root_index])
    end
    scales.base.context.counts.envelope_point_calls += 1
    scales.base.context.counts.envelope_point_seconds +=
        (time_ns() - started) * 1e-9
    return out
end

function (scales::OMRFLazyNodeGroupedPersonHCVBoxScales)(out, state,
        flow::PDMPSamplers.AnyBoomerang, left, right)
    started = time_ns()
    midpoint = (left + right) / 2
    midpoint_sine, midpoint_cosine = sincos(midpoint)
    half_width = (right - left) / 2
    _cache_omrf_node_grouped_cell!(scales.base, state, flow,
        midpoint_sine, midpoint_cosine, half_width)
    @inbounds for root_index in eachindex(scales.roots, out)
        component = scales.roots[root_index]
        node = scales.base.component_nodes[component]
        out[root_index] = _omrf_node_grouped_cell_bound(
            scales.base, component, node, midpoint_sine, midpoint_cosine,
            half_width)
    end
    scales.base.context.counts.envelope_cell_calls += 1
    scales.base.context.counts.envelope_cell_seconds +=
        (time_ns() - started) * 1e-9
    return out
end

(scales::OMRFLazyNodeGroupedPersonHCVBoxScales)(out, state,
    flow::PDMPSamplers.PreconditionedDynamics, args...) =
    scales(out, state, flow.dynamics, args...)

function PDMPSamplers.screening_residual_bound(
        envelope::PDMPSamplers.GroupedResidualEnvelope{
            <:OMRFLazyNodeGroupedPersonHCVBoxScales}, state, flow, t)
    scales = envelope.component_scales!
    scales(envelope.scales, state, flow, t)
    cumulative = 0.0
    @inbounds for root_index in eachindex(envelope.scales, scales.root_totals)
        cumulative += envelope.scales[root_index] *
            scales.root_totals[root_index]
        envelope.cumulative_masses[root_index] = cumulative
    end
    return cumulative
end

function PDMPSamplers.prepare_residual_sampling!(
        envelope::PDMPSamplers.GroupedResidualEnvelope{
            <:OMRFLazyNodeGroupedPersonHCVBoxScales}, state, flow, t)
    scales = envelope.component_scales!
    # `screening_residual_bound` was evaluated at the same candidate state
    # immediately before this hook. Keep a defensive fallback for direct API
    # use and future event-process changes.
    if !scales.cache_valid
        _cache_omrf_lazy_point!(scales, state, flow, t)
    end
    return nothing
end

function PDMPSamplers.draw_subset!(rng::Random.AbstractRNG,
        cv::PDMPSamplers.SubsampledControlVariate{<:Any,<:Any,
            <:PDMPSamplers.GroupedResidualEnvelope{
                <:OMRFLazyNodeGroupedPersonHCVBoxScales}},
        D::Real, B::Real)
    envelope = cv.envelope
    scales = envelope.component_scales!
    counts = scales.base.context.counts
    started = time_ns()
    N = PDMPSamplers.n_observations(envelope)
    m = cv.m
    total = D + B

    distinguished = 0
    hierarchy_rejected = false
    if !iszero(B) && (iszero(D) || rand(rng) * total > D)
        u = rand(rng) * B
        root_index = min(searchsortedlast(
            envelope.cumulative_masses, u) + 1,
            lastindex(envelope.cumulative_masses))
        component = scales.roots[root_index]
        parent_scale = envelope.scales[root_index]
        while !scales.leaves[component]
            left, right = scales.children[component]
            left_scale = _omrf_lazy_component_point_bound(scales, left)
            right_scale = _omrf_lazy_component_point_bound(scales, right)
            counts.lazy_branch_evaluations += 2
            left_mass = left_scale * length(scales.members[left])
            child_mass = left_mass +
                right_scale * length(scales.members[right])
            parent_mass = parent_scale * length(scales.members[component])
            child_mass <= parent_mass * (1 + 64eps(Float64)) ||
                throw(ErrorException(
                    "lazy node-support child mass exceeds its parent certificate"))
            if iszero(child_mass) || rand(rng) * parent_mass > child_mass
                counts.lazy_branch_rejections += 1
                hierarchy_rejected = true
                break
            end
            if rand(rng) * child_mass <= left_mass
                component = left
                parent_scale = left_scale
            else
                component = right
                parent_scale = right_scale
            end
        end
        if scales.leaves[component] && !iszero(parent_scale)
            distinguished = rand(rng, scales.members[component])
        end
    end

    if hierarchy_rejected
        PDMPSamplers._draw_base_subset!(rng, cv, cv.subset_design)
        oracle = cv.residual_oracle
        oracle isa PDMPSamplers.WithResidualStats && (oracle = oracle.f)
        oracle isa OMRFAnalyticSubsamplingOracle &&
            (oracle.lazy_residual_subset_bound = 0.0)
        counts.lazy_sampling_seconds += (time_ns() - started) * 1e-9
        return 0.0
    elseif iszero(distinguished)
        PDMPSamplers._draw_base_subset!(rng, cv, cv.subset_design)
    else
        PDMPSamplers._draw_base_subset_conditional!(
            rng, cv, cv.subset_design, distinguished)
    end

    subset_weight = 0.0
    @inbounds for observation in cv.subset
        bound = 0.0
        for partition in axes(scales.leaf_groups, 1)
            component = scales.leaf_groups[partition, observation]
            value = _omrf_lazy_component_point_bound(scales, component)
            counts.lazy_leaf_evaluations += 1
            bound += value
        end
        subset_weight += bound
        scales.person_leaf_bounds[observation] = bound
    end
    oracle = cv.residual_oracle
    if oracle isa PDMPSamplers.WithResidualStats
        oracle = oracle.f
    end
    oracle isa OMRFAnalyticSubsamplingOracle &&
        (oracle.lazy_residual_subset_bound = (N / m) * subset_weight)
    counts.lazy_sampling_seconds += (time_ns() - started) * 1e-9
    return D + (N / m) * subset_weight
end

function _cache_omrf_batched_harmonic_fields!(
        scales::OMRFBatchedFactorHCVScales, state,
        flow::PDMPSamplers.AnyBoomerang)
    context = scales.context
    fill!(scales.interaction_offsets, 0.0)
    fill!(scales.interaction_cosines, 0.0)
    fill!(scales.interaction_sines, 0.0)
    P = length(context.seen)
    edge = 1
    @inbounds for left in 1:(P - 1), right in (left + 1):P
        index = context.interaction_indices[edge]
        left_anchor = _omrf_factor_anchor(scales, index, left)
        right_anchor = _omrf_factor_anchor(scales, index, right)
        if _omrf_coordinate_free(state, index)
            cosine = state.ξ.x[index] - flow.μ[index]
            sine = state.ξ.θ[index]
            # `X * interaction_offsets` stores the field of `node` in column
            # `node`, so row `neighbour`, column `node` must use that node's
            # independently selected source.
            scales.interaction_offsets[right, left] = flow.μ[index] - left_anchor
            scales.interaction_offsets[left, right] = flow.μ[index] - right_anchor
            scales.interaction_cosines[left, right] = cosine
            scales.interaction_cosines[right, left] = cosine
            scales.interaction_sines[left, right] = sine
            scales.interaction_sines[right, left] = sine
        else
            scales.interaction_offsets[right, left] =
                state.ξ.x[index] - left_anchor
            scales.interaction_offsets[left, right] =
                state.ξ.x[index] - right_anchor
        end
        edge += 1
    end
    X = context.X_float
    mul!(scales.field_offsets, X, scales.interaction_offsets)
    mul!(scales.field_cosines, X, scales.interaction_cosines)
    mul!(scales.field_sines, X, scales.interaction_sines)
    return nothing
end

function (scales::OMRFBatchedFactorHCVScales)(out, state,
        flow::PDMPSamplers.AnyBoomerang, t)
    started = time_ns()
    _cache_omrf_batched_harmonic_fields!(scales, state, flow)
    sine, cosine = sincos(t)
    context = scales.context
    N, P = size(context.X)
    context.counts.event_bound_observation_visits += N
    destination = 1
    @inbounds for node in 1:P
        q = context.seen[node] - 1
        start = context.threshold_starts[node]
        for n in 1:N
            field_displacement = scales.field_offsets[n, node] +
                scales.field_cosines[n, node] * cosine +
                scales.field_sines[n, node] * sine
            field_velocity = -scales.field_cosines[n, node] * sine +
                scales.field_sines[n, node] * cosine
            minimum_displacement = 0.0
            maximum_displacement = 0.0
            minimum_velocity = 0.0
            maximum_velocity = 0.0
            for u in 1:q
                index = context.threshold_indices[start + u - 1]
                anchor = _omrf_factor_anchor(scales, index, node)
                if _omrf_coordinate_free(state, index)
                    centered = state.ξ.x[index] - flow.μ[index]
                    displacement = flow.μ[index] - anchor +
                        centered * cosine + state.ξ.θ[index] * sine
                    velocity = -centered * sine + state.ξ.θ[index] * cosine
                else
                    displacement = state.ξ.x[index] - anchor
                    velocity = 0.0
                end
                logit_displacement = displacement + u * field_displacement
                logit_velocity = velocity + u * field_velocity
                minimum_displacement = min(minimum_displacement, logit_displacement)
                maximum_displacement = max(maximum_displacement, logit_displacement)
                minimum_velocity = min(minimum_velocity, logit_velocity)
                maximum_velocity = max(maximum_velocity, logit_velocity)
            end
            displacement_range = maximum_displacement - minimum_displacement
            velocity_range = maximum_velocity - minimum_velocity
            out[destination] = context.hcv_active ?
                velocity_range * displacement_range^2 / 8 :
                velocity_range * displacement_range / 4
            destination += 1
        end
    end
    context.counts.envelope_point_calls += 1
    context.counts.envelope_point_seconds += (time_ns() - started) * 1e-9
    return out
end

function (scales::OMRFBatchedFactorHCVScales)(out, state,
        flow::PDMPSamplers.AnyBoomerang, left, right)
    started = time_ns()
    _cache_omrf_batched_harmonic_fields!(scales, state, flow)
    midpoint_sine, midpoint_cosine = sincos((left + right) / 2)
    half_width = (right - left) / 2
    context = scales.context
    N, P = size(context.X)
    context.counts.event_bound_observation_visits += N
    destination = 1
    @inbounds for node in 1:P
        q = context.seen[node] - 1
        start = context.threshold_starts[node]
        for n in 1:N
            displacement_range = 0.0
            velocity_range = 0.0
            for u in 1:q
                u_index = context.threshold_indices[start + u - 1]
                u_anchor = _omrf_factor_anchor(scales, u_index, node)
                if _omrf_coordinate_free(state, u_index)
                    u_offset = flow.μ[u_index] - u_anchor
                    u_cosine = state.ξ.x[u_index] - flow.μ[u_index]
                    u_sine = state.ξ.θ[u_index]
                else
                    u_offset = state.ξ.x[u_index] - u_anchor
                    u_cosine = 0.0
                    u_sine = 0.0
                end
                for v in 0:(u - 1)
                    if iszero(v)
                        v_offset = 0.0
                        v_cosine = 0.0
                        v_sine = 0.0
                    else
                        v_index = context.threshold_indices[start + v - 1]
                        v_anchor = _omrf_factor_anchor(scales, v_index, node)
                        if _omrf_coordinate_free(state, v_index)
                            v_offset = flow.μ[v_index] - v_anchor
                            v_cosine = state.ξ.x[v_index] - flow.μ[v_index]
                            v_sine = state.ξ.θ[v_index]
                        else
                            v_offset = state.ξ.x[v_index] - v_anchor
                            v_cosine = 0.0
                            v_sine = 0.0
                        end
                    end
                    multiplier = u - v
                    offset = u_offset - v_offset +
                        multiplier * scales.field_offsets[n, node]
                    cosine = u_cosine - v_cosine +
                        multiplier * scales.field_cosines[n, node]
                    sine = u_sine - v_sine +
                        multiplier * scales.field_sines[n, node]
                    displacement_range = max(displacement_range,
                        _omrf_harmonic_abs_cell_bound(offset, cosine, sine,
                            midpoint_sine, midpoint_cosine, half_width))
                    velocity_range = max(velocity_range,
                        _omrf_harmonic_abs_cell_bound(0.0, sine, -cosine,
                            midpoint_sine, midpoint_cosine, half_width))
                end
            end
            out[destination] = context.hcv_active ?
                velocity_range * displacement_range^2 / 8 :
                velocity_range * displacement_range / 4
            destination += 1
        end
    end
    context.counts.envelope_cell_calls += 1
    context.counts.envelope_cell_seconds += (time_ns() - started) * 1e-9
    return out
end

(scales::OMRFBatchedFactorHCVScales)(out, state,
    flow::PDMPSamplers.PreconditionedDynamics, args...) =
    scales(out, state, flow.dynamics, args...)

(scales::OMRFPatternLocalScales)(out, state,
    flow::PDMPSamplers.PreconditionedDynamics, args...) =
    scales(out, state, flow.dynamics, args...)

function _omrf_linear_node_geometry(state, anchor, support, t)
    displacement2 = 0.0
    velocity2 = 0.0
    @inbounds for j in support
        displacement = state.ξ.x[j] + t * state.ξ.θ[j] - anchor[j]
        displacement2 += displacement * displacement
        if _omrf_coordinate_free(state, j)
            velocity2 += state.ξ.θ[j]^2
        end
    end
    return sqrt(displacement2), sqrt(velocity2)
end

function _omrf_boomerang_node_geometry(flow, state, anchor, support, t)
    s, c = sincos(t)
    displacement2 = 0.0
    velocity2 = 0.0
    @inbounds for j in support
        if _omrf_coordinate_free(state, j)
            centered = state.ξ.x[j] - flow.μ[j]
            displacement = flow.μ[j] - anchor[j] +
                centered * c + state.ξ.θ[j] * s
            velocity = -centered * s + state.ξ.θ[j] * c
            displacement2 += displacement * displacement
            velocity2 += velocity * velocity
        else
            displacement = state.ξ.x[j] - anchor[j]
            displacement2 += displacement * displacement
        end
    end
    return sqrt(displacement2), sqrt(velocity2)
end

function _omrf_linear_node_cell_geometry(state, anchor, support, left, right)
    left_displacement2 = 0.0
    right_displacement2 = 0.0
    velocity2 = 0.0
    @inbounds for j in support
        left_displacement = state.ξ.x[j] + left * state.ξ.θ[j] - anchor[j]
        right_displacement = state.ξ.x[j] + right * state.ξ.θ[j] - anchor[j]
        left_displacement2 += left_displacement^2
        right_displacement2 += right_displacement^2
        if _omrf_coordinate_free(state, j)
            velocity2 += state.ξ.θ[j]^2
        end
    end
    return sqrt(max(left_displacement2, right_displacement2)), sqrt(velocity2)
end

function _omrf_boomerang_node_cell_geometry(
        flow, state, anchor, support, left, right)
    displacement2 = 0.0
    velocity2 = 0.0
    @inbounds for j in support
        if _omrf_coordinate_free(state, j)
            centered = state.ξ.x[j] - flow.μ[j]
            displacement = PDMPSamplers._harmonic_abs_max(
                flow.μ[j] - anchor[j], centered, state.ξ.θ[j], left, right)
            velocity = PDMPSamplers._harmonic_abs_max(
                0.0, state.ξ.θ[j], -centered, left, right)
            displacement2 += displacement^2
            velocity2 += velocity^2
        else
            displacement2 += (state.ξ.x[j] - anchor[j])^2
        end
    end
    return sqrt(displacement2), sqrt(velocity2)
end

_omrf_node_geometry(state, flow::Union{ZigZag,BouncyParticle}, anchor,
    support, t) = _omrf_linear_node_geometry(state, anchor, support, t)
_omrf_node_geometry(state, flow::PDMPSamplers.AnyBoomerang, anchor, support, t) =
    _omrf_boomerang_node_geometry(flow, state, anchor, support, t)
_omrf_node_geometry(state, flow::PDMPSamplers.PreconditionedDynamics, anchor, support, t) =
    PDMPSamplers.trajectory_geometry_bounds(flow, state, anchor, t)
_omrf_node_geometry(state, flow::ContinuousDynamics, anchor, support, t) =
    throw(ArgumentError("no OMRF node-local trajectory geometry for $(typeof(flow))"))

_omrf_node_cell_geometry(state, flow::Union{ZigZag,BouncyParticle}, anchor,
    support, left, right) = _omrf_linear_node_cell_geometry(
        state, anchor, support, left, right)
_omrf_node_cell_geometry(state, flow::PDMPSamplers.AnyBoomerang, anchor, support, left, right) =
    _omrf_boomerang_node_cell_geometry(
        flow, state, anchor, support, left, right)
_omrf_node_cell_geometry(state, flow::PDMPSamplers.PreconditionedDynamics, anchor,
    support, left, right) = PDMPSamplers.trajectory_geometry_cell_bounds(
        flow, state, anchor, left, right)
_omrf_node_cell_geometry(state, flow::ContinuousDynamics, anchor,
    support, left, right) = throw(ArgumentError(
        "no OMRF node-local trajectory geometry for $(typeof(flow))"))

function (scales::OMRFNodeLocalScales)(out, state, flow, t)
    @inbounds for node in eachindex(out, scales.supports)
        displacement, velocity = _omrf_node_geometry(
            state, flow, scales.anchor, scales.supports[node], t)
        out[node] = displacement * velocity
    end
    return out
end

function (scales::OMRFNodeLocalScales)(out, state, flow, left, right)
    @inbounds for node in eachindex(out, scales.supports)
        displacement, velocity = _omrf_node_cell_geometry(
            state, flow, scales.anchor, scales.supports[node], left, right)
        out[node] = displacement * velocity
    end
    return out
end

@inline function _fill_omrf_damped_hcv_scales!(out, node, P,
        displacement, velocity, α)
    out[node] = velocity * displacement * (1 - α)
    out[P + node] = velocity * displacement^2 * α
    return nothing
end

function (scales::OMRFDampedHCVNodeLocalScales)(out, state, flow, t)
    P = length(scales.supports)
    if !scales.context.hcv_active
        @inbounds for node in 1:P
            displacement, velocity = _omrf_node_geometry(
                state, flow, scales.anchor, scales.supports[node], t)
            out[node] = velocity * displacement
            out[P + node] = 0.0
        end
        return out
    end
    global_displacement, _ = PDMPSamplers.trajectory_geometry_bounds(
        flow, state, scales.anchor, t)
    damping = scales.context.hcv_damping
    α = damping / (damping + global_displacement^2)
    @inbounds for node in 1:P
        displacement, velocity = _omrf_node_geometry(
            state, flow, scales.anchor, scales.supports[node], t)
        _fill_omrf_damped_hcv_scales!(out, node, P,
            displacement, velocity, α)
    end
    return out
end

function _omrf_linear_displacement_lower_bound(state, anchor, left, right)
    velocity2 = 0.0
    displacement_velocity = 0.0
    @inbounds for j in eachindex(anchor)
        displacement = state.ξ.x[j] - anchor[j]
        velocity = state.ξ.θ[j]
        velocity2 += velocity^2
        displacement_velocity += displacement * velocity
    end
    minimizer = iszero(velocity2) ? left :
        clamp(-displacement_velocity / velocity2, left, right)
    displacement2 = 0.0
    @inbounds for j in eachindex(anchor)
        displacement = state.ξ.x[j] + minimizer * state.ξ.θ[j] - anchor[j]
        displacement2 += displacement^2
    end
    return sqrt(displacement2)
end

_omrf_displacement_lower_bound(state, flow::Union{ZigZag,BouncyParticle},
    anchor, left, right) =
    _omrf_linear_displacement_lower_bound(state, anchor, left, right)
_omrf_displacement_lower_bound(state, flow::ContinuousDynamics,
    anchor, left, right) = 0.0

function (scales::OMRFDampedHCVNodeLocalScales)(out, state, flow, left, right)
    P = length(scales.supports)
    if !scales.context.hcv_active
        @inbounds for node in 1:P
            displacement, velocity = _omrf_node_cell_geometry(
                state, flow, scales.anchor, scales.supports[node], left, right)
            out[node] = velocity * displacement
            out[P + node] = 0.0
        end
        return out
    end
    global_displacement, _ = PDMPSamplers.trajectory_geometry_cell_bounds(
        flow, state, scales.anchor, left, right)
    one_minus_α_bound = global_displacement^2 /
        (scales.context.hcv_damping + global_displacement^2)
    displacement_lower = _omrf_displacement_lower_bound(
        state, flow, scales.anchor, left, right)
    α_bound = scales.context.hcv_damping /
        (scales.context.hcv_damping + displacement_lower^2)
    @inbounds for node in 1:P
        displacement, velocity = _omrf_node_cell_geometry(
            state, flow, scales.anchor, scales.supports[node], left, right)
        # The two damping factors are bounded separately over the cell; their
        # extrema need not occur at the same time. Linear flows have an exact
        # minimum-distance calculation, while other flows conservatively use
        # the lower bound zero.
        out[node] = velocity * displacement * one_minus_α_bound
        out[P + node] = velocity * displacement^2 * α_bound
    end
    return out
end

function _omrf_node_supports(context::OMRFResidualContext)
    P = length(context.seen)
    supports = [Int[] for _ in 1:P]
    @inbounds for node in 1:P
        q = context.seen[node] - 1
        start = context.threshold_starts[node]
        for u in 1:q
            push!(supports[node], context.threshold_indices[start + u - 1])
        end
        for edge in context.incident_edges[node]
            push!(supports[node], context.interaction_indices[edge])
        end
        sort!(supports[node])
    end
    return supports
end

function _build_stan_residual_envelope(envelope_spec, anchor,
        oracle::OMRFAnalyticSubsamplingOracle)
    return _build_omrf_residual_envelope(
        envelope_spec, anchor, oracle.residual_context)
end

function _build_omrf_residual_envelope(envelope_spec, anchor,
        residual_context::OMRFResidualContext)
    bound_type = _haskey(envelope_spec, :bound_type) ?
        String(_rget(envelope_spec, :bound_type)) : ""
    if bound_type in ("node_grouped_person_hcv_box",
            "cell_cached_node_grouped_person_hcv_box",
            "lazy_node_grouped_person_hcv_box")
        residual_context.factorization === :person || throw(ArgumentError(
            "node-grouped HCV boxes require person factorization"))
        residual_context.use_hcv && isinf(residual_context.hcv_damping) ||
            throw(ArgumentError(
                "node-grouped HCV boxes require the undamped HCV"))
        groups = Matrix{Int}(_rget(envelope_spec, :node_box_membership))
        size(groups) == size(transpose(residual_context.X)) ||
            throw(DimensionMismatch(
                "node-grouped membership must contain one group per node/person"))
        component_nodes = _as_int_vector(
            _rget(envelope_spec, :node_box_nodes))
        node_component_ranges = [begin
            first_component = findfirst(==(node), component_nodes)
            last_component = findlast(==(node), component_nodes)
            first_component === nothing ? (1:0) :
                (first_component:last_component)
        end for node in eachindex(residual_context.seen)]
        @inbounds for node in eachindex(node_component_ranges)
            all(==(node), view(component_nodes,
                node_component_ranges[node])) || throw(ArgumentError(
                "node-grouped HCV components must be contiguous by node"))
        end
        minima = _as_float_matrix(_rget(envelope_spec, :node_box_minima))
        maxima = _as_float_matrix(_rget(envelope_spec, :node_box_maxima))
        n_components = length(component_nodes)
        size(minima) == size(maxima) ==
            (n_components, size(residual_context.X, 2)) ||
            throw(DimensionMismatch(
                "node-grouped HCV extrema must span every component/node"))
        all(group -> 1 <= group <= n_components, groups) ||
            throw(ArgumentError(
                "node-grouped HCV membership is out of range"))
        all(minima .<= maxima) || throw(ArgumentError(
            "node-grouped HCV minima must not exceed maxima"))
        pair_counts = [q * (q + 1) ÷ 2 for q in residual_context.seen .- 1]
        pair_offsets = cumsum(vcat(1, pair_counts[1:(end - 1)]))
        edge_counts = length.(residual_context.incident_edges)
        edge_offsets_by_node = cumsum(vcat(0, edge_counts[1:(end - 1)]))
        n_pairs = sum(pair_counts)
        n_edges = sum(edge_counts)
        base_scales = OMRFNodeGroupedPersonHCVBoxScales(
            anchor, nothing, Int[], component_nodes, node_component_ranges,
            minima, maxima,
            residual_context, pair_offsets, edge_offsets_by_node,
            zeros(n_pairs), zeros(n_pairs), zeros(n_edges), zeros(n_edges),
            zeros(n_pairs), zeros(n_pairs), zeros(n_pairs),
            zeros(n_edges), zeros(n_edges), zeros(n_edges),
            zeros(n_edges), zeros(n_edges))
        if bound_type == "node_grouped_person_hcv_box"
            return PDMPSamplers.GroupedResidualEnvelope(
                groups, n_components, base_scales;
                component_cell_scales! = base_scales)
        elseif bound_type == "cell_cached_node_grouped_person_hcv_box"
            cached_scales = OMRFCellCachedNodeGroupedPersonHCVBoxScales(
                base_scales, Float64[], Float64[],
                zeros(n_components, 0), 0, nothing, nothing, 0.0, NaN)
            return PDMPSamplers.GroupedResidualEnvelope(
                groups, n_components, cached_scales;
                component_cell_scales! = cached_scales)
        end
        parents = _as_int_vector(_rget(envelope_spec, :node_box_parents))
        roots = _as_int_vector(_rget(envelope_spec, :node_box_roots))
        leaf_components = _as_int_vector(
            _rget(envelope_spec, :node_box_leaves))
        length(parents) == n_components || throw(DimensionMismatch(
            "lazy node-grouped parents must span every component"))
        all(parent -> 0 <= parent <= n_components, parents) ||
            throw(ArgumentError("lazy node-grouped parent is out of range"))
        children_lists = [Int[] for _ in 1:n_components]
        @inbounds for child in eachindex(parents)
            parent = parents[child]
            iszero(parent) || push!(children_lists[parent], child)
        end
        leaves = falses(n_components)
        leaves[leaf_components] .= true
        children = Vector{NTuple{2,Int}}(undef, n_components)
        @inbounds for component in 1:n_components
            if leaves[component]
                isempty(children_lists[component]) || throw(ArgumentError(
                    "lazy node-grouped leaf cannot own children"))
                children[component] = (0, 0)
            else
                length(children_lists[component]) == 2 || throw(ArgumentError(
                    "lazy node-grouped internal component must own two children"))
                children[component] =
                    (children_lists[component][1], children_lists[component][2])
            end
        end
        all(root -> 1 <= root <= n_components, roots) ||
            throw(ArgumentError("lazy node-grouped root is out of range"))
        members = [Int[] for _ in 1:n_components]
        @inbounds for person in axes(groups, 2), node in axes(groups, 1)
            component = groups[node, person]
            push!(members[component], person)
            parent = parents[component]
            while !iszero(parent)
                push!(members[parent], person)
                parent = parents[parent]
            end
        end
        root_totals = Float64[length(members[root]) for root in roots]
        lazy_scales = OMRFLazyNodeGroupedPersonHCVBoxScales(
            base_scales, parents, children, roots, leaves, members,
            root_totals, zeros(length(roots)), zeros(length(roots)),
            zeros(n_components), groups, nothing, nothing, 0.0, false,
            zeros(size(groups, 2)))
        # Map every person/node leaf to its ancestor in the exposed root
        # frontier. This lets the generic envelope own the correct number of
        # observations and exact root multiplicities, while leaf membership
        # and deferred bounds remain in `lazy_scales`.
        root_lookup = zeros(Int, n_components)
        root_lookup[roots] .= eachindex(roots)
        root_groups = similar(groups)
        @inbounds for person in axes(groups, 2), node in axes(groups, 1)
            component = groups[node, person]
            while iszero(root_lookup[component])
                component = parents[component]
                iszero(component) && throw(ArgumentError(
                    "lazy node-grouped leaf has no exposed root ancestor"))
            end
            root_groups[node, person] = root_lookup[component]
        end
        return PDMPSamplers.GroupedResidualEnvelope(
            root_groups, length(roots), lazy_scales;
            component_cell_scales! = lazy_scales)
    end
    if bound_type == "grouped_person_hcv_box"
        residual_context.factorization === :person || throw(ArgumentError(
            "grouped-person HCV boxes require person factorization"))
        residual_context.use_hcv && isinf(residual_context.hcv_damping) ||
            throw(ArgumentError(
                "grouped-person HCV boxes require the undamped HCV"))
        membership = _as_int_vector(_rget(envelope_spec, :box_membership))
        length(membership) == size(residual_context.X, 1) ||
            throw(DimensionMismatch(
                "grouped-person HCV membership must contain one group per person"))
        minima = _as_float_matrix(_rget(envelope_spec, :box_minima))
        maxima = _as_float_matrix(_rget(envelope_spec, :box_maxima))
        size(minima) == size(maxima) || throw(DimensionMismatch(
            "grouped-person HCV box extrema must have matching dimensions"))
        size(minima, 2) == size(residual_context.X, 2) ||
            throw(DimensionMismatch(
                "grouped-person HCV boxes must span every OMRF node"))
        group_count = size(minima, 1)
        all(group -> 1 <= group <= group_count, membership) ||
            throw(ArgumentError("grouped-person HCV membership is out of range"))
        all(minima .<= maxima) || throw(ArgumentError(
            "grouped-person HCV box minima must not exceed maxima"))
        groups = reshape(membership, 1, :)
        scales = OMRFGroupedPersonHCVBoxScales(anchor, nothing, Int[],
            minima, maxima, residual_context)
        return PDMPSamplers.GroupedResidualEnvelope(
            groups, group_count, scales; component_cell_scales! = scales)
    end
    if _haskey(envelope_spec, :node_weights)
        bound_type = isempty(bound_type) ? "node_local_norm_low_rank" : bound_type
        if bound_type == "pattern_local_hcv_range" &&
                _haskey(envelope_spec, :structural_groups) &&
                residual_context.use_hcv
            isinf(residual_context.hcv_damping) || throw(ArgumentError(
                "pattern-local HCV currently requires exact undamped damping Inf"))
            if residual_context.factorization === :person_node
                N, P = size(residual_context.X)
                groups = reshape(collect(1:(N * P)), 1, N * P)
                scales = OMRFBatchedFactorHCVScales(
                    anchor, nothing, Int[], residual_context,
                    zeros(P, P), zeros(P, P), zeros(P, P),
                    zeros(N, P), zeros(N, P), zeros(N, P))
                return PDMPSamplers.GroupedResidualEnvelope(
                    groups, N * P, scales;
                    component_cell_scales! = scales)
            end
            groups = Matrix{Int}(_rget(envelope_spec, :structural_groups))
            n_components = Int(_rget(envelope_spec, :n_structural_components))
            component_nodes = _as_int_vector(
                _rget(envelope_spec, :component_nodes))
            component_covariates = Matrix{Int}(
                _rget(envelope_spec, :component_covariates))
            length(component_nodes) == n_components ||
                throw(DimensionMismatch(
                    "OMRF structural components must match envelope rows"))
            size(component_covariates) ==
                (length(component_nodes), size(residual_context.X, 2)) ||
                throw(DimensionMismatch(
                    "OMRF structural covariates have the wrong dimensions"))
            workspace_length = maximum(residual_context.seen)
            node_component_ranges = [begin
                first_component = findfirst(==(node), component_nodes)
                last_component = findlast(==(node), component_nodes)
                first_component === nothing ? (1:0) :
                    (first_component:last_component)
            end for node in eachindex(residual_context.seen)]
            scales = OMRFPatternLocalHCVScales(anchor, nothing, Int[],
                component_nodes, node_component_ranges,
                component_covariates, residual_context,
                zeros(workspace_length), zeros(workspace_length),
                zeros(workspace_length), zeros(length(anchor)),
                zeros(length(anchor)), zeros(length(anchor)))
            return PDMPSamplers.GroupedResidualEnvelope(
                groups, n_components, scales;
                component_cell_scales! = scales)
        end
        if bound_type == "pattern_local_range" &&
                _haskey(envelope_spec, :structural_groups) &&
                !residual_context.use_hcv
            groups = Matrix{Int}(_rget(envelope_spec, :structural_groups))
            n_components = Int(_rget(envelope_spec, :n_structural_components))
            component_nodes = _as_int_vector(
                _rget(envelope_spec, :component_nodes))
            component_covariates = Matrix{Int}(
                _rget(envelope_spec, :component_covariates))
            length(component_nodes) == n_components ||
                throw(DimensionMismatch(
                    "OMRF structural components must match envelope rows"))
            size(component_covariates) ==
                (length(component_nodes), size(residual_context.X, 2)) ||
                throw(DimensionMismatch(
                    "OMRF structural covariates have the wrong dimensions"))
            workspace_length = maximum(residual_context.seen)
            scales = OMRFPatternLocalScales(anchor, component_nodes,
                component_covariates, residual_context,
                zeros(workspace_length), zeros(workspace_length),
                zeros(workspace_length), zeros(length(anchor)),
                zeros(length(anchor)), zeros(length(anchor)))
            return PDMPSamplers.GroupedResidualEnvelope(
                groups, n_components, scales;
                component_cell_scales! = scales)
        end
        if bound_type == "covariate_local_expansion" &&
                _haskey(envelope_spec, :covariate_weights) &&
                !residual_context.use_hcv
            weights = _as_float_matrix(
                _rget(envelope_spec, :covariate_weights))
            expected_rows = sum((length(edges) + 1)^2
                for edges in residual_context.incident_edges)
            size(weights) == (expected_rows, size(residual_context.X, 1)) ||
                throw(DimensionMismatch(
                    "OMRF covariate-local weights have the wrong dimensions"))
            max_degree = maximum(length, residual_context.incident_edges;
                init=0)
            scales = OMRFCovariateLocalScales(anchor, residual_context,
                zeros(max_degree), zeros(max_degree), zeros(length(anchor)),
                zeros(length(anchor)), zeros(length(anchor)))
            if residual_context.factorization === :person_node
                component_blocks = _as_int_vector(
                    _rget(envelope_spec, :covariate_component_blocks))
                return PDMPSamplers.BlockSeparableResidualEnvelope(
                    weights, component_blocks, size(residual_context.X, 2),
                    scales; component_cell_scales! = scales)
            end
            return SeparableResidualEnvelope(weights, scales;
                component_cell_scales! = scales)
        end
        if _haskey(envelope_spec, :norm_weights) &&
                !residual_context.use_hcv
            weights = _as_float_matrix(_rget(envelope_spec, :norm_weights))
            size(weights, 1) == 3 * size(residual_context.X, 2) ||
                throw(DimensionMismatch(
                    "OMRF norm-local weights must have three rows per node"))
            scales = OMRFNormLocalScales(anchor, residual_context)
            return SeparableResidualEnvelope(weights, scales;
                component_cell_scales! = scales)
        end
        weights = _as_float_matrix(_rget(envelope_spec, :node_weights))
        supports = _omrf_node_supports(residual_context)
        if residual_context.use_hcv
            remainder = _as_float_matrix(
                _rget(envelope_spec, :hcv_remainder_weights))
            size(remainder) == size(weights) || throw(DimensionMismatch(
                "OMRF HCV remainder weights must match node weights"))
            scales = OMRFDampedHCVNodeLocalScales(
                anchor, supports, residual_context)
            return SeparableResidualEnvelope(vcat(weights, remainder), scales;
                component_cell_scales! = scales)
        end
        scales = OMRFNodeLocalScales(
            anchor, supports)
        return SeparableResidualEnvelope(weights, scales;
            component_cell_scales! = scales)
    end
    weights = _as_float_matrix(_rget(envelope_spec, :weights))
    return TrajectoryResidualEnvelope(weights, anchor; growth_rates=growth)
end

function _build_stan_residual_envelope(envelope_spec, anchor,
        oracle::StanSubsamplingOracle)
    weights = _as_float_matrix(_rget(envelope_spec, :weights))
    growth = _as_float_vector(_rget(envelope_spec, :growth_rates))
    return TrajectoryResidualEnvelope(weights, anchor; growth_rates=growth)
end

"""Atomically installable OMRF control-variate, cache, and envelope state."""
struct OMRFPreparedAnchorState
    anchor::Vector{Float64}
    full_anchor::Vector{Float64}
    prior_anchor::Vector{Float64}
    residual_context::OMRFResidualContext
    envelope::PDMPSamplers.AbstractResidualEnvelope
end

mutable struct OMRFAnchorEntry
    state::OMRFPreparedAnchorState
    age::Int
end

"""Chain-local LRU bank of fully prepared OMRF anchor states."""
mutable struct OMRFAnchorManager
    oracle::OMRFAnalyticSubsamplingOracle
    envelope_spec::Any
    unc_names::Vector{String}
    entries::Vector{OMRFAnchorEntry}
    active_idx::Int
    capacity::Int
    proposal::Vector{Float64}
    distance_scales::Vector{Float64}
    preparations::Int
    preparation_seconds::Float64
    activations::Int
    main_activations::Int
    main_refreshes::Int
    main_selections::Int
    main_distance_sum::Float64
    max_main_distance::Float64
    main_refresh_distance::Float64
    hcv_active::Bool
    recycled_state::Union{Nothing,OMRFPreparedAnchorState}
    recycled_preparations::Int
    warmup_finish_preparations::Int
    warmup_finish_native_full_data_scans::Int
    warmup_finish_event_bound_observation_visits::Int
    selection_seconds::Float64
    selection_component_evaluations::Int
    node_source_check_distance::Float64
    node_source_last_position::Vector{Float64}
    node_source_checks::Int
    node_source_skips::Int
    node_source_changes::Int
end

function _prepare_omrf_anchor_envelope(manager::OMRFAnchorManager,
        requested, residual_context)
    spec = manager.envelope_spec
    bound_type = _haskey(spec, :bound_type) ?
        String(_rget(spec, :bound_type)) : ""
    if bound_type == "covariate_local_expansion" &&
            !residual_context.use_hcv && !isempty(manager.entries)
        template = manager.entries[1].state.envelope
        if template isa Union{SeparableResidualEnvelope,
                PDMPSamplers.BlockSeparableResidualEnvelope}
            max_degree = maximum(length,
                residual_context.incident_edges; init=0)
            scales = OMRFCovariateLocalScales(requested, residual_context,
                zeros(max_degree), zeros(max_degree), zeros(length(requested)),
                zeros(length(requested)), zeros(length(requested)))
            # Weights, totals, and alias tables depend only on the observed
            # covariates. Anchor-bank entries share them; only trajectory
            # callbacks and their workspaces are anchor-specific.
            if template isa PDMPSamplers.BlockSeparableResidualEnvelope
                return PDMPSamplers.BlockSeparableResidualEnvelope(
                    template.weights, template.component_blocks,
                    template.block_components, template.n_blocks, scales,
                    scales, template.totals, template.alias_tables,
                    zeros(length(template.scales)),
                    zeros(length(template.cell_scales)),
                    zeros(length(template.cumulative_masses)))
            end
            return SeparableResidualEnvelope(template.weights, scales, scales,
                template.totals, template.alias_tables,
                zeros(length(template.scales)),
                zeros(length(template.cell_scales)),
                zeros(length(template.cumulative_masses)))
        end
    end
    return _build_omrf_residual_envelope(spec, requested, residual_context)
end

function _prepare_omrf_anchor(manager::OMRFAnchorManager,
        anchor::AbstractVector)
    started = time_ns()
    oracle = manager.oracle
    recycled = manager.recycled_state
    if recycled !== nothing && oracle.analytic_prior !== nothing &&
            !manager.hcv_active && !recycled.residual_context.use_hcv &&
            recycled.envelope isa Union{SeparableResidualEnvelope,
                PDMPSamplers.BlockSeparableResidualEnvelope} &&
            recycled.envelope.component_scales! isa OMRFCovariateLocalScales
        # The recycled state is not installed in the bank, so it can be
        # mutated without exposing a partially prepared anchor. Its envelope
        # callbacks retain the same anchor/context array identities.
        manager.recycled_state = nothing
        copyto!(recycled.anchor, anchor)
        scales = recycled.envelope.component_scales!
        scales.anchor === recycled.anchor || copyto!(scales.anchor, anchor)
        recycled.residual_context.hcv_active = false
        _omrf_prior_gradient!(
            recycled.prior_anchor, oracle.analytic_prior, recycled.anchor)
        _cache_omrf_anchor_and_full_gradient!(recycled.residual_context,
            recycled.full_anchor, oracle.analytic_prior, recycled.anchor)
        oracle.ctx.counts.full_gradient += 1
        manager.preparations += 1
        manager.recycled_preparations += 1
        manager.preparation_seconds += (time_ns() - started) / 1e9
        return recycled
    end
    envelope_spec = manager.envelope_spec
    requested = collect(Float64, anchor)
    residual_context = _omrf_residual_context(envelope_spec,
        manager.unc_names, requested, oracle.ctx.counts)
    residual_context.hcv_active = manager.hcv_active &&
        residual_context.use_hcv
    full_anchor = zeros(length(requested))
    prior_anchor = zeros(length(requested))
    _clear_gradient!(oracle.ctx, :full, full_anchor, requested)
    if oracle.analytic_prior === nothing
        _clear_gradient!(oracle.ctx, :prior, prior_anchor, requested)
    else
        _omrf_prior_gradient!(prior_anchor, oracle.analytic_prior, requested)
    end
    envelope = _prepare_omrf_anchor_envelope(
        manager, requested, residual_context)
    manager.preparations += 1
    manager.preparation_seconds += (time_ns() - started) / 1e9
    return OMRFPreparedAnchorState(
        requested, full_anchor, prior_anchor, residual_context, envelope)
end

function _install_omrf_anchor!(manager::OMRFAnchorManager,
        state::OMRFPreparedAnchorState)
    oracle = manager.oracle
    oracle.anchor = state.anchor
    oracle.full_anchor = state.full_anchor
    oracle.prior_anchor = state.prior_anchor
    oracle.residual_context = state.residual_context
    oracle.prior_x_valid = false
    return state.envelope
end

function _insert_omrf_anchor!(manager::OMRFAnchorManager,
        state::OMRFPreparedAnchorState)
    idx = if length(manager.entries) < manager.capacity
        push!(manager.entries, OMRFAnchorEntry(state, 0))
        length(manager.entries)
    else
        replace_idx = argmax(map(entry -> entry.age, manager.entries))
        replaced = manager.entries[replace_idx].state
        manager.entries[replace_idx] = OMRFAnchorEntry(state, 0)
        manager.recycled_state = replaced
        replace_idx
    end
    return idx
end

function _new_omrf_anchor_manager(oracle::OMRFAnalyticSubsamplingOracle,
        envelope_spec, unc_names, initial_envelope, capacity::Integer)
    capacity >= 1 || throw(ArgumentError("OMRF anchor-bank capacity must be positive"))
    manager = OMRFAnchorManager(oracle, envelope_spec, String.(unc_names),
        OMRFAnchorEntry[], 0, Int(capacity), zeros(length(oracle.full_anchor)),
        ones(length(oracle.full_anchor)),
        0, 0.0, 0, 0, 0, 0, 0.0, 0.0, Inf,
        oracle.residual_context.hcv_active, nothing, 0, 0,
        oracle.ctx.counts.native_full_data_scans,
        oracle.ctx.counts.event_bound_observation_visits, 0.0, 0,
        0.0, copy(oracle.anchor), 0, 0, 0)
    initial = OMRFPreparedAnchorState(oracle.anchor, oracle.full_anchor,
        oracle.prior_anchor, oracle.residual_context, initial_envelope)
    manager.active_idx = _insert_omrf_anchor!(manager, initial)
    return manager
end

function _omrf_warmup_candidates(trace::PDMPSamplers.AbstractPDMPTrace,
        position::AbstractVector, capacity::Integer;
        retained_fraction::Real=0.5)
    pdmp_trace = trace isa PDMPSamplers.PDMPTrace ? trace :
        PDMPSamplers.PDMPTrace(trace)
    n_events = length(pdmp_trace.times)
    n_events >= 2 || return reshape(collect(Float64, position), :, 1)
    0 < retained_fraction <= 1 || throw(ArgumentError(
        "warmup candidate retained fraction must lie in (0, 1]"))
    n_grid = min(512, max(8 * Int(capacity), Int(capacity)))
    final_time = last(pdmp_trace.times)
    initial_time = first(pdmp_trace.times)
    pool_start = final_time - retained_fraction * (final_time - initial_time)
    grid = range(pool_start, final_time; length=n_grid)
    candidates = Matrix{Float64}(undef, length(position), n_grid + 1)
    copyto!(view(candidates, :, 1), position)
    for (offset, time) in enumerate(grid)
        column = offset + 1
        idx = clamp(searchsortedlast(pdmp_trace.times, time), 1, n_events)
        copyto!(view(candidates, :, column), view(pdmp_trace.positions, :, idx))
        if idx < n_events && time > pdmp_trace.times[idx]
            velocity = collect(Float64, view(pdmp_trace.velocities, :, idx))
            point = PDMPSamplers.SkeletonPoint(view(candidates, :, column), velocity)
            PDMPSamplers.move_forward_time!(
                point, time - pdmp_trace.times[idx], pdmp_trace.flow)
        end
    end
    return candidates
end

function _omrf_representative_anchors(trace::PDMPSamplers.AbstractPDMPTrace,
        position::AbstractVector, capacity::Integer)
    # The early warmup transient is useful for adapting the flow but harmful
    # as a frozen retained-phase reference. Use the latter half of the same
    # online warmup path; no optimization or extra model evaluation is added.
    candidates = _omrf_warmup_candidates(
        trace, position, capacity; retained_fraction=0.5)
    d, n = size(candidates)
    means = vec(sum(candidates; dims=2)) ./ n
    scales = zeros(d)
    if n > 1
        @inbounds for column in axes(candidates, 2), j in axes(candidates, 1)
            delta = candidates[j, column] - means[j]
            scales[j] += delta * delta
        end
        scales .= sqrt.(scales ./ (n - 1))
    end
    positive = filter(scale -> scale > sqrt(eps(Float64)), scales)
    scale_floor = isempty(positive) ? 1.0 : 0.1 * Statistics.median(positive)
    @. scales = max(scales, scale_floor)

    k = min(Int(capacity), n)
    selected = Matrix{Float64}(undef, d, k)
    # The first anchor exactly covers the warmup/main phase boundary. Greedy
    # farthest-point additions then cover the time-uniform warmup occupation
    # measure rather than the much more event-dense portions of the skeleton.
    copyto!(view(selected, :, 1), view(candidates, :, 1))
    minimum_distance2 = fill(Inf, n)
    for selection in 1:k
        if selection > 1
            idx = argmax(minimum_distance2)
            copyto!(view(selected, :, selection), view(candidates, :, idx))
        end
        @inbounds for column in axes(candidates, 2)
            distance2 = 0.0
            for j in axes(candidates, 1)
                delta = (candidates[j, column] - selected[j, selection]) /
                    scales[j]
                distance2 += delta * delta
            end
            minimum_distance2[column] = min(
                minimum_distance2[column], distance2)
        end
    end
    return selected, scales
end

"""
Learn warmup-bank source labels per node against the actual retained-person
certificate.

The score is a time-uniform average over the latter half of the online warmup
trajectory. For every sampled warmup state it sums the exact observed-pattern
HCV certificate over people. All `N * P * K` work is performed once at the
warmup/main boundary and is included in the warmup observation-visit counter.
Retained sampling only looks up these labels by node-local state geometry; it
never evaluates a person certificate or visits data rows.
"""
function _select_omrf_warmup_pattern_sources!(manager::OMRFAnchorManager,
        cv::SubsampledControlVariate, trace::PDMPSamplers.AbstractPDMPTrace,
        state, flow; retained_fraction::Real=0.5)
    started = time_ns()
    bank = something(manager.oracle.node_local_hcv)
    pdmp_trace = trace isa PDMPSamplers.PDMPTrace ? trace :
        PDMPSamplers.PDMPTrace(trace)
    N, P = size(bank.contexts[1].X)
    K = size(bank.candidates, 2)
    n_events = length(pdmp_trace.times)
    n_events >= 2 || return _select_omrf_node_local_sources!(
        manager, cv, state, flow)
    0 < retained_fraction <= 1 || throw(ArgumentError(
        "warmup scoring retained fraction must lie in (0, 1]"))
    n_grid = min(128, max(8K, 32))
    final_time = last(pdmp_trace.times)
    initial_time = first(pdmp_trace.times)
    grid_start = final_time - retained_fraction *
        (final_time - initial_time)
    scores = zeros(P, K)
    state_scores = zeros(P, K)
    selector_positions = Matrix{Float64}(
        undef, length(state.ξ.x), n_grid)
    selector_velocities = similar(selector_positions)
    selector_sources = Matrix{Int}(undef, P, n_grid)
    position = zeros(length(state.ξ.x))
    velocity = similar(position)
    scoring_state = PDMPSamplers.PDMPState(
        0.0, PDMPSamplers.SkeletonPoint(position, velocity))
    envelope_scales = cv.envelope.component_scales!
    node_scales = envelope_scales isa
        OMRFCellCachedNodeGroupedPersonHCVBoxScales ?
        envelope_scales.base :
        (envelope_scales isa OMRFNodeGroupedPersonHCVBoxScales ?
            envelope_scales : nothing)
    @inbounds for (grid_index, time) in enumerate(
            range(grid_start, final_time; length=n_grid))
        idx = clamp(searchsortedlast(pdmp_trace.times, time), 1, n_events)
        copyto!(position, view(pdmp_trace.positions, :, idx))
        copyto!(velocity, view(pdmp_trace.velocities, :, idx))
        if idx < n_events && time > pdmp_trace.times[idx]
            PDMPSamplers.move_forward_time!(
                scoring_state.ξ, time - pdmp_trace.times[idx], flow)
        end
        copyto!(view(selector_positions, :, grid_index), position)
        copyto!(view(selector_velocities, :, grid_index), velocity)
        fill!(state_scores, 0.0)
        for source in 1:K
            context = bank.contexts[source]
            anchor = view(bank.candidates, :, source)
            for node in 1:P
                if node_scales === nothing
                    for person in 1:N
                        value = _omrf_node_hcv_bound_at(
                            context, anchor, scoring_state, flow, person, node)
                        scores[node, source] += value
                        state_scores[node, source] += value
                    end
                else
                    # Reuse the node/category trajectory geometry across all
                    # observed rows, just as retained selected-person screens
                    # do. This is exact but avoids repeating coordinate work.
                    _cache_omrf_node_grouped_point_node!(
                        node_scales, scoring_state, flow, node, anchor)
                    for person in 1:N
                        value = _omrf_node_grouped_person_point_bound(
                            node_scales, person, node)
                        scores[node, source] += value
                        state_scores[node, source] += value
                    end
                end
            end
        end
        for node in 1:P
            selector_sources[node, grid_index] =
                argmin(view(state_scores, node, :))
        end
    end
    requested = bank.requested_sources
    @inbounds for node in 1:P
        requested[node] = argmin(view(scores, node, :))
    end
    changed = _install_omrf_node_local_sources!(manager, cv, requested)
    bank.selector_positions = selector_positions
    bank.selector_velocities = selector_velocities
    bank.selector_sources = selector_sources
    manager.node_source_checks += 1
    changed && (manager.node_source_changes += 1)
    evaluations = n_grid * N * P * K
    manager.selection_component_evaluations += evaluations
    manager.oracle.ctx.counts.event_bound_observation_visits += evaluations
    manager.selection_seconds += (time_ns() - started) * 1e-9
    return changed
end

function _select_omrf_warmup_mapped_sources!(manager::OMRFAnchorManager,
        cv::SubsampledControlVariate, state)
    bank = something(manager.oracle.node_local_hcv)
    x = state.ξ.x
    velocity = state.ξ.θ
    n_reference = size(bank.selector_positions, 2)
    iszero(n_reference) && return false
    requested = bank.requested_sources
    context = bank.contexts[1]
    @inbounds for node in eachindex(requested)
        nearest = 1
        nearest_distance = Inf
        start = context.threshold_starts[node]
        q = context.seen[node] - 1
        for reference in 1:n_reference
            distance = 0.0
            for offset in 0:(q - 1)
                j = context.threshold_indices[start + offset]
                delta = (x[j] - bank.selector_positions[j, reference]) /
                    manager.distance_scales[j]
                velocity_delta = (velocity[j] -
                    bank.selector_velocities[j, reference]) /
                    manager.distance_scales[j]
                distance += delta * delta + velocity_delta * velocity_delta
            end
            for edge in context.incident_edges[node]
                j = context.interaction_indices[edge]
                delta = (x[j] - bank.selector_positions[j, reference]) /
                    manager.distance_scales[j]
                velocity_delta = (velocity[j] -
                    bank.selector_velocities[j, reference]) /
                    manager.distance_scales[j]
                distance += delta * delta + velocity_delta * velocity_delta
            end
            if distance < nearest_distance
                nearest = reference
                nearest_distance = distance
            end
        end
        requested[node] = bank.selector_sources[node, nearest]
    end
    changed = _install_omrf_node_local_sources!(manager, cv, requested)
    manager.node_source_checks += 1
    if changed
        manager.node_source_changes += 1
        manager.main_activations += 1
    end
    manager.main_selections += 1
    return changed
end

"""
Choose the exact aggregate observed-pattern certificate minimizer at one
retained event boundary using only warmup-compiled pattern representatives.

Pattern counts are sufficient for this certificate because a node factor
depends on a person only through its incident-neighbour covariates. No person
row, gradient, probability cache, or full-data object is touched here.
"""
function _select_omrf_pattern_count_sources!(manager::OMRFAnchorManager,
        cv::SubsampledControlVariate, state, flow)
    started = time_ns()
    bank = something(manager.oracle.node_local_hcv)
    P = length(bank.sources)
    K = size(bank.candidates, 2)
    requested = bank.requested_sources
    envelope_scales = cv.envelope.component_scales!
    node_scales = envelope_scales isa
        OMRFCellCachedNodeGroupedPersonHCVBoxScales ?
        envelope_scales.base :
        (envelope_scales isa OMRFNodeGroupedPersonHCVBoxScales ?
            envelope_scales : nothing)
    @inbounds for node in 1:P
        best_source = 1
        best_score = Inf
        persons = bank.node_pattern_persons[node]
        counts = bank.node_pattern_counts[node]
        for source in 1:K
            context = bank.contexts[source]
            anchor = view(bank.candidates, :, source)
            score = 0.0
            node_scales === nothing ||
                _cache_omrf_node_grouped_point_node!(
                    node_scales, state, flow, node, anchor)
            for pattern in eachindex(persons)
                bound = node_scales === nothing ? _omrf_node_hcv_bound_at(
                    context, anchor, state, flow, persons[pattern], node) :
                    _omrf_node_grouped_person_point_bound(
                        node_scales, persons[pattern], node)
                score += counts[pattern] * bound
            end
            if score < best_score
                best_source = source
                best_score = score
            end
        end
        requested[node] = best_source
    end
    changed = _install_omrf_node_local_sources!(manager, cv, requested)
    manager.node_source_checks += 1
    if changed
        manager.node_source_changes += 1
        manager.main_activations += 1
    end
    manager.main_selections += 1
    manager.selection_component_evaluations += K * sum(
        length, bank.node_pattern_persons)
    manager.selection_seconds += (time_ns() - started) * 1e-9
    return changed
end

function _finish_omrf_warmup!(manager::OMRFAnchorManager,
        cv::SubsampledControlVariate, position, trace_manager=nothing;
        state=nothing, flow=nothing)
    previous_anchor = copy(manager.oracle.anchor)
    frozen_bank = isinf(manager.main_refresh_distance) && manager.capacity > 1
    trace = trace_manager === nothing ? nothing :
        PDMPSamplers.get_warmup_trace(trace_manager)
    if frozen_bank && trace !== nothing &&
            PDMPSamplers._has_integrable_segment(trace)
        representatives, scales = _omrf_representative_anchors(
            trace, position, manager.capacity)
        prepared = [_prepare_omrf_anchor(manager, view(representatives, :, idx))
            for idx in axes(representatives, 2)]
        empty!(manager.entries)
        manager.active_idx = 0
        manager.recycled_state = nothing
        copyto!(manager.distance_scales, scales)
        for state in prepared
            _insert_omrf_anchor!(manager, state)
        end
        _select_omrf_anchor!(manager, cv, position; phase=:warmup)
    elseif iszero(manager.preparations)
        prepared = _prepare_omrf_anchor(manager, position)
        manager.active_idx = _insert_omrf_anchor!(manager, prepared)
        PDMPSamplers.refresh_anchor!(cv, prepared.anchor)
    else
        # A scheduled warmup anchor is already an exact control variate. Reuse
        # the closest prepared entry instead of paying for a second anchor that
        # immediately supersedes it at the phase boundary.
        _select_omrf_anchor!(manager, cv, position; phase=:warmup)
    end
    anchor_changed = manager.oracle.anchor != previous_anchor
    hcv_changed = manager.oracle.residual_context.use_hcv &&
        !manager.hcv_active
    if hcv_changed
        manager.hcv_active = true
        for entry in manager.entries
            entry.state.residual_context.hcv_active = true
        end
        active = manager.entries[manager.active_idx].state
        PDMPSamplers.refresh_anchor!(cv, active.anchor)
    end
    active_scales = cv.envelope.component_scales!
    node_local_hcv = manager.oracle.residual_context.factorization in
            (:person, :person_node) && manager.hcv_active &&
        manager.capacity > 1 && active_scales isa Union{
            OMRFBatchedFactorHCVScales, OMRFGroupedPersonHCVBoxScales,
            OMRFNodeGroupedPersonHCVBoxScales,
            OMRFCellCachedNodeGroupedPersonHCVBoxScales,
            OMRFLazyNodeGroupedPersonHCVBoxScales,
            OMRFPatternLocalHCVScales}
    if node_local_hcv
        state === nothing && throw(ArgumentError(
            "node-local HCV warmup completion requires the PDMP state"))
        flow === nothing && throw(ArgumentError(
            "node-local HCV warmup completion requires the PDMP flow"))
        manager.oracle.node_local_hcv = _prepare_omrf_node_local_hcv_bank(
            manager.entries, manager.oracle.ctx.counts)
        if trace !== nothing && PDMPSamplers._has_integrable_segment(trace)
            _select_omrf_warmup_pattern_sources!(
                manager, cv, trace, state, flow)
        else
            _select_omrf_node_local_sources!(manager, cv, state, flow)
        end
        copyto!(manager.node_source_last_position, state.ξ.x)
    end
    # This phase-boundary snapshot makes the no-full-data-main invariant
    # directly testable.  A preparation includes the exact complete anchor
    # sum, whether obtained through BridgeStan or the native OMRF algebra.
    manager.warmup_finish_preparations = manager.preparations
    manager.warmup_finish_native_full_data_scans =
        manager.oracle.ctx.counts.native_full_data_scans
    manager.warmup_finish_event_bound_observation_visits =
        manager.oracle.ctx.counts.event_bound_observation_visits
    return anchor_changed || hcv_changed
end

@inline function _omrf_anchor_distance2(manager::OMRFAnchorManager, x,
        state::OMRFPreparedAnchorState)
    distance = 0.0
    @inbounds for j in eachindex(x, state.anchor)
        delta = (x[j] - state.anchor[j]) / manager.distance_scales[j]
        distance += delta * delta
    end
    return distance
end

function _activate_omrf_anchor!(manager::OMRFAnchorManager, requested)
    entry = manager.entries[manager.active_idx]
    entry.state.anchor == requested || throw(ArgumentError(
        "active OMRF anchor entry does not match the requested anchor"))
    manager.activations += 1
    return _install_omrf_anchor!(manager, entry.state)
end

function _select_omrf_anchor!(manager::OMRFAnchorManager,
        cv::SubsampledControlVariate, x; phase::Symbol=:unknown,
        state=nothing, flow=nothing)
    if phase === :main && manager.oracle.node_local_hcv !== nothing
        state === nothing && throw(ArgumentError(
            "node-local HCV selection requires the PDMP state"))
        flow === nothing && throw(ArgumentError(
            "node-local HCV selection requires the PDMP flow"))
        if !_omrf_node_source_moved_enough(manager, state.ξ.x)
            manager.node_source_skips += 1
            return false
        end
        copyto!(manager.node_source_last_position, state.ξ.x)
        manager.node_source_checks += 1
        return _select_omrf_node_local_sources!(manager, cv, state, flow)
    end
    previous = manager.active_idx
    best_idx = firstindex(manager.entries)
    best_distance = Inf
    best_score = Inf
    envelope_selection = phase === :main &&
        isinf(manager.main_refresh_distance) && manager.capacity > 1 &&
        state !== nothing && !(manager.entries[1].state.envelope.
            component_scales! isa OMRFPatternLocalHCVScales)
    selection_start = time_ns()
    @inbounds for idx in eachindex(manager.entries)
        entry = manager.entries[idx]
        entry.age += 1
        distance = _omrf_anchor_distance2(manager, x, entry.state)
        score = envelope_selection ? PDMPSamplers.total_residual_bound(
            entry.state.envelope, state, flow, 0.0) : distance
        if score < best_score
            best_idx = idx
            best_distance = distance
            best_score = score
        end
    end
    manager.selection_seconds += (time_ns() - selection_start) * 1e-9
    if envelope_selection
        manager.selection_component_evaluations += sum(
            length(entry.state.envelope.totals) for entry in manager.entries)
    end
    manager.entries[best_idx].age = 0
    manager.active_idx = best_idx
    if phase === :main
        distance = sqrt(best_distance)
        manager.main_selections += 1
        manager.main_distance_sum += distance
        manager.max_main_distance = max(manager.max_main_distance, distance)
        if distance > manager.main_refresh_distance
            prepared = _prepare_omrf_anchor(manager, x)
            best_idx = _insert_omrf_anchor!(manager, prepared)
            manager.active_idx = best_idx
            PDMPSamplers.refresh_anchor!(cv, prepared.anchor)
            manager.main_activations += 1
            manager.main_refreshes += 1
            return true
        end
    end
    if best_idx != previous
        PDMPSamplers.refresh_anchor!(cv,
            manager.entries[best_idx].state.anchor)
        phase === :main && (manager.main_activations += 1)
        return true
    end
    return false
end

function _add_omrf_anchor!(manager::OMRFAnchorManager,
        cv::SubsampledControlVariate, trace)
    # A frozen multi-anchor bank is constructed once from representative
    # warmup states at the phase boundary. Successive cumulative means cluster
    # near one point and waste the fixed full-data preparation budget.
    if isinf(manager.main_refresh_distance) && manager.capacity > 1
        return manager.active_idx
    end
    Statistics.mean!(manager.proposal, trace)
    prepared = _prepare_omrf_anchor(manager, manager.proposal)
    previous = manager.active_idx
    idx = _insert_omrf_anchor!(manager, prepared)
    if idx == previous
        manager.active_idx = idx
        PDMPSamplers.refresh_anchor!(cv, prepared.anchor)
    end
    return idx
end

mutable struct StanSubsamplingSlabDeterministic{O,S}
    oracle::O
    slab::S
    active::BitVector
    all_active::BitVector
    all_buf::Vector{Float64}
    active_buf::Vector{Float64}
end

function StanSubsamplingSlabDeterministic(oracle, slab, d)
    p = length(PDMPSamplers.beta_indices(slab))
    return StanSubsamplingSlabDeterministic(oracle, slab, trues(p), trues(p),
        zeros(d), zeros(d))
end

function (target::StanSubsamplingSlabDeterministic)(out, x)
    target.oracle(out, x)
    PDMPSamplers.active_prior_grad!(target.slab, target.all_buf, x,
        target.all_active)
    PDMPSamplers.active_prior_grad!(target.slab, target.active_buf, x,
        target.active)
    @. out = out - target.all_buf + target.active_buf
    return out
end

function PDMPSamplers.set_active_set!(target::StanSubsamplingSlabDeterministic,
        free::BitVector)
    indices = PDMPSamplers.beta_indices(target.slab)
    @inbounds for j in eachindex(indices)
        target.active[j] = free[indices[j]]
    end
    return nothing
end

function _new_stan_subsampling_context(lib_path::String, full_data::String,
        prior_data::String, m::Int)
    counts = StanSubsamplingCallCounts()
    construct = function(data_path)
        model = lock(_BRIDGESTAN_CALL_LOCK) do
            BridgeStan.StanModel(lib_path, data_path; warn=false)
        end
        counts.model_constructions += 1
        counts.data_constructions += 1
        return model
    end
    sm_full = construct(full_data)
    sm_prior = construct(prior_data)
    sm_selected = construct(full_data)
    names, prior_names, selected_names = lock(_BRIDGESTAN_CALL_LOCK) do
        (BridgeStan.param_unc_names(sm_full),
         BridgeStan.param_unc_names(sm_prior),
         BridgeStan.param_unc_names(sm_selected))
    end
    prior_names == names || throw(ArgumentError(
        "full and prior-only data must produce identical unconstrained parameter names"))
    selected_names == names || throw(ArgumentError(
        "selected and full contexts must have identical unconstrained parameter names"))
    set_fn, clear_fn = _resolve_subset_hooks(sm_selected)
    ctx = StanSubsamplingContext(PDMPModel(sm_full), PDMPModel(sm_prior),
        PDMPModel(sm_selected), set_fn, clear_fn, Vector{Int32}(undef, m), m,
        counts)
    _clear_stan_subset!(ctx)
    return ctx, String.(names)
end

function _build_stan_subsampling_model(ctx::StanSubsamplingContext,
        unc_names, subsampling, anchor, slab_prior, can_stick;
        anchor_capacity::Integer=0)
    N = Int(_rget(subsampling, :n_observations))
    m = Int(_rget(subsampling, :subsample_size))
    envelope_spec = _rget(subsampling, :residual_envelope)
    weights = _as_float_matrix(_rget(envelope_spec, :weights))
    size(weights, 2) == N || throw(DimensionMismatch(
        "residual envelope must have N observation columns"))
    d = ctx.full.d
    length(anchor) == d || throw(DimensionMismatch("subsampling anchor has the wrong dimension"))
    is_omrf = _haskey(envelope_spec, :type) &&
        String(_rget(envelope_spec, :type)) == "omrf"
    backend = is_omrf && _haskey(envelope_spec, :backend) ?
        String(_rget(envelope_spec, :backend)) : "stan"
    backend in ("analytic", "stan") || throw(ArgumentError(
        "unknown OMRF residual backend: $backend"))
    oracle = is_omrf && backend == "analytic" ?
        OMRFAnalyticSubsamplingOracle(ctx, envelope_spec, unc_names, anchor) :
        StanSubsamplingOracle(ctx, anchor, N)
    if oracle isa OMRFAnalyticSubsamplingOracle
        residual_context = oracle.residual_context
        people, nodes = size(residual_context.X)
        expected = residual_context.factorization === :person ?
            people : people * nodes
        N == expected || throw(DimensionMismatch(
            "OMRF $(residual_context.factorization) factorization requires " *
            "$expected likelihood contributions; received $N"))
    end
    # Analytic OMRF oracles own their anchor storage. Point the initial
    # envelope at that owned vector as well, so later state recycling cannot
    # mutate a caller-owned initial-position vector.
    envelope_anchor = oracle isa OMRFAnalyticSubsamplingOracle ?
        oracle.anchor : anchor
    envelope = _build_stan_residual_envelope(
        envelope_spec, envelope_anchor, oracle)
    staged_hcv = oracle isa OMRFAnalyticSubsamplingOracle &&
        oracle.residual_context.use_hcv && !oracle.residual_context.hcv_active
    manager = (anchor_capacity > 0 || staged_hcv) &&
        oracle isa OMRFAnalyticSubsamplingOracle ?
        _new_omrf_anchor_manager(
            oracle, envelope_spec, unc_names, envelope, max(anchor_capacity, 1)) : nothing
    deterministic = if isnothing(slab_prior)
        oracle
    else
        provider = build_slab_provider(slab_prior, unc_names, can_stick, d)
        StanSubsamplingSlabDeterministic(oracle, provider, d)
    end
    refresh_callback = manager === nothing ? nothing :
        (requested -> _activate_omrf_anchor!(manager, requested))
    factor_sampling = is_omrf && _haskey(envelope_spec, :factor_sampling) ?
        Symbol(String(_rget(envelope_spec, :factor_sampling))) : :size_biased
    subset_design = if factor_sampling === :node_stratified
        oracle isa OMRFAnalyticSubsamplingOracle || throw(ArgumentError(
            "node-stratified subsets require the analytic OMRF provider"))
        people, nodes = size(oracle.residual_context.X)
        oracle.residual_context.factorization === :person_node || throw(ArgumentError(
            "node-stratified subsets require person-node factorization"))
        PDMPSamplers.balanced_stratified_subsampling_design(
            people * nodes, nodes, m)
    elseif factor_sampling === :size_biased
        PDMPSamplers.UniformSubsamplingDesign()
    else
        throw(ArgumentError("unknown OMRF factor sampling design: $factor_sampling"))
    end
    cv = SubsampledControlVariate(deterministic, oracle, envelope, anchor, m;
        refresh_anchor! = refresh_callback, subset_design)
    model = PDMPModel(d, cv)
    if manager === nothing
        return model, ctx, unc_names
    end
    adapter = PDMPSamplers.SubsamplingAnchorBankAdapter(
        (grad, state, flow, phase) -> _select_omrf_anchor!(manager, grad,
            state.ξ.x; phase, state, flow),
        (grad, trace) -> _add_omrf_anchor!(manager, grad, trace),
        (grad, state, flow, trace_manager, stats) ->
            _finish_omrf_warmup!(manager, grad, state.ξ.x, trace_manager;
                state, flow),
        Inf, 0.0)
    return model, ctx, unc_names, adapter, manager
end

function prepare_stan_subsampling(lib_path::String, full_data::String,
        prior_data::String, subsampling, n_chains::Integer)
    n_chains >= 1 || throw(ArgumentError("n_chains must be positive"))
    m = Int(_rget(subsampling, :subsample_size))
    contexts = StanSubsamplingContext[]
    unc_names = String[]
    for _ in 1:n_chains
        ctx, names = _new_stan_subsampling_context(
            lib_path, full_data, prior_data, m)
        isempty(unc_names) || names == unc_names || throw(ArgumentError(
            "subsampling contexts have inconsistent unconstrained parameter names"))
        push!(contexts, ctx)
        unc_names = names
    end
    d = first(contexts).full.d
    return (; contexts, unc_names, d)
end

function r_stan_subsampling_diagnostics(lib_path::String, full_data::String,
        prior_data::String, subsampling, position, subset, velocity,
        flow_type::String, flow_mean, flow_cov, diagnostic_phase::String)
    N = Int(_rget(subsampling, :n_observations))
    m = Int(_rget(subsampling, :subsample_size))
    subset_vec = _as_int_vector(subset)
    length(subset_vec) == m || throw(DimensionMismatch(
        "diagnostic subset must contain exactly the configured subsample size"))
    all(i -> 1 <= i <= N, subset_vec) || throw(ArgumentError(
        "diagnostic subset indices must lie in 1:N"))
    allunique(subset_vec) || throw(ArgumentError(
        "diagnostic subset indices must be unique"))

    x = _as_float_vector(position)
    subsampling_anchor = _rget(subsampling, :anchor)
    anchor = isnothing(subsampling_anchor) ? copy(x) : _as_float_vector(subsampling_anchor)
    ctx, unc_names = _new_stan_subsampling_context(lib_path, full_data, prior_data, m)
    d = ctx.full.d
    length(x) == d || throw(DimensionMismatch("diagnostic position has the wrong dimension"))
    length(anchor) == d || throw(DimensionMismatch("diagnostic anchor has the wrong dimension"))
    envelope_spec = _rget(subsampling, :residual_envelope)
    is_omrf = _haskey(envelope_spec, :type) &&
        String(_rget(envelope_spec, :type)) == "omrf"
    backend = is_omrf && _haskey(envelope_spec, :backend) ?
        String(_rget(envelope_spec, :backend)) : "stan"
    factorization = is_omrf && _haskey(envelope_spec, :factorization) ?
        String(_rget(envelope_spec, :factorization)) : "person"
    factorization == "person_node" && throw(ArgumentError(
        "stan_subsampling_diagnostics currently supports person-level OMRF " *
        "contributions; person-node sampling is validated by the factor-local contract tests"))
    diagnostic_phase in ("warmup", "main") || throw(ArgumentError(
        "diagnostic phase must be warmup or main"))
    oracle = is_omrf && backend == "analytic" ?
        OMRFAnalyticSubsamplingOracle(ctx, envelope_spec, unc_names, anchor) :
        StanSubsamplingOracle(ctx, anchor, N)
    if oracle isa OMRFAnalyticSubsamplingOracle &&
            oracle.residual_context.use_hcv &&
            _haskey(envelope_spec, :hcv_after_warmup) &&
            Bool(_rget(envelope_spec, :hcv_after_warmup))
        oracle.residual_context.hcv_active = diagnostic_phase == "main"
    end
    deterministic = zeros(d)
    residual = zeros(d)
    full = zeros(d)
    oracle(deterministic, x)
    oracle(residual, x, subset_vec, anchor)
    prior_gradient = copy(_prior_gradient_at_x!(oracle, x))
    if oracle isa OMRFAnalyticSubsamplingOracle
        selected_gradient = zeros(d)
        selected_anchor_gradient = zeros(d)
        _selected_gradient!(ctx, selected_gradient, x, subset_vec)
        _selected_gradient!(ctx, selected_anchor_gradient, anchor, subset_vec)
        selected_likelihood = selected_gradient - prior_gradient
        selected_anchor_likelihood = selected_anchor_gradient - oracle.prior_anchor
    else
        selected_gradient = copy(oracle.selected_x)
        selected_likelihood = oracle.selected_x - oracle.prior_x
        selected_anchor_likelihood = oracle.selected_anchor - oracle.prior_anchor
    end
    _clear_gradient!(ctx, :full, full, x)
    subsampled_gradient = deterministic + (N / m) * residual
    anchor_deterministic = zeros(d)
    anchor_residual = zeros(d)
    oracle(anchor_deterministic, anchor)
    oracle(anchor_residual, anchor, subset_vec, anchor)
    anchor_subsampling = anchor_deterministic + (N / m) * anchor_residual
    anchor_closure_error = maximum(abs, anchor_subsampling - oracle.full_anchor)

    envelope = _build_stan_residual_envelope(envelope_spec, anchor, oracle)
    diagnostic_velocity = _as_float_vector(velocity)
    length(diagnostic_velocity) == d || throw(DimensionMismatch(
        "diagnostic velocity has the wrong dimension"))
    flow = build_flow(flow_type,
        _to_precision(_as_flow_cov(flow_cov, d), d),
        _as_flow_mean(flow_mean, d))
    state = PDMPState(0.0, SkeletonPoint(copy(x), diagnostic_velocity))
    PDMPSamplers.initialize_flow_state!(state, flow)
    scales = PDMPSamplers.component_scales!(
        envelope.scales, envelope, state, flow, 0.0)
    envelope_rate = (N / m) * sum(i ->
        PDMPSamplers.observation_residual_bound(envelope, i), subset_vec)
    scaled_residual = (N / m) * residual
    residual_rate = PDMPSamplers.λ(state, scaled_residual, flow) +
        PDMPSamplers.λ(state, -scaled_residual, flow)

    rate_decomposition = nothing
    if m == 1 && oracle isa OMRFAnalyticSubsamplingOracle
        per_observation_residual_rate = zeros(N)
        per_observation_envelope = zeros(N)
        per_observation_first_order_envelope = zeros(N)
        per_observation_remainder_envelope = zeros(N)
        per_observation_stochastic_rate = zeros(N)
        diagnostic_residual = zeros(d)
        diagnostic_gradient = zeros(d)
        deterministic_rate = PDMPSamplers.λ(state, deterministic, flow)
        for observation in 1:N
            oracle(diagnostic_residual, x, (observation,), anchor)
            per_observation_residual_rate[observation] =
                PDMPSamplers.λ(state, diagnostic_residual, flow) +
                PDMPSamplers.λ(state, -diagnostic_residual, flow)
            per_observation_envelope[observation] =
                PDMPSamplers.observation_residual_bound(
                    envelope, observation)
            if envelope.component_scales! isa OMRFDampedHCVNodeLocalScales
                P = length(oracle.residual_context.seen)
                per_observation_first_order_envelope[observation] = dot(
                    view(envelope.weights, 1:P, observation),
                    view(scales, 1:P))
                per_observation_remainder_envelope[observation] = dot(
                    view(envelope.weights, (P + 1):(2P), observation),
                    view(scales, (P + 1):(2P)))
            end
            @. diagnostic_gradient = deterministic + N * diagnostic_residual
            per_observation_stochastic_rate[observation] =
                PDMPSamplers.λ(state, diagnostic_gradient, flow)
        end
        exact_residual_total = sum(per_observation_residual_rate)
        configured_residual_total = sum(per_observation_envelope)
        first_order_envelope_total =
            sum(per_observation_first_order_envelope)
        remainder_envelope_total =
            sum(per_observation_remainder_envelope)
        expected_stochastic_rate = sum(per_observation_stochastic_rate) / N
        exact_proposal_rate = deterministic_rate + exact_residual_total
        configured_proposal_rate = deterministic_rate + configured_residual_total
        full_rate = PDMPSamplers.λ(state, full, flow)
        rate_decomposition = Dict{String,Any}(
            "deterministic_rate" => deterministic_rate,
            "full_rate" => full_rate,
            "expected_stochastic_rate" => expected_stochastic_rate,
            "exact_residual_total" => exact_residual_total,
            "configured_residual_total" => configured_residual_total,
            "first_order_envelope_total" => first_order_envelope_total,
            "remainder_envelope_total" => remainder_envelope_total,
            "exact_proposal_rate" => exact_proposal_rate,
            "configured_proposal_rate" => configured_proposal_rate,
            "envelope_inflation" => configured_residual_total /
                max(exact_residual_total, eps(Float64)),
            "intrinsic_exact_acceptance" => expected_stochastic_rate /
                max(exact_proposal_rate, eps(Float64)),
            "configured_acceptance" => expected_stochastic_rate /
                max(configured_proposal_rate, eps(Float64)),
            "stochastic_event_inflation" => expected_stochastic_rate /
                max(full_rate, eps(Float64)),
            "per_observation_residual_rate" => per_observation_residual_rate,
            "per_observation_envelope" => per_observation_envelope,
            "per_observation_first_order_envelope" =>
                per_observation_first_order_envelope,
            "per_observation_remainder_envelope" =>
                per_observation_remainder_envelope,
            "per_observation_stochastic_rate" => per_observation_stochastic_rate,
        )
    end

    return Dict{String,Any}(
        "parameter_names" => unc_names,
        "full_gradient" => full,
        "prior_gradient" => prior_gradient,
        "selected_gradient" => selected_gradient,
        "selected_likelihood_gradient" => selected_likelihood,
        "selected_anchor_likelihood_gradient" => selected_anchor_likelihood,
        "deterministic_gradient" => deterministic,
        "residual_gradient" => residual,
        "subsampled_gradient" => subsampled_gradient,
        "anchor_full_gradient" => copy(oracle.full_anchor),
        "anchor_closure_error" => anchor_closure_error,
        "residual_rate" => residual_rate,
        "envelope_rate" => envelope_rate,
        "rate_decomposition" => rate_decomposition,
        "selected_scale" => N / m,
        "model_constructions" => ctx.counts.model_constructions,
        "data_constructions" => ctx.counts.data_constructions,
    )
end

function r_pdmp_stan_subsampling(prepared, subsampling, x0, flow_type::String,
        algorithm_type::String, flow_mean, flow_cov;
        theta0=nothing,
        c0::Float64=1e-2, grid_n::Int=30, grid_t_max::Float64=2.0,
        use_fd_hvp::Bool=true,
        curvature_backend::String="finite_difference",
        post_warmup_simplify::Bool=false,
        grid_bound::String="constant",
        grid_curvature_bound::Union{Nothing,Float64}=nothing,
        linear_area_threshold::Float64=0.95,
        linear_min_area_gain::Float64=0.0,
        lazy_low_tightness_threshold::Float64=0.1,
        lazy_max_low_tightness_rejections::Int=3,
        lazy_max_rejections::Int=0,
        t0::Float64=0.0, T::Float64=10000.0, t_warmup::Float64=0.0,
        warmup_adaptation_interval::Union{Nothing,Float64}=nothing,
        sticky::Bool=false, can_stick=nothing, model_prior=nothing,
        parameter_prior=nothing, slab_prior=nothing,
        show_progress::Bool=true, n_chains::Int=1, threaded::Bool=false,
        seed::Union{Integer,Nothing}=nothing,
        adaptive_scheme::String="diagonal",
        warmup_mode::String="subsampled",
        n_anchor_updates::Int=0,
        use_anchor_bank::Bool=false,
        anchor_bank_capacity::Int=8,
        main_anchor_refresh_distance::Float64=Inf,
        node_source_check_distance::Float64=0.0,
        require_no_main_full_data::Bool=false,
        support_boundary_mode::String="error",
        support_boundary_max_bisection_steps::Int=60,
        support_boundary_time_rtol::Float64=1e-8,
        support_boundary_time_atol::Float64=1e-10,
        support_boundary_clip_fraction::Float64=1 - 1e-10,
        support_boundary_max_refresh_attempts::Int=20,
        support_boundary_refresh_probe_time::Float64=1e-4,
        support_boundary_min_safe_time::Float64=1e-12)
    algorithm_type in ("GridThinningStrategy", "ThinningStrategy") ||
        throw(ArgumentError("subsampling custom-Stan sampling requires GridThinningStrategy or ThinningStrategy"))
    warmup_mode in ("subsampled", "full") || throw(ArgumentError(
        "warmup_mode must be \"subsampled\" or \"full\""))
    warmup_mode == "full" && t_warmup <= 0 && throw(ArgumentError(
        "full-gradient warmup requires positive t_warmup"))
    x0_vec = _as_float_vector(x0)
    subsampling_anchor = _rget(subsampling, :anchor)
    anchor = isnothing(subsampling_anchor) ? copy(x0_vec) :
        _as_float_vector(subsampling_anchor)
    can_stick_vec = _as_bool_vector(can_stick)
    models = PDMPModel[]
    anchor_adapters = Any[]
    anchor_managers = Union{Nothing,OMRFAnchorManager}[]
    n_anchor_updates >= 0 || throw(ArgumentError(
        "n_anchor_updates must be nonnegative"))
    anchor_bank_capacity >= 1 || throw(ArgumentError(
        "anchor_bank_capacity must be positive"))
    (isinf(main_anchor_refresh_distance) ||
     isfinite(main_anchor_refresh_distance) && main_anchor_refresh_distance > 0) ||
        throw(ArgumentError(
            "main_anchor_refresh_distance must be positive or Inf"))
    isfinite(node_source_check_distance) && node_source_check_distance >= 0 ||
        throw(ArgumentError(
            "node_source_check_distance must be finite and nonnegative"))
    anchor_capacity = warmup_mode == "full" ? max(1, use_anchor_bank ?
        anchor_bank_capacity : (n_anchor_updates > 0 ? 1 : 0)) :
        (use_anchor_bank ? anchor_bank_capacity : (n_anchor_updates > 0 ? 1 : 0))
    contexts = prepared.contexts
    length(contexts) == n_chains || throw(DimensionMismatch(
        "prepared subsampling context count must equal n_chains"))
    unc_names = prepared.unc_names
    for ctx in contexts
        built = _build_stan_subsampling_model(ctx, unc_names,
            subsampling, anchor, slab_prior, can_stick_vec; anchor_capacity)
        if length(built) == 5
            model, _, _, adapter, manager = built
            adapter.update_dt = n_anchor_updates > 0 ?
                t_warmup / n_anchor_updates : Inf
            adapter.last_update = t0
            manager.main_refresh_distance = main_anchor_refresh_distance
            manager.node_source_check_distance = node_source_check_distance
            push!(models, model)
            push!(anchor_adapters, adapter)
            push!(anchor_managers, manager)
        else
            model, _, _ = built
            push!(models, model)
            push!(anchor_managers, nothing)
        end
    end
    construction_snapshots = [(
        model=ctx.counts.model_constructions,
        data=ctx.counts.data_constructions,
        full_gradient=ctx.counts.full_gradient,
        native_full_data_scans=ctx.counts.native_full_data_scans,
        event_bound_observation_visits=
            ctx.counts.event_bound_observation_visits,
    ) for ctx in contexts]
    d = first(models).d
    flow = build_flow(flow_type, _to_precision(_as_flow_cov(flow_cov, d), d),
        _as_flow_mean(flow_mean, d); adaptive_scheme)
    alg0 = build_algorithm(algorithm_type; c0, d, grid_n, grid_t_max,
        use_fd_hvp, curvature_backend,
        post_warmup_simplify, grid_bound, grid_curvature_bound,
        linear_area_threshold, linear_min_area_gain,
        lazy_low_tightness_threshold, lazy_max_low_tightness_rejections,
        lazy_max_rejections)
    alg = isnothing(slab_prior) ?
        wrap_sticky(alg0, sticky, model_prior,
            isnothing(parameter_prior) ? nothing : _as_float_vector(parameter_prior),
            can_stick_vec) :
        wrap_dependent_sticky(alg0, sticky, model_prior, slab_prior,
            can_stick_vec, flow, unc_names)
    warmup_models = warmup_mode == "full" ?
        [build_omrf_full_model(_rget(subsampling, :residual_envelope),
            unc_names, d) for _ in 1:n_chains] : nothing
    warmup_alg0 = warmup_mode == "full" ? build_algorithm(algorithm_type;
        c0, d, grid_n, grid_t_max, use_fd_hvp=false,
        curvature_backend="exact", post_warmup_simplify,
        grid_bound, grid_curvature_bound, linear_area_threshold,
        linear_min_area_gain, lazy_low_tightness_threshold,
        lazy_max_low_tightness_rejections, lazy_max_rejections) : nothing
    warmup_alg = if warmup_mode == "full"
        isnothing(slab_prior) ?
            wrap_sticky(warmup_alg0, sticky, model_prior,
                isnothing(parameter_prior) ? nothing :
                    _as_float_vector(parameter_prior), can_stick_vec) :
            wrap_dependent_sticky(warmup_alg0, sticky, model_prior,
                slab_prior, can_stick_vec, flow, unc_names)
    else
        nothing
    end
    sbopts = SupportBoundaryOptions(;
        detect_boundaries = support_boundary_mode != "error",
        mode = Symbol(support_boundary_mode),
        max_bisection_steps = support_boundary_max_bisection_steps,
        time_rtol = support_boundary_time_rtol,
        time_atol = support_boundary_time_atol,
        clip_fraction = support_boundary_clip_fraction,
        max_refresh_attempts = support_boundary_max_refresh_attempts,
        refresh_probe_time = support_boundary_refresh_probe_time,
        min_safe_time = support_boundary_min_safe_time)
    initial = isnothing(theta0) ? x0_vec :
        SkeletonPoint(x0_vec, _as_float_vector(theta0))
    chains = if all(isnothing, anchor_managers)
        adapter = PDMPSamplers.default_warmup_adapter(
            flow, first(models).grad, t_warmup, t0;
            warmup_adaptation_interval,
            can_stick=sticky ? can_stick_vec : nothing)
        pdmp_sample(initial, flow, models, alg, t0, T, t_warmup;
            progress=show_progress, threaded, seed,
            adapter,
            warmup_models, warmup_algorithm=warmup_alg,
            warmup_stop=build_warmup_stop(t_warmup),
            support_boundary_options=sbopts,
            statistic_counter=PDMPSamplers.DevelStatisticCounter)
    else
        run_chain = function(i)
            seed_i = isnothing(seed) ? nothing : seed + i - 1
            flow_i = deepcopy(flow)
            dynamics_adapter = PDMPSamplers.default_warmup_adapter(
                flow_i, models[i].grad, t_warmup, t0;
                warmup_adaptation_interval,
                can_stick=sticky ? can_stick_vec : nothing)
            adapter_i = dynamics_adapter isa PDMPSamplers.NoAdaptation ?
                anchor_adapters[i] : PDMPSamplers.SequenceAdapter(
                    (dynamics_adapter, anchor_adapters[i]))
            chain = pdmp_sample(initial, flow_i, [models[i]], alg,
                t0, T, t_warmup; progress=show_progress && i == 1,
                adapter=adapter_i, seed=seed_i,
                warmup_models=isnothing(warmup_models) ? nothing :
                    [warmup_models[i]], warmup_algorithm=warmup_alg,
                warmup_stop=build_warmup_stop(t_warmup),
                support_boundary_options=sbopts,
                statistic_counter=PDMPSamplers.DevelStatisticCounter)
            return chain.traces[1], chain.stats[1]
        end
        chain_results = if threaded && n_chains > 1
            fetch.([Threads.@spawn run_chain(i) for i in 1:n_chains])
        else
            [run_chain(i) for i in 1:n_chains]
        end
        PDMPChains(first.(chain_results), last.(chain_results))
    end
    result = _pack_result(chains)
    subsampling_context_counters = [Dict(
        "full_gradient_calls" => ctx.counts.full_gradient,
        "prior_gradient_calls" => ctx.counts.prior_gradient,
        "selected_gradient_calls" => ctx.counts.selected_gradient,
        "persons_evaluated" => ctx.counts.persons_evaluated +
            ctx.counts.analytic_persons,
        "selected_persons_evaluated" => ctx.counts.persons_evaluated,
        "analytic_residual_calls" => ctx.counts.analytic_residual,
        "analytic_persons_evaluated" => ctx.counts.analytic_persons,
        "analytic_factors_evaluated" => ctx.counts.analytic_factors,
        "node_conditionals_evaluated" => ctx.counts.analytic_node_conditionals,
        "native_full_gradient_seconds" =>
            ctx.counts.native_full_gradient_seconds,
        "analytic_residual_seconds" => ctx.counts.analytic_residual_seconds,
        "envelope_point_calls" => ctx.counts.envelope_point_calls,
        "envelope_point_seconds" => ctx.counts.envelope_point_seconds,
        "envelope_cell_calls" => ctx.counts.envelope_cell_calls,
        "envelope_cell_seconds" => ctx.counts.envelope_cell_seconds,
        "proposal_deterministic_bound_sum" =>
            ctx.counts.proposal_deterministic_bound_sum,
        "proposal_aggregate_residual_bound_sum" =>
            ctx.counts.proposal_aggregate_residual_bound_sum,
        "proposal_subset_bound_sum" => ctx.counts.proposal_subset_bound_sum,
        "proposal_deterministic_actual_sum" =>
            ctx.counts.proposal_deterministic_actual_sum,
        "proposal_exact_residual_rate_sum" =>
            ctx.counts.proposal_exact_residual_rate_sum,
        "proposal_actual_rate_sum" => ctx.counts.proposal_actual_rate_sum,
        "proposal_count" => ctx.counts.proposal_count,
        "proposal_mean_acceptance_probability" =>
            ctx.counts.proposal_acceptance_probability_sum /
            max(ctx.counts.proposal_count, 1),
        "proposal_mean_certificate_ratio" =>
            ctx.counts.proposal_certificate_ratio_sum /
            max(ctx.counts.proposal_count, 1),
        "proposal_mean_intrinsic_ratio" =>
            ctx.counts.proposal_intrinsic_ratio_sum /
            max(ctx.counts.proposal_count, 1),
        "lazy_branch_evaluations" => ctx.counts.lazy_branch_evaluations,
        "lazy_branch_rejections" => ctx.counts.lazy_branch_rejections,
        "lazy_leaf_evaluations" => ctx.counts.lazy_leaf_evaluations,
        "lazy_sampling_seconds" => ctx.counts.lazy_sampling_seconds,
        "cell_cached_point_screen_seconds" =>
            ctx.counts.cell_cached_point_screen_seconds,
        "proposal_envelope_acceptance_ratio" =>
            ctx.counts.proposal_actual_rate_sum /
            max(ctx.counts.proposal_subset_bound_sum, eps(Float64)),
        "proposal_certificate_efficiency" =>
            (ctx.counts.proposal_deterministic_actual_sum +
             ctx.counts.proposal_exact_residual_rate_sum) /
            max(ctx.counts.proposal_subset_bound_sum, eps(Float64)),
        "proposal_intrinsic_acceptance_ratio" =>
            ctx.counts.proposal_actual_rate_sum /
            max(ctx.counts.proposal_deterministic_actual_sum +
                ctx.counts.proposal_exact_residual_rate_sum, eps(Float64)),
        "accepted_proposal_bound_sum" => ctx.counts.accepted_proposal_bound_sum,
        "accepted_actual_rate_sum" => ctx.counts.accepted_actual_rate_sum,
        "omrf_factorization" => model.grad.residual_oracle isa
            OMRFAnalyticSubsamplingOracle ? String(
                model.grad.residual_oracle.residual_context.factorization) : "person",
        "anchor_cache_gradient_calls" => ctx.counts.anchor_cache_gradient,
        "anchor_cache_persons_evaluated" => ctx.counts.anchor_cache_persons,
        "anchor_likelihood_cache_bytes" => _anchor_cache_bytes(
            model.grad.residual_oracle),
        "model_constructions" => ctx.counts.model_constructions,
        "data_constructions" => ctx.counts.data_constructions,
        "native_full_data_scans" => ctx.counts.native_full_data_scans,
        "event_bound_observation_visits" =>
            ctx.counts.event_bound_observation_visits,
        "initialization_model_constructions" => snapshot.model,
        "initialization_data_constructions" => snapshot.data,
        "initialization_full_gradient_calls" => snapshot.full_gradient,
        "sampling_model_constructions" => ctx.counts.model_constructions - snapshot.model,
        "sampling_data_constructions" => ctx.counts.data_constructions - snapshot.data,
        "sampling_full_gradient_calls" => ctx.counts.full_gradient - snapshot.full_gradient,
        "sampling_native_full_data_scans" =>
            ctx.counts.native_full_data_scans - snapshot.native_full_data_scans,
        "anchor_preparations" => manager === nothing ? 0 : manager.preparations,
        "anchor_warmup_finish_preparations" => manager === nothing ? 0 :
            manager.warmup_finish_preparations,
        "anchor_main_preparations" => manager === nothing ? 0 :
            manager.preparations - manager.warmup_finish_preparations,
        "main_native_full_data_scans" => manager === nothing ?
            ctx.counts.native_full_data_scans - snapshot.native_full_data_scans :
            ctx.counts.native_full_data_scans -
                manager.warmup_finish_native_full_data_scans,
        "main_event_bound_observation_visits" => manager === nothing ?
            ctx.counts.event_bound_observation_visits -
                snapshot.event_bound_observation_visits :
            ctx.counts.event_bound_observation_visits -
                manager.warmup_finish_event_bound_observation_visits,
        "anchor_recycled_preparations" => manager === nothing ? 0 :
            manager.recycled_preparations,
        "anchor_preparation_seconds" => manager === nothing ? 0.0 : manager.preparation_seconds,
        "anchor_activations" => manager === nothing ? 0 : manager.activations,
        "anchor_main_activations" => manager === nothing ? 0 : manager.main_activations,
        "anchor_main_refreshes" => manager === nothing ? 0 : manager.main_refreshes,
        "anchor_main_selections" => manager === nothing ? 0 : manager.main_selections,
        "anchor_selection_seconds" => manager === nothing ? 0.0 :
            manager.selection_seconds,
        "anchor_selection_component_evaluations" => manager === nothing ? 0 :
            manager.selection_component_evaluations,
        "node_source_check_distance" => manager === nothing ? 0.0 :
            manager.node_source_check_distance,
        "node_source_checks" => manager === nothing ? 0 :
            manager.node_source_checks,
        "node_source_skips" => manager === nothing ? 0 :
            manager.node_source_skips,
        "node_source_changes" => manager === nothing ? 0 :
            manager.node_source_changes,
        "anchor_mean_main_distance" => manager === nothing ||
            iszero(manager.main_selections) ? 0.0 :
            manager.main_distance_sum / manager.main_selections,
        "anchor_max_main_distance" => manager === nothing ? 0.0 :
            manager.max_main_distance,
        "anchor_bank_entries" => manager === nothing ? 0 : length(manager.entries),
        "active_anchor" => manager === nothing ? Float64[] : copy(
            manager.entries[manager.active_idx].state.anchor),
        "anchor_positions" => manager === nothing ? zeros(0, 0) :
            reduce(hcat, (entry.state.anchor for entry in manager.entries)),
        "anchor_bank_bytes" => manager === nothing ? 0 :
            Base.summarysize((manager.entries, manager.recycled_state)),
        "analytic_hcv" => model.grad.residual_oracle isa OMRFAnalyticSubsamplingOracle &&
            model.grad.residual_oracle.residual_context.use_hcv,
        "analytic_prior" => model.grad.residual_oracle isa
            OMRFAnalyticSubsamplingOracle &&
            model.grad.residual_oracle.analytic_prior !== nothing,
        "hcv_damping" => model.grad.residual_oracle isa OMRFAnalyticSubsamplingOracle &&
            model.grad.residual_oracle.residual_context.use_hcv ?
            model.grad.residual_oracle.residual_context.hcv_damping : 0.0,
    ) for (model, ctx, snapshot, manager) in
        zip(models, contexts, construction_snapshots, anchor_managers)]
    result["subsampling_context_counters"] = subsampling_context_counters
    require_no_main_full_data && _enforce_no_main_full_data!(
        result["stats"]["main_full_gradient_calls"],
        subsampling_context_counters)
    result["subsampling"] = true
    return result
end

function _enforce_no_main_full_data!(main_full_gradient_calls,
        subsampling_context_counters)
    length(main_full_gradient_calls) == length(subsampling_context_counters) ||
        throw(DimensionMismatch(
            "main full-gradient and subsampling counter lengths differ"))
    violations = String[]
    for (chain, counters) in enumerate(subsampling_context_counters)
        full_gradients = main_full_gradient_calls[chain]
        anchor_preparations = counters["anchor_main_preparations"]
        native_scans = counters["main_native_full_data_scans"]
        event_bound_visits = get(
            counters, "main_event_bound_observation_visits", 0)
        if !iszero(full_gradients) || !iszero(anchor_preparations) ||
                !iszero(native_scans) || !iszero(event_bound_visits)
            push!(violations, "chain $chain: main_full_gradient_calls=" *
                "$full_gradients, anchor_main_preparations=" *
                "$anchor_preparations, main_native_full_data_scans=" *
                "$native_scans, main_event_bound_observation_visits=" *
                "$event_bound_visits")
        end
    end
    isempty(violations) || throw(ErrorException(
        "retained sampling performed forbidden full-data work (" *
        join(violations, "; ") * ")"))
    return nothing
end

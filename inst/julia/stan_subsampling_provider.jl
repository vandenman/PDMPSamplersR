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
end
StanSubsamplingCallCounts() = StanSubsamplingCallCounts(
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0, 0.0, 0.0, 0.0)

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
    seen::Vector{Int}
    threshold_indices::Vector{Int}
    interaction_indices::Vector{Int}
    threshold_starts::Vector{Int}
    incident_edges::Vector{Vector{Int}}
    incident_neighbours::Vector{Vector{Int}}
    anchor_probabilities::Matrix{Float64}
    anchor_expected_scores::Matrix{Float64}
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

function _omrf_analytic_prior(envelope_spec, context::OMRFResidualContext)
    _haskey(envelope_spec, :analytic_prior) || return nothing
    spec = _rget(envelope_spec, :analytic_prior)
    interaction_type = Symbol(String(_rget(spec, :type)))
    interaction_type in (:cauchy, :gaussian) || throw(ArgumentError(
        "unsupported analytic OMRF interaction prior: $interaction_type"))
    α = Float64(_rget(spec, :threshold_alpha))
    β = Float64(_rget(spec, :threshold_beta))
    scale = Float64(_rget(spec, :interaction_scale))
    all(isfinite, (α, β, scale)) && α > 0 && β > 0 && scale > 0 ||
        throw(ArgumentError("analytic OMRF prior parameters must be finite and positive"))
    return OMRFIndependentAnalyticPrior(
        context.threshold_indices, context.interaction_indices,
        α, β, interaction_type, scale)
end

@inline function _omrf_logistic(x::Float64)
    if x >= 0
        return inv(1 + exp(-x))
    end
    ex = exp(x)
    return ex / (1 + ex)
end

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
    prior::OMRFIndependentAnalyticPrior
end

function Base.copy(provider::OMRFFullGradientProvider)
    return deepcopy(provider)
end

function (provider::OMRFFullGradientProvider)(out, x)
    context = provider.context
    _omrf_prior_gradient!(out, provider.prior, x)
    N, P = size(context.X)
    @inbounds for n in 1:N, node in 1:P
        expected = _omrf_node_probabilities!(context, x, n, node)
        observed = context.X[n, node]
        q = context.seen[node] - 1
        start = context.threshold_starts[node]
        for u in 1:q
            out[context.threshold_indices[start + u - 1]] +=
                context.probability_buffer[u] - (observed == u)
        end
        for k in eachindex(context.incident_edges[node])
            edge = context.incident_edges[node][k]
            neighbour = context.incident_neighbours[node][k]
            out[context.interaction_indices[edge]] +=
                context.X[n, neighbour] * (expected - observed)
        end
    end
    return out
end

"""Refresh an OMRF anchor cache and its exact full gradient in one pass."""
function _cache_omrf_anchor_and_full_gradient!(context::OMRFResidualContext,
        full, prior::OMRFIndependentAnalyticPrior, anchor)
    _omrf_prior_gradient!(full, prior, anchor)
    N, P = size(context.X)
    @inbounds for n in 1:N, node in 1:P
        expected = _omrf_node_probabilities!(context, anchor, n, node)
        context.anchor_expected_scores[node, n] = expected
        observed = context.X[n, node]
        q = context.seen[node] - 1
        start = context.threshold_starts[node]
        for u in 1:q
            probability = context.probability_buffer[u]
            context.anchor_probabilities[start + u - 1, n] = probability
            full[context.threshold_indices[start + u - 1]] +=
                probability - (observed == u)
        end
        for k in eachindex(context.incident_edges[node])
            edge = context.incident_edges[node][k]
            neighbour = context.incident_neighbours[node][k]
            full[context.interaction_indices[edge]] +=
                context.X[n, neighbour] * (expected - observed)
        end
    end
    return full
end

"""Exact Hessian-vector product paired with `OMRFFullGradientProvider`."""
mutable struct OMRFFullHVPProvider <: Function
    context::OMRFResidualContext
    prior::OMRFIndependentAnalyticPrior
end

function Base.copy(provider::OMRFFullHVPProvider)
    return deepcopy(provider)
end

function (provider::OMRFFullHVPProvider)(out, x, direction)
    context = provider.context
    fill!(out, 0.0)
    prior = provider.prior
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
    return out
end

"""
    build_omrf_full_model(spec, unconstrained_names, d)

Build a non-factorized `PDMPModel` whose complete OMRF negative gradient and
Hessian-vector product are evaluated in native Julia.  The analytic prior and
the threshold/interaction blocks must cover all `d` unconstrained coordinates;
this intentionally rejects models with unimplemented nuisance parameters.
"""
function build_omrf_full_model(spec, unconstrained_names, d::Integer)
    d = Int(d)
    anchor = zeros(d)
    gradient_context = _omrf_residual_context(
        spec, unconstrained_names, anchor, StanSubsamplingCallCounts())
    covered = sort!(vcat(copy(gradient_context.threshold_indices),
        copy(gradient_context.interaction_indices)))
    covered == collect(1:d) || throw(ArgumentError(
        "native OMRF full gradients require threshold and interaction blocks " *
        "to cover every unconstrained coordinate; covered $covered of 1:$d"))
    gradient_prior = _omrf_analytic_prior(spec, gradient_context)
    gradient_prior === nothing && throw(ArgumentError(
        "native OMRF full gradients require a supported analytic independent prior"))
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
    isfinite(hcv_damping) && hcv_damping > 0 || throw(ArgumentError(
        "OMRF HCV damping must be finite and positive"))
    d = length(anchor)
    factorization = _haskey(envelope_spec, :factorization) ?
        Symbol(String(_rget(envelope_spec, :factorization))) : :person
    factorization in (:person, :person_node) || throw(ArgumentError(
        "unknown OMRF factorization: $factorization"))
    context = OMRFResidualContext(X, seen, threshold_indices,
        interaction_indices, threshold_starts, incident_edges,
        incident_neighbours,
        Matrix{Float64}(undef, no_thresholds, N),
        Matrix{Float64}(undef, P, N), zeros(maximum(seen) - 1),
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
    fill!(out, 0.0)
    @inbounds for index in context.touched_positions
        context.touched_mask[index] = false
    end
    empty!(context.touched_positions)
    N, P = size(context.X)
    α = if context.hcv_active
        displacement2 = dot(context.displacement, context.displacement)
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
    analytic_prior::Union{Nothing,OMRFIndependentAnalyticPrior}
end

function OMRFAnalyticSubsamplingOracle(ctx::StanSubsamplingContext,
        envelope_spec, unc_names, anchor::Vector{Float64})
    d = length(anchor)
    residual_context = _omrf_residual_context(
        envelope_spec, unc_names, anchor, ctx.counts)
    analytic_prior = _omrf_analytic_prior(envelope_spec, residual_context)
    oracle = OMRFAnalyticSubsamplingOracle(ctx, copy(anchor), zeros(d), zeros(d),
        zeros(d), zeros(d), false, residual_context, analytic_prior)
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
        α = context.hcv_damping / (context.hcv_damping + displacement2)
        mul!(context.hcv_product, something(context.likelihood_hessian),
            context.displacement)
        @. out += α * context.hcv_product
    end
    return out
end

function (oracle::OMRFAnalyticSubsamplingOracle)(out, x, subset, anchor)
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

@inline function _omrf_coordinate_at(state, flow::PDMPSamplers.AnyBoomerang,
        anchor, j, t)
    if !_omrf_coordinate_free(state, j)
        return state.ξ.x[j] - anchor[j], 0.0
    end
    s, c = sincos(t)
    centered = state.ξ.x[j] - flow.μ[j]
    displacement = flow.μ[j] - anchor[j] + centered * c + state.ξ.θ[j] * s
    velocity = -centered * s + state.ξ.θ[j] * c
    return displacement, velocity
end

@inline _omrf_coordinate_at(state, flow::PDMPSamplers.PreconditionedDynamics,
    anchor, j, t) = _omrf_coordinate_at(state, flow.dynamics, anchor, j, t)

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
        scales::Union{OMRFPatternLocalScales,OMRFCovariateLocalScales},
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
        scales::Union{OMRFPatternLocalScales,OMRFCovariateLocalScales},
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

function PDMPSamplers.subsampling_residual_subset_bound(
        oracle::OMRFAnalyticSubsamplingOracle, state, flow,
        D, M, subset, scale)
    context = oracle.residual_context
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

function _omrf_harmonic_coefficients(state, flow, anchor, j)
    if !_omrf_coordinate_free(state, j)
        return state.ξ.x[j] - anchor[j], 0.0, 0.0
    end
    return flow.μ[j] - anchor[j], state.ξ.x[j] - flow.μ[j], state.ξ.θ[j]
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
        scales::Union{OMRFPatternLocalScales,OMRFCovariateLocalScales},
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
    if _haskey(envelope_spec, :node_weights)
        bound_type = _haskey(envelope_spec, :bound_type) ?
            String(_rget(envelope_spec, :bound_type)) : "node_local_norm_low_rank"
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
end

function _prepare_omrf_anchor_envelope(manager::OMRFAnchorManager,
        requested, residual_context)
    spec = manager.envelope_spec
    bound_type = _haskey(spec, :bound_type) ?
        String(_rget(spec, :bound_type)) : ""
    if bound_type == "covariate_local_expansion" &&
            !residual_context.use_hcv && !isempty(manager.entries)
        template = manager.entries[1].state.envelope
        if template isa SeparableResidualEnvelope
            max_degree = maximum(length,
                residual_context.incident_edges; init=0)
            scales = OMRFCovariateLocalScales(requested, residual_context,
                zeros(max_degree), zeros(max_degree), zeros(length(requested)),
                zeros(length(requested)), zeros(length(requested)))
            # Weights, totals, and alias tables depend only on the observed
            # covariates. Anchor-bank entries share them; only trajectory
            # callbacks and their workspaces are anchor-specific.
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
            recycled.envelope isa SeparableResidualEnvelope &&
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
        0, 0.0, 0, 0, 0, 0, 0.0, 0.0, Inf,
        oracle.residual_context.hcv_active, nothing, 0)
    initial = OMRFPreparedAnchorState(oracle.anchor, oracle.full_anchor,
        oracle.prior_anchor, oracle.residual_context, initial_envelope)
    manager.active_idx = _insert_omrf_anchor!(manager, initial)
    return manager
end

function _finish_omrf_warmup!(manager::OMRFAnchorManager,
        cv::SubsampledControlVariate, position)
    previous = manager.active_idx
    if iszero(manager.preparations)
        prepared = _prepare_omrf_anchor(manager, position)
        manager.active_idx = _insert_omrf_anchor!(manager, prepared)
        PDMPSamplers.refresh_anchor!(cv, prepared.anchor)
    else
        # A scheduled warmup anchor is already an exact control variate. Reuse
        # the closest prepared entry instead of paying for a second anchor that
        # immediately supersedes it at the phase boundary.
        _select_omrf_anchor!(manager, cv, position; phase=:warmup)
    end
    anchor_changed = manager.active_idx != previous
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
    return anchor_changed || hcv_changed
end

@inline function _omrf_anchor_distance2(x, state::OMRFPreparedAnchorState)
    distance = 0.0
    @inbounds for j in eachindex(x, state.anchor)
        distance += (x[j] - state.anchor[j])^2
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
        cv::SubsampledControlVariate, x; phase::Symbol=:unknown)
    previous = manager.active_idx
    best_idx = firstindex(manager.entries)
    best_distance = Inf
    @inbounds for idx in eachindex(manager.entries)
        entry = manager.entries[idx]
        entry.age += 1
        distance = _omrf_anchor_distance2(x, entry.state)
        if distance < best_distance
            best_idx = idx
            best_distance = distance
        end
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
            return nothing
        end
    end
    if best_idx != previous
        PDMPSamplers.refresh_anchor!(cv,
            manager.entries[best_idx].state.anchor)
        phase === :main && (manager.main_activations += 1)
    end
    return nothing
end

function _add_omrf_anchor!(manager::OMRFAnchorManager,
        cv::SubsampledControlVariate, trace)
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
        (grad, x, phase) -> _select_omrf_anchor!(manager, grad, x; phase),
        (grad, trace) -> _add_omrf_anchor!(manager, grad, trace),
        (grad, state, args...) ->
            _finish_omrf_warmup!(manager, grad, state.ξ.x),
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
    warmup_mode == "full" && !isnothing(slab_prior) && throw(ArgumentError(
        "full-gradient subsampling warmup does not yet support dependent slab priors"))
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
    warmup_alg = warmup_mode == "full" ? wrap_sticky(warmup_alg0, sticky,
        model_prior, isnothing(parameter_prior) ? nothing :
            _as_float_vector(parameter_prior), can_stick_vec) : nothing
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
            warmup_adaptation_interval)
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
                warmup_adaptation_interval)
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
    result["subsampling_context_counters"] = [Dict(
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
        "initialization_model_constructions" => snapshot.model,
        "initialization_data_constructions" => snapshot.data,
        "initialization_full_gradient_calls" => snapshot.full_gradient,
        "sampling_model_constructions" => ctx.counts.model_constructions - snapshot.model,
        "sampling_data_constructions" => ctx.counts.data_constructions - snapshot.data,
        "sampling_full_gradient_calls" => ctx.counts.full_gradient - snapshot.full_gradient,
        "anchor_preparations" => manager === nothing ? 0 : manager.preparations,
        "anchor_recycled_preparations" => manager === nothing ? 0 :
            manager.recycled_preparations,
        "anchor_preparation_seconds" => manager === nothing ? 0.0 : manager.preparation_seconds,
        "anchor_activations" => manager === nothing ? 0 : manager.activations,
        "anchor_main_activations" => manager === nothing ? 0 : manager.main_activations,
        "anchor_main_refreshes" => manager === nothing ? 0 : manager.main_refreshes,
        "anchor_main_selections" => manager === nothing ? 0 : manager.main_selections,
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
    result["subsampling"] = true
    return result
end

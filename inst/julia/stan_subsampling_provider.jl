# Persistent selected-gradient backend for public custom-Stan subsampling sampling.


mutable struct StanSubsamplingCallCounts
    full_gradient::Int
    prior_gradient::Int
    selected_gradient::Int
    persons_evaluated::Int
    anchor_cache_gradient::Int
    anchor_cache_persons::Int
    analytic_residual::Int
    selected_refinement_calls::Int
    analytic_persons::Int
    analytic_node_conditionals::Int
    model_constructions::Int
    data_constructions::Int
    native_full_data_scans::Int
    analytic_residual_seconds::Float64
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
    compact_anchor_likelihood_seconds::Float64
    direct_factor_bound_scans::Int
    direct_factors_evaluated::Int
    direct_cell_certificate_seconds::Float64
    direct_cached_cell_hits::Int
    direct_selected_factor_exact_evaluations::Int
    direct_path_preparations::Int
    direct_path_factors_prepared::Int
    direct_path_preparation_seconds::Float64
    direct_adjacent_endpoint_reuses::Int
    direct_left_endpoint_evaluations::Int
    direct_right_endpoint_evaluations::Int
    direct_residual_table_constructions::Int
    direct_internal_cell_uses::Int
    deterministic_hvp_calls::Int
    deterministic_hvp_seconds::Float64
    last_mark_residual::Bool
    deterministic_source_proposals::Int
    residual_source_proposals::Int
    deterministic_source_acceptances::Int
    residual_source_acceptances::Int
    shared_current_pattern_evaluations::Int
    shared_current_cell_evaluations::Int
    shared_current_rebuilt_nodes::Int
    shared_current_rebuild_seconds::Float64
    shared_current_dense_fallback_cells::Int
    joint_cell_width_bins::Vector{Int}
    direct_path_anchor_changes::Int
    direct_path_mask_changes::Int
    direct_path_flow_changes::Int
    direct_path_state_changes::Int
    direct_path_forced_rebuilds::Int
end
StanSubsamplingCallCounts() = StanSubsamplingCallCounts(
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0, 0.0, 0.0, 0.0, 0.0, 0, 0, 0.0, 0, 0, 0, 0, 0.0, 0, 0, 0, 0, 0,
    0, 0.0, false, 0, 0, 0, 0, 0, 0, 0, 0.0, 0, zeros(Int, 10), 0, 0, 0, 0, 0)

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

# ---------------------------------------------------------------------------
# Reuse of products with the anchor likelihood Hessian H.
#
# Global HCV evaluates H * w for many vectors w that repeat or move along the
# current flight.  `OMRFHMemo` returns the stored product when w is
# bit-identical to its key (exact reuse, e.g. H * v for a fixed velocity).
# `OMRFHAffine` keeps H * base and the products of up to two directions and
# answers any displacement in that affine span in O(d); a linear flight needs
# one direction and survives one velocity change, a harmonic flight needs two.
# Entries are tied to the Hessian generation, which `_omrf_likelihood_hessian!`
# advances.  `PDMPSAMPLERS_HCV_CACHE=0` disables reuse for comparison.
# ---------------------------------------------------------------------------
mutable struct OMRFHMemo
    generation::Int
    key::Vector{Float64}
    product::Vector{Float64}
    updates::Int
    changed::Vector{Int}
end
OMRFHMemo(d::Int) = OMRFHMemo(-1, zeros(d), zeros(d), 0, sizehint!(Int[], d))

mutable struct OMRFHAffine
    generation::Int
    base::Vector{Float64}
    pbase::Vector{Float64}
    u::Matrix{Float64}
    pu::Matrix{Float64}
    umax::Vector{Float64}
    ndir::Int
    r::Vector{Float64}
end
OMRFHAffine(d::Int) = OMRFHAffine(-1, zeros(d), zeros(d), zeros(d, 2),
    zeros(d, 2), zeros(2), 0, zeros(d))

mutable struct OMRFHCVProductCache
    enabled::Bool
    generation::Int
    velocity::OMRFHMemo      # H * (free part of v)
    state::OMRFHMemo         # H * (x - anchor)
    offset::OMRFHMemo        # H * (x - mu) on Boomerang flights
    reference::OMRFHMemo     # H * (mu - anchor) on Boomerang flights
    path_offset::OMRFHMemo   # H * (x - anchor) at the start of linear paths
    displacement::OMRFHAffine    # x - anchor along linear flights
    harmonic::OMRFHAffine        # x - mu and v along a Boomerang orbit
    direct::Int
    reused::Int
end
OMRFHCVProductCache(d::Int) = OMRFHCVProductCache(
    get(ENV, "PDMPSAMPLERS_HCV_CACHE", "1") != "0", 0,
    OMRFHMemo(d), OMRFHMemo(d), OMRFHMemo(d), OMRFHMemo(d), OMRFHMemo(d),
    OMRFHAffine(d), OMRFHAffine(d), 0, 0)

_omrf_hcv_cache_invalidate!(cache::OMRFHCVProductCache) =
    (cache.generation += 1; cache)

# H * w from the memo: reused when w equals the key, updated column by column
# when w differs in a few coordinates (a Zig-Zag flip, or a coordinate that
# sticks or unsticks; H is symmetric so columns are rows), otherwise
# recomputed.  A full product is forced after 64 updates to bound rounding
# drift.
function _omrf_hmemo(cache::OMRFHCVProductCache, memo::OMRFHMemo, H, w)
    if cache.enabled && memo.generation == cache.generation
        changed = memo.changed
        empty!(changed)
        limit = max(4, length(w) >> 4)
         for j in eachindex(w)
            if memo.key[j] != w[j]
                push!(changed, j)
                length(changed) > limit && break
            end
        end
        if isempty(changed)
            cache.reused += 1
            return memo.product
        end
        if length(changed) <= limit && memo.updates < 64
             for j in changed
                step = w[j] - memo.key[j]
                for i in axes(H, 1)
                    memo.product[i] += step * H[i, j]
                end
                memo.key[j] = w[j]
            end
            memo.updates += 1
            cache.reused += 1
            return memo.product
        end
    end
    mul!(memo.product, H, w)
    copyto!(memo.key, w)
    memo.generation = cache.enabled ? cache.generation : -1
    memo.updates = 0
    cache.direct += 1
    return memo.product
end

# Index of the largest |u[i, column]|, without allocating.
function _omrf_pivot(u, column)
    best = -1.0
    k = 1
    @inbounds for i in axes(u, 1)
        a = abs(u[i, column])
        if a > best
            best = a
            k = i
        end
    end
    return k
end

# H * (x - anchor) at the start of a linear path.  The previous path started at
# x_prev and moved with the velocity v stored in `cache.velocity`, so
# x - anchor = (x_prev - anchor) + tau * v and the product follows from
# H v in O(d).  tau is solved at the largest |v| and the identity is checked
# coordinatewise to a few ulps; otherwise, and after 64 updates, the product is
# recomputed.  Call before the velocity memo is updated for the new path.
function _omrf_path_offset!(out, cache::OMRFHCVProductCache, H, delta)
    memo = cache.path_offset
    velocity = cache.velocity
    if cache.enabled && memo.generation == cache.generation &&
            velocity.generation == cache.generation && memo.updates < 64
        v = velocity.key
        k = 1
        best = -1.0
        @inbounds for i in eachindex(v)
            a = abs(v[i])
            if a > best
                best = a
                k = i
            end
        end
        if best > 0
            tau = (delta[k] - memo.key[k]) / v[k]
            scale = maximum(abs, delta) + maximum(abs, memo.key) + abs(tau) * best
            tol = 32 * eps(Float64) * scale
            ok = isfinite(tau)
            @inbounds for i in eachindex(v)
                ok || break
                ok = abs(delta[i] - memo.key[i] - tau * v[i]) <= tol
            end
            if ok
                @inbounds for i in eachindex(out)
                    memo.product[i] += tau * velocity.product[i]
                end
                copyto!(memo.key, delta)
                memo.updates += 1
                copyto!(out, memo.product)
                cache.reused += 1
                return out
            end
        end
    end
    mul!(memo.product, H, delta)
    copyto!(memo.key, delta)
    memo.generation = cache.enabled ? cache.generation : -1
    memo.updates = 0
    copyto!(out, memo.product)
    cache.direct += 1
    return out
end

function _omrf_haffine_reset!(out, cache, A::OMRFHAffine, H, delta)
    mul!(A.pbase, H, delta)
    copyto!(A.base, delta)
    A.generation = cache.generation
    A.ndir = 0
    copyto!(out, A.pbase)
    cache.direct += 1
    return out
end

# out = H * delta.  When delta - base lies in the span of the stored directions
# (checked coordinatewise to a few ulps of the magnitudes involved) the
# product is assembled from stored products in O(d); otherwise the new offset
# becomes a stored direction, or the cache restarts at delta.
function _omrf_haffine!(out, cache::OMRFHCVProductCache, H, delta,
        A::OMRFHAffine=cache.displacement)
    if !cache.enabled
        mul!(out, H, delta)
        cache.direct += 1
        return out
    end
    A.generation == cache.generation ||
        return _omrf_haffine_reset!(out, cache, A, H, delta)
    r = A.r
    @. r = delta - A.base
    rmax = maximum(abs, r)
    if iszero(rmax)
        copyto!(out, A.pbase)
        cache.reused += 1
        return out
    end
    s = 0.0
    t = 0.0
    solved = false
    if A.ndir == 1
        k = _omrf_pivot(A.u, 1)
        s = r[k] / A.u[k, 1]
        solved = isfinite(s)
    elseif A.ndir == 2
        k1 = _omrf_pivot(A.u, 1)
        ratio = A.u[k1, 2] / A.u[k1, 1]
        best = 0.0
        k2 = 0
        @inbounds for i in axes(A.u, 1)
            e = abs(A.u[i, 2] - ratio * A.u[i, 1])
            if e > best
                best = e
                k2 = i
            end
        end
        if k2 > 0
            a11 = A.u[k1, 1]; a12 = A.u[k1, 2]
            a21 = A.u[k2, 1]; a22 = A.u[k2, 2]
            det = a11 * a22 - a12 * a21
            s = (r[k1] * a22 - a12 * r[k2]) / det
            t = (a11 * r[k2] - a21 * r[k1]) / det
            solved = isfinite(s) && isfinite(t)
        end
    end
    if solved
        scale = maximum(abs, delta) + maximum(abs, A.base) +
            abs(s) * A.umax[1] + abs(t) * A.umax[2]
        tol = 32 * eps(Float64) * scale
        ok = true
        @inbounds for i in eachindex(r)
            if abs(r[i] - s * A.u[i, 1] - t * A.u[i, 2]) > tol
                ok = false
                break
            end
        end
        if ok
            @inbounds for i in eachindex(out)
                out[i] = A.pbase[i] + s * A.pu[i, 1] + t * A.pu[i, 2]
            end
            cache.reused += 1
            return out
        end
    end
    A.ndir < 2 || return _omrf_haffine_reset!(out, cache, A, H, delta)
    A.ndir += 1
    column = A.ndir
    copyto!(view(A.u, :, column), r)
    mul!(view(A.pu, :, column), H, r)
    A.umax[column] = rmax
    out .= A.pbase .+ view(A.pu, :, column)
    cache.direct += 1
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
    anchor_reference_probabilities::Matrix{Float64}
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
    displacement::Vector{Float64}
    hcv_product::Vector{Float64}
    node_supports::Vector{Vector{Int}}
    touched_positions::Vector{Int}
    touched_mask::BitVector
    counts::StanSubsamplingCallCounts
    active_free::BitVector
    active_generation::Int
    hcv_cache::OMRFHCVProductCache
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
    # Two variants of the node-shared slab, sd = scale * global * sqrt(node_j node_k):
    #   half-normal global and node scales (`global_scale_index > 0`,
    #   `lognormal_nodes = false`; inst/stan/omrf/omrf_node_shared.stan), or
    #   no global scale and lognormal(0, node_scale_sd) node scales
    #   (`global_scale_index == 0`, `lognormal_nodes = true`).
    lognormal_nodes::Bool
end

const OMRFAnalyticPrior = Union{OMRFIndependentAnalyticPrior,
    OMRFNodeScaleAnalyticPrior}

# A field that is present and not R's NULL.
_has_value(x, key::Symbol) = _haskey(x, key) && _rget(x, key) !== nothing

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
        lognormal_nodes = _has_value(spec, :node_scale_prior) &&
            String(_rget(spec, :node_scale_prior)) == "lognormal"
        lognormal_nodes == (context.global_scale_index == 0) || throw(ArgumentError(
            "node-scale OMRF prior: lognormal node scales require no global scale, " *
            "half-normal node scales require one"))
        global_scale_sd = lognormal_nodes ? 1.0 : Float64(_rget(spec, :global_scale_sd))
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
            zeros(length(context.node_scale_indices)), lognormal_nodes)
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
    global_scale = global_index > 0 ? exp(x[global_index]) : 1.0
    global_sum = global_index > 0 ? global_scale^2 / prior.global_scale_sd^2 - 1.0 : 0.0
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
    global_index > 0 && (out[global_index] = global_sum)
    @inbounds for node in eachindex(node_sums)
        scale = node_scales[node]
        node_index = prior.node_scale_indices[node]
        # -d log p / d log(scale): lognormal, or half-normal with its Jacobian.
        out[node_index] = node_sums[node] + (prior.lognormal_nodes ?
            x[node_index] / prior.node_scale_sd^2 :
            scale^2 / prior.node_scale_sd^2 - 1.0)
    end
    return out
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
            reference_probability = denominator
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
            if cache_anchor
                context.anchor_expected_scores[node, n] = expected
                context.anchor_reference_probabilities[node, n] =
                    reference_probability / denominator
            end
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


@inline function _omrf_signed_interval_upper(glo, ghi, vlo, vhi)
    upper = 0.0
    scale = 0.0
    @inbounds for i in eachindex(glo, ghi, vlo, vhi)
        _, hi = _omrf_imul(vlo[i], vhi[i], glo[i], ghi[i])
        upper += hi
        scale += abs(hi)
    end
    return upper, scale
end


"""Refresh an OMRF anchor cache and its exact full gradient in one pass."""
function _cache_omrf_anchor_and_full_gradient!(context::OMRFResidualContext,
        full, prior::OMRFAnalyticPrior, anchor)
    _omrf_prior_gradient!(full, prior, anchor)
    _omrf_add_likelihood_gradient!(full, context, anchor; cache_anchor=true)
    context.counts.native_full_data_scans += 1
    return full
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
    global_scale = global_index > 0 ? exp(x[global_index]) : 1.0
    global_direction = global_index > 0 ? direction[global_index] : 0.0
    global_index > 0 && (out[global_index] = 2 * global_scale^2 /
        prior.global_scale_sd^2 * global_direction)

    node_scales = prior.node_scales
    @inbounds for node in eachindex(node_scales)
        node_scales[node] = exp(x[prior.node_scale_indices[node]])
        node_index = prior.node_scale_indices[node]
        out[node_index] = (prior.lognormal_nodes ? 1.0 :
            2 * node_scales[node]^2) / prior.node_scale_sd^2 * direction[node_index]
    end

    edge = 1
    @inbounds for left in 1:(length(node_scales) - 1), right in (left + 1):length(node_scales)
        interaction_index = prior.interaction_indices[edge]
        interaction = x[interaction_index]
        scale2 = prior.interaction_scale^2 * global_scale^2 *
            node_scales[left] * node_scales[right]
        inv_scale2 = inv(scale2)
        ratio = interaction^2 * inv_scale2
        # The edge variance is proportional to
        # `global_scale^2 * node_scale[left] * node_scale[right]`.  Its
        # log-variance directional derivative therefore has coefficient two
        # only for the global log scale; the incident node log scales enter
        # once each.
        log_variance_direction = 2 * global_direction +
            direction[prior.node_scale_indices[left]] +
            direction[prior.node_scale_indices[right]]
        out[interaction_index] = inv_scale2 * direction[interaction_index] -
            interaction * inv_scale2 * log_variance_direction
        dc = -2 * interaction * inv_scale2 * direction[interaction_index] +
            ratio * log_variance_direction
        global_index > 0 && (out[global_index] += dc)
        out[prior.node_scale_indices[left]] += 0.5 * dc
        out[prior.node_scale_indices[right]] += 0.5 * dc
        edge += 1
    end
    return out
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

# Direct interval certificate for the position-constant first-order OMRF
# control variate.  Every primitive is rounded outwards by one representable
# Float64 value.  This is deliberately restricted to the ordinary CV: an HCV
# has an additional position-dependent deterministic likelihood term.
# x -+ (|x| 2^-52 + 2^-1074) lies beyond prevfloat/nextfloat(x) for every finite
# x, so this is a (marginally wider) outward rounding without the branches of
# prevfloat/nextfloat.
@inline _omrf_down(x::Float64) =
    ifelse(isfinite(x), x - muladd(abs(x), 0x1p-52, 0x1p-1074), x)
@inline _omrf_up(x::Float64) =
    ifelse(isfinite(x), x + muladd(abs(x), 0x1p-52, 0x1p-1074), x)
@inline _omrf_iadd(al, ah, bl, bh) =
    (_omrf_down(al + bl), _omrf_up(ah + bh))
@inline _omrf_iscale(al, ah, c) = c >= 0 ?
    (_omrf_down(c * al), _omrf_up(c * ah)) :
    (_omrf_down(c * ah), _omrf_up(c * al))
@inline function _omrf_imul(al, ah, bl, bh)
    p1, p2, p3, p4 = al * bl, al * bh, ah * bl, ah * bh
    return _omrf_down(min(p1, p2, p3, p4)),
        _omrf_up(max(p1, p2, p3, p4))
end
@inline function _omrf_isquare(al, ah)
    lo = al <= 0 <= ah ? 0.0 : min(al * al, ah * ah)
    hi = max(al * al, ah * ah)
    return _omrf_down(lo), _omrf_up(hi)
end
@inline function _omrf_iexp(al, ah)
    lo, hi = exp(al), exp(ah)
    isfinite(lo) && isfinite(hi) || throw(OverflowError(
        "OMRF deterministic interval exponential overflow; subdivide the cell"))
    return _omrf_down(lo), _omrf_up(hi)
end

@inline function _omrf_harmonic_range(offset::Float64, cosine::Float64,
        sine::Float64, left::Float64, right::Float64)
    left <= right || throw(ArgumentError("cell endpoints are reversed"))
    function endpoint(t)
        cl, ch = _omrf_down(cos(t)), _omrf_up(cos(t))
        sl, sh = _omrf_down(sin(t)), _omrf_up(sin(t))
        al, ah = _omrf_iscale(cl, ch, cosine)
        bl, bh = _omrf_iscale(sl, sh, sine)
        yl, yh = _omrf_iadd(al, ah, bl, bh)
        return _omrf_iadd(yl, yh, offset, offset)
    end
    ll, lh = endpoint(left); rl, rh = endpoint(right)
    lo = min(ll, rl); hi = max(lh, rh)
    if !iszero(cosine) || !iszero(sine)
        phase = atan(sine, cosine)
        # Search one critical point beyond either rounded quotient and admit a
        # small ULP neighbourhood at the cell boundary.  At an admitted
        # critical point use offset +/- hypot directly; evaluating sin/cos at
        # the rounded critical time could otherwise underestimate an extremum.
        first = floor(Int, (left - phase) / π) - 1
        last = ceil(Int, (right - phase) / π) + 1
        radius = hypot(cosine, sine)
        tolerance = 16eps(max(abs(left), abs(right), 1.0))
        for k in first:last
            critical = phase + k * π
            if left - tolerance <= critical <= right + tolerance
                y = iseven(k) ? offset + radius : offset - radius
                lo = min(lo, _omrf_down(y)); hi = max(hi, _omrf_up(y))
            end
        end
    end
    return _omrf_down(lo), _omrf_up(hi)
end

@inline _omrf_actual_flow(flow) = flow isa PDMPSamplers.PreconditionedDynamics ?
    flow.dynamics : flow
@inline _omrf_is_free(state, i) = !hasproperty(state, :free) || state.free[i]

function _omrf_fill_path_intervals!(xlo, xhi, vlo, vhi, state, flow,
        left::Float64, right::Float64)
    actual = _omrf_actual_flow(flow)
    @inbounds for i in eachindex(xlo)
        if !_omrf_is_free(state, i)
            xlo[i] = _omrf_down(state.ξ.x[i]); xhi[i] = _omrf_up(state.ξ.x[i])
            vlo[i] = -0.0; vhi[i] = 0.0
        elseif actual isa PDMPSamplers.AnyBoomerang
            y = state.ξ.x[i] - actual.μ[i]
            xlo[i], xhi[i] = _omrf_harmonic_range(
                actual.μ[i], y, state.ξ.θ[i], left, right)
            vlo[i], vhi[i] = _omrf_harmonic_range(
                0.0, state.ξ.θ[i], -y, left, right)
        else
            ll, lh = _omrf_iscale(state.ξ.θ[i], state.ξ.θ[i], left)
            rl, rh = _omrf_iscale(state.ξ.θ[i], state.ξ.θ[i], right)
            ll, lh = _omrf_iadd(ll, lh, state.ξ.x[i], state.ξ.x[i])
            rl, rh = _omrf_iadd(rl, rh, state.ξ.x[i], state.ξ.x[i])
            xlo[i] = min(ll, rl); xhi[i] = max(lh, rh)
            vlo[i] = _omrf_down(state.ξ.θ[i]); vhi[i] = _omrf_up(state.ξ.θ[i])
        end
    end
    return nothing
end

function _omrf_prior_intervals!(glo, ghi, prior::OMRFIndependentAnalyticPrior,
        xlo, xhi, active=nothing)
    fill!(glo, 0.0); fill!(ghi, 0.0)
    α, β = prior.threshold_alpha, prior.threshold_beta
    @inbounds for i in prior.threshold_indices
        gl = (α + β) * _omrf_logistic(xlo[i]) - α
        gh = (α + β) * _omrf_logistic(xhi[i]) - α
        glo[i] = _omrf_down(gl); ghi[i] = _omrf_up(gh)
    end
    s = prior.interaction_scale
    s2 = s * s
    @inbounds for i in prior.interaction_indices
        active !== nothing && !active[i] && continue
        if prior.interaction_type === :gaussian
            glo[i] = _omrf_down(xlo[i] / s2)
            ghi[i] = _omrf_up(xhi[i] / s2)
        else
            f(x) = 2x / (s2 + x * x)
            lo = min(f(xlo[i]), f(xhi[i])); hi = max(f(xlo[i]), f(xhi[i]))
            xlo[i] <= -s <= xhi[i] && (lo = min(lo, -inv(s)))
            xlo[i] <=  s <= xhi[i] && (hi = max(hi,  inv(s)))
            glo[i] = _omrf_down(lo); ghi[i] = _omrf_up(hi)
        end
    end
    return nothing
end

function _omrf_prior_intervals!(glo, ghi, prior::OMRFNodeScaleAnalyticPrior,
        xlo, xhi, active)
    fill!(glo, 0.0); fill!(ghi, 0.0)
    α, β = prior.threshold_alpha, prior.threshold_beta
    @inbounds for i in prior.threshold_indices
        glo[i] = _omrf_down((α + β) * _omrf_logistic(xlo[i]) - α)
        ghi[i] = _omrf_up((α + β) * _omrf_logistic(xhi[i]) - α)
    end
    gi = prior.global_scale_index
    if gi > 0
        hglo, hghi = _omrf_iexp(2xlo[gi] - 2log(prior.global_scale_sd),
                                 2xhi[gi] - 2log(prior.global_scale_sd))
        glo[gi], ghi[gi] = _omrf_down(hglo - 1), _omrf_up(hghi - 1)
    end
    @inbounds for (node, ni) in pairs(prior.node_scale_indices)
        if prior.lognormal_nodes
            s2 = prior.node_scale_sd^2
            glo[ni], ghi[ni] = _omrf_down(xlo[ni] / s2), _omrf_up(xhi[ni] / s2)
        else
            hlo, hhi = _omrf_iexp(2xlo[ni] - 2log(prior.node_scale_sd),
                                   2xhi[ni] - 2log(prior.node_scale_sd))
            glo[ni], ghi[ni] = _omrf_down(hlo - 1), _omrf_up(hhi - 1)
        end
    end
    edge = 1
    @inbounds for u in 1:(length(prior.node_scale_indices) - 1),
            v in (u + 1):length(prior.node_scale_indices)
        bi = prior.interaction_indices[edge]
        if active === nothing || active[bi]
            eta_lo = log(prior.interaction_scale) + (gi > 0 ? xlo[gi] : 0.0) +
                0.5xlo[prior.node_scale_indices[u]] +
                0.5xlo[prior.node_scale_indices[v]]
            eta_hi = log(prior.interaction_scale) + (gi > 0 ? xhi[gi] : 0.0) +
                0.5xhi[prior.node_scale_indices[u]] +
                0.5xhi[prior.node_scale_indices[v]]
            wlo, whi = _omrf_iexp(-2eta_hi, -2eta_lo)
            glo[bi], ghi[bi] = _omrf_imul(xlo[bi], xhi[bi], wlo, whi)
            b2lo, b2hi = _omrf_isquare(xlo[bi], xhi[bi])
            rlo, rhi = _omrf_imul(b2lo, b2hi, wlo, whi)
            elo, ehi = _omrf_down(1 - rhi), _omrf_up(1 - rlo)
            gi > 0 && ((glo[gi], ghi[gi]) = _omrf_iadd(glo[gi], ghi[gi], elo, ehi))
            for node in (u, v)
                ni = prior.node_scale_indices[node]
                sl, sh = _omrf_iscale(elo, ehi, 0.5)
                glo[ni], ghi[ni] = _omrf_iadd(glo[ni], ghi[ni], sl, sh)
            end
        end
        edge += 1
    end
    return nothing
end

function _omrf_interval_rate_bound(likelihood_baseline, prior_baseline,
        glo, ghi, xlo, xhi, vlo, vhi,
        state, flow)
    actual = _omrf_actual_flow(flow)
    zigzag = actual isa ZigZag
    signed_lo = 0.0; signed_hi = 0.0; total = 0.0
    indices = likelihood_baseline === nothing ? eachindex(glo) :
        eachindex(likelihood_baseline)
    @inbounds for i in indices
        c = likelihood_baseline === nothing ? 0.0 : likelihood_baseline[i] -
            (prior_baseline === nothing ? 0.0 : prior_baseline[i])
        gl, gh = _omrf_iadd(glo[i], ghi[i], c, c)
        if actual isa PDMPSamplers.AnyBoomerang
            gamma = actual.Γ isa Diagonal ? actual.Γ.diag[i] : throw(ArgumentError(
                "direct OMRF deterministic certificate requires diagonal Boomerang precision"))
            rl, rh = _omrf_iscale(xlo[i] - actual.μ[i],
                                   xhi[i] - actual.μ[i], gamma)
            gl, gh = _omrf_iadd(gl, gh, -rh, -rl)
        end
        pl, ph = _omrf_imul(vlo[i], vhi[i], gl, gh)
        if zigzag
            total = _omrf_up(total + max(0.0, ph))
        else
            signed_lo = _omrf_down(signed_lo + pl)
            signed_hi = _omrf_up(signed_hi + ph)
        end
    end
    result = zigzag ? total : max(0.0, signed_hi)
    isfinite(result) || throw(OverflowError(
        "OMRF deterministic cell rate is not finite; subdivide or reject the state"))
    return _omrf_up(result)
end

function _omrf_direct_deterministic_cell_bound(oracle, active, state, flow,
        left::Real, right::Real)
    context = oracle.residual_context
    xlo, xhi = oracle.interval_xlo, oracle.interval_xhi
    vlo, vhi = oracle.interval_vlo, oracle.interval_vhi
    glo, ghi = oracle.interval_glo, oracle.interval_ghi
    _omrf_fill_path_intervals!(xlo, xhi, vlo, vhi, state, flow,
        Float64(left), Float64(right))
    prior = something(oracle.analytic_prior)
    _omrf_prior_intervals!(glo, ghi, prior, xlo, xhi, active)
    if context.hcv_active
        H = something(context.likelihood_hessian)
        anchor = oracle.anchor
        if _omrf_actual_flow(flow) isa BouncyParticle
            # On a fixed linear sticky stratum the HCV likelihood term is
            # affine and H is a likelihood Hessian, hence positive
            # semidefinite.  Its signed BPS rate has derivative v'Hv >= 0 and
            # is maximized at the right endpoint.  Evaluate it there exactly
            # and combine it with the separately certified prior interval.
            prior_upper, _ = _omrf_signed_interval_upper(
                glo, ghi, vlo, vhi)
            t = Float64(right)
            @inbounds for j in eachindex(anchor)
                context.displacement[j] = _omrf_is_free(state, j) ?
                    muladd(t, state.ξ.θ[j], state.ξ.x[j]) - anchor[j] :
                    state.ξ.x[j] - anchor[j]
            end
            # H is symmetric, so v_f' H d = (H v_f)' d; H v_f is fixed along
            # the flight and is reused until the velocity or stratum changes.
            @inbounds for j in eachindex(anchor)
                context.hcv_product[j] = _omrf_is_free(state, j) ?
                    state.ξ.θ[j] : 0.0
            end
            Hv = _omrf_hmemo(context.hcv_cache, context.hcv_cache.velocity,
                H, context.hcv_product)
            likelihood_rate = 0.0
            likelihood_absolute = 0.0
            @inbounds for i in eachindex(anchor)
                term = Hv[i] * context.displacement[i]
                if _omrf_is_free(state, i)
                    anchor_term = state.ξ.θ[i] *
                        (oracle.full_anchor[i] - oracle.prior_anchor[i])
                    term += anchor_term
                    likelihood_absolute += abs(anchor_term)
                end
                likelihood_rate += term
                likelihood_absolute += abs(Hv[i] * context.displacement[i])
            end
            # Reassociating the product changes rounding; cover it.
            likelihood_rate += 1e-10 * likelihood_absolute
            result = max(0.0, likelihood_rate + prior_upper)
            isfinite(result) || throw(OverflowError(
                "OMRF HCV BPS deterministic rate is not finite"))
            return nextfloat(result)
        end
        if _omrf_actual_flow(flow) isa ZigZag
            # Each HCV likelihood-gradient coordinate is affine along the
            # linear Zig-Zag stratum.  Preserve that shared-time dependence
            # exactly instead of interval-multiplying H by independent
            # coordinate boxes, which discards cancellations and badly
            # inflates the sum of coordinatewise positive rates.
            @inbounds for j in eachindex(anchor)
                context.displacement[j] = state.ξ.x[j] - anchor[j]
            end
            copyto!(xlo, _omrf_hmemo(context.hcv_cache,
                context.hcv_cache.state, H, context.displacement))
            @inbounds for j in eachindex(anchor)
                context.displacement[j] = _omrf_is_free(state, j) ?
                    state.ξ.θ[j] : 0.0
            end
            copyto!(xhi, _omrf_hmemo(context.hcv_cache,
                context.hcv_cache.velocity, H, context.displacement))
            l = Float64(left)
            r = Float64(right)
            @inbounds for i in eachindex(anchor)
                baseline = oracle.full_anchor[i] - oracle.prior_anchor[i]
                gl = baseline + xlo[i] + l * xhi[i]
                gr = baseline + xlo[i] + r * xhi[i]
                glo[i] = _omrf_down(glo[i] + min(gl, gr))
                ghi[i] = _omrf_up(ghi[i] + max(gl, gr))
            end
            return _omrf_interval_rate_bound(nothing, nothing,
                glo, ghi, xlo, xhi, vlo, vhi, state, flow)
        end
        if _omrf_actual_flow(flow) isa PDMPSamplers.AnyBoomerang
            actual = _omrf_actual_flow(flow)
            actual.Γ isa Diagonal || throw(ArgumentError(
                "direct HCV Boomerang certificate requires diagonal precision"))
            prior_upper, _ = _omrf_signed_interval_upper(
                glo, ghi, vlo, vhi)

            # On the free stratum write d(t)=c+u*cos(t)+w*sin(t) and
            # v(t)=-u*sin(t)+w*cos(t).  The HCV likelihood minus Gaussian
            # reference gradient has coefficients g0, gc and gs.  Their
            # scalar product is a degree-two trigonometric polynomial.  Bound
            # its first and second harmonics separately using exact closed
            # harmonic intervals; this is much tighter than multiplying
            # coordinate boxes and needs only three dense HCV products.
            u = xlo
            w = xhi
            g0 = glo
            gc = ghi
            gs = vlo
            @inbounds for j in eachindex(anchor)
                if _omrf_is_free(state, j)
                    u[j] = state.ξ.x[j] - actual.μ[j]
                    w[j] = state.ξ.θ[j]
                    context.displacement[j] = actual.μ[j] - anchor[j]
                else
                    u[j] = 0.0
                    w[j] = 0.0
                    context.displacement[j] = state.ξ.x[j] - anchor[j]
                end
            end
            cache = context.hcv_cache
            Hd = _omrf_hmemo(cache, cache.reference, H, context.displacement)
            @inbounds for i in eachindex(anchor)
                g0[i] = oracle.full_anchor[i] - oracle.prior_anchor[i] + Hd[i]
            end
            Hu = _omrf_hmemo(cache, cache.offset, H, u)
            @inbounds for i in eachindex(anchor)
                gc[i] = Hu[i] - actual.Γ.diag[i] * u[i]
            end
            Hw = _omrf_hmemo(cache, cache.velocity, H, w)
            @inbounds for i in eachindex(anchor)
                gs[i] = Hw[i] - actual.Γ.diag[i] * w[i]
            end
            A1 = dot(w, g0)
            B1 = -dot(u, g0)
            ugs = dot(u, gs)
            wgc = dot(w, gc)
            C = (-ugs + wgc) / 2
            A2 = (ugs + wgc) / 2
            B2 = (-dot(u, gc) + dot(w, gs)) / 2
            _, first_upper = _omrf_harmonic_cell_interval(
                0.0, A1, B1, Float64(left), Float64(right))
            _, second_upper = _omrf_harmonic_cell_interval(
                0.0, A2, B2, 2Float64(left), 2Float64(right))
            result = max(0.0, prior_upper + C + first_upper + second_upper)
            isfinite(result) || throw(OverflowError(
                "OMRF HCV Boomerang deterministic rate is not finite"))
            return nextfloat(result)
        end
        # The deterministic HCV likelihood gradient is exactly
        #   full(anchor) - prior(anchor) + H(anchor) * (x-anchor).
        # Add an interval matrix-vector product to the exact analytic-prior
        # intervals.  This is valid for linear, harmonic and sticky strata and
        # retains each coordinate separately for Zig-Zag positive-part rates.
        @inbounds for i in eachindex(glo, ghi)
            lo = oracle.full_anchor[i] - oracle.prior_anchor[i]
            hi = lo
            for j in eachindex(anchor)
                hij = H[i, j]
                iszero(hij) && continue
                dl = _omrf_down(xlo[j] - anchor[j])
                dh = _omrf_up(xhi[j] - anchor[j])
                tl, th = _omrf_iscale(dl, dh, hij)
                lo = _omrf_down(lo + tl)
                hi = _omrf_up(hi + th)
            end
            glo[i] = _omrf_down(glo[i] + lo)
            ghi[i] = _omrf_up(ghi[i] + hi)
        end
        return _omrf_interval_rate_bound(nothing, nothing,
            glo, ghi, xlo, xhi, vlo, vhi, state, flow)
    end
    return _omrf_interval_rate_bound(oracle.full_anchor, oracle.prior_anchor,
        glo, ghi, xlo, xhi,
        vlo, vhi, state, flow)
end

"""Signed deterministic HCV/prior upper bound on a linear BPS cell.

The likelihood HCV contribution is affine with nonnegative directional
derivative `v'Hv`, so its maximum is at the right endpoint. The active-stratum
analytic prior is enclosed coordinatewise on the whole cell. Frozen
coordinates have zero physical velocity and their prior contribution is
excluded through the parameter-indexed active mask.
"""
function _omrf_joint_signed_deterministic_upper(oracle, active, state, flow,
        left::Float64, right::Float64)
    _omrf_actual_flow(flow) isa BouncyParticle || throw(ArgumentError(
        "joint signed deterministic bound supports linear BPS only"))
    context = oracle.residual_context
    context.hcv_active || throw(ArgumentError(
        "joint signed bound requires the installed global HCV"))
    _omrf_fill_path_intervals!(oracle.interval_xlo, oracle.interval_xhi,
        oracle.interval_vlo, oracle.interval_vhi, state, flow, left, right)
    _omrf_prior_intervals!(oracle.interval_glo, oracle.interval_ghi,
        something(oracle.analytic_prior), oracle.interval_xlo,
        oracle.interval_xhi, active)
    prior_upper, _ = _omrf_signed_interval_upper(oracle.interval_glo,
        oracle.interval_ghi, oracle.interval_vlo, oracle.interval_vhi)
    H = something(context.likelihood_hessian)
    anchor = oracle.anchor
    # H (x + right v_f - anchor) = H (x - anchor) + right H v_f.  Both
    # products are fixed on the flight segment and reused across its cells.
    cache = context.hcv_cache
    @. context.displacement = state.ξ.x - anchor
    Hd = _omrf_hmemo(cache, cache.state, H, context.displacement)
    @inbounds for j in eachindex(anchor)
        context.displacement[j] = _omrf_is_free(state, j) ? state.ξ.θ[j] : 0.0
    end
    Hv = _omrf_hmemo(cache, cache.velocity, H, context.displacement)
    likelihood_rate = 0.0
    likelihood_absolute = 0.0
    @inbounds for i in eachindex(anchor)
        active[i] || continue
        baseline = oracle.full_anchor[i] - oracle.prior_anchor[i]
        term = state.ξ.θ[i] * (baseline + Hd[i] + right * Hv[i])
        likelihood_rate += term
        likelihood_absolute += abs(state.ξ.θ[i]) *
            (abs(baseline) + abs(Hd[i]) + right * abs(Hv[i]))
    end
    n = length(anchor)
    gamma = n * eps(Float64) / (1 - n * eps(Float64))
    result = nextfloat(nextfloat(likelihood_rate +
        gamma * likelihood_absolute) + prior_upper)
    isfinite(result) || throw(OverflowError(
        "joint signed deterministic bound is not finite"))
    return result
end


function _omrf_observed_threshold_counts(X::Matrix{Int},
        threshold_starts::Vector{Int}, no_thresholds::Int)
    counts = zeros(Float64, no_thresholds)
    @inbounds for node in axes(X, 2), person in axes(X, 1)
        observed = X[person, node]
        if observed > 0
            counts[threshold_starts[node] + observed - 1] += 1.0
        end
    end
    return counts
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
            if _has_value(prior_spec, :global_scale)
                global_scale_index = only(_resolve_unc_spec(
                    _rget(prior_spec, :global_scale), unc_names,
                    "OMRF global scale"))
            end
            node_scale_indices = _resolve_unc_spec(
                _rget(prior_spec, :node_scale), unc_names,
                "OMRF node scales")
            length(node_scale_indices) == P || throw(DimensionMismatch(
                "resolved OMRF node-scale block has $(length(node_scale_indices)) " *
                "coordinates; expected $P"))
            isempty(intersect(node_scale_indices,
                [filter(>(0), [global_scale_index]); threshold_indices; interaction_indices])) ||
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
    # The HCV is undamped; the specification still records `hcv_damping = Inf`.
    use_hcv && _haskey(envelope_spec, :hcv_damping) &&
        !isinf(Float64(_rget(envelope_spec, :hcv_damping))) && throw(
            ArgumentError("only the undamped OMRF HCV (hcv_damping = Inf) is supported"))
    hcv_scope = use_hcv && _haskey(envelope_spec, :hcv_scope) ?
        String(_rget(envelope_spec, :hcv_scope)) : "global"
    hcv_scope == "global" || throw(ArgumentError(
        "only the global OMRF HCV scope is supported"))
    d = length(anchor)
    factorization = _haskey(envelope_spec, :factorization) ?
        String(_rget(envelope_spec, :factorization)) : "person"
    factorization == "person" || throw(ArgumentError(
        "only the person OMRF factorization is supported"))
    _haskey(envelope_spec, :zigzag_sign_aware) &&
        Bool(_rget(envelope_spec, :zigzag_sign_aware)) && throw(ArgumentError(
            "sign-aware Zig-Zag OMRF bounds are no longer supported"))
    X_float = Matrix{Float64}(X)
    observed_threshold_counts = _omrf_observed_threshold_counts(
        X, threshold_starts, no_thresholds)
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
        Matrix{Float64}(undef, P, N), Matrix{Float64}(undef, P, N),
        zeros(N, P), zeros(P, P),
        zeros(N, P), zeros(P, P), observed_threshold_counts,
        observed_interaction_sums, zeros(maximum(seen) - 1),
        zeros(maximum(seen) - 1), nothing, use_hcv,
        use_hcv && !hcv_after_warmup,
        zeros(d), zeros(d), Vector{Vector{Int}}(),
        sizehint!(Int[], d), falses(d), counts, trues(d), 0,
        OMRFHCVProductCache(d))
    context.node_supports = _omrf_node_supports(context)
    _cache_omrf_anchor!(context, anchor)
    if use_hcv
        context.likelihood_hessian = _omrf_likelihood_hessian(context, anchor)
    end
    return context
end

function _omrf_likelihood_hessian(context::OMRFResidualContext, anchor)
    hessian = zeros(length(anchor), length(anchor))
    _omrf_likelihood_hessian!(hessian, context)
    return hessian
end

"""Fill the aggregate anchor likelihood Hessian from cached probabilities."""
function _omrf_likelihood_hessian!(hessian::Matrix{Float64},
        context::OMRFResidualContext)
    fill!(hessian, 0.0)
    _omrf_hcv_cache_invalidate!(context.hcv_cache)
    N, P = size(context.X)
    @inbounds for n in 1:N, node in 1:P
        q = context.seen[node] - 1
        iszero(q) && continue
        start = context.threshold_starts[node]
        expected = context.anchor_expected_scores[node, n]
        expected2 = 0.0
        for u in 1:q
            probability = context.anchor_probabilities[start + u - 1, n]
            expected2 += u * u * probability
        end
        variance = expected2 - expected * expected
        for u in 1:q
            row = context.threshold_indices[start + u - 1]
            probability_u = context.anchor_probabilities[start + u - 1, n]
            for v in 1:q
                col = context.threshold_indices[start + v - 1]
                probability_v =
                    context.anchor_probabilities[start + v - 1, n]
                hessian[row, col] += probability_u *
                    ((u == v ? 1.0 : 0.0) - probability_v)
            end
            covariance_score = probability_u * (u - expected)
            for k in eachindex(context.incident_edges[node])
                edge_index = context.interaction_indices[
                    context.incident_edges[node][k]]
                value = covariance_score * context.X[n,
                    context.incident_neighbours[node][k]]
                hessian[row, edge_index] += value
                hessian[edge_index, row] += value
            end
        end
        for k in eachindex(context.incident_edges[node])
            row = context.interaction_indices[context.incident_edges[node][k]]
            xk = context.X[n, context.incident_neighbours[node][k]]
            for l in eachindex(context.incident_edges[node])
                col = context.interaction_indices[
                    context.incident_edges[node][l]]
                xl = context.X[n, context.incident_neighbours[node][l]]
                hessian[row, col] += variance * xk * xl
            end
        end
    end
    context.counts.native_full_data_scans += 1
    return hessian
end


@inline function _omrf_node_probabilities_with_reference!(context::OMRFResidualContext,
        x, n::Int, node::Int)
    q = context.seen[node] - 1
    iszero(q) && return (0.0, 1.0)
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
    reference_weight = exp(-max_eta)
    denominator = reference_weight
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
    return expected, reference_weight / denominator
end

@inline function _omrf_node_probabilities!(context::OMRFResidualContext,
        x, n::Int, node::Int)
    expected, _ = _omrf_node_probabilities_with_reference!(
        context, x, n, node)
    return expected
end

function _cache_omrf_anchor!(context::OMRFResidualContext, anchor)
    N, P = size(context.X)
    @inbounds for n in 1:N, node in 1:P
        expected, reference_probability =
            _omrf_node_probabilities_with_reference!(context, anchor, n, node)
        context.anchor_expected_scores[node, n] = expected
        context.anchor_reference_probabilities[node, n] = reference_probability
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
    # Weight of the (undamped) Hessian correction.
    α = context.hcv_active ? 1.0 : 0.0
    # A person touches the support of every node; record it once.
    if !isempty(subset)
        @inbounds for node in 1:P
            _record_omrf_touched_support!(context, node)
        end
    end
    @inbounds for n in subset
        1 <= n <= N || throw(BoundsError(context.X, (n, :)))
        for node in 1:P
            _omrf_node_residual!(context, out, x, n, node, α)
        end
    end
    context.counts.analytic_persons += length(subset)
    context.counts.analytic_node_conditionals += P * length(subset)
    context.counts.analytic_residual += 1
    context.counts.analytic_residual_seconds += (time_ns() - started) * 1e-9
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
    lazy_residual_subset_bound::Float64
    cached_deterministic_signed_rate::Float64
    exact_selected_refinement::Bool
    cached_residual::Vector{Float64}
    cached_residual_valid::Bool
    active_free::BitVector
    active_generation::Int
    interval_xlo::Vector{Float64}
    interval_xhi::Vector{Float64}
    interval_vlo::Vector{Float64}
    interval_vhi::Vector{Float64}
    interval_glo::Vector{Float64}
    interval_ghi::Vector{Float64}
end

@inline _omrf_base_oracle(oracle) = oracle
@inline _omrf_base_oracle(oracle::PDMPSamplers.WithStats) =
    _omrf_base_oracle(oracle.f)
@inline _omrf_base_oracle(oracle::PDMPSamplers.WithResidualStats) =
    _omrf_base_oracle(oracle.f)

function OMRFAnalyticSubsamplingOracle(ctx::StanSubsamplingContext,
        envelope_spec, unc_names, anchor::Vector{Float64})
    d = length(anchor)
    residual_context = _omrf_residual_context(
        envelope_spec, unc_names, anchor, ctx.counts)
    analytic_prior = _omrf_analytic_prior(envelope_spec, residual_context)
    # R lists arrive through JuliaCall with Symbol keys, whereas research
    # fixtures commonly use String keys.  Resolve both forms consistently;
    # `get(..., "name", default)` silently disabled the requested refinement
    # for the ordinary R-origin path.
    exact_selected_refinement = _haskey(
        envelope_spec, :exact_selected_refinement) && Bool(_rget(
            envelope_spec, :exact_selected_refinement))
    oracle = OMRFAnalyticSubsamplingOracle(ctx, copy(anchor), zeros(d), zeros(d),
        zeros(d), zeros(d), false, residual_context, analytic_prior,
        NaN, NaN, exact_selected_refinement, zeros(d), false, trues(d), 0,
        zeros(d), zeros(d), zeros(d), zeros(d), zeros(d), zeros(d))
    _clear_gradient!(ctx, :full, oracle.full_anchor, anchor)
    if analytic_prior === nothing
        _clear_gradient!(ctx, :prior, oracle.prior_anchor, anchor)
    else
        _omrf_prior_gradient!(oracle.prior_anchor, analytic_prior, anchor)
    end
    return oracle
end

function PDMPSamplers.set_active_set!(
        oracle::OMRFAnalyticSubsamplingOracle, free::BitVector)
    length(free) == length(oracle.active_free) || throw(DimensionMismatch(
        "active OMRF set has the wrong dimension"))
    changed = false
    @inbounds for index in eachindex(free, oracle.active_free)
        changed |= free[index] != oracle.active_free[index]
    end
    changed || return nothing
    copyto!(oracle.active_free, free)
    oracle.active_generation += 1
    context = oracle.residual_context
    copyto!(context.active_free, free)
    context.active_generation = oracle.active_generation
    return nothing
end

function PDMPSamplers.deterministic_rate_cell_bound(
        oracle::OMRFAnalyticSubsamplingOracle,
        state::PDMPSamplers.AbstractPDMPState,
        flow::PDMPSamplers.ContinuousDynamics, left::Real, right::Real)
    oracle.analytic_prior === nothing && return nothing
    active = oracle.analytic_prior isa OMRFNodeScaleAnalyticPrior ?
        oracle.active_free : nothing
    return _omrf_direct_deterministic_cell_bound(
        oracle, active, state, flow, left, right)
end
PDMPSamplers.has_direct_deterministic_rate_cell_bound(
    oracle::OMRFAnalyticSubsamplingOracle) =
    oracle.analytic_prior !== nothing

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
        _omrf_haffine!(context.hcv_product, context.hcv_cache,
            something(context.likelihood_hessian), context.displacement)
        @. out += context.hcv_product
    end
    return out
end

function (oracle::OMRFAnalyticSubsamplingOracle)(out, x, subset, anchor)
    context = oracle.residual_context
    context.counts.direct_selected_factor_exact_evaluations += length(subset)
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
        if counts.last_mark_residual
            counts.residual_source_acceptances += 1
        else
            counts.deterministic_source_acceptances += 1
        end
    end
    return nothing
end

function PDMPSamplers.record_subsampling_mark_source!(
        oracle::OMRFAnalyticSubsamplingOracle, residual_source::Bool)
    counts = oracle.residual_context.counts
    counts.last_mark_residual = residual_source
    if residual_source
        counts.residual_source_proposals += 1
    else
        counts.deterministic_source_proposals += 1
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
Certified directional categorical HCV remainder bound.

The local Taylor bound is `V * Delta^2 / 8`.  Independently, exponential
tilting across a log-weight range `Delta` changes a categorical law by at
most `tanh(Delta/4)` in total variation, while Hoeffding's covariance
inequality gives `|z'C(a)delta| <= V * Delta / 4`.  Taking the smaller bound
preserves the certificate while preventing a quadratic blow-up far from the
anchor.
"""
@inline function _omrf_categorical_hcv_bound(
        velocity_range::Real, displacement_range::Real)
    V = max(0.0, Float64(velocity_range))
    Delta = max(0.0, Float64(displacement_range))
    local_bound = Delta * Delta / 8
    # Unique positive solution of
    #   Delta^2/8 = tanh(Delta/4) + Delta/4.
    # Below it the quadratic branch is tighter, so the transcendental call is
    # unnecessary.  Testing against the coarser `1 + Delta/4` is invalid: it
    # switches too late and makes the resulting minimum non-monotone, which
    # can make a closed-cell certificate smaller than an interior point.
    Delta <= 3.592978281031538 && return V * local_bound
    # Exponential tilting by a log-weight whose range is Delta changes a
    # categorical probability vector by at most tanh(Delta/4) in total
    # variation.  Combining this with the covariance term V*Delta/4 tightens
    # the far-anchor enclosure.  `tanh` is evaluated only when the global
    # branch beats the inexpensive local quadratic branch.
    return V * (tanh(Delta / 4) + Delta / 4)
end

@inline function _omrf_range_with_reference(values)
    lower = 0.0
    upper = 0.0
    @inbounds for value in values
        lower = min(lower, value)
        upper = max(upper, value)
    end
    return upper - lower
end

"""Direct categorical directional Taylor remainder and its certified bound."""
function _omrf_categorical_hcv_remainder(anchor_logits, delta, z)
    q = length(anchor_logits)
    length(delta) == q == length(z) || throw(DimensionMismatch(
        "categorical HCV vectors must have equal lengths"))
    maximum_logit = 0.0
    displaced_maximum = 0.0
    @inbounds for u in 1:q
        maximum_logit = max(maximum_logit, anchor_logits[u])
        displaced_maximum = max(displaced_maximum,
            anchor_logits[u] + delta[u])
    end
    anchor_denominator = exp(-maximum_logit)
    displaced_denominator = exp(-displaced_maximum)
    anchor_mean_delta = 0.0
    anchor_mean_z = 0.0
    displaced_mean_z = 0.0
    @inbounds for u in 1:q
        anchor_weight = exp(anchor_logits[u] - maximum_logit)
        displaced_weight = exp(
            anchor_logits[u] + delta[u] - displaced_maximum)
        anchor_denominator += anchor_weight
        displaced_denominator += displaced_weight
        anchor_mean_delta += anchor_weight * delta[u]
        anchor_mean_z += anchor_weight * z[u]
        displaced_mean_z += displaced_weight * z[u]
    end
    anchor_mean_delta /= anchor_denominator
    anchor_mean_z /= anchor_denominator
    displaced_mean_z /= displaced_denominator
    # Include the reference category explicitly: its z and delta values are
    # zero, but its centered-product contribution is generally nonzero.
    covariance = exp(-maximum_logit) / anchor_denominator *
        anchor_mean_z * anchor_mean_delta
    @inbounds for u in 1:q
        probability = exp(anchor_logits[u] - maximum_logit) /
            anchor_denominator
        covariance += probability *
            (z[u] - anchor_mean_z) * (delta[u] - anchor_mean_delta)
    end
    remainder = displaced_mean_z - anchor_mean_z - covariance
    bound = _omrf_categorical_hcv_bound(
        _omrf_range_with_reference(z), _omrf_range_with_reference(delta))
    return remainder, bound
end

function _omrf_direct_categorical_cell_bound(
        anchor_logits, delta_left, z, width::Real)
    width >= 0 || throw(ArgumentError("cell width must be nonnegative"))
    left, _ = _omrf_categorical_hcv_remainder(
        anchor_logits, delta_left, z)
    right_delta = similar(delta_left)
    @inbounds for index in eachindex(right_delta, delta_left, z)
        right_delta[index] = delta_left[index] + width * z[index]
    end
    right, _ = _omrf_categorical_hcv_remainder(
        anchor_logits, right_delta, z)
    range_z = _omrf_range_with_reference(z)
    delta_range = max(_omrf_range_with_reference(delta_left),
        _omrf_range_with_reference(right_delta))
    return min(max(abs(left), abs(right)) + range_z^3 * width^2 / 32,
        _omrf_categorical_hcv_bound(range_z, delta_range))
end


"""
Research-only direct person-node HCV certificate for BPS, Zig-Zag and
Boomerang cells.

Along a linear path the exact directional remainder `h(t)` satisfies
`h''(t) = D^3 logsumexp(eta(t))[z,z,z]`. With the reference-category value
zero included, this is a categorical third central moment. Hoeffding's
covariance inequality gives `abs(h''(t)) <= R(z)^3/4`; the interpolation
error on `[l,r]` is at most `sup(abs(h''))*(r-l)^2/8`. Therefore
`max(abs(h(l)),abs(h(r))) + R(z)^3*(r-l)^2/32` is certified on the cell.

Each cache column stores cumulative bounds for all person-node factors in one
active grid cell. It is both the factor-bound store (by differencing adjacent
entries) and the allocation-free distinguished-factor sampling table.

For the ordinary first-order CV mode, `h(t)=z'(p(t)-p(a))` and
`h'(t)=z'(diag(p(t))-p(t)p(t)')z >= 0`. Hence `h` is nondecreasing and
`max(abs(h(t)))` on a cell is attained at an endpoint. The shared fill uses
`nextfloat(max(abs(h(left)), abs(h(right))))`; its curvature coefficient is
exactly zero. This mode stores no likelihood Hessian.

For Zig-Zag, the cached cell mass is the sum of coordinatewise categorical
HCV bounds. The velocity range is bounded by the sum of absolute threshold
velocities plus the category range times the absolute incident-interaction
field. Because the categorical displacement range is convex along a linear
cell and the capped bound is monotone, its maximum is attained at an endpoint.

For Boomerang, each factor stores its harmonic offset/cosine/sine triples.
Pairwise analytic harmonic maxima give the displacement and velocity ranges
over a closed cell, after which the same capped categorical HCV inequality is
used. The cumulative factor table is built once for the cell and reused by
every proposal drawn from that cell.
"""

mutable struct OMRFSharedCurrentWorkspace
    representatives::Matrix{Int}
    person_group::Matrix{Int}
    multiplicities::Matrix{Int}
    group_count::Vector{Int}
    hash_keys::Vector{UInt}
    hash_values::Vector{Int}
    hash_epochs::Vector{Int}
    epoch::Int
    previous_free::BitVector
    previous_frozen_nonzero::BitVector
    dirty_nodes::BitVector
    active_incident_slots::Matrix{Int}
    active_incident_count::Vector{Int}
    initialized::Bool
    use_fallback::Bool
    field0::Matrix{Float64}
    field0_upper::Matrix{Float64}
    fieldv::Matrix{Float64}
    fieldv_upper::Matrix{Float64}
    left_moment::Matrix{Float64}
    right_moment::Matrix{Float64}
    left_error::Matrix{Float64}
    right_error::Matrix{Float64}
    last_boundary::Float64
    last_valid::Bool
    rebuild_count::Int
    rebuild_seconds::Float64
    pattern_evaluations::Int
end

function _new_omrf_shared_current_workspace(context::OMRFResidualContext)
    N, P = size(context.X)
    slots = nextpow(2, max(4, 2N))
    E = length(context.interaction_indices)
    return OMRFSharedCurrentWorkspace(zeros(Int, N, P), zeros(Int, P, N),
        zeros(Int, N, P), zeros(Int, P), zeros(UInt, slots),
        zeros(Int, slots), zeros(Int, slots), 0, falses(length(context.active_free)),
        falses(E), falses(P), zeros(Int, P - 1, P), zeros(Int, P),
        false, false, zeros(N, P), zeros(N, P),
        zeros(N, P), zeros(N, P), zeros(N, P), zeros(N, P),
        zeros(N, P), zeros(N, P),
        0.0, false, 0, 0.0, 0)
end

mutable struct OMRFDirectPersonNodeHCVEnvelope <:
        PDMPSamplers.AbstractResidualEnvelope
    n_persons::Int
    n_nodes::Int
    anchor::Vector{Float64}
    context::OMRFResidualContext
    person_signed_roof::Bool
    signed_only_person_roof::Bool
    component_scales!::Nothing
    component_cell_scales!::Nothing
    totals::Vector{Float64}
    scales::Vector{Float64}
    cell_scales::Vector{Float64}
    cumulative_masses::Vector{Float64}
    max_cells::Int
    cell_left::Vector{Float64}
    cell_right::Vector{Float64}
    cell_mass::Vector{Float64}
    cell_cumulative::Matrix{Float64}
    factor_workspace::Vector{Float64}
    person_cell_cumulative::Matrix{Float64}
    person_workspace::Vector{Float64}
    person_left_signed::Vector{Float64}
    person_right_signed::Vector{Float64}
    person_left_abs::Vector{Float64}
    person_right_abs::Vector{Float64}
    person_curvature::Vector{Float64}
    joint_group_size::Int
    joint_groups::Matrix{Int}
    joint_group_offsets::Vector{Int}
    joint_group_lengths::Vector{Int}
    joint_cell_caps::Matrix{Float64}
    joint_cell_cumulative::Matrix{Float64}
    joint_cell_mass::Vector{Float64}
    shared_current::Union{Nothing,OMRFSharedCurrentWorkspace}
    cell_count::Int
    active_cell::Int
    active_time::Float64
    path_kind::UInt8
    origin_time::Float64
    path_x::Vector{Float64}
    path_theta::Vector{Float64}
    path_anchor::Vector{Float64}
    hcv_path_offset::Vector{Float64}
    hcv_path_cosine::Vector{Float64}
    hcv_path_sine::Vector{Float64}
    path_active_generation::Int
    path_flow_token::UInt
    qmax::Int
    path_threshold_delta0::Matrix{Float64}
    path_threshold_z::Matrix{Float64}
    path_threshold_sine::Matrix{Float64}
    interaction_delta::Matrix{Float64}
    interaction_velocity::Matrix{Float64}
    interaction_abs_velocity::Matrix{Float64}
    field_delta::Matrix{Float64}
    field_velocity::Matrix{Float64}
    field_abs_velocity::Matrix{Float64}
    path_delta0::Matrix{Float64}
    path_z::Matrix{Float64}
    path_harmonic_sine::Matrix{Float64}
    path_anchor_reference_probability::Vector{Float64}
    path_anchor_mean_z::Vector{Float64}
    path_anchor_cov_z_delta0::Vector{Float64}
    path_anchor_var_z::Vector{Float64}
    path_curvature_coefficient::Vector{Float64}
    path_velocity_range::Vector{Float64}
    path_zigzag_velocity_sum::Vector{Float64}
    last_boundary_remainders::Vector{Float64}
    last_boundary::Float64
    last_boundary_valid::Bool
    path_valid::Bool
end

function _new_omrf_direct_hcv_envelope(anchor,
        context::OMRFResidualContext; max_cells::Int=128,
        person_signed_roof::Bool=false,
        signed_only_person_roof::Bool=false,
        joint_group_size::Int=0,
        shared_current::Bool=false)
    max_cells >= 1 || throw(ArgumentError(
        "direct HCV cache must hold at least one cell"))
    N, P = size(context.X)
    n_factors = N * P
    qmax = maximum(context.seen) - 1
    signed_only_person_roof && !person_signed_roof && throw(ArgumentError(
        "signed-only person roofs require per-cell signed HCV"))
    joint_group_size in (0, 4, 16, 32, 64) || throw(ArgumentError(
        "research joint signed route supports groups of 4, 16, 32 or 64 only"))
    joint_group_size > 0 &&
        !signed_only_person_roof &&
        throw(ArgumentError(
            "joint groups require a signed-only person roof"))
    joint_group_size == 4 && N % 4 != 0 && throw(ArgumentError(
        "fixed groups of four require a divisible population"))
    joint_group_size >= 16 &&
        !(shared_current && context.use_hcv) && throw(ArgumentError(
            "grouped joint signed bounds require shared-current evaluation and the HCV; got shared=$shared_current use_hcv=$(context.use_hcv)"))
    joint_group_size >= 16 && N < 2joint_group_size && throw(ArgumentError(
        "grouped joint signed bounds need at least two full groups (N >= 2m)"))
    shared_current && joint_group_size == 0 && throw(ArgumentError(
        "shared-current categorical evaluation requires joint groups"))
    shared_current && any(!=(4), context.seen) && throw(ArgumentError(
        "shared-current categorical evaluation currently requires four categories at every node"))
    shared_current && !context.hcv_active && throw(ArgumentError(
        "shared-current JS4 requires the active HCV"))
    n_groups = joint_group_size == 0 ? 0 : cld(N, joint_group_size)
    permutation = joint_group_size == 0 ? Int[] :
        randperm(Random.Xoshiro(230923), N)
    groups = zeros(Int, joint_group_size, n_groups)
    group_offsets = Vector{Int}(undef, n_groups)
    group_lengths = Vector{Int}(undef, n_groups)
    @inbounds for b in 1:n_groups
        offset = (b - 1) * joint_group_size + 1
        length_b = min(joint_group_size, N - offset + 1)
        group_offsets[b] = offset
        group_lengths[b] = length_b
        for slot in 1:length_b
            groups[slot, b] = permutation[offset + slot - 1]
        end
    end
    return OMRFDirectPersonNodeHCVEnvelope(
        N, P, anchor, context, person_signed_roof, signed_only_person_roof,
        nothing, nothing,
        ones(1), zeros(1), zeros(1),
        zeros(1), max_cells, zeros(max_cells), zeros(max_cells),
        zeros(max_cells), zeros(n_factors, max_cells), zeros(n_factors),
        zeros(N, max_cells), zeros(N), zeros(N), zeros(N), zeros(N),
        zeros(N), zeros(N),
        joint_group_size, groups, group_offsets, group_lengths,
        zeros(n_groups, max_cells),
        zeros(n_groups, max_cells), zeros(max_cells),
        shared_current ? _new_omrf_shared_current_workspace(context) : nothing,
        0, 0, 0.0, UInt8(0), 0.0,
        zeros(length(anchor)), zeros(length(anchor)),
        zeros(length(anchor)), zeros(length(anchor)),
        zeros(length(anchor)), zeros(length(anchor)),
        -1, UInt(0), qmax,
        zeros(qmax, P), zeros(qmax, P), zeros(qmax, P),
        zeros(P, P), zeros(P, P), zeros(P, P),
        zeros(N, P), zeros(N, P), zeros(N, P),
        zeros(qmax, n_factors), zeros(qmax, n_factors),
        zeros(qmax, n_factors), zeros(n_factors), zeros(n_factors),
        zeros(n_factors), zeros(n_factors),
        zeros(n_factors), zeros(n_factors), zeros(n_factors),
        zeros(n_factors), 0.0, false, false)
end

PDMPSamplers.n_observations(envelope::OMRFDirectPersonNodeHCVEnvelope) =
    envelope.n_persons
PDMPSamplers.trajectory_scale_anchor(
    envelope::OMRFDirectPersonNodeHCVEnvelope) = envelope.anchor

@inline function _omrf_shared_equal(workspace::OMRFSharedCurrentWorkspace,
        context, node::Int, person::Int,
        representative::Int, free, x)
    neighbours = context.incident_neighbours[node]
    @inbounds for slot in 1:workspace.active_incident_count[node]
        k = workspace.active_incident_slots[slot, node]
        context.X[person, neighbours[k]] ==
            context.X[representative, neighbours[k]] || return false
    end
    return true
end

function _omrf_shared_rebuild_node!(workspace::OMRFSharedCurrentWorkspace,
        context, node::Int, free, x)
    workspace.epoch = workspace.epoch == typemax(Int) ? 1 :
        workspace.epoch + 1
    workspace.epoch == 1 && fill!(workspace.hash_epochs, 0)
    epoch = workspace.epoch
    edges = context.incident_edges[node]
    neighbours = context.incident_neighbours[node]
    N = size(context.X, 1)
    slots = length(workspace.hash_keys)
    groups = 0
    active_count = 0
    @inbounds for k in eachindex(edges)
        index = context.interaction_indices[edges[k]]
        if free[index] || !iszero(x[index])
            active_count += 1
            workspace.active_incident_slots[active_count, node] = k
        end
    end
    workspace.active_incident_count[node] = active_count
    @inbounds for person in 1:N
        key = UInt(0x7b8a9cdef0123456)
        for slot in 1:active_count
            k = workspace.active_incident_slots[slot, node]
            key = hash(context.X[person, neighbours[k]], key)
        end
        slot = Int(key & UInt(slots - 1)) + 1
        group = 0
        while workspace.hash_epochs[slot] == epoch
            if workspace.hash_keys[slot] == key
                candidate = workspace.hash_values[slot]
                representative = workspace.representatives[candidate, node]
                if _omrf_shared_equal(workspace, context, node, person,
                        representative, free, x)
                    group = candidate
                    break
                end
            end
            slot = slot == slots ? 1 : slot + 1
        end
        if group == 0
            groups += 1
            group = groups
            workspace.representatives[group, node] = person
            workspace.multiplicities[group, node] = 0
            workspace.hash_epochs[slot] = epoch
            workspace.hash_keys[slot] = key
            workspace.hash_values[slot] = group
        end
        workspace.multiplicities[group, node] += 1
        workspace.person_group[node, person] = group
    end
    workspace.group_count[node] = groups
    return nothing
end

function _omrf_shared_sync_patterns!(workspace::OMRFSharedCurrentWorkspace,
        context, x)
    started = time_ns()
    free = context.active_free
    fill!(workspace.dirty_nodes, !workspace.initialized)
    @inbounds for node in eachindex(context.incident_edges)
        for edge in context.incident_edges[node]
            index = context.interaction_indices[edge]
            frozen_nonzero = !free[index] && !iszero(x[index])
            if workspace.previous_free[index] != free[index] ||
                    workspace.previous_frozen_nonzero[edge] != frozen_nonzero
                workspace.dirty_nodes[node] = true
            end
        end
    end
    @inbounds for node in eachindex(workspace.dirty_nodes)
        if workspace.dirty_nodes[node]
            _omrf_shared_rebuild_node!(workspace, context, node, free, x)
            workspace.rebuild_count += 1
            context.counts.shared_current_rebuilt_nodes += 1
        end
    end
    copyto!(workspace.previous_free, free)
    @inbounds for edge in eachindex(context.interaction_indices)
        index = context.interaction_indices[edge]
        workspace.previous_frozen_nonzero[edge] =
            !free[index] && !iszero(x[index])
    end
    workspace.initialized = true
    duration = (time_ns() - started) * 1e-9
    workspace.rebuild_seconds += duration
    context.counts.shared_current_rebuild_seconds += duration
    return nothing
end

function _omrf_shared_prepare_fields!(workspace::OMRFSharedCurrentWorkspace,
        context, x, v)
    @inbounds for node in eachindex(workspace.group_count)
        edges = context.incident_edges[node]
        neighbours = context.incident_neighbours[node]
        for group in 1:workspace.group_count[node]
            person = workspace.representatives[group, node]
            f0lo = 0.0; f0hi = 0.0
            fvlo = 0.0; fvhi = 0.0
            for slot in 1:workspace.active_incident_count[node]
                k = workspace.active_incident_slots[slot, node]
                index = context.interaction_indices[edges[k]]
                observed = context.X[person, neighbours[k]]
                # Exact zero products contribute no interval width for finite
                # four-category observations. The active list still contains
                # free-zero edges, so later velocity/position changes retain
                # the correct fixed current-pattern membership.
                if !iszero(x[index])
                    tl, th = _omrf_iscale(x[index], x[index], observed)
                    f0lo, f0hi = _omrf_iadd(f0lo, f0hi, tl, th)
                end
                if !iszero(v[index])
                    tl, th = _omrf_iscale(v[index], v[index], observed)
                    fvlo, fvhi = _omrf_iadd(fvlo, fvhi, tl, th)
                end
            end
            workspace.field0[group, node] = f0lo
            workspace.field0_upper[group, node] = f0hi
            workspace.fieldv[group, node] = fvlo
            workspace.fieldv_upper[group, node] = fvhi
        end
    end
    workspace.last_valid = false
    return nothing
end

@inline function _omrf_shared_logit_interval(context, node::Int, u::Int,
        f0lo::Float64, f0hi::Float64, fvlo::Float64, fvhi::Float64,
        x, v, t::Float64)
    index = context.threshold_indices[context.threshold_starts[node] + u - 1]
    tl, th = _omrf_iscale(v[index], v[index], t)
    xl, xh = _omrf_iadd(x[index], x[index], tl, th)
    tl, th = _omrf_iscale(fvlo, fvhi, t)
    fl, fh = _omrf_iadd(f0lo, f0hi, tl, th)
    tl, th = _omrf_iscale(fl, fh, u)
    ll, lh = _omrf_iadd(xl, xh, tl, th)
    tl, th = _omrf_iscale(fvlo, fvhi, u)
    zl, zh = _omrf_iadd(v[index], v[index], tl, th)
    return ll, lh, zl, zh
end

@inline function _omrf_shared_moment(context, node::Int,
        f0lo::Float64, f0hi::Float64, fvlo::Float64, fvhi::Float64,
        x, v, t::Float64)
    q = context.seen[node] - 1
    f0 = f0lo + 0.5 * (f0hi - f0lo)
    fv = fvlo + 0.5 * (fvhi - fvlo)
    maximum_logit = 0.0
    max_logit_error = 0.0
    max_z_error = 0.0
    max_abs_z = 0.0
    @inbounds for u in 1:q
        index = context.threshold_indices[context.threshold_starts[node] + u - 1]
        logit = x[index] + t * v[index] + u * (f0 + t * fv)
        z = v[index] + u * fv
        ll, lh, zl, zh = _omrf_shared_logit_interval(context, node, u,
            f0lo, f0hi, fvlo, fvhi, x, v, t)
        max_logit_error = max(max_logit_error,
            abs(logit - ll), abs(lh - logit))
        max_z_error = max(max_z_error, abs(z - zl), abs(zh - z))
        max_abs_z = max(max_abs_z, abs(z))
        maximum_logit = max(maximum_logit, logit)
    end
    denominator = exp(-maximum_logit)
    numerator = 0.0
    @inbounds for u in 1:q
        index = context.threshold_indices[context.threshold_starts[node] + u - 1]
        logit = x[index] + t * v[index] + u * (f0 + t * fv)
        z = v[index] + u * fv
        weight = exp(logit - maximum_logit)
        denominator += weight
        numerator += weight * z
    end
    moment = numerator / denominator
    # If all logits move by at most delta, every softmax probability ratio
    # moves by at most exp(2delta). Thus the expectation changes by at most
    # expm1(2delta)*max(abs(z)); derivative uncertainty adds max_z_error.
    # The interval above includes all field, path and logit rounding. The
    # gamma_32 term encloses the four-category exp/sum/divide evaluation;
    # this is operation-count based, not calibrated to saved discrepancies.
    sensitivity = _omrf_up(expm1(_omrf_up(2max_logit_error)))
    isfinite(sensitivity) || throw(OverflowError(
        "shared-current softmax sensitivity overflow"))
    gamma32 = 32eps(Float64) / (1 - 32eps(Float64))
    error = _omrf_up(max_z_error +
        _omrf_up(sensitivity * max_abs_z) +
        gamma32 * (1 + max_abs_z))
    return moment, error
end

function _omrf_shared_fill_signed_cell!(envelope, workspace, state,
        left::Float64, right::Float64, cell::Int)
    context = envelope.context
    context.counts.shared_current_cell_evaluations += 1
    adjacent = workspace.last_valid && workspace.last_boundary == left
    @inbounds for node in 1:envelope.n_nodes
        for group in 1:workspace.group_count[node]
            f0 = workspace.field0[group, node]
            f0hi = workspace.field0_upper[group, node]
            fv = workspace.fieldv[group, node]
            fvhi = workspace.fieldv_upper[group, node]
            if adjacent
                workspace.left_moment[group, node] =
                    workspace.right_moment[group, node]
                workspace.left_error[group, node] =
                    workspace.right_error[group, node]
            else
                moment, error = _omrf_shared_moment(context, node, f0,
                    f0hi, fv, fvhi, state.ξ.x, state.ξ.θ, left)
                workspace.left_moment[group, node] = moment
                workspace.left_error[group, node] = error
                workspace.pattern_evaluations += 1
                context.counts.shared_current_pattern_evaluations += 1
            end
            moment, error = _omrf_shared_moment(context, node, f0,
                f0hi, fv, fvhi, state.ξ.x, state.ξ.θ, right)
            workspace.right_moment[group, node] = moment
            workspace.right_error[group, node] = error
            workspace.pattern_evaluations += 1
            context.counts.shared_current_pattern_evaluations += 1
        end
    end
    fill!(envelope.person_left_signed, 0.0)
    fill!(envelope.person_right_signed, 0.0)
    fill!(envelope.person_left_abs, 0.0)
    fill!(envelope.person_right_abs, 0.0)
    fill!(envelope.person_curvature, 0.0)
    factor = 1
    @inbounds for node in 1:envelope.n_nodes
        for person in 1:envelope.n_persons
            group = workspace.person_group[node, person]
            left_error = workspace.left_error[group, node]
            right_error = workspace.right_error[group, node]
            left_h = workspace.left_moment[group, node] -
                envelope.path_anchor_mean_z[factor] -
                envelope.path_anchor_cov_z_delta0[factor] -
                left * envelope.path_anchor_var_z[factor]
            right_h = workspace.right_moment[group, node] -
                envelope.path_anchor_mean_z[factor] -
                envelope.path_anchor_cov_z_delta0[factor] -
                right * envelope.path_anchor_var_z[factor]
            gamma8 = 8eps(Float64) / (1 - 8eps(Float64))
            left_error = _omrf_up(left_error + gamma8 *
                (1 + abs(workspace.left_moment[group, node]) +
                 abs(envelope.path_anchor_mean_z[factor]) +
                 abs(envelope.path_anchor_cov_z_delta0[factor]) +
                 abs(left * envelope.path_anchor_var_z[factor])))
            right_error = _omrf_up(right_error + gamma8 *
                (1 + abs(workspace.right_moment[group, node]) +
                 abs(envelope.path_anchor_mean_z[factor]) +
                 abs(envelope.path_anchor_cov_z_delta0[factor]) +
                 abs(right * envelope.path_anchor_var_z[factor])))
            envelope.last_boundary_remainders[factor] = right_h
            envelope.person_left_signed[person] += left_h + left_error
            envelope.person_right_signed[person] += right_h + right_error
            envelope.person_left_abs[person] += abs(left_h) + left_error
            envelope.person_right_abs[person] += abs(right_h) + right_error
            envelope.person_curvature[person] +=
                envelope.path_curvature_coefficient[factor]
            factor += 1
        end
    end
    envelope.cell_left[cell] = left
    envelope.cell_right[cell] = right
    width2 = nextfloat((right - left)^2)
    rounding = envelope.n_nodes * eps(Float64) /
        (1 - envelope.n_nodes * eps(Float64))
    cumulative = 0.0
    @inbounds for person in 1:envelope.n_persons
        # The signed sums above include a one-sided error allowance for JS4.
        # This auxiliary person mass remains a valid *absolute* fallback by
        # using the sum of absolute factor endpoints, including that error.
        endpoint_bound = max(
            nextfloat(envelope.person_left_abs[person] * (1 + rounding)),
            nextfloat(envelope.person_right_abs[person] * (1 + rounding)))
        curvature_bound = nextfloat(envelope.person_curvature[person] *
            (1 + rounding))
        bound = nextfloat(endpoint_bound +
            nextfloat(curvature_bound * width2))
        envelope.person_workspace[person] = bound
        cumulative = nextfloat(cumulative + bound)
        envelope.person_cell_cumulative[person, cell] = cumulative
    end
    envelope.cell_mass[cell] = cumulative
    envelope.last_boundary = right
    envelope.last_boundary_valid = true
    workspace.last_boundary = right
    workspace.last_valid = true
    counts = context.counts
    counts.direct_factor_bound_scans += 1
    counts.direct_factors_evaluated += factor - 1
    counts.direct_residual_table_constructions += 1
    counts.direct_left_endpoint_evaluations += adjacent ? 0 : factor - 1
    counts.direct_right_endpoint_evaluations += factor - 1
    counts.direct_adjacent_endpoint_reuses += adjacent ? factor - 1 : 0
    return cumulative
end

# The direct envelope exposes its certified Boomerang cell callback by Julia
# dispatch rather than through the legacy callback field.  Tell the generic
# grid validator about that explicit implementation; this does not weaken the
# validation for any other envelope.
PDMPSamplers.validate_subsampling_grid_envelope(
    ::PDMPSamplers.AnyBoomerang,
    ::OMRFDirectPersonNodeHCVEnvelope) = nothing

function _omrf_direct_prepare_linear_fields!(
        envelope::OMRFDirectPersonNodeHCVEnvelope, state, flow)
    context = envelope.context
    actual_flow = _omrf_actual_flow(flow)
    zigzag = actual_flow isa ZigZag
    harmonic = actual_flow isa PDMPSamplers.AnyBoomerang
    envelope.path_kind = harmonic ? UInt8(2) : zigzag ? UInt8(1) : UInt8(0)
    fill!(envelope.interaction_delta, 0.0)
    fill!(envelope.interaction_velocity, 0.0)
    (zigzag || harmonic) && fill!(envelope.interaction_abs_velocity, 0.0)
    edge = 1
    @inbounds for left in 1:(envelope.n_nodes - 1)
        for right in (left + 1):envelope.n_nodes
        index = context.interaction_indices[edge]
        delta, velocity, sine = if harmonic
            _omrf_harmonic_coefficients(
                state, actual_flow, envelope.anchor, index)
        else
            d, v = _omrf_coordinate_at(
                state, actual_flow, envelope.anchor, index, 0.0)
            (d, v, 0.0)
        end
            envelope.interaction_delta[left, right] = delta
            envelope.interaction_delta[right, left] = delta
            envelope.interaction_velocity[left, right] = velocity
            envelope.interaction_velocity[right, left] = velocity
            if zigzag || harmonic
                third = harmonic ? sine : abs(velocity)
                envelope.interaction_abs_velocity[left, right] = third
                envelope.interaction_abs_velocity[right, left] = third
            end
            edge += 1
        end
    end
    mul!(envelope.field_delta, context.X_float,
        envelope.interaction_delta)
    mul!(envelope.field_velocity, context.X_float,
        envelope.interaction_velocity)
    (zigzag || harmonic) && mul!(envelope.field_abs_velocity,
        context.X_float, envelope.interaction_abs_velocity)
    return nothing
end


@inline function _omrf_direct_hcv_remainder_at(
        envelope::OMRFDirectPersonNodeHCVEnvelope, state, flow,
        person::Int, node::Int, t::Float64)
    context = envelope.context
    q = context.seen[node] - 1
    iszero(q) && return 0.0, 0.0
    field_delta = 0.0
    field_velocity = 0.0
    @inbounds for k in eachindex(context.incident_edges[node])
        covariate = context.X[person, context.incident_neighbours[node][k]]
        iszero(covariate) && continue
        index = context.interaction_indices[context.incident_edges[node][k]]
        delta, velocity = _omrf_coordinate_at(
            state, flow, envelope.anchor, index, t)
        field_delta += covariate * delta
        field_velocity += covariate * velocity
    end
    start = context.threshold_starts[node]
    mean_z = 0.0
    mean_delta = 0.0
    mean_z_delta = 0.0
    probability_sum = 0.0
    minimum_delta = 0.0
    maximum_delta = 0.0
    minimum_z = 0.0
    maximum_z = 0.0
    @inbounds for u in 1:q
        index = context.threshold_indices[start + u - 1]
        threshold_delta, threshold_velocity = _omrf_coordinate_at(
            state, flow, envelope.anchor, index, t)
        delta = threshold_delta + u * field_delta
        z = threshold_velocity + u * field_velocity
        probability = context.anchor_probabilities[start + u - 1, person]
        probability_sum += probability
        mean_z += probability * z
        mean_delta += probability * delta
        mean_z_delta += probability * z * delta
        minimum_delta = min(minimum_delta, delta)
        maximum_delta = max(maximum_delta, delta)
        minimum_z = min(minimum_z, z)
        maximum_z = max(maximum_z, z)
    end
    reference_probability = context.anchor_reference_probabilities[node, person]
    denominator = reference_probability * exp(-maximum_delta)
    displaced_mean_z = 0.0
    @inbounds for u in 1:q
        index = context.threshold_indices[start + u - 1]
        threshold_delta, threshold_velocity = _omrf_coordinate_at(
            state, flow, envelope.anchor, index, t)
        delta = threshold_delta + u * field_delta
        z = threshold_velocity + u * field_velocity
        weight = context.anchor_probabilities[start + u - 1, person] *
            exp(delta - maximum_delta)
        denominator += weight
        displaced_mean_z += weight * z
    end
    displaced_mean_z /= denominator
    covariance = mean_z_delta - mean_z * mean_delta
    return displaced_mean_z - mean_z - covariance,
        maximum_z - minimum_z, maximum_delta - minimum_delta
end

@inline function _omrf_direct_hcv_factor_cell_bound(
        envelope::OMRFDirectPersonNodeHCVEnvelope, state, flow,
        person::Int, node::Int, left::Float64, right::Float64)
    left_h, range_z, left_delta_range = _omrf_direct_hcv_remainder_at(
        envelope, state, flow, person, node, left)
    right_h, _, right_delta_range = _omrf_direct_hcv_remainder_at(
        envelope, state, flow, person, node, right)
    width = right - left
    global_cap = range_z * (1 + max(left_delta_range,
        right_delta_range) / 4)
    return nextfloat(min(max(abs(left_h), abs(right_h)) +
        range_z^3 * width^2 / 32, global_cap))
end

function _omrf_direct_path_matches(
        envelope::OMRFDirectPersonNodeHCVEnvelope, state,
        flow_token::UInt=envelope.path_flow_token)
    envelope.path_valid || return false
    state.t[] == envelope.origin_time || return false
    envelope.path_active_generation == envelope.context.active_generation ||
        return false
    envelope.path_flow_token == flow_token || return false
    @inbounds for index in eachindex(envelope.path_x)
        envelope.path_x[index] == state.ξ.x[index] || return false
        envelope.path_theta[index] == state.ξ.θ[index] || return false
        envelope.path_anchor[index] == envelope.anchor[index] || return false
    end
    return true
end

function _omrf_direct_prepare_path!(
        envelope::OMRFDirectPersonNodeHCVEnvelope, state, flow,
        flow_token::UInt=_omrf_direct_flow_token(flow))
    context = envelope.context
    counts = context.counts
    anchor_changed = false
    state_changed = state.t[] != envelope.origin_time
    @inbounds for index in eachindex(envelope.path_x)
        anchor_changed |= envelope.path_anchor[index] != envelope.anchor[index]
        state_changed |= envelope.path_x[index] != state.ξ.x[index] ||
            envelope.path_theta[index] != state.ξ.θ[index]
    end
    if anchor_changed
        counts.direct_path_anchor_changes += 1
    elseif envelope.path_active_generation != context.active_generation
        counts.direct_path_mask_changes += 1
    elseif envelope.path_flow_token != flow_token
        counts.direct_path_flow_changes += 1
    elseif state_changed
        counts.direct_path_state_changes += 1
    else
        counts.direct_path_forced_rebuilds += 1
    end
    copyto!(envelope.path_x, state.ξ.x)
    copyto!(envelope.path_theta, state.ξ.θ)
    copyto!(envelope.path_anchor, envelope.anchor)
    envelope.origin_time = state.t[]
    envelope.path_active_generation = envelope.context.active_generation
    envelope.path_flow_token = flow_token
    N = envelope.n_persons
    actual_flow = _omrf_actual_flow(flow)
    harmonic = actual_flow isa PDMPSamplers.AnyBoomerang
    if context.hcv_active
        H = something(context.likelihood_hessian)
        if harmonic
            @inbounds for index in eachindex(envelope.anchor)
                if _omrf_is_free(state, index)
                    context.displacement[index] =
                        actual_flow.μ[index] - envelope.anchor[index]
                    envelope.hcv_path_cosine[index] =
                        state.ξ.x[index] - actual_flow.μ[index]
                    envelope.hcv_path_sine[index] = state.ξ.θ[index]
                else
                    context.displacement[index] =
                        state.ξ.x[index] - envelope.anchor[index]
                    envelope.hcv_path_cosine[index] = 0.0
                    envelope.hcv_path_sine[index] = 0.0
                end
            end
            cache = context.hcv_cache
            copyto!(envelope.hcv_path_offset,
                _omrf_hmemo(cache, cache.reference, H, context.displacement))
            # x - mu and v stay in a fixed two-dimensional span along an orbit.
            _omrf_haffine!(context.hcv_product, cache, H,
                envelope.hcv_path_cosine, cache.harmonic)
            copyto!(envelope.hcv_path_cosine, context.hcv_product)
            _omrf_haffine!(context.hcv_product, cache, H,
                envelope.hcv_path_sine, cache.harmonic)
            copyto!(envelope.hcv_path_sine, context.hcv_product)
        else
            @inbounds for index in eachindex(envelope.anchor)
                context.displacement[index] =
                    state.ξ.x[index] - envelope.anchor[index]
                envelope.hcv_path_sine[index] = _omrf_is_free(state, index) ?
                    state.ξ.θ[index] : 0.0
            end
            cache = context.hcv_cache
            _omrf_path_offset!(envelope.hcv_path_offset, cache, H,
                context.displacement)
            copyto!(envelope.hcv_path_sine,
                _omrf_hmemo(cache, cache.velocity, H, envelope.hcv_path_sine))
            fill!(envelope.hcv_path_cosine, 0.0)
        end
    end
    # Cache linear signed/absolute fields or harmonic offset/cosine/sine
    # fields with matrix products once per physical path instead of repeating
    # an incident-edge traversal for every factor.
    _omrf_direct_prepare_linear_fields!(envelope, state, flow)
    destination = 1
    @inbounds for node in 1:envelope.n_nodes
        q = context.seen[node] - 1
        start = context.threshold_starts[node]
        for u in 1:q
            index = context.threshold_indices[start + u - 1]
            threshold_delta0, threshold_velocity, threshold_sine = if harmonic
                _omrf_harmonic_coefficients(
                    state, actual_flow, envelope.anchor, index)
            else
                d, v = _omrf_coordinate_at(
                    state, actual_flow, envelope.anchor, index, 0.0)
                (d, v, 0.0)
            end
            envelope.path_threshold_delta0[u, node] = threshold_delta0
            envelope.path_threshold_z[u, node] = threshold_velocity
            envelope.path_threshold_sine[u, node] = threshold_sine
        end
        for person in 1:N
            field_delta0 = envelope.field_delta[person, node]
            field_velocity = envelope.field_velocity[person, node]
            field_sine = harmonic ?
                envelope.field_abs_velocity[person, node] : 0.0
            probability_sum = 0.0
            mean_z = 0.0
            mean_delta0 = 0.0
            mean_z_delta0 = 0.0
            mean_z2 = 0.0
            zigzag_velocity_sum = 0.0
            minimum_z = 0.0
            maximum_z = 0.0
            for u in 1:q
                threshold_delta0 = envelope.path_threshold_delta0[u, node]
                threshold_velocity = envelope.path_threshold_z[u, node]
                threshold_sine = envelope.path_threshold_sine[u, node]
                delta0 = threshold_delta0 + u * field_delta0
                z = threshold_velocity + u * field_velocity
                sine = threshold_sine + u * field_sine
                envelope.path_delta0[u, destination] = delta0
                envelope.path_z[u, destination] = z
                envelope.path_harmonic_sine[u, destination] = sine
                probability = context.anchor_probabilities[
                    start + u - 1, person]
                probability_sum += probability
                mean_z += probability * z
                mean_delta0 += probability * delta0
                mean_z_delta0 += probability * z * delta0
                mean_z2 += probability * z * z
                zigzag_velocity_sum += abs(threshold_velocity)
                minimum_z = min(minimum_z, z)
                maximum_z = max(maximum_z, z)
            end
            for u in (q + 1):envelope.qmax
                envelope.path_delta0[u, destination] = 0.0
                envelope.path_z[u, destination] = 0.0
                envelope.path_harmonic_sine[u, destination] = 0.0
            end
            # The reference category has delta=z=0; its probability therefore
            # does not contribute to these raw moments.
            envelope.path_anchor_reference_probability[destination] =
                context.anchor_reference_probabilities[node, person]
            envelope.path_anchor_mean_z[destination] = mean_z
            envelope.path_anchor_cov_z_delta0[destination] =
                mean_z_delta0 - mean_z * mean_delta0
            envelope.path_anchor_var_z[destination] = mean_z2 - mean_z^2
            range_z = maximum_z - minimum_z
            envelope.path_curvature_coefficient[destination] =
                range_z^3 / 32
            envelope.path_velocity_range[destination] = range_z
            envelope.path_zigzag_velocity_sum[destination] =
                zigzag_velocity_sum +
                q * envelope.field_abs_velocity[person, node]
            destination += 1
        end
    end
    if envelope.shared_current !== nothing
        workspace = envelope.shared_current::OMRFSharedCurrentWorkspace
        _omrf_shared_sync_patterns!(workspace, context, state.ξ.x)
        # The dense initial D2 fixture showed no useful arithmetic reuse.
        # Freeze this conservative cutoff before any route comparison.
        total_patterns = 0
        @inbounds for count in workspace.group_count
            total_patterns += count
        end
        workspace.use_fallback = total_patterns >= 4000
        workspace.last_valid = false
        if !workspace.use_fallback
            _omrf_shared_prepare_fields!(workspace, context, state.ξ.x,
                state.ξ.θ)
        end
    end
    envelope.cell_count = 0
    envelope.active_cell = 0
    envelope.last_boundary_valid = false
    envelope.path_valid = true
    counts = context.counts
    counts.direct_path_preparations += 1
    counts.direct_path_factors_prepared += N * envelope.n_nodes
    return envelope
end

function PDMPSamplers.deterministic_gradient!(out,
        cv::PDMPSamplers.SubsampledControlVariate{F,O,E}, state, flow) where
        {F,O,E<:OMRFDirectPersonNodeHCVEnvelope}
    provider = cv.deterministic_gradient!
    if provider isa PDMPSamplers.WithStats
        PDMPSamplers._inc_counter_∇f_calls(provider.stats)
        PDMPSamplers._inc_counter_deterministic_gradient_calls(provider.stats)
    end
    # A dependent slab (node-scale prior) wraps the oracle in a slab correction;
    # the fast path below then adds the same correction, so that the
    # precomputed control-variate path is used for every slab.
    slab = _omrf_slab_wrapper(provider)
    oracle = slab === nothing ? _omrf_base_oracle(provider) :
        _omrf_base_oracle(slab.oracle)
    oracle isa OMRFAnalyticSubsamplingOracle ||
        return cv.deterministic_gradient!(out, state.ξ.x)
    context = oracle.residual_context
    envelope = cv.envelope
    if !context.hcv_active || !envelope.path_valid ||
            envelope.path_active_generation != context.active_generation
        return cv.deterministic_gradient!(out, state.ξ.x)
    end
    _omrf_deterministic_gradient!(out, oracle,
        something(oracle.analytic_prior), state.ξ.x)
    _omrf_add_path_terms!(out, envelope, state, flow)
    if slab !== nothing
        x = state.ξ.x
        PDMPSamplers.active_slab_potential_gradient!(slab.slab, slab.all_buf, x, slab.all_active)
        PDMPSamplers.active_slab_potential_gradient!(slab.slab, slab.active_buf, x, slab.active)
        @. out = out - slab.all_buf + slab.active_buf
    end
    _OMRF_CHECK_FAST_DETERMINISTIC[] && _omrf_check_fast_deterministic(out, cv, state)
    return out
end

# Verification switch: compare the fast deterministic gradient with the
# generic (full Hessian-product) evaluation at every call.
const _OMRF_CHECK_FAST_DETERMINISTIC = Ref(
    get(ENV, "PDMPSAMPLERSR_CHECK_FAST_DETERMINISTIC", "0") == "1")
const _OMRF_FAST_DETERMINISTIC_MAXERR = Ref(0.0)
function _omrf_check_fast_deterministic(out, cv, state)
    reference = similar(out)
    cv.deterministic_gradient!(reference, state.ξ.x)
    err = maximum(abs.(reference .- out)) / max(1.0, maximum(abs.(reference)))
    _OMRF_FAST_DETERMINISTIC_MAXERR[] = max(_OMRF_FAST_DETERMINISTIC_MAXERR[], err)
    return nothing
end

# The control-variate Hessian term along the current path, from the vectors
# precomputed by the envelope at the start of the cell.
function _omrf_add_path_terms!(out, envelope, state, flow)
    elapsed = Float64(state.t[] - envelope.origin_time)
    actual = _omrf_actual_flow(flow)
    if actual isa PDMPSamplers.AnyBoomerang
        sine, cosine = sincos(elapsed)
        @inbounds @simd for index in eachindex(out)
            out[index] += envelope.hcv_path_offset[index] +
                cosine * envelope.hcv_path_cosine[index] +
                sine * envelope.hcv_path_sine[index]
        end
    else
        @inbounds @simd for index in eachindex(out)
            out[index] += envelope.hcv_path_offset[index] +
                elapsed * envelope.hcv_path_sine[index]
        end
    end
    return out
end

function PDMPSamplers.subsampling_deterministic_signed_rate(
        oracle::OMRFAnalyticSubsamplingOracle,
        cv::PDMPSamplers.SubsampledControlVariate{F,O,E}, state, flow) where
        {F,O,E<:OMRFDirectPersonNodeHCVEnvelope}
    actual = _omrf_actual_flow(flow)
    actual isa Union{BouncyParticle,PDMPSamplers.AnyBoomerang} ||
        return PDMPSamplers.DeferredSubsamplingRate(0.0, false)
    context = oracle.residual_context
    context.hcv_active && oracle.exact_selected_refinement ||
        return PDMPSamplers.DeferredSubsamplingRate(0.0, false)
    scratch = oracle.interval_glo
    PDMPSamplers.deterministic_gradient!(scratch, cv, state, flow)
    if actual isa PDMPSamplers.AnyBoomerang
        actual.Γ isa Diagonal ||
            return PDMPSamplers.DeferredSubsamplingRate(0.0, false)
        @inbounds @simd for index in eachindex(scratch)
            scratch[index] -= actual.Γ.diag[index] *
                (state.ξ.x[index] - actual.μ[index])
        end
    end
    signed = dot(state.ξ.θ, scratch)
    oracle.cached_deterministic_signed_rate = signed
    return PDMPSamplers.DeferredSubsamplingRate(signed, true)
end

function PDMPSamplers.subsampling_deferred_candidate_rate(
        oracle::OMRFAnalyticSubsamplingOracle, state, residual, scale,
        flow::Union{BouncyParticle,PDMPSamplers.AnyBoomerang,
            PDMPSamplers.PreconditionedDynamics}, signed_rate, subset)
    return max(0.0, signed_rate + scale * dot(state.ξ.θ, residual))
end

@inline function _omrf_direct_cached_zigzag_bound(
        envelope::OMRFDirectPersonNodeHCVEnvelope,
        factor::Int, node::Int, left::Float64, right::Float64)
    velocity_sum = envelope.path_zigzag_velocity_sum[factor]
    left_delta = _omrf_direct_cached_delta_range(
        envelope, factor, node, left)
    right_delta = _omrf_direct_cached_delta_range(
        envelope, factor, node, right)
    return nextfloat(max(
        _omrf_categorical_hcv_bound(velocity_sum, left_delta),
        _omrf_categorical_hcv_bound(velocity_sum, right_delta)))
end

@inline function _omrf_direct_harmonic_abs_cell_bound(
        offset::Float64, cosine::Float64, sine::Float64,
        midpoint_sine::Float64, midpoint_cosine::Float64,
        half_width::Float64)
    amplitude = sqrt(muladd(cosine, cosine, sine * sine))
    midpoint = abs(offset + cosine * midpoint_cosine +
        sine * midpoint_sine)
    return nextfloat(min(abs(offset) + amplitude,
        midpoint + half_width * amplitude))
end

@inline function _omrf_direct_cached_harmonic_bound(
        envelope::OMRFDirectPersonNodeHCVEnvelope,
        factor::Int, node::Int, midpoint_sine::Float64,
        midpoint_cosine::Float64, half_width::Float64)
    q = envelope.context.seen[node] - 1
    iszero(q) && return 0.0
    displacement_range = 0.0
    velocity_range = 0.0
    # Include the reference category as the identically zero coefficient
    # triple.  Pairwise harmonic extrema give the exact range envelope over
    # the complete closed cell; the categorical cap then certifies the HCV
    # directional remainder for this factor.
    @inbounds for u in 1:q, v in 0:(u - 1)
        offset = envelope.path_delta0[u, factor]
        cosine = envelope.path_z[u, factor]
        sine = envelope.path_harmonic_sine[u, factor]
        if v > 0
            offset -= envelope.path_delta0[v, factor]
            cosine -= envelope.path_z[v, factor]
            sine -= envelope.path_harmonic_sine[v, factor]
        end
        # These coefficients are ordinary finite posterior coordinates.  The
        # generic `hypot` helper pays for overflow/underflow rescaling that is
        # unnecessary in this hot factor loop.  The fused square norm is the
        # same harmonic amplitude for the model's numerical range.
        displacement_bound = _omrf_direct_harmonic_abs_cell_bound(
            offset, cosine, sine, midpoint_sine, midpoint_cosine, half_width)
        velocity_bound = _omrf_direct_harmonic_abs_cell_bound(
            0.0, sine, -cosine, midpoint_sine, midpoint_cosine, half_width)
        displacement_range = max(displacement_range, displacement_bound)
        velocity_range = max(velocity_range, velocity_bound)
    end
    return nextfloat(_omrf_categorical_hcv_bound(
        velocity_range, displacement_range))
end

@inline function _omrf_direct_cached_remainder_at(
        envelope::OMRFDirectPersonNodeHCVEnvelope,
        factor::Int, person::Int, node::Int, t::Float64)
    context = envelope.context
    q = context.seen[node] - 1
    iszero(q) && return 0.0
    start = context.threshold_starts[node]
    maximum_delta = 0.0
    @inbounds for u in 1:q
        delta = envelope.path_delta0[u, factor] +
            t * envelope.path_z[u, factor]
        maximum_delta = max(maximum_delta, delta)
    end
    denominator = 0.0
    displaced_mean_z = 0.0
    @inbounds for u in 1:q
        probability = context.anchor_probabilities[start + u - 1, person]
        delta = envelope.path_delta0[u, factor] +
            t * envelope.path_z[u, factor]
        weight = probability * exp(delta - maximum_delta)
        denominator += weight
        displaced_mean_z += weight * envelope.path_z[u, factor]
    end
    denominator += envelope.path_anchor_reference_probability[factor] *
        exp(-maximum_delta)
    displaced_mean_z /= denominator
    residual = displaced_mean_z - envelope.path_anchor_mean_z[factor]
    residual -= envelope.path_anchor_cov_z_delta0[factor] +
        t * envelope.path_anchor_var_z[factor]
    return residual
end

@inline function _omrf_direct_cached_delta_range(
        envelope::OMRFDirectPersonNodeHCVEnvelope,
        factor::Int, node::Int, t::Float64)
    q = envelope.context.seen[node] - 1
    minimum_delta = 0.0
    maximum_delta = 0.0
    @inbounds for u in 1:q
        delta = envelope.path_delta0[u, factor] +
            t * envelope.path_z[u, factor]
        minimum_delta = min(minimum_delta, delta)
        maximum_delta = max(maximum_delta, delta)
    end
    return maximum_delta - minimum_delta
end

function _omrf_direct_reset_path!(
        envelope::OMRFDirectPersonNodeHCVEnvelope, state, flow,
        flow_token::UInt=_omrf_direct_flow_token(flow))
    started = time_ns()
    _omrf_direct_prepare_path!(envelope, state, flow, flow_token)
    envelope.context.counts.direct_path_preparation_seconds +=
        (time_ns() - started) * 1e-9
    return envelope
end

@inline function _omrf_direct_find_cell(
        envelope::OMRFDirectPersonNodeHCVEnvelope,
        left::Float64, right::Float64)
    @inbounds for cell in 1:envelope.cell_count
        envelope.cell_left[cell] == left &&
            envelope.cell_right[cell] == right && return cell
    end
    return 0
end


function _omrf_direct_finish_person_cell!(
        envelope::OMRFDirectPersonNodeHCVEnvelope,
        cell::Int, signed_available::Bool)
    fill!(envelope.person_workspace, 0.0)
    N = envelope.n_persons
    @inbounds for node in 1:envelope.n_nodes, person in 1:N
        envelope.person_workspace[person] +=
            envelope.factor_workspace[(node - 1) * N + person]
    end
    cumulative = 0.0
    rounding = envelope.n_nodes * eps(Float64) /
        (1 - envelope.n_nodes * eps(Float64))
    @inbounds for person in 1:N
        factor_bound = envelope.person_workspace[person]
        if signed_available && envelope.person_signed_roof
            # Each factor coefficient already contains K_nj/8. Summing the
            # signed endpoints before taking absolute values preserves node
            # cancellation, while the curvature sum remains conservative.
            width2 = nextfloat((envelope.cell_right[cell] -
                envelope.cell_left[cell])^2)
            endpoint_bound = max(
                nextfloat(abs(envelope.person_left_signed[person]) +
                    rounding * envelope.person_left_abs[person]),
                nextfloat(abs(envelope.person_right_signed[person]) +
                    rounding * envelope.person_right_abs[person]))
            curvature_bound = nextfloat(envelope.person_curvature[person] *
                (1 + rounding))
            signed_bound = nextfloat(endpoint_bound +
                nextfloat(curvature_bound * width2))
            factor_bound = min(factor_bound, nextfloat(signed_bound))
        end
        cumulative += factor_bound
        envelope.person_cell_cumulative[person, cell] = cumulative
    end
    envelope.cell_mass[cell] = cumulative
    return cumulative
end

"""Harmonic (Boomerang) cell of the cached direct envelope.

Computes the same factor bounds as `_omrf_direct_cached_harmonic_bound`, but
pair by pair over all persons of a node, so that the inner loop runs over
persons and vectorises. Every factor sees the same category pairs in the same
order, so the running maxima, and hence all outputs, are bit-identical to the
factor-by-factor loop. `person_left_signed` and `person_right_signed` hold the
per-person running maxima: they are used only by linear (BPS) signed roofs.
"""
function _omrf_direct_hcv_fill_harmonic_cell!(
        envelope::OMRFDirectPersonNodeHCVEnvelope,
        left::Float64, right::Float64, cell::Int, adjacent::Bool)
    midpoint_sine, midpoint_cosine = sincos((left + right) / 2)
    half_width = (right - left) / 2
    delta0 = envelope.path_delta0
    velocity = envelope.path_z
    sine_coefficient = envelope.path_harmonic_sine
    factor_bound = envelope.factor_workspace
    displacement_range = envelope.person_left_signed
    velocity_range = envelope.person_right_signed
    n_persons = envelope.n_persons
    @inbounds for node in 1:envelope.n_nodes
        q = envelope.context.seen[node] - 1
        first = (node - 1) * n_persons
        if iszero(q)
            for person in 1:n_persons
                factor_bound[first + person] = 0.0
            end
            continue
        end
        fill!(displacement_range, 0.0)
        fill!(velocity_range, 0.0)
        # The reference category is the identically zero coefficient triple.
        for u in 1:q, v in 0:(u - 1)
            @simd ivdep for person in 1:n_persons
                factor = first + person
                offset = delta0[u, factor]
                cosine = velocity[u, factor]
                sine = sine_coefficient[u, factor]
                if v > 0
                    offset -= delta0[v, factor]
                    cosine -= velocity[v, factor]
                    sine -= sine_coefficient[v, factor]
                end
                displacement_bound = _omrf_direct_harmonic_abs_cell_bound(
                    offset, cosine, sine, midpoint_sine, midpoint_cosine,
                    half_width)
                velocity_bound = _omrf_direct_harmonic_abs_cell_bound(
                    0.0, sine, -cosine, midpoint_sine, midpoint_cosine,
                    half_width)
                displacement_range[person] = max(displacement_range[person],
                    displacement_bound)
                velocity_range[person] = max(velocity_range[person],
                    velocity_bound)
            end
        end
        for person in 1:n_persons
            factor_bound[first + person] = nextfloat(_omrf_categorical_hcv_bound(
                velocity_range[person], displacement_range[person]))
        end
    end
    cumulative = 0.0
    @inbounds for factor in eachindex(factor_bound)
        cumulative += factor_bound[factor]
        envelope.cell_cumulative[factor, cell] = cumulative
    end
    envelope.cell_left[cell] = left
    envelope.cell_right[cell] = right
    envelope.cell_mass[cell] = cumulative
    cumulative = _omrf_direct_finish_person_cell!(envelope, cell, false)
    envelope.last_boundary = right
    envelope.last_boundary_valid = true
    n_factors = envelope.n_persons * envelope.n_nodes
    counts = envelope.context.counts
    counts.direct_factor_bound_scans += 1
    counts.direct_factors_evaluated += n_factors
    counts.direct_residual_table_constructions += 1
    counts.direct_left_endpoint_evaluations += adjacent ? 0 : n_factors
    counts.direct_right_endpoint_evaluations += n_factors
    counts.direct_adjacent_endpoint_reuses += adjacent ? n_factors : 0
    return cumulative
end

function _omrf_direct_hcv_fill_cell_cached!(
        envelope::OMRFDirectPersonNodeHCVEnvelope,
        left::Float64, right::Float64, cell::Int,
        reuse_adjacent::Bool, zigzag::Bool=false, harmonic::Bool=false)
    adjacent = reuse_adjacent && envelope.last_boundary_valid &&
        envelope.last_boundary == left
    harmonic && return _omrf_direct_hcv_fill_harmonic_cell!(
        envelope, left, right, cell, adjacent)
    cumulative = 0.0
    factor = 1
    width2 = (right - left)^2
    midpoint_sine, midpoint_cosine = harmonic ?
        sincos((left + right) / 2) : (0.0, 1.0)
    half_width = harmonic ? (right - left) / 2 : 0.0
    signed_linear = envelope.person_signed_roof && !zigzag && !harmonic
    if signed_linear
        fill!(envelope.person_left_signed, 0.0)
        fill!(envelope.person_right_signed, 0.0)
        fill!(envelope.person_left_abs, 0.0)
        fill!(envelope.person_right_abs, 0.0)
        fill!(envelope.person_curvature, 0.0)
    end
    @inbounds for node in 1:envelope.n_nodes
        for person in 1:envelope.n_persons
            if harmonic
                bound = _omrf_direct_cached_harmonic_bound(
                    envelope, factor, node, midpoint_sine,
                    midpoint_cosine, half_width)
            elseif zigzag
                bound = _omrf_direct_cached_zigzag_bound(
                    envelope, factor, node, left, right)
            else
                if adjacent
                    left_h = envelope.last_boundary_remainders[factor]
                else
                    left_h = _omrf_direct_cached_remainder_at(
                        envelope, factor, person, node, left)
                end
                right_h = _omrf_direct_cached_remainder_at(
                    envelope, factor, person, node, right)
                envelope.last_boundary_remainders[factor] = right_h
                if signed_linear
                    envelope.person_left_signed[person] += left_h
                    envelope.person_right_signed[person] += right_h
                    envelope.person_left_abs[person] += abs(left_h)
                    envelope.person_right_abs[person] += abs(right_h)
                    envelope.person_curvature[person] +=
                        envelope.path_curvature_coefficient[factor]
                end
                interpolation_bound = max(abs(left_h), abs(right_h)) +
                    envelope.path_curvature_coefficient[factor] * width2
                delta_range = max(
                    _omrf_direct_cached_delta_range(
                        envelope, factor, node, left),
                    _omrf_direct_cached_delta_range(
                        envelope, factor, node, right))
                velocity_range = envelope.path_velocity_range[factor]
                global_cap = _omrf_categorical_hcv_bound(
                    velocity_range, delta_range)
                bound = nextfloat(min(interpolation_bound, global_cap))
            end
            envelope.factor_workspace[factor] = bound
            cumulative += bound
            envelope.cell_cumulative[factor, cell] = cumulative
            factor += 1
        end
    end
    envelope.cell_left[cell] = left
    envelope.cell_right[cell] = right
    envelope.cell_mass[cell] = cumulative
    cumulative = _omrf_direct_finish_person_cell!(
        envelope, cell, signed_linear)
    envelope.last_boundary = right
    envelope.last_boundary_valid = true
    counts = envelope.context.counts
    counts.direct_factor_bound_scans += 1
    counts.direct_factors_evaluated += envelope.n_persons * envelope.n_nodes
    counts.direct_residual_table_constructions += 1
    counts.direct_left_endpoint_evaluations += adjacent ? 0 : factor - 1
    counts.direct_right_endpoint_evaluations += factor - 1
    counts.direct_adjacent_endpoint_reuses += adjacent ? factor - 1 : 0
    return cumulative
end

"""Certified person endpoints without constructing unused factor fallback masses.

The signed Taylor enclosure is identical to the S1a branch of
`_omrf_direct_finish_person_cell!`. Only BPS uses this research policy.
Selected-person point refinement still evaluates its own valid factor caps.
"""
function _omrf_direct_hcv_fill_signed_only_cell!(
        envelope::OMRFDirectPersonNodeHCVEnvelope,
        left::Float64, right::Float64, cell::Int)
    adjacent = envelope.last_boundary_valid &&
        envelope.last_boundary == left
    fill!(envelope.person_left_signed, 0.0)
    fill!(envelope.person_right_signed, 0.0)
    fill!(envelope.person_left_abs, 0.0)
    fill!(envelope.person_right_abs, 0.0)
    fill!(envelope.person_curvature, 0.0)
    factor = 1
    @inbounds for node in 1:envelope.n_nodes
        for person in 1:envelope.n_persons
            left_h = adjacent ? envelope.last_boundary_remainders[factor] :
                _omrf_direct_cached_remainder_at(
                    envelope, factor, person, node, left)
            right_h = _omrf_direct_cached_remainder_at(
                envelope, factor, person, node, right)
            envelope.last_boundary_remainders[factor] = right_h
            envelope.person_left_signed[person] += left_h
            envelope.person_right_signed[person] += right_h
            envelope.person_left_abs[person] += abs(left_h)
            envelope.person_right_abs[person] += abs(right_h)
            envelope.person_curvature[person] +=
                envelope.path_curvature_coefficient[factor]
            factor += 1
        end
    end
    envelope.cell_left[cell] = left
    envelope.cell_right[cell] = right
    width2 = nextfloat((right - left)^2)
    rounding = envelope.n_nodes * eps(Float64) /
        (1 - envelope.n_nodes * eps(Float64))
    cumulative = 0.0
    @inbounds for person in 1:envelope.n_persons
        endpoint_bound = max(
            nextfloat(abs(envelope.person_left_signed[person]) +
                rounding * envelope.person_left_abs[person]),
            nextfloat(abs(envelope.person_right_signed[person]) +
                rounding * envelope.person_right_abs[person]))
        curvature_bound = nextfloat(envelope.person_curvature[person] *
            (1 + rounding))
        signed_bound = nextfloat(endpoint_bound +
            nextfloat(curvature_bound * width2))
        envelope.person_workspace[person] = signed_bound
        cumulative = nextfloat(cumulative + signed_bound)
        envelope.person_cell_cumulative[person, cell] = cumulative
    end
    envelope.cell_mass[cell] = cumulative
    envelope.last_boundary = right
    envelope.last_boundary_valid = true
    counts = envelope.context.counts
    counts.direct_factor_bound_scans += 1
    counts.direct_factors_evaluated += factor - 1
    counts.direct_residual_table_constructions += 1
    counts.direct_left_endpoint_evaluations += adjacent ? 0 : factor - 1
    counts.direct_right_endpoint_evaluations += factor - 1
    counts.direct_adjacent_endpoint_reuses += adjacent ? factor - 1 : 0
    return cumulative
end

function _omrf_direct_hcv_fill_cell!(
        envelope::OMRFDirectPersonNodeHCVEnvelope, state, flow,
        left::Float64, right::Float64, cell::Int)
    if envelope.signed_only_person_roof
        _omrf_actual_flow(flow) isa BouncyParticle || throw(ArgumentError(
            "signed-only person roofs currently support linear BPS"))
        if envelope.shared_current !== nothing
            workspace = envelope.shared_current::OMRFSharedCurrentWorkspace
            if !workspace.use_fallback
                return _omrf_shared_fill_signed_cell!(envelope, workspace,
                    state, left, right, cell)
            end
            envelope.context.counts.shared_current_dense_fallback_cells += 1
        end
        return _omrf_direct_hcv_fill_signed_only_cell!(
            envelope, left, right, cell)
    end
    return _omrf_direct_hcv_fill_cell_cached!(
        envelope, left, right, cell, true,
        _omrf_actual_flow(flow) isa ZigZag,
        _omrf_actual_flow(flow) isa PDMPSamplers.AnyBoomerang)
end

@inline _omrf_direct_flow_token(::BouncyParticle) = UInt(0x01)
@inline _omrf_direct_flow_token(::ZigZag) = UInt(0x02)
@inline _omrf_direct_flow_token(::PDMPSamplers.AnyBoomerang) = UInt(0x03)
@inline _omrf_direct_flow_token(flow) = UInt(hash(typeof(flow)))
@inline function _omrf_direct_flow_token(
        flow::PDMPSamplers.PreconditionedDynamics)
    generation = hasproperty(flow.metric, :generation) ?
        flow.metric.generation : UInt(0)
    return xor(UInt(objectid(flow.metric)), generation)
end

function _omrf_direct_component_cell_scales!(out,
        envelope::OMRFDirectPersonNodeHCVEnvelope, state,
        actual_flow::Union{BouncyParticle,ZigZag,PDMPSamplers.AnyBoomerang},
        left, right, flow_token::UInt)
    left = Float64(left)
    right = Float64(right)
    right >= left || throw(ArgumentError("direct HCV cell is reversed"))
    _omrf_direct_path_matches(envelope, state, flow_token) ||
        _omrf_direct_reset_path!(
            envelope, state, actual_flow, flow_token)
    envelope.context.counts.direct_internal_cell_uses += 1
    cell = _omrf_direct_find_cell(envelope, left, right)
    if iszero(cell)
        envelope.cell_count < envelope.max_cells || throw(ArgumentError(
            "direct HCV active grid exceeds the preallocated cell cache"))
        cell = envelope.cell_count + 1
        started = time_ns()
        _omrf_direct_hcv_fill_cell!(
            envelope, state, actual_flow, left, right, cell)
        envelope.context.counts.direct_cell_certificate_seconds +=
            (time_ns() - started) * 1e-9
        envelope.cell_count = cell
    else
        envelope.context.counts.direct_cached_cell_hits += 1
    end
    envelope.active_cell = cell
    envelope.scales[1] = envelope.cell_mass[cell]
    envelope.cell_scales[1] = envelope.cell_mass[cell]
    envelope.cumulative_masses[1] = envelope.cell_mass[cell]
    out[1] = envelope.cell_mass[cell]
    return out
end

function PDMPSamplers.component_cell_scales!(out,
        envelope::OMRFDirectPersonNodeHCVEnvelope, state,
        actual_flow::BouncyParticle, left, right)
    return _omrf_direct_component_cell_scales!(out, envelope, state,
        actual_flow, left, right, _omrf_direct_flow_token(actual_flow))
end

function PDMPSamplers.component_cell_scales!(out,
        envelope::OMRFDirectPersonNodeHCVEnvelope, state,
        actual_flow::ZigZag, left, right)
    return _omrf_direct_component_cell_scales!(out, envelope, state,
        actual_flow, left, right, _omrf_direct_flow_token(actual_flow))
end

function PDMPSamplers.component_cell_scales!(out,
        envelope::OMRFDirectPersonNodeHCVEnvelope, state,
        actual_flow::PDMPSamplers.AnyBoomerang, left, right)
    return _omrf_direct_component_cell_scales!(out, envelope, state,
        actual_flow, left, right, _omrf_direct_flow_token(actual_flow))
end

function PDMPSamplers.component_cell_scales!(out,
        envelope::OMRFDirectPersonNodeHCVEnvelope, state,
        flow::PDMPSamplers.PreconditionedDynamics, left, right)
    flow.dynamics isa Union{BouncyParticle,ZigZag,
            PDMPSamplers.AnyBoomerang} || throw(ArgumentError(
        "direct person-node HCV cells require BPS, Zig-Zag, or Boomerang"))
    return _omrf_direct_component_cell_scales!(out, envelope, state,
        flow.dynamics, left, right, _omrf_direct_flow_token(flow))
end

@inline function _omrf_direct_activate_time!(
        envelope::OMRFDirectPersonNodeHCVEnvelope, state, t)::Float64
    # Aggregate screening receives the live cell-origin state and a future
    # relative time. Selected-subset preparation receives an already moved
    # candidate and `t == 0`. Distinguish these without scanning factors.
    local_time = _omrf_direct_path_matches(envelope, state) ? Float64(t) :
        Float64(state.t[] - envelope.origin_time + t)
    tolerance = 16eps(max(1.0, abs(local_time)))
    @inbounds for cell in 1:envelope.cell_count
        if envelope.cell_left[cell] - tolerance <= local_time <=
                envelope.cell_right[cell] + tolerance
            envelope.active_cell = cell
            envelope.active_time = local_time
            mass = envelope.cell_mass[cell]
            envelope.scales[1] = mass
            envelope.cumulative_masses[1] = mass
            envelope.context.counts.direct_cached_cell_hits += 1
            return mass
        end
    end
    throw(ArgumentError(
        "direct HCV proposal time is outside the cached active grid"))
end

PDMPSamplers.joint_signed_group_enabled(
    envelope::OMRFDirectPersonNodeHCVEnvelope) =
    envelope.joint_group_size > 0

function PDMPSamplers.joint_signed_cell_mass!(cv,
        envelope::OMRFDirectPersonNodeHCVEnvelope, state, flow,
        left::Float64, right::Float64)
    envelope.joint_group_size > 0 || throw(ArgumentError(
        "joint cell requested for an ordinary person envelope"))
    if envelope.shared_current !== nothing
        oracle_check = _omrf_base_oracle(cv.deterministic_gradient!)
        oracle_check.analytic_prior isa OMRFIndependentAnalyticPrior ||
            throw(ArgumentError(
                "shared-current joint groups require an independent analytic prior"))
        _omrf_actual_flow(flow) isa BouncyParticle || throw(ArgumentError(
            "shared-current joint groups require linear BPS"))
    end
    width = right - left
    width_bin = width <= 0.001 ? 1 :
        width <= 0.005 ? 2 :
        width <= 0.01 ? 3 :
        width <= 0.02 ? 4 :
        width <= 0.05 ? 5 :
        width <= 0.1 ? 6 :
        width <= 0.2 ? 7 :
        width <= 0.5 ? 8 :
        width <= 1.0 ? 9 : 10
    envelope.context.counts.joint_cell_width_bins[width_bin] += 1
    PDMPSamplers.component_cell_scales!(envelope.cell_scales,
        envelope, state, flow, left, right)
    cell = envelope.active_cell
    oracle = _omrf_base_oracle(cv.deterministic_gradient!)
    Ud = _omrf_joint_signed_deterministic_upper(oracle,
        oracle.active_free, state, flow, left, right)
    B = size(envelope.joint_groups, 2)
    width2 = nextfloat((right - left)^2)
    equal_groups = envelope.joint_group_size == 4
    n_terms = envelope.joint_group_size * envelope.n_nodes
    gamma = n_terms * eps(Float64) / (1 - n_terms * eps(Float64))
    total = 0.0
    @inbounds for b in 1:B
        left_signed = 0.0; right_signed = 0.0
        left_abs = 0.0; right_abs = 0.0; curvature = 0.0
        m_b = envelope.joint_group_lengths[b]
        for slot in 1:m_b
            person = envelope.joint_groups[slot, b]
            left_signed += envelope.person_left_signed[person]
            right_signed += envelope.person_right_signed[person]
            left_abs += envelope.person_left_abs[person]
            right_abs += envelope.person_right_abs[person]
            curvature += envelope.person_curvature[person]
        end
        if !equal_groups
            n_terms_b = m_b * envelope.n_nodes
            gamma = n_terms_b * eps(Float64) /
                (1 - n_terms_b * eps(Float64))
        end
        upper_left = equal_groups ?
            nextfloat(left_signed + gamma * left_abs) :
            nextfloat(left_signed + nextfloat(gamma * left_abs))
        upper_right = equal_groups ?
            nextfloat(right_signed + gamma * right_abs) :
            nextfloat(right_signed + nextfloat(gamma * right_abs))
        curvature = nextfloat(curvature * (1 + gamma))
        Ub = nextfloat(max(upper_left, upper_right) +
            nextfloat(curvature * width2))
        cap = if equal_groups
            max(0.0, nextfloat(Ud + B * Ub))
        else
            # The final JS16 mark has four people. Bracket N/m_b before
            # multiplication, choosing the endpoint by the signed bound.
            scale = envelope.n_persons / m_b
            scaled = Ub >= 0 ? nextfloat(nextfloat(scale) * Ub) :
                nextfloat(prevfloat(scale) * Ub)
            max(0.0, nextfloat(Ud + scaled))
        end
        envelope.joint_cell_caps[b, cell] = cap
        # Equal JS4 retains its original cumulative CDF exactly. JS16 uses
        # integer group sizes as weights: m_b*C_b, then divides by N once.
        total = equal_groups ? nextfloat(total + cap) :
            nextfloat(total + nextfloat(m_b * cap))
        envelope.joint_cell_cumulative[b, cell] = total
    end
    mass = equal_groups ? nextfloat(total / B) :
        nextfloat(total / envelope.n_persons)
    envelope.joint_cell_mass[cell] = mass
    return mass
end

function PDMPSamplers.joint_signed_active_mass!(
        envelope::OMRFDirectPersonNodeHCVEnvelope, state, t)
    _omrf_direct_activate_time!(envelope, state, t)
    return @inbounds envelope.joint_cell_mass[envelope.active_cell]
end

function PDMPSamplers.joint_signed_draw_mark!(rng::Random.AbstractRNG,
        envelope::OMRFDirectPersonNodeHCVEnvelope, subset::Vector{Int})
    cell = envelope.active_cell
    B = size(envelope.joint_groups, 2)
    total = @inbounds envelope.joint_cell_cumulative[B, cell]
    total > 0 || return 0.0
    u = rand(rng) * total
    lo = 1; hi = B
    @inbounds while lo < hi
        mid = (lo + hi) >>> 1
        if envelope.joint_cell_cumulative[mid, cell] <= u
            lo = mid + 1
        else
            hi = mid
        end
    end
    m_b = @inbounds envelope.joint_group_lengths[lo]
    resize!(subset, m_b)
    @inbounds for slot in 1:m_b
        subset[slot] = envelope.joint_groups[slot, lo]
    end
    return @inbounds envelope.joint_cell_caps[lo, cell]
end

function PDMPSamplers.component_scales!(out,
        envelope::OMRFDirectPersonNodeHCVEnvelope, state, flow, t)
    out[1] = _omrf_direct_activate_time!(envelope, state, t)
    return out
end

@inline function PDMPSamplers.screening_residual_bound(
    envelope::OMRFDirectPersonNodeHCVEnvelope, state, flow, t)::Float64
    return _omrf_direct_activate_time!(envelope, state, t)
end

function PDMPSamplers.prepare_residual_sampling!(
        envelope::OMRFDirectPersonNodeHCVEnvelope, state, flow, t)
    _omrf_direct_activate_time!(envelope, state, t)
    return nothing
end

function PDMPSamplers.observation_residual_bound(
        envelope::OMRFDirectPersonNodeHCVEnvelope, observation::Integer)
    envelope.signed_only_person_roof && throw(ArgumentError(
        "signed-only person roofs do not construct factor cell bounds"))
    # The internal table is always factor-indexed even when the public
    # sampling population is persons for the hybrid route.
    n_factors = envelope.n_persons * envelope.n_nodes
    1 <= observation <= n_factors ||
        throw(BoundsError(1:n_factors, observation))
    cell = envelope.active_cell
    iszero(cell) && throw(ArgumentError(
        "direct HCV factor bound requested before cell activation"))
    @inbounds if observation == 1
        return envelope.cell_cumulative[1, cell]
    end
    upper = envelope.cell_cumulative[observation, cell]
    lower = envelope.cell_cumulative[observation - 1, cell]
    # Adjacent-prefix subtraction can lose several ulps when the selected
    # factor is small relative to the accumulated mass. Add the standard
    # floating summation/subtraction error scale so the recovered value is an
    # upper bound on the bit-identical scalar value kept during construction.
    return (upper - lower) + 4eps(max(abs(upper), abs(lower), 1.0))
end


function PDMPSamplers.draw_component(rng::Random.AbstractRNG,
        envelope::OMRFDirectPersonNodeHCVEnvelope,
        B::Real)::Int
    cell = envelope.active_cell
    iszero(cell) && throw(ArgumentError(
        "direct HCV distinguished person requested before cell activation"))
    target = rand(rng) * B
    left = 1
    right = envelope.n_persons
    @inbounds while left < right
        middle = (left + right) >>> 1
        if envelope.person_cell_cumulative[middle, cell] < target
            left = middle + 1
        else
            right = middle
        end
    end
    return left
end

@inline function _omrf_direct_person_cell_bound(
        envelope::OMRFDirectPersonNodeHCVEnvelope, person::Int)
    envelope.signed_only_person_roof &&
        return @inbounds envelope.person_workspace[person]
    cell = envelope.active_cell
    @inbounds return person == 1 ? envelope.person_cell_cumulative[1, cell] :
        envelope.person_cell_cumulative[person, cell] -
        envelope.person_cell_cumulative[person - 1, cell]
end

PDMPSamplers.draw_distinguished(rng::Random.AbstractRNG,
    ::OMRFDirectPersonNodeHCVEnvelope, observation::Int) = observation

@inline function PDMPSamplers.begin_subsampling_mark!(
        oracle::OMRFAnalyticSubsamplingOracle)
    oracle.cached_residual_valid = false
    oracle.cached_deterministic_signed_rate = NaN
    return nothing
end

@inline function _omrf_direct_factor_point_bound(
        envelope::OMRFDirectPersonNodeHCVEnvelope, factor::Int)
    context = envelope.context
    N = envelope.n_persons
    node0, person0 = divrem(factor - 1, N)
    node = node0 + 1
    q = context.seen[node] - 1
    iszero(q) && return 0.0
    t = envelope.active_time
    minimum_displacement = 0.0
    maximum_displacement = 0.0
    minimum_velocity = 0.0
    maximum_velocity = 0.0
    if envelope.path_kind == UInt8(2)
        sine, cosine = sincos(t)
        @inbounds for u in 1:q
            displacement = envelope.path_delta0[u, factor] +
                envelope.path_z[u, factor] * cosine +
                envelope.path_harmonic_sine[u, factor] * sine
            velocity = -envelope.path_z[u, factor] * sine +
                envelope.path_harmonic_sine[u, factor] * cosine
            minimum_displacement = min(minimum_displacement, displacement)
            maximum_displacement = max(maximum_displacement, displacement)
            minimum_velocity = min(minimum_velocity, velocity)
            maximum_velocity = max(maximum_velocity, velocity)
        end
    else
        @inbounds for u in 1:q
            displacement = envelope.path_delta0[u, factor] +
                t * envelope.path_z[u, factor]
            minimum_displacement = min(minimum_displacement, displacement)
            maximum_displacement = max(maximum_displacement, displacement)
        end
    end
    displacement_range = maximum_displacement - minimum_displacement
    velocity_range = if envelope.path_kind == UInt8(2)
        maximum_velocity - minimum_velocity
    elseif envelope.path_kind == UInt8(1)
        envelope.path_zigzag_velocity_sum[factor]
    else
        envelope.path_velocity_range[factor]
    end
    return context.hcv_active ? _omrf_categorical_hcv_bound(
        velocity_range, displacement_range) :
        velocity_range * displacement_range / 4
end

"""Size-biased direct-table draw with a cheap selected-factor point screen.

The cell table still defines the exact proposal distribution.  The point
bound only performs an additional valid thinning step before construction of
the dense deterministic HCV gradient; it therefore does not alter the
accepted marked-event law.
"""
function PDMPSamplers.draw_subset!(rng::Random.AbstractRNG,
        cv::PDMPSamplers.SubsampledControlVariate{F,O,E},
        D::Real, B::Real) where {F,O,
            E<:OMRFDirectPersonNodeHCVEnvelope}
    envelope = cv.envelope
    oracle = cv.residual_oracle
    oracle isa PDMPSamplers.WithResidualStats && (oracle = oracle.f)
    m = cv.m
    isone(m) || throw(ArgumentError(
        "factor-roof person HCV currently requires one selected person"))
    total = D + B
    residual_source = !(iszero(B) ||
        (!iszero(D) && rand(rng) * total <= D))
    if residual_source
        cv.subset[1] = PDMPSamplers.draw_component(rng, envelope, B)
    else
        PDMPSamplers.draw_base_subset!(rng, cv, cv.subset_design)
    end
    PDMPSamplers.record_subsampling_mark_source!(
        cv.residual_oracle, residual_source)
    person = cv.subset[1]
    cell_weight = _omrf_direct_person_cell_bound(envelope, person)
    point_weight = 0.0
    if envelope.signed_only_person_roof
        @inbounds for node in 1:envelope.n_nodes
            factor = (node - 1) * envelope.n_persons + person
            point_weight = nextfloat(point_weight +
                _omrf_direct_factor_point_bound(envelope, factor))
        end
        point_weight = min(cell_weight, point_weight)
    else
        @inbounds for node in 1:envelope.n_nodes
            factor = (node - 1) * envelope.n_persons + person
            factor_cell_weight = PDMPSamplers.observation_residual_bound(
                envelope, factor)
            point_weight += min(factor_cell_weight,
                _omrf_direct_factor_point_bound(envelope, factor))
        end
    end
    scale = envelope.n_persons
    oracle isa OMRFAnalyticSubsamplingOracle &&
        (oracle.lazy_residual_subset_bound = scale * point_weight)
    return D + scale * cell_weight
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


@inline function _omrf_exact_selected_residual_abs_rate(
        state, residual, flow::PDMPSamplers.PreconditionedDynamics,
        touched_positions)
    return _omrf_exact_selected_residual_abs_rate(
        state, residual, flow.dynamics, touched_positions)
end

@inline function _omrf_exact_selected_residual_abs_rate(
        state, residual, ::ZigZag, touched_positions)
    result = 0.0
    @inbounds for index in touched_positions
        result += abs(state.ξ.θ[index] * residual[index])
    end
    return result
end

@inline function _omrf_exact_selected_residual_abs_rate(
        state, residual,
        ::Union{BouncyParticle,PDMPSamplers.AnyBoomerang},
        touched_positions)
    result = 0.0
    @inbounds for index in touched_positions
        result += state.ξ.θ[index] * residual[index]
    end
    return abs(result)
end

function _omrf_prepare_exact_selected_residual!(
        oracle::OMRFAnalyticSubsamplingOracle, state, subset)
    oracle.ctx.counts.selected_refinement_calls += 1
    fill!(oracle.cached_residual, 0.0)
    oracle(oracle.cached_residual, state.ξ.x, subset, oracle.anchor)
    # The generic event loop consumes this value before drawing another mark
    # or changing an anchor.  draw_subset! invalidates it at mark ownership,
    # and anchor/source installation invalidates it separately.  This avoids
    # evaluating the same selected residual once for refinement and again for
    # the final thinning test.
    oracle.cached_residual_valid = true
    return nothing
end

function PDMPSamplers.subsampling_subset_bound(
        oracle::OMRFAnalyticSubsamplingOracle, state, gradient, flow,
        D, M, subset, scale)
    oracle.exact_selected_refinement || return PDMPSamplers.signed_subset_bound(
        state, gradient, flow, D, M)
    supported_flow = flow isa Union{ZigZag,BouncyParticle,
        PDMPSamplers.AnyBoomerang,PDMPSamplers.PreconditionedDynamics}
    supported_flow || return PDMPSamplers.signed_subset_bound(
        state, gradient, flow, D, M)
    oracle.cached_residual_valid ||
        _omrf_prepare_exact_selected_residual!(oracle, state, subset)
    residual_rate_bound = scale *
        _omrf_exact_selected_residual_abs_rate(
            state, oracle.cached_residual, flow,
            oracle.residual_context.touched_positions)
    deterministic_actual = gradient === nothing ?
        max(0.0, oracle.cached_deterministic_signed_rate) :
        λ(state, gradient, flow)
    return min(M, deterministic_actual + residual_rate_bound)
end

function PDMPSamplers.subsampling_cached_residual!(
        oracle::OMRFAnalyticSubsamplingOracle, out, state, subset, anchor)
    oracle.cached_residual_valid || return false
    copyto!(out, oracle.cached_residual)
    oracle.cached_residual_valid = false
    return true
end

function PDMPSamplers.subsampling_residual_subset_bound(
        oracle::OMRFAnalyticSubsamplingOracle, state, gradient, flow,
        D, M, subset, scale)
    # The direct HCV envelope stores a tighter bound for the selected subset
    # while drawing it; without one, the drawn roof is kept.
    if isfinite(oracle.lazy_residual_subset_bound)
        bound = D + oracle.lazy_residual_subset_bound
        oracle.lazy_residual_subset_bound = NaN
        return min(M, bound)
    end
    return M
end

function PDMPSamplers.subsampling_residual_subset_bound(
        oracle::OMRFAnalyticSubsamplingOracle, state, flow,
        D, M, subset, scale)
    return PDMPSamplers.subsampling_residual_subset_bound(
        oracle, state, nothing, flow, D, M, subset, scale)
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


@inline function _omrf_harmonic_phase_in_cell(phase, left, right)
    period = 2pi
    occurrence = phase + ceil((left - phase) / period) * period
    return occurrence <= right
end

"""Exact closed interval for `c + a*cos(t) + b*sin(t)`."""
@inline function _omrf_harmonic_cell_interval(c, a, b, left, right)
    left_value = c + a * cos(left) + b * sin(left)
    right_value = c + a * cos(right) + b * sin(right)
    lower = min(left_value, right_value)
    upper = max(left_value, right_value)
    amplitude = hypot(a, b)
    if !iszero(amplitude)
        maximum_phase = atan(b, a)
        _omrf_harmonic_phase_in_cell(maximum_phase, left, right) &&
            (upper = c + amplitude)
        _omrf_harmonic_phase_in_cell(maximum_phase + pi, left, right) &&
            (lower = c - amplitude)
    end
    return prevfloat(lower), nextfloat(upper)
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
    bound_type == "direct_factor_roof_person_hcv" || throw(ArgumentError(
        "unsupported OMRF residual geometry \"$bound_type\"; only " *
        "direct_factor_roof_person_hcv is supported"))
    residual_context.use_hcv && residual_context.hcv_active || throw(
        ArgumentError("direct HCV geometry requires the active HCV"))
    cell_scope = _haskey(envelope_spec, :direct_cv_cell_scope) ?
        String(_rget(envelope_spec, :direct_cv_cell_scope)) : "cell"
    cell_scope == "cell" || throw(ArgumentError(
        "direct HCV cell scope must be cell"))
    return _omrf_direct_hcv_envelope_from_spec(
        envelope_spec, anchor, residual_context)
end

"""Direct whole-person HCV envelope with the roof options of `spec`."""
function _omrf_direct_hcv_envelope_from_spec(spec, anchor,
        residual_context::OMRFResidualContext)
    return _new_omrf_direct_hcv_envelope(anchor, residual_context;
        person_signed_roof=_haskey(spec, :person_signed_roof) &&
            Bool(_rget(spec, :person_signed_roof)),
        signed_only_person_roof=_haskey(spec, :signed_only_person_roof) &&
            Bool(_rget(spec, :signed_only_person_roof)),
        joint_group_size=_haskey(spec, :joint_signed_group_size) ?
            Int(_rget(spec, :joint_signed_group_size)) : 0,
        shared_current=_haskey(spec, :shared_current_categorical) &&
            Bool(_rget(spec, :shared_current_categorical)))
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

mutable struct OMRFWindowDiagnostics
    width::Float64
    origin::Float64
    next_boundary::Float64
    last_time::Float64
    last_wall_ns::UInt64
    last_reflections::Int
    last_roof_proposals::Int
    last_subset_evaluations::Int
    last_proposal_count::Int
    last_certificate_sum::Float64
    last_intrinsic_sum::Float64
    last_allocated_bytes::UInt64
    anchor_distance_sum::Float64
    anchor_distance_count::Int
    diagnostic_seconds::Float64
    rows::Vector{Dict{String,Any}}
    pending_post::Any
    row_pushed_this_call::Bool
    last_pipeline_values::Dict{Symbol,Float64}
end

function OMRFWindowDiagnostics(width::Real)
    width > 0 || throw(ArgumentError("OMRF diagnostic window width must be positive"))
    return OMRFWindowDiagnostics(Float64(width), NaN, NaN, NaN, UInt(0),
        0, 0, 0, 0, 0.0, 0.0, UInt(0), 0.0, 0,
        0.0, Dict{String,Any}[], nothing, false, Dict{Symbol,Float64}())
end

@inline _omrf_allocated_bytes() = _OMRF_RESEARCH_DIAGNOSTICS ? begin
    gc = Base.gc_num()
    UInt(gc.total_allocd + gc.allocd)
end : UInt(0)

"""Measure the residual funnel at one retained BPS trajectory state.

This is deliberately opt-in through the existing diagnostic-window path. It
reuses the exact factor residual oracle for the diagnostic only; it is not
part of ordinary sampling. The signed quantity is the existing person-node
pattern certificate, while the compact quantity is the pointwise compact
aggregate bound.
"""
function _omrf_window_funnel(manager, state, flow)
    actual_flow = flow isa PDMPSamplers.PreconditionedDynamics ?
        flow.dynamics : flow
    actual_flow isa Union{BouncyParticle,PDMPSamplers.AnyBoomerang} ||
        return nothing
    manager.active_idx > 0 || return nothing
    prepared = manager.entries[manager.active_idx].state
    envelope = prepared.envelope
    oracle = manager.oracle
    context = oracle.residual_context
    compact = PDMPSamplers.total_residual_bound(
        envelope, state, flow, 0.0)
    signed = 0.0
    people, nodes = size(context.X)
    @inbounds for node in 1:nodes, person in 1:people
        signed += _omrf_person_node_pattern_bound_at(
            context, oracle.anchor, state, flow, person, node, 0.0)
    end
    deterministic_gradient = zeros(Float64, length(state.ξ.x))
    oracle(deterministic_gradient, state.ξ.x)
    corrected_deterministic = copy(deterministic_gradient)
    if actual_flow isa PDMPSamplers.AnyBoomerang
        z = state.ξ.x .- actual_flow.μ
        mul!(corrected_deterministic, actual_flow.Γ, z,
            -1.0, 1.0)
    end
    deterministic_rate = PDMPSamplers.λ(
        state, corrected_deterministic, flow)
    residual = zeros(Float64, length(state.ξ.x))
    full_gradient = copy(deterministic_gradient)
    one_factor = [0]
    exact_sum = 0.0
    stochastic_sum = 0.0
    @inbounds for factor in 1:(people * nodes)
        one_factor[1] = factor
        oracle(residual, state.ξ.x, one_factor, oracle.anchor)
        exact_sum += _omrf_exact_selected_residual_abs_rate(
            state, residual, flow, context.touched_positions)
        full_gradient .+= residual
        stochastic_gradient = deterministic_gradient .+
            (people * nodes) .* residual
        if actual_flow isa PDMPSamplers.AnyBoomerang
            z = state.ξ.x .- actual_flow.μ
            mul!(stochastic_gradient, actual_flow.Γ, z,
                -1.0, 1.0)
        end
        stochastic_sum += PDMPSamplers.λ(
            state, stochastic_gradient, flow)
    end
    if actual_flow isa PDMPSamplers.AnyBoomerang
        z = state.ξ.x .- actual_flow.μ
        mul!(full_gradient, actual_flow.Γ, z, -1.0, 1.0)
    end
    full_rate = PDMPSamplers.λ(state, full_gradient, flow)
    return (compact=compact, signed=signed, exact_sum=exact_sum,
        exact_mean=exact_sum / (people * nodes),
        stochastic_rate_mean=stochastic_sum / (people * nodes),
        deterministic_rate=deterministic_rate,
        deterministic_grid_roof=deterministic_rate,
        actual_stochastic_rate=stochastic_sum / (people * nodes),
        full_rate=full_rate)
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
    warmup_finish_preparation_seconds::Float64
    warmup_finish_selection_seconds::Float64
    warmup_finish_native_full_data_scans::Int
    selection_seconds::Float64
    selection_component_evaluations::Int
    window_diagnostics::Any
    selection_policy::Symbol
    pending_boundary_idx::Int
    pending_boundary_time::Float64
    boundary_generation::UInt64
    pending_boundary_generation::UInt64
    selection_boundaries::Int
    # Position-policy scores, see `_omrf_position_score`.
    position_center::Vector{Float64}
    position_weights::Matrix{Float64}
    position_offsets::Vector{Float64}
    position_metric::Symbol
    position_valid::Bool
end

@inline function _omrf_window_counter(stats, name::Symbol)
    return Int(getproperty(stats, name))
end

const _OMRF_PIPELINE_WINDOW_FIELDS = (
    :clock_candidates, :grid_bound_evaluations,
    :pointwise_screen_evaluations, :pointwise_screen_passes,
    :subset_draws, :selected_subset_bound_evaluations,
    :selected_subset_bound_passes, :residual_oracle_evaluations,
    :final_thinning_acceptances, :reflection_events, :sticky_events,
    :grid_bound_seconds, :pointwise_screen_seconds, :subset_draw_seconds,
    :selected_subset_refinement_seconds, :residual_oracle_seconds,
    :final_thinning_seconds, :reflection_update_seconds,
    :sticky_update_seconds,
)

function _omrf_window_values(manager::OMRFAnchorManager, stats)
    names = _OMRF_PIPELINE_WINDOW_FIELDS
    values = Dict{Symbol,Float64}(
        name => Float64(getproperty(stats, name)) for name in names)
    values[:anchor_selections] = Float64(manager.main_selections)
    values[:anchor_switches] = Float64(manager.main_activations)
    values[:anchor_refreshes] = Float64(manager.main_refreshes)
    values[:native_full_data_scans] = Float64(
        manager.oracle.ctx.counts.native_full_data_scans -
        manager.warmup_finish_native_full_data_scans)
    return values
end

function _omrf_nearest_anchor_distance(manager::OMRFAnchorManager, x)
    isempty(manager.entries) && return NaN
    return sqrt(minimum(_omrf_anchor_distance2(manager, x, entry.state)
        for entry in manager.entries))
end

function _omrf_window_post!(manager::OMRFAnchorManager, state, flow)
    window = manager.window_diagnostics
    window === nothing && return nothing
    # Funnel evaluation scans every person-node factor and is intended only
    # for bounded research diagnostics.  Evaluate the post-event side only
    # when the matching pre-event call actually emitted a window row.
    window.row_pushed_this_call || return nothing
    diagnostic_started = time_ns()
    funnel = _omrf_window_funnel(manager, state, flow)
    window.diagnostic_seconds += (time_ns() - diagnostic_started) * 1e-9
    distance = manager.active_idx > 0 ? sqrt(_omrf_anchor_distance2(
        manager, state.ξ.x, manager.entries[manager.active_idx].state)) : NaN
    snapshot = (funnel=funnel, distance=distance)
    if window.row_pushed_this_call && !isempty(window.rows)
        row = window.rows[end]
        row["post_anchor_distance"] = distance
        row["post_compact_residual_bound"] = funnel === nothing ? NaN : funnel.compact
        row["post_signed_predictor_bound"] = funnel === nothing ? NaN : funnel.signed
        row["post_stochastic_bps_rate_mean"] = funnel === nothing ? NaN : funnel.stochastic_rate_mean
        row["post_exact_factor_rate_sum"] = funnel === nothing ? NaN : funnel.exact_sum
        row["post_exact_factor_rate_mean"] = funnel === nothing ? NaN : funnel.exact_mean
        row["post_full_gradient_bps_rate"] = funnel === nothing ? NaN : funnel.full_rate
        row["post_deterministic_rate"] = funnel === nothing ? NaN : funnel.deterministic_rate
        row["post_deterministic_grid_roof"] = funnel === nothing ? NaN : funnel.deterministic_grid_roof
        row["post_actual_stochastic_rate"] = funnel === nothing ? NaN : funnel.actual_stochastic_rate
    else
        window.pending_post = snapshot
    end
    return nothing
end

function _record_omrf_window!(manager::OMRFAnchorManager, state, flow, stats;
        stage::Symbol=:pre)
    window = manager.window_diagnostics
    window === nothing && return nothing
    stage === :post && return _omrf_window_post!(manager, state, flow)
    window.row_pushed_this_call = false
    reflections = _omrf_window_counter(stats, :reflections_events)
    roofs = _omrf_window_counter(stats, :subsampling_cell_roof_proposals)
    subsets = _omrf_window_counter(stats, :subsampling_subset_evaluations)
    proposal_count = manager.oracle.ctx.counts.proposal_count
    certificate_sum = manager.oracle.ctx.counts.proposal_certificate_ratio_sum
    intrinsic_sum = manager.oracle.ctx.counts.proposal_intrinsic_ratio_sum
    time_now = Float64(state.t[])
    if isnan(window.origin)
        # Trace times are absolute (warmup plus retained phase).  Establish a
        # local retained-phase origin at the first main-phase callback and
        # initialize cumulative baselines there; finalization extends this
        # local origin by the requested retained horizon.
        window.origin = time_now
        window.last_time = time_now
        window.next_boundary = time_now + window.width
        window.last_wall_ns = UInt(time_ns())
        window.last_reflections = reflections
        window.last_roof_proposals = roofs
        window.last_subset_evaluations = subsets
        window.last_proposal_count = proposal_count
        window.last_certificate_sum = certificate_sum
        window.last_allocated_bytes = _omrf_allocated_bytes()
        empty!(window.last_pipeline_values)
        merge!(window.last_pipeline_values, _omrf_window_values(manager, stats))
    end
    distance = sqrt(_omrf_anchor_distance2(manager, state.ξ.x,
        manager.entries[manager.active_idx].state))
    window.anchor_distance_sum += distance
    window.anchor_distance_count += 1
    # Do not pay an O(NP) diagnostic scan at every event.  A funnel is needed
    # only when this call crosses a requested diagnostic-window boundary.
    pre_funnel = nothing
    if time_now >= window.next_boundary
        diagnostic_started = time_ns()
        pre_funnel = _omrf_window_funnel(manager, state, flow)
        window.diagnostic_seconds += (time_ns() - diagnostic_started) * 1e-9
    end
    pre_distance = distance
    while time_now >= window.next_boundary
        elapsed = window.next_boundary - window.last_time
        elapsed = max(elapsed, eps(Float64))
        wall_now = UInt(time_ns())
        wall_seconds = max(0.0, (wall_now - window.last_wall_ns) * 1e-9)
        reflection_delta = reflections - window.last_reflections
        roof_delta = roofs - window.last_roof_proposals
        subset_delta = subsets - window.last_subset_evaluations
        proposal_delta = proposal_count - window.last_proposal_count
        certificate_delta = certificate_sum - window.last_certificate_sum
        intrinsic_delta = intrinsic_sum - window.last_intrinsic_sum
        allocated_now = _omrf_allocated_bytes()
        allocated_delta = allocated_now - window.last_allocated_bytes
        current_pipeline = _omrf_window_values(manager, stats)
        pipeline_delta = Dict{Symbol,Float64}(
            name => current_pipeline[name] -
                get(window.last_pipeline_values, name, 0.0)
            for name in keys(current_pipeline))
        funnel = pre_funnel
        pending_post = window.pending_post
        push!(window.rows, Dict{String,Any}(
            "window_start" => window.last_time - window.origin,
            "window_end" => window.next_boundary - window.origin,
            "anchor_mean_distance" => window.anchor_distance_count == 0 ?
                NaN : window.anchor_distance_sum / window.anchor_distance_count,
            "reflection_count" => reflection_delta,
            "reflection_rate" => reflection_delta / elapsed,
            "roof_proposals" => roof_delta,
            "proposal_rate" => roof_delta / elapsed,
            "subset_evaluations" => subset_delta,
            "subsets_per_reflection" => reflection_delta == 0 ?
                NaN : subset_delta / reflection_delta,
            "proposals_per_reflection" => reflection_delta == 0 ?
                NaN : roof_delta / reflection_delta,
            "proposal_count" => proposal_delta,
            "certificate_ratio" => proposal_delta == 0 ?
                NaN : certificate_delta / proposal_delta,
            "intrinsic_ratio" => proposal_delta == 0 ?
                NaN : intrinsic_delta / proposal_delta,
            "wall_seconds_per_pdmp_time" => wall_seconds / elapsed,
            "main_phase_wall_seconds" => wall_seconds,
            "anchor_distance_active" => pre_distance,
            "nearest_prepared_anchor_distance" =>
                _omrf_nearest_anchor_distance(manager, state.ξ.x),
            "allocated_bytes" => allocated_delta,
            "allocated_bytes_per_roof_proposal" => roof_delta == 0 ?
                NaN : allocated_delta / roof_delta,
            "allocated_bytes_per_subset_evaluation" => subset_delta == 0 ?
                NaN : allocated_delta / subset_delta,
            "anchor_distance" => pre_distance,
            "compact_residual_bound" => funnel === nothing ? NaN : funnel.compact,
            "signed_predictor_bound" => funnel === nothing ? NaN : funnel.signed,
            "stochastic_bps_rate_mean" => funnel === nothing ? NaN : funnel.stochastic_rate_mean,
            "exact_factor_rate_sum" => funnel === nothing ? NaN : funnel.exact_sum,
            "exact_factor_rate_mean" => funnel === nothing ? NaN : funnel.exact_mean,
            "full_gradient_bps_rate" => funnel === nothing ? NaN : funnel.full_rate,
            "deterministic_rate" => funnel === nothing ? NaN : funnel.deterministic_rate,
            "deterministic_grid_roof" => funnel === nothing ? NaN : funnel.deterministic_grid_roof,
            "actual_stochastic_rate" => funnel === nothing ? NaN : funnel.actual_stochastic_rate,
            "pre_anchor_distance" => pre_distance,
            "pre_compact_residual_bound" => funnel === nothing ? NaN : funnel.compact,
            "pre_signed_predictor_bound" => funnel === nothing ? NaN : funnel.signed,
            "pre_stochastic_bps_rate_mean" => funnel === nothing ? NaN : funnel.stochastic_rate_mean,
            "pre_exact_factor_rate_sum" => funnel === nothing ? NaN : funnel.exact_sum,
            "pre_exact_factor_rate_mean" => funnel === nothing ? NaN : funnel.exact_mean,
            "pre_full_gradient_bps_rate" => funnel === nothing ? NaN : funnel.full_rate,
            "post_anchor_distance" => pending_post === nothing ? NaN : pending_post.distance,
            "post_compact_residual_bound" => pending_post === nothing || pending_post.funnel === nothing ? NaN : pending_post.funnel.compact,
            "post_signed_predictor_bound" => pending_post === nothing || pending_post.funnel === nothing ? NaN : pending_post.funnel.signed,
            "post_stochastic_bps_rate_mean" => pending_post === nothing || pending_post.funnel === nothing ? NaN : pending_post.funnel.stochastic_rate_mean,
            "post_exact_factor_rate_sum" => pending_post === nothing || pending_post.funnel === nothing ? NaN : pending_post.funnel.exact_sum,
            "post_exact_factor_rate_mean" => pending_post === nothing || pending_post.funnel === nothing ? NaN : pending_post.funnel.exact_mean,
            "post_full_gradient_bps_rate" => pending_post === nothing || pending_post.funnel === nothing ? NaN : pending_post.funnel.full_rate,
            "post_deterministic_rate" => pending_post === nothing || pending_post.funnel === nothing ? NaN : pending_post.funnel.deterministic_rate,
            "post_deterministic_grid_roof" => pending_post === nothing || pending_post.funnel === nothing ? NaN : pending_post.funnel.deterministic_grid_roof,
            "post_actual_stochastic_rate" => pending_post === nothing || pending_post.funnel === nothing ? NaN : pending_post.funnel.actual_stochastic_rate))
        for name in keys(current_pipeline)
            row = window.rows[end]
            row[string(name)] = pipeline_delta[name]
            row[string(Symbol("cumulative_", name))] = current_pipeline[name]
        end
        window.pending_post = nothing
        window.row_pushed_this_call = true
        window.last_time = window.next_boundary
        window.next_boundary += window.width
        window.last_wall_ns = wall_now
        window.last_reflections = reflections
        window.last_roof_proposals = roofs
        window.last_subset_evaluations = subsets
        window.last_proposal_count = proposal_count
        window.last_certificate_sum = certificate_sum
        window.last_intrinsic_sum = intrinsic_sum
        window.last_allocated_bytes = allocated_now
        empty!(window.last_pipeline_values)
        merge!(window.last_pipeline_values, current_pipeline)
        window.anchor_distance_sum = 0.0
        window.anchor_distance_count = 0
    end
    return nothing
end

function _finish_omrf_windows!(manager::OMRFAnchorManager, stats, final_time, flow)
    window = manager.window_diagnostics
    window === nothing && return nothing
    if !isnan(window.origin)
        window.next_boundary = max(window.next_boundary, Float64(final_time))
        _record_omrf_window!(manager,
            PDMPState(Float64(final_time), SkeletonPoint(
                copy(manager.oracle.anchor),
                zeros(Float64, length(manager.oracle.anchor)))),
            flow, stats)
    end
    return nothing
end

function _omrf_window_matrix(manager::OMRFAnchorManager)
    window = manager.window_diagnostics
    window === nothing && return zeros(Float64, 0, 17)
    keys = ["window_start", "window_end", "anchor_mean_distance",
        "reflection_count", "reflection_rate", "roof_proposals",
        "proposal_rate", "subset_evaluations", "subsets_per_reflection",
        "proposals_per_reflection", "proposal_count", "certificate_ratio",
        "intrinsic_ratio", "wall_seconds_per_pdmp_time",
        "main_phase_wall_seconds", "anchor_distance_active",
        "nearest_prepared_anchor_distance"]
    window_fields = (_OMRF_PIPELINE_WINDOW_FIELDS...,
        :anchor_selections, :anchor_switches, :anchor_refreshes,
        :native_full_data_scans)
    append!(keys, [string(name) for name in window_fields])
    append!(keys, [string(Symbol("cumulative_", name)) for name in
        window_fields])
    append!(keys, ["allocated_bytes", "allocated_bytes_per_roof_proposal",
        "allocated_bytes_per_subset_evaluation", "compact_residual_bound",
        "signed_predictor_bound", "exact_factor_rate_sum",
        "exact_factor_rate_mean", "full_gradient_bps_rate",
        "anchor_distance", "stochastic_bps_rate_mean",
        "pre_anchor_distance", "pre_compact_residual_bound",
        "pre_signed_predictor_bound", "pre_stochastic_bps_rate_mean",
        "pre_exact_factor_rate_sum", "pre_exact_factor_rate_mean",
        "pre_full_gradient_bps_rate", "post_anchor_distance",
        "post_compact_residual_bound", "post_signed_predictor_bound",
        "post_stochastic_bps_rate_mean", "post_exact_factor_rate_sum",
        "post_exact_factor_rate_mean", "post_full_gradient_bps_rate"])
    keys = (keys..., "deterministic_rate", "deterministic_grid_roof",
        "actual_stochastic_rate", "post_deterministic_rate",
        "post_deterministic_grid_roof", "post_actual_stochastic_rate")
    out = Matrix{Float64}(undef, length(window.rows), length(keys))
    @inbounds for i in axes(out, 1), j in axes(out, 2)
        out[i, j] = Float64(window.rows[i][keys[j]])
    end
    return out
end

_omrf_window_matrix(::Nothing) = zeros(Float64, 0, 0)

function _prepare_omrf_anchor_envelope(manager::OMRFAnchorManager,
        requested, residual_context)
    residual_context.use_hcv || throw(ArgumentError(
        "the OMRF anchor envelope requires the Hessian control variate"))
    return _omrf_direct_hcv_envelope_from_spec(
        manager.envelope_spec, requested, residual_context)
end

function _prepare_omrf_anchor(manager::OMRFAnchorManager,
        anchor::AbstractVector)
    started = time_ns()
    oracle = manager.oracle
    recycled = manager.recycled_state
    if recycled !== nothing && oracle.analytic_prior !== nothing &&
            manager.hcv_active && recycled.residual_context.use_hcv &&
            recycled.envelope isa OMRFDirectPersonNodeHCVEnvelope
        # The recycled state is not installed in the bank, so it can be
        # mutated without exposing a partially prepared anchor. Its envelope
        # callbacks retain the same anchor/context array identities.
        manager.recycled_state = nothing
        copyto!(recycled.anchor, anchor)
        anchor_storage = recycled.envelope.anchor
        anchor_storage === recycled.anchor || copyto!(anchor_storage, anchor)
        recycled.residual_context.hcv_active = manager.hcv_active &&
            recycled.residual_context.use_hcv
        likelihood_started = time_ns()
        _omrf_prior_gradient!(
            recycled.prior_anchor, oracle.analytic_prior, recycled.anchor)
        _cache_omrf_anchor_and_full_gradient!(recycled.residual_context,
            recycled.full_anchor, oracle.analytic_prior, recycled.anchor)
        hessian = recycled.residual_context.likelihood_hessian::Matrix{Float64}
        _omrf_likelihood_hessian!(hessian, recycled.residual_context)
        recycled.envelope.cell_count = 0
        recycled.envelope.active_cell = 0
        recycled.envelope.path_valid = false
        oracle.ctx.counts.compact_anchor_likelihood_seconds +=
            (time_ns() - likelihood_started) * 1e-9
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
    likelihood_started = time_ns()
    _clear_gradient!(oracle.ctx, :full, full_anchor, requested)
    if oracle.analytic_prior === nothing
        _clear_gradient!(oracle.ctx, :prior, prior_anchor, requested)
    else
        _omrf_prior_gradient!(prior_anchor, oracle.analytic_prior, requested)
    end
    oracle.ctx.counts.compact_anchor_likelihood_seconds +=
        (time_ns() - likelihood_started) * 1e-9
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
    copyto!(state.residual_context.active_free, oracle.active_free)
    state.residual_context.active_generation = oracle.active_generation
    oracle.prior_x_valid = false
    return state.envelope
end

function _insert_omrf_anchor!(manager::OMRFAnchorManager,
        state::OMRFPreparedAnchorState)
    manager.position_valid = false
    idx = if length(manager.entries) < manager.capacity
        push!(manager.entries, OMRFAnchorEntry(state, 0))
        length(manager.entries)
    else
        replace_idx = firstindex(manager.entries)
        oldest_age = manager.entries[replace_idx].age
        @inbounds for candidate_idx in (replace_idx + 1):lastindex(manager.entries)
            candidate_age = manager.entries[candidate_idx].age
            if candidate_age > oldest_age
                replace_idx = candidate_idx
                oldest_age = candidate_age
            end
        end
        entry = manager.entries[replace_idx]
        replaced = entry.state
        entry.state = state
        entry.age = 0
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
        oracle.residual_context.hcv_active, nothing, 0, 0, 0.0, 0.0,
        oracle.ctx.counts.native_full_data_scans, 0.0, 0, nothing,
        :adaptive, 0, Inf, UInt64(0), UInt64(0), 0,
        Float64[], zeros(0, 0), Float64[], :scaled, false)
    initial = OMRFPreparedAnchorState(oracle.anchor, oracle.full_anchor,
        oracle.prior_anchor, oracle.residual_context, initial_envelope)
    manager.active_idx = _insert_omrf_anchor!(manager, initial)
    return manager
end


function _omrf_warmup_candidates(trace::PDMPSamplers.AbstractPDMPTrace,
        position::AbstractVector, capacity::Integer;
        retained_fraction::Real=0.5, n_grid::Union{Nothing,Integer}=nothing)
    0 < retained_fraction <= 1 || throw(ArgumentError(
        "warmup candidate retained fraction must lie in (0, 1]"))
    n_grid = isnothing(n_grid) ?
        min(512, max(8 * Int(capacity), Int(capacity))) : Int(n_grid)
    n_grid >= 1 || throw(ArgumentError(
        "warmup candidate grid size must be positive"))
    if trace isa PDMPSamplers.StreamingPDMPTrace
        grid_positions = PDMPSamplers.streaming_dense_time_uniform_states(
            trace, retained_fraction, n_grid)
        candidates = Matrix{Float64}(undef, length(position), n_grid + 1)
        copyto!(view(candidates, :, 1), position)
        copyto!(view(candidates, :, 2:(n_grid + 1)), grid_positions)
        return candidates
    end
    pdmp_trace = trace isa PDMPSamplers.PDMPTrace ? trace :
        PDMPSamplers.PDMPTrace(trace)
    n_events = length(pdmp_trace.times)
    n_events >= 2 || return reshape(collect(Float64, position), :, 1)
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


function _finish_omrf_warmup!(manager::OMRFAnchorManager,
        cv::SubsampledControlVariate, position, trace_manager=nothing)
    previous_anchor = copy(manager.oracle.anchor)
    # The position policy owns a fixed bank explicitly.  The legacy `Inf`
    # refresh sentinel remains meaningful only for the adaptive comparator.
    frozen_bank = manager.selection_policy === :position ||
        (isinf(manager.main_refresh_distance) && manager.capacity > 1)
    trace = trace_manager === nothing ? nothing :
        PDMPSamplers.get_warmup_trace(trace_manager)
    if frozen_bank && trace !== nothing &&
            PDMPSamplers.has_integrable_segment(trace)
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
        _omrf_prepare_position_metric!(manager, :scaled)
        _select_omrf_anchor!(manager, cv, position; phase=:warmup)
    elseif iszero(manager.preparations)
        prepared = _prepare_omrf_anchor(manager, position)
        manager.active_idx = _insert_omrf_anchor!(manager, prepared)
        PDMPSamplers.refresh_anchor_owned!(cv, prepared.anchor)
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
        PDMPSamplers.refresh_anchor_owned!(cv, active.anchor)
    end
    # `with_stats` reconstructs the retained control variate with a copied
    # envelope, while the manager and the wrapped deterministic/residual
    # oracles retain their shared prepared-anchor state.  If warmup selection
    # keeps the manager's current entry, `_select_omrf_anchor!` legitimately
    # returns without calling the refresh callback; the copied live envelope
    # can then enter the first retained event at its pre-finish anchor.  Make
    # the phase boundary authoritative even when the selected entry index did
    # not change.
    active_state = manager.entries[manager.active_idx].state
    envelope_anchor = PDMPSamplers.trajectory_scale_anchor(cv.envelope)
    needs_sync = cv.anchor != active_state.anchor ||
        manager.oracle.anchor != active_state.anchor ||
        (envelope_anchor !== nothing && envelope_anchor != active_state.anchor)
    needs_sync && PDMPSamplers.refresh_anchor_owned!(cv, active_state.anchor)
    # This phase-boundary snapshot makes the no-full-data-main invariant
    # directly testable.  A preparation includes the exact complete anchor
    # sum, whether obtained through BridgeStan or the native OMRF algebra.
    manager.warmup_finish_preparations = manager.preparations
    manager.warmup_finish_preparation_seconds = manager.preparation_seconds
    manager.warmup_finish_selection_seconds = manager.selection_seconds
    manager.warmup_finish_native_full_data_scans =
        manager.oracle.ctx.counts.native_full_data_scans
    return anchor_changed || hcv_changed || needs_sync
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


# Position-policy anchor scores.  For a shared positive definite metric M and
# y = x - c (c the first anchor),
#   (x - a_k)' M (x - a_k) = y'My - 2 y'w_k + s_k,
# with w_k = M (a_k - c) and s_k = (a_k - c)' M (a_k - c).  The common y'My
# term does not change the nearest anchor, so the score s_k - 2 y'w_k makes
# every pairwise comparison affine in x: Voronoi boundaries are hyperplanes
# for any M.  `distance_scales` gives M = diag(scales)^-2; the "hessian"
# metric uses the mean anchor likelihood Hessian plus the identity.
@inline function _omrf_position_score(manager::OMRFAnchorManager, x, k::Int)
    weights = manager.position_weights
    c = manager.position_center
    acc = 0.0
    @inbounds for j in eachindex(x, c)
        acc += (x[j] - c[j]) * weights[j, k]
    end
    return manager.position_offsets[k] - 2acc
end

# Time derivative of the score along a linear path with velocity v.
@inline function _omrf_position_score_rate(manager::OMRFAnchorManager, v,
        k::Int)
    weights = manager.position_weights
    acc = 0.0
    @inbounds for j in eachindex(v)
        acc += v[j] * weights[j, k]
    end
    return -2acc
end

function _omrf_prepare_position_metric!(manager::OMRFAnchorManager,
        metric::Symbol)
    K = length(manager.entries)
    d = length(manager.distance_scales)
    c = copy(manager.entries[1].state.anchor)
    offsets = zeros(d, K)
    for k in 1:K
        @. offsets[:, k] = manager.entries[k].state.anchor - c
    end
    weights = if metric === :scaled
        offsets ./ manager.distance_scales .^ 2
    elseif metric === :hessian
        M = Matrix{Float64}(I, d, d)
        for entry in manager.entries
            H = entry.state.residual_context.likelihood_hessian
            H === nothing && throw(ArgumentError(
                "the hessian anchor metric requires global HCV anchors"))
            M .+= H ./ K
        end
        M * offsets
    else
        throw(ArgumentError("unknown position-anchor metric: $metric"))
    end
    manager.position_center = c
    manager.position_weights = weights
    manager.position_offsets = [dot(view(offsets, :, k), view(weights, :, k))
        for k in 1:K]
    manager.position_metric = metric
    manager.position_valid = true
    return manager
end

function _omrf_position_anchor_index(manager::OMRFAnchorManager, x)
    manager.position_valid ||
        _omrf_prepare_position_metric!(manager, manager.position_metric)
    best = firstindex(manager.entries)
    best_distance = _omrf_position_score(manager, x, best)
    @inbounds for idx in eachindex(manager.entries)
        idx == best && continue
        distance = _omrf_position_score(manager, x, idx)
        tolerance = 64eps(Float64) * max(1.0, abs(distance),
            abs(best_distance))
        # Numerically coincident Voronoi distances use deterministic
        # lowest-index ties; this is a selection rule, not physical travel.
        if distance < best_distance - tolerance
            best = idx
            best_distance = distance
        end
    end
    return best
end

"""Rebuild a fixed position bank from its versioned mathematical spec.

Only values are restored.  Prepared gradients, envelopes, callbacks, and
scratch caches are reconstructed against the current immutable provider.
"""
function _restore_omrf_position_bank!(manager::OMRFAnchorManager,
        cv::SubsampledControlVariate, spec, position)
    Int(_rget(spec, :version)) == 1 || throw(ArgumentError(
        "unsupported position-anchor bank specification version"))
    policy = String(_rget(spec, :selection_policy))
    policy == "position" || throw(ArgumentError(
        "restored anchor bank is not a supported position policy"))
    Symbol(policy) == manager.selection_policy || throw(ArgumentError(
        "restored anchor policy does not match the requested live policy"))
    String(_rget(spec, :tie_policy)) == "lowest_index_outgoing" ||
        throw(ArgumentError("unsupported position-anchor tie policy"))
    saved_names = String.(_rget(spec, :parameter_order))
    saved_names == manager.unc_names || throw(ArgumentError(
        "position-anchor parameter order does not match the live provider"))
    anchors = Matrix{Float64}(_rget(spec, :ordered_anchor_positions))
    scales = _as_float_vector(_rget(spec, :distance_scales))
    size(anchors, 1) == length(manager.unc_names) || throw(DimensionMismatch(
        "restored anchor positions have the wrong parameter dimension"))
    size(anchors, 2) <= manager.capacity || throw(DimensionMismatch(
        "restored anchor bank exceeds the requested capacity"))
    length(scales) == size(anchors, 1) || throw(DimensionMismatch(
        "restored distance scales have the wrong dimension"))
    all(x -> isfinite(x) && x > 0, scales) || throw(ArgumentError(
        "restored distance scales must be positive and finite"))
    prepared = [_prepare_omrf_anchor(manager, view(anchors, :, idx))
        for idx in axes(anchors, 2)]
    empty!(manager.entries)
    manager.active_idx = 0
    manager.recycled_state = nothing
    copyto!(manager.distance_scales, scales)
    for prepared_state in prepared
        _insert_omrf_anchor!(manager, prepared_state)
    end
    metric = _haskey(spec, :distance_metric) ?
        Symbol(String(_rget(spec, :distance_metric))) : :scaled
    _omrf_prepare_position_metric!(manager, metric)
    manager.active_idx = _omrf_position_anchor_index(manager, position)
    PDMPSamplers.refresh_anchor_owned!(cv,
        manager.entries[manager.active_idx].state.anchor)
    manager.pending_boundary_idx = 0
    manager.pending_boundary_time = Inf
    manager.pending_boundary_generation = UInt64(0)
    return true
end

# Anchor nearest just after the current time: scores are compared at x + h v,
# where h is a few ulps of the current time.  A boundary crossing closer than
# h cannot be reached by advancing time, so it is resolved here instead;
# otherwise the selector could switch anchors and immediately switch back
# without time advancing.
_omrf_position_lookahead(t) = 16 * eps(max(1.0, abs(Float64(t))))

function _omrf_outgoing_position_anchor_index(manager::OMRFAnchorManager,
        x, v, h::Float64=0.0)
    manager.position_valid ||
        _omrf_prepare_position_metric!(manager, manager.position_metric)
    minimum_score = Inf
    @inbounds for idx in eachindex(manager.entries)
        minimum_score = min(minimum_score, _omrf_position_score(manager, x, idx) +
            h * _omrf_position_score_rate(manager, v, idx))
    end
    best = firstindex(manager.entries)
    best_slope = Inf
    @inbounds for idx in eachindex(manager.entries)
        score = _omrf_position_score(manager, x, idx) +
            h * _omrf_position_score_rate(manager, v, idx)
        tolerance = 64eps(Float64) * max(1.0, abs(score),
            abs(minimum_score))
        abs(score - minimum_score) <= tolerance || continue
        slope = _omrf_position_score_rate(manager, v, idx)
        if slope < best_slope || (slope == best_slope && idx < best)
            best = idx
            best_slope = slope
        end
    end
    return best
end

"""Exact next Voronoi-boundary time for a frozen weighted anchor bank.

For a linear path `x + t v`, each competitor difference is affine.  The
pending target records the outgoing cell at a transverse crossing, allowing a
zero-time computational boundary without moving the physical state.
"""
function _omrf_position_anchor_horizon_impl!(manager::OMRFAnchorManager,
        state, flow, phase::Symbol)
    manager.boundary_generation += UInt64(1)
    manager.pending_boundary_idx = 0
    manager.pending_boundary_time = Inf
    manager.pending_boundary_generation = manager.boundary_generation
    x = state.ξ.x
    v = state.ξ.θ
    current = manager.active_idx
    h = _omrf_position_lookahead(state.t[])
    outgoing = _omrf_outgoing_position_anchor_index(manager, x, v, h)
    if outgoing != current
        manager.pending_boundary_idx = outgoing
        manager.pending_boundary_time = state.t[]
        return 0.0
    end
    current_score = _omrf_position_score(manager, x, current)
    current_rate = _omrf_position_score_rate(manager, v, current)
    earliest = Inf
    target = 0
    @inbounds for idx in eachindex(manager.entries)
        idx == current && continue
        difference = _omrf_position_score(manager, x, idx) - current_score
        slope = _omrf_position_score_rate(manager, v, idx) - current_rate
        if slope < 0.0
            crossing = max(h, difference / -slope)
            if crossing < earliest || (crossing == earliest && idx < target)
                earliest = crossing
                target = idx
            end
        end
    end
    manager.pending_boundary_idx = target
    manager.pending_boundary_time = isfinite(earliest) ? state.t[] + earliest : Inf
    return earliest
end

# Score of anchor `anchor_index` along a sticky harmonic path and its first
# three time derivatives (the common y'My term is omitted, see
# `_omrf_position_score`).
@inline function _omrf_harmonic_anchor_jet(manager::OMRFAnchorManager,
        state, flow::PDMPSamplers.AnyBoomerang, anchor_index::Int,
        t::Float64)
    weights = manager.position_weights
    c = manager.position_center
    cosine, sine = cos(t), sin(t)
    value = 0.0
    first = 0.0
    second = 0.0
    third = 0.0
    @inbounds for j in eachindex(c, state.ξ.x, state.ξ.θ, flow.μ)
        weight = weights[j, anchor_index]
        if _omrf_is_free(state, j)
            u = state.ξ.x[j] - flow.μ[j]
            w = state.ξ.θ[j]
            y = flow.μ[j] + u * cosine + w * sine - c[j]
            y1 = -u * sine + w * cosine
            y2 = -u * cosine - w * sine
            value += y * weight
            first += y1 * weight
            second += y2 * weight
            third -= y1 * weight
        else
            value += (state.ξ.x[j] - c[j]) * weight
        end
    end
    return manager.position_offsets[anchor_index] - 2value, -2first,
        -2second, -2third
end

@inline function _omrf_harmonic_jet_precedes(candidate, incumbent,
        candidate_index::Int, incumbent_index::Int)
    @inbounds for order in 1:4
        left, right = candidate[order], incumbent[order]
        tolerance = 128eps(Float64) * max(1.0, abs(left), abs(right))
        left < right - tolerance && return true
        left > right + tolerance && return false
    end
    return candidate_index < incumbent_index
end

function _omrf_harmonic_outgoing_position_anchor_index(
        manager::OMRFAnchorManager, state,
        flow::PDMPSamplers.AnyBoomerang, t::Float64)
    manager.position_valid ||
        _omrf_prepare_position_metric!(manager, manager.position_metric)
    best = firstindex(manager.entries)
    best_jet = _omrf_harmonic_anchor_jet(manager, state, flow, best, t)
    @inbounds for idx in eachindex(manager.entries)
        idx == best && continue
        candidate = _omrf_harmonic_anchor_jet(manager, state, flow, idx, t)
        if _omrf_harmonic_jet_precedes(candidate, best_jet, idx, best)
            best = idx
            best_jet = candidate
        end
    end
    return best
end

# Score difference (competitor minus current) along a harmonic path as
# `constant + cosine*cos(t) + sine*sin(t)`.
@inline function _omrf_harmonic_anchor_difference_coefficients(
        manager::OMRFAnchorManager, state,
        flow::PDMPSamplers.AnyBoomerang, current::Int, competitor::Int)
    weights = manager.position_weights
    c = manager.position_center
    constant = 0.0
    cosine = 0.0
    sine = 0.0
    @inbounds for j in eachindex(c, state.ξ.x, state.ξ.θ, flow.μ)
        delta = weights[j, competitor] - weights[j, current]
        if _omrf_is_free(state, j)
            constant -= 2(flow.μ[j] - c[j]) * delta
            cosine -= 2(state.ξ.x[j] - flow.μ[j]) * delta
            sine -= 2state.ξ.θ[j] * delta
        else
            constant -= 2(state.ξ.x[j] - c[j]) * delta
        end
    end
    constant += manager.position_offsets[competitor] -
        manager.position_offsets[current]
    return constant, cosine, sine
end

@inline function _omrf_positive_periodic_root(base::Float64)
    root = mod(base, 2pi)
    tolerance = 128eps(Float64) * max(1.0, abs(base))
    return root <= tolerance ? root + 2pi : root
end

"""Exact next weighted-Voronoi boundary along a sticky harmonic path.

Each pairwise distance difference is `C + A*cos(t) + B*sin(t)`. Roots are
only deadlines: a root is retained when the deterministic outgoing selector
shows that the minimizer genuinely changes immediately after it.
"""
function _omrf_harmonic_position_anchor_horizon_impl!(
        manager::OMRFAnchorManager, state,
        flow::PDMPSamplers.AnyBoomerang, phase::Symbol)
    manager.boundary_generation += UInt64(1)
    manager.pending_boundary_idx = 0
    manager.pending_boundary_time = Inf
    manager.pending_boundary_generation = manager.boundary_generation
    h = _omrf_position_lookahead(state.t[])
    current = manager.active_idx
    outgoing = _omrf_harmonic_outgoing_position_anchor_index(
        manager, state, flow, 0.0)
    if outgoing != current
        manager.pending_boundary_idx = outgoing
        manager.pending_boundary_time = state.t[]
        return 0.0
    end
    earliest = Inf
    target = 0
    @inbounds for competitor in eachindex(manager.entries)
        competitor == current && continue
        constant, cosine, sine =
            _omrf_harmonic_anchor_difference_coefficients(
                manager, state, flow, current, competitor)
        radius = hypot(cosine, sine)
        tolerance = 128eps(Float64) * max(1.0, abs(constant), radius)
        radius <= tolerance && continue
        abs(constant) > radius + tolerance && continue
        phase_angle = atan(sine, cosine)
        offset = acos(clamp(-constant / radius, -1.0, 1.0))
        for base in (phase_angle - offset, phase_angle + offset)
            crossing = _omrf_positive_periodic_root(base)
            crossing = max(crossing, h)
            crossing > 2pi + tolerance && continue
            selected = _omrf_harmonic_outgoing_position_anchor_index(
                manager, state, flow, crossing)
            selected == current && continue
            if crossing < earliest ||
                    (crossing == earliest && selected < target)
                earliest = crossing
                target = selected
            end
        end
    end
    manager.pending_boundary_idx = target
    manager.pending_boundary_time = isfinite(earliest) ?
        state.t[] + earliest : Inf
    return earliest
end

function _omrf_position_anchor_horizon!(manager::OMRFAnchorManager,
        state, flow, phase::Symbol)
    manager.selection_policy === :position || return Inf
    phase === :main || return Inf
    started = time_ns()
    actual = _omrf_actual_flow(flow)
    horizon = if actual isa Union{BouncyParticle,ZigZag}
        _omrf_position_anchor_horizon_impl!(manager, state, flow, phase)
    elseif actual isa PDMPSamplers.AnyBoomerang
        _omrf_harmonic_position_anchor_horizon_impl!(
            manager, state, actual, phase)
    else
        throw(ArgumentError(
            "position-selected anchor banks require linear or harmonic dynamics"))
    end
    manager.selection_seconds += (time_ns() - started) * 1e-9
    manager.selection_component_evaluations += length(manager.entries)
    return horizon
end


"""Attach observational envelope counters to the authoritative live oracle.

`with_stats` copies the envelope while sharing the prepared oracle.  This
attachment must happen before the first retained proposal (and after every
subsequent wrapper reconstruction), rather than waiting for an anchor change.
It changes counter ownership only; certificate values and caches are untouched.
"""
function _attach_omrf_live_envelope_context!(manager::OMRFAnchorManager,
        cv::SubsampledControlVariate)
    authoritative = manager.oracle.residual_context
    cv.envelope isa OMRFDirectPersonNodeHCVEnvelope || return false
    cv.envelope.context === authoritative && return false
    cv.envelope.context = authoritative
    return true
end

function _select_omrf_position_anchor!(manager::OMRFAnchorManager,
        cv::SubsampledControlVariate, x; phase::Symbol=:unknown,
        state=nothing, boundary_won::Bool=false)
    # `with_stats` copies the envelope but deliberately shares the live
    # oracle.  Reattach the cell certificate's context on first selection so
    # observational counters (and active-set generation) belong to that
    # authoritative oracle even when the bank index itself is unchanged.
    _attach_omrf_live_envelope_context!(manager, cv)
    previous = manager.active_idx
    owns_boundary = boundary_won && phase === :main &&
        manager.pending_boundary_idx > 0 && state !== nothing &&
        manager.pending_boundary_generation == manager.boundary_generation &&
        state.t[] == manager.pending_boundary_time
    started = time_ns()
    selected = if owns_boundary
        idx = manager.pending_boundary_idx
        manager.selection_boundaries += 1
        idx
    else
        _omrf_position_anchor_index(manager, x)
    end
    manager.selection_seconds += (time_ns() - started) * 1e-9
    manager.selection_component_evaluations += length(manager.entries)
    # Every completed step consumes or invalidates the prediction.  The next
    # search derives a fresh deadline from the actual outgoing state.
    manager.pending_boundary_idx = 0
    manager.pending_boundary_time = Inf
    manager.pending_boundary_generation = UInt64(0)
    manager.active_idx = selected
    phase === :main && (manager.main_selections += 1)
    if selected != previous
        PDMPSamplers.refresh_anchor_owned!(cv,
            manager.entries[selected].state.anchor)
        phase === :main && (manager.main_activations += 1)
        return true
    end
    return false
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
        state=nothing, flow=nothing, boundary_won::Bool=false,
        event_type::Symbol=:unknown)
    if manager.selection_policy === :position
        return _select_omrf_position_anchor!(manager, cv, x;
            phase, state, boundary_won)
    end
    previous = manager.active_idx
    best_idx = firstindex(manager.entries)
    best_distance = Inf
    selection_start = time_ns()
    @inbounds for idx in eachindex(manager.entries)
        entry = manager.entries[idx]
        entry.age += 1
        distance = _omrf_anchor_distance2(manager, x, entry.state)
        if distance < best_distance
            best_idx = idx
            best_distance = distance
        end
    end
    manager.selection_seconds += (time_ns() - selection_start) * 1e-9
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
            PDMPSamplers.refresh_anchor_owned!(cv, prepared.anchor)
            manager.main_activations += 1
            manager.main_refreshes += 1
            return true
        end
    end
    if best_idx != previous
        PDMPSamplers.refresh_anchor_owned!(cv,
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
        PDMPSamplers.refresh_anchor_owned!(cv, prepared.anchor)
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

# The slab wrapper inside a statistics wrapper, or nothing (used by the fast
# OMRF deterministic gradient).
@inline _omrf_slab_wrapper(provider) = nothing
@inline _omrf_slab_wrapper(provider::PDMPSamplers.WithStats) =
    _omrf_slab_wrapper(provider.f)
@inline _omrf_slab_wrapper(provider::PDMPSamplers.WithResidualStats) =
    _omrf_slab_wrapper(provider.f)
@inline _omrf_slab_wrapper(provider::StanSubsamplingSlabDeterministic) = provider


function (target::StanSubsamplingSlabDeterministic)(out, x)
    target.oracle(out, x)
    PDMPSamplers.active_slab_potential_gradient!(target.slab, target.all_buf, x,
        target.all_active)
    PDMPSamplers.active_slab_potential_gradient!(target.slab, target.active_buf, x,
        target.active)
    @. out = out - target.all_buf + target.active_buf
    return out
end

# The ordinary OMRF control variate has a position-constant likelihood term,
# so its deterministic HVP is the analytic OMRF prior HVP plus the same
# all-active-to-current-active slab correction used by the gradient.  Keeping
# this composition analytic avoids finite-difference errors in Zig--Zag grid
# endpoint derivatives.
function (target::StanSubsamplingSlabDeterministic{
        <:OMRFAnalyticSubsamplingOracle,<:LogLinearGaussianScaleSlab})(
        out, x, direction)
    started = _OMRF_RESEARCH_DIAGNOSTICS ? time_ns() : UInt64(0)
    oracle = target.oracle
    _omrf_prior_hvp!(out, something(oracle.analytic_prior), x, direction)
    PDMPSamplers.active_prior_hvp!(target.slab, target.all_buf, x,
        direction, target.all_active)
    PDMPSamplers.active_prior_hvp!(target.slab, target.active_buf, x,
        direction, target.active)
    @. out = out - target.all_buf + target.active_buf
    if _OMRF_RESEARCH_DIAGNOSTICS
        target.oracle.residual_context.counts.deterministic_hvp_calls += 1
        target.oracle.residual_context.counts.deterministic_hvp_seconds +=
            (time_ns() - started) * 1e-9
    end
    return out
end

function PDMPSamplers.set_active_set!(target::StanSubsamplingSlabDeterministic,
        free::BitVector)
    indices = PDMPSamplers.beta_indices(target.slab)
    @inbounds for j in eachindex(indices)
        target.active[j] = free[indices[j]]
    end
    # The direct fixed-stratum deterministic certificate is owned by the
    # analytic oracle.  Keep its full-coordinate active mask synchronized
    # with the slab wrapper; otherwise it bounds the all-active prior while
    # the evaluated target subtracts inactive slab contributions.
    PDMPSamplers.set_active_set!(target.oracle, free)
    return nothing
end

function PDMPSamplers.deterministic_rate_cell_bound(
        target::StanSubsamplingSlabDeterministic{
            <:OMRFAnalyticSubsamplingOracle,<:LogLinearGaussianScaleSlab},
        state::PDMPSamplers.AbstractPDMPState,
        flow::PDMPSamplers.ContinuousDynamics, left::Real, right::Real)
    return _omrf_direct_deterministic_cell_bound(
        target.oracle, target.oracle.active_free, state, flow, left, right)
end
PDMPSamplers.has_direct_deterministic_rate_cell_bound(
    target::StanSubsamplingSlabDeterministic{
        <:OMRFAnalyticSubsamplingOracle,<:LogLinearGaussianScaleSlab}) =
    true

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

function _require_omrf_analytic_backend(envelope_spec)
    backend = _haskey(envelope_spec, :backend) ?
        String(_rget(envelope_spec, :backend)) : "analytic"
    backend == "analytic" || throw(ArgumentError(
        "OMRF residual envelopes require the analytic backend"))
    return true
end

function _build_stan_subsampling_model(ctx::StanSubsamplingContext,
        unc_names, subsampling, anchor, slab_prior, can_stick;
        anchor_capacity::Integer=0, manage_during_warmup::Bool=true)
    N = Int(_rget(subsampling, :n_observations))
    m = Int(_rget(subsampling, :subsample_size))
    envelope_spec = _rget(subsampling, :residual_envelope)
    is_omrf = _haskey(envelope_spec, :type) &&
        String(_rget(envelope_spec, :type)) == "omrf"
    is_omrf && _require_omrf_analytic_backend(envelope_spec)
    # Generic envelopes retain the matrix-column contract. Analytic OMRF
    # validates its contribution population from X below, without
    # materializing an unused generic matrix.
    if !is_omrf
        weights = _as_float_matrix(_rget(envelope_spec, :weights))
        size(weights, 2) == N || throw(DimensionMismatch(
            "residual envelope must have N observation columns"))
    end
    d = ctx.full.d
    length(anchor) == d || throw(DimensionMismatch("subsampling anchor has the wrong dimension"))
    oracle = is_omrf ?
        OMRFAnalyticSubsamplingOracle(ctx, envelope_spec, unc_names, anchor) :
        StanSubsamplingOracle(ctx, anchor, N)
    if oracle isa OMRFAnalyticSubsamplingOracle
        residual_context = oracle.residual_context
        people = size(residual_context.X, 1)
        N == people || throw(DimensionMismatch(
            "the OMRF person factorization requires $people likelihood " *
            "contributions; received $N"))
    end
    # Analytic OMRF oracles own their anchor storage. Point the initial
    # envelope at that owned vector as well, so later state recycling cannot
    # mutate a caller-owned initial-position vector.
    envelope_anchor = oracle isa OMRFAnalyticSubsamplingOracle ?
        oracle.anchor : anchor
    envelope = _build_stan_residual_envelope(
        envelope_spec, envelope_anchor, oracle)
    if envelope isa OMRFDirectPersonNodeHCVEnvelope &&
            envelope.joint_group_size >= 16
        m == envelope.joint_group_size || throw(ArgumentError(
            "grouped joint signed bounds require subsample_size equal to the group size"))
        oracle isa OMRFAnalyticSubsamplingOracle &&
            oracle.analytic_prior isa OMRFIndependentAnalyticPrior ||
            throw(ArgumentError(
                "research JS16 requires an independent analytic prior"))
    end
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
    factor_sampling === :size_biased || throw(ArgumentError(
        "unknown OMRF factor sampling design: $factor_sampling"))
    subset_design = PDMPSamplers.UniformSubsamplingDesign()
    deterministic_hvp = if oracle isa OMRFAnalyticSubsamplingOracle &&
            oracle.analytic_prior !== nothing
        if slab_prior === nothing
            nothing
        elseif deterministic isa StanSubsamplingSlabDeterministic &&
                deterministic.slab isa LogLinearGaussianScaleSlab
            deterministic
        else
            nothing
        end
    else
        nothing
    end
    cv = SubsampledControlVariate(deterministic, oracle, envelope, anchor, m;
        deterministic_hvp! = deterministic_hvp,
        refresh_anchor! = refresh_callback, subset_design)
    model = PDMPModel(d, cv)
    if manager === nothing
        return model, ctx, unc_names
    end
    adapter = PDMPSamplers.SubsamplingAnchorBankAdapter(
        (live_grad, state, flow, phase, stats, event_type) -> begin
            phase === :warmup && !manage_during_warmup && return false
            # `with_stats` reconstructs the retained control variate. Anchor
            # selection must therefore update the chain-local live gradient,
            # not the pre-wrapping `cv` captured while the model was built.
            # Updating the captured object leaves the live envelope and its
            # factor certificate at an old anchor while the shared analytic
            # oracle advances to the newly installed bank entry.
            active_cv = live_grad isa PDMPSamplers.SubsampledControlVariate ?
                live_grad : cv
            if phase === :main
                _record_omrf_window!(manager, state, flow, stats; stage=:pre)
            end
            changed = _select_omrf_anchor!(manager, active_cv, state.ξ.x;
                phase, state, flow, boundary_won=event_type === :anchor_selection_boundary,
                event_type)
            if phase === :main
                _record_omrf_window!(manager, state, flow, stats; stage=:post)
            end
            return changed
        end,
        (warmup_grad, trace) -> begin
            active_cv = warmup_grad isa PDMPSamplers.SubsampledControlVariate ?
                warmup_grad : cv
            _add_omrf_anchor!(manager, active_cv, trace)
        end,
        (warmup_grad, state, flow, trace_manager, stats) -> begin
            active_cv = warmup_grad isa PDMPSamplers.SubsampledControlVariate ?
                warmup_grad : cv
            _finish_omrf_warmup!(manager, active_cv, state.ξ.x, trace_manager)
        end,
        (live_grad, state, flow, phase, stats) -> begin
            active_cv = live_grad isa PDMPSamplers.SubsampledControlVariate ?
                live_grad : cv
            # The horizon callback is invoked before event search, including
            # the first retained proposal after `with_stats` or restoration.
            _attach_omrf_live_envelope_context!(manager, active_cv)
            _omrf_position_anchor_horizon!(manager, state, flow, phase)
        end,
        Inf, 0.0, false)
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

# Residual-rate bound of one sampled observation (a person for the direct
# OMRF envelope, whose generic observation bounds are indexed by factor).
_diagnostic_observation_bound(envelope, observation) =
    PDMPSamplers.observation_residual_bound(envelope, observation)
_diagnostic_observation_bound(envelope::OMRFDirectPersonNodeHCVEnvelope,
    person) = _omrf_direct_person_cell_bound(envelope, person)

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
    is_omrf && _require_omrf_analytic_backend(envelope_spec)
    diagnostic_phase in ("warmup", "main") || throw(ArgumentError(
        "diagnostic phase must be warmup or main"))
    oracle = is_omrf ?
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
    # The direct OMRF envelope certifies grid cells; a zero-width cell at the
    # diagnostic state gives its pointwise bounds.
    envelope isa OMRFDirectPersonNodeHCVEnvelope &&
        PDMPSamplers.component_cell_scales!(
            envelope.cell_scales, envelope, state, flow, 0.0, 0.0)
    scales = PDMPSamplers.component_scales!(
        envelope.scales, envelope, state, flow, 0.0)
    envelope_rate = (N / m) * sum(i ->
        _diagnostic_observation_bound(envelope, i), subset_vec)
    scaled_residual = (N / m) * residual
    residual_rate = PDMPSamplers.λ(state, scaled_residual, flow) +
        PDMPSamplers.λ(state, -scaled_residual, flow)

    rate_decomposition = nothing
    if m == 1 && oracle isa OMRFAnalyticSubsamplingOracle
        per_observation_residual_rate = zeros(N)
        per_observation_envelope = zeros(N)
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
                _diagnostic_observation_bound(envelope, observation)
            @. diagnostic_gradient = deterministic + N * diagnostic_residual
            per_observation_stochastic_rate[observation] =
                PDMPSamplers.λ(state, diagnostic_gradient, flow)
        end
        exact_residual_total = sum(per_observation_residual_rate)
        configured_residual_total = sum(per_observation_envelope)
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
        free0=nothing,
        stored_velocity0=nothing,
        c0::Float64=1e-2, grid_n::Int=30, grid_t_max::Float64=2.0,
        curvature_backend::String="finite_difference",
        post_warmup_simplify::Bool=false,
        grid_bound::String="value_quadratic",
        grid_curvature_bound::Union{Nothing,Float64}=nothing,
        lazy_low_tightness_threshold::Float64=0.1,
        lazy_max_low_tightness_rejections::Int=3,
        lazy_max_rejections::Int=0,
        t0::Float64=0.0, T::Float64=10000.0, t_warmup::Float64=0.0,
        warmup_adaptation_interval::Union{Nothing,Float64}=nothing,
        boomerang_variance_floor_fraction::Float64=0.0,
        sticky::Bool=false, can_stick=nothing, model_prior=nothing,
        parameter_prior=nothing, slab_prior=nothing,
        show_progress::Bool=true, n_chains::Int=1, threaded::Bool=false,
        seed::Union{Integer,Nothing}=nothing,
        adaptive_scheme::String="diagonal",
        refresh_rate::Float64=1.0,
        n_anchor_updates::Int=0,
        use_anchor_bank::Bool=false,
        anchor_bank_capacity::Int=8,
        anchor_selection_policy::String="adaptive",
        position_anchor_bank_spec=nothing,
        main_anchor_refresh_distance::Float64=Inf,
        require_no_main_full_data::Bool=false,
        support_boundary_mode::String="error",
        support_boundary_max_bisection_steps::Int=60,
        support_boundary_time_rtol::Float64=1e-8,
        support_boundary_time_atol::Float64=1e-10,
        support_boundary_clip_fraction::Float64=1 - 1e-10,
        support_boundary_max_refresh_attempts::Int=20,
        support_boundary_refresh_probe_time::Float64=1e-4,
        support_boundary_min_safe_time::Float64=1e-12,
        diagnostic_window_width::Float64=0.0,
        streaming_trace_directory::Union{Nothing,String}=nothing,
        streaming_trace_buffer_events::Int=4096)
    bridge_profile_start_ns = time_ns()
    algorithm_type in ("GridThinningStrategy", "ThinningStrategy") ||
        throw(ArgumentError("subsampling custom-Stan sampling requires GridThinningStrategy or ThinningStrategy"))
    x0_vec = _as_float_vector(x0)
    envelope_specification = _rget(subsampling, :residual_envelope)
    hcv_bank_learning = _haskey(envelope_specification,
            :hcv_bank_learning) ?
        Symbol(String(_rget(envelope_specification, :hcv_bank_learning))) :
        :geometric
    hcv_bank_learning === :geometric || throw(ArgumentError(
        "only geometric HCV bank learning is supported"))
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
    anchor_selection_policy in ("adaptive", "position") ||
        throw(ArgumentError(
            "anchor_selection_policy must be \"adaptive\" or \"position\""))
    anchor_selection_policy == "position" && !use_anchor_bank &&
        throw(ArgumentError(
            "position-based anchor selection requires use_anchor_bank=true"))
    anchor_selection_policy == "position" &&
        !(flow_type in ("BouncyParticle", "PreconditionedBPS",
            "ZigZag", "PreconditionedZigZag", "Boomerang",
            "AdaptiveBoomerang")) &&
        throw(ArgumentError(
            "position anchor selection requires linear Zig-Zag/BPS or harmonic Boomerang dynamics"))
    !isnothing(position_anchor_bank_spec) &&
        anchor_selection_policy != "position" && throw(ArgumentError(
            "a position anchor bank specification requires position selection"))
    !isnothing(position_anchor_bank_spec) && t_warmup > 0 &&
        throw(ArgumentError(
            "restored position anchor banks require zero warmup"))
    (isinf(main_anchor_refresh_distance) ||
     isfinite(main_anchor_refresh_distance) && main_anchor_refresh_distance > 0) ||
        throw(ArgumentError(
            "main_anchor_refresh_distance must be positive or Inf"))
    anchor_capacity =
        use_anchor_bank ? anchor_bank_capacity : (n_anchor_updates > 0 ? 1 : 0)
    contexts = prepared.contexts
    length(contexts) == n_chains || throw(DimensionMismatch(
        "prepared subsampling context count must equal n_chains"))
    unc_names = prepared.unc_names
    isfinite(boomerang_variance_floor_fraction) &&
        boomerang_variance_floor_fraction >= 0 || throw(ArgumentError(
            "boomerang_variance_floor_fraction must be finite and nonnegative"))
    for ctx in contexts
        built = _build_stan_subsampling_model(ctx, unc_names,
            subsampling, anchor, slab_prior, can_stick_vec; anchor_capacity,
            manage_during_warmup=true)
        if length(built) == 5
            model, _, _, adapter, manager = built
            adapter.update_dt = n_anchor_updates > 0 ?
                t_warmup / n_anchor_updates : Inf
            adapter.last_update = t0
            manager.main_refresh_distance = main_anchor_refresh_distance
            manager.selection_policy = Symbol(anchor_selection_policy)
            if !isnothing(position_anchor_bank_spec)
                chain_spec = _haskey(position_anchor_bank_spec, :version) ?
                    position_anchor_bank_spec :
                    position_anchor_bank_spec[length(position_anchor_bank_spec) == 1 ?
                        1 : length(anchor_managers) + 1]
                _restore_omrf_position_bank!(manager, model.grad, chain_spec,
                    x0_vec)
            end
            if diagnostic_window_width > 0
                manager.window_diagnostics = OMRFWindowDiagnostics(
                    diagnostic_window_width)
            end
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
    ) for ctx in contexts]
    d = first(models).d
    flow = build_flow(flow_type, _to_precision(_as_flow_cov(flow_cov, d), d),
        _as_flow_mean(flow_mean, d); adaptive_scheme, refresh_rate)
    boom_options = PDMPSamplers.BoomerangAdaptationOptions(
        initial_variance_floor_fraction=boomerang_variance_floor_fraction)
    boom_floor_mask = _boomerang_logscale_floor_mask(unc_names, d)
    alg0 = build_algorithm(algorithm_type; c0, d, grid_n, grid_t_max,
        curvature_backend,
        post_warmup_simplify, grid_bound, grid_curvature_bound,
        lazy_low_tightness_threshold, lazy_max_low_tightness_rejections,
        lazy_max_rejections)
    parameter_prior_vec = isnothing(parameter_prior) ? nothing :
        _as_float_vector(parameter_prior)
    alg = _build_sticky_event_algorithm(alg0, sticky, model_prior,
        parameter_prior_vec, slab_prior, can_stick_vec, flow, unc_names)
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
    bridge_profile_before_sample_ns = time_ns()
    chains = if all(isnothing, anchor_managers)
        adaptation_grad = models[1].grad
        adapter = PDMPSamplers.default_warmup_adapter(
            flow, adaptation_grad, t_warmup, t0;
            warmup_adaptation_interval,
            can_stick=sticky ? can_stick_vec : nothing,
            boomerang_options=boom_options,
            boomerang_variance_floor_mask=boom_floor_mask)
        pdmp_sample(initial, flow, models, alg, t0, T, t_warmup;
            progress=show_progress, threaded, seed,
            adapter,
            warmup_stop=build_warmup_stop(t_warmup),
            support_boundary_options=sbopts,
            statistic_counter=PDMPSamplers.DevelStatisticCounter,
            initial_free=isnothing(free0) ? nothing : _as_bool_vector(free0),
            initial_stored_velocity=isnothing(stored_velocity0) ? nothing :
                _as_float_vector(stored_velocity0),
            trace_storage=isnothing(streaming_trace_directory) ? nothing :
                _bridge_trace_storage(streaming_trace_directory,
                    streaming_trace_buffer_events, T))
    else
        run_chain = function(i)
            seed_i = isnothing(seed) ? nothing : seed + i - 1
            flow_i = deepcopy(flow)
            adaptation_grad = models[i].grad
            dynamics_adapter = PDMPSamplers.default_warmup_adapter(
                flow_i, adaptation_grad, t_warmup, t0;
                warmup_adaptation_interval,
                can_stick=sticky ? can_stick_vec : nothing,
                boomerang_options=boom_options,
                boomerang_variance_floor_mask=boom_floor_mask)
            adapter_i = dynamics_adapter isa PDMPSamplers.NoAdaptation ?
                anchor_adapters[i] : PDMPSamplers.SequenceAdapter(
                    (dynamics_adapter, anchor_adapters[i]))
            chain = pdmp_sample(initial, flow_i, [models[i]], alg,
                t0, T, t_warmup; progress=show_progress && i == 1,
                adapter=adapter_i, seed=seed_i,
                warmup_stop=build_warmup_stop(t_warmup),
                support_boundary_options=sbopts,
                statistic_counter=PDMPSamplers.DevelStatisticCounter,
                initial_free=isnothing(free0) ? nothing : _as_bool_vector(free0),
                initial_stored_velocity=isnothing(stored_velocity0) ? nothing :
                    _as_float_vector(stored_velocity0),
                trace_storage=isnothing(streaming_trace_directory) ? nothing :
                    _bridge_trace_storage(joinpath(
                        streaming_trace_directory,
                        "chain_" * lpad(string(i), 4, '0')),
                        streaming_trace_buffer_events, T))
            return chain.traces[1], chain.stats[1], chain.initial_states[1],
                chain.retained_initial_states[1], chain.terminal_states[1]
        end
        chain_results = if threaded && n_chains > 1
            fetch.([Threads.@spawn run_chain(i) for i in 1:n_chains])
        else
            [run_chain(i) for i in 1:n_chains]
        end
        PDMPChains(getindex.(chain_results, 1), getindex.(chain_results, 2),
            getindex.(chain_results, 3), getindex.(chain_results, 4),
            getindex.(chain_results, 5))
    end
    if diagnostic_window_width > 0
        for (manager, trace, stats) in zip(
                anchor_managers, chains.traces, chains.stats)
            manager === nothing && continue
            # PDMPTrace stores event times in a flat `times` vector.  Keep
            # finalization compatible with the current Julia trace layout;
            # this is reached only when opt-in window diagnostics are active.
            # Window diagnostics cover the requested local retained interval,
            # not merely the last event time.  T is the total horizon passed
            # by the R interface, so remove the warmup portion here.
            final_time = isnan(manager.window_diagnostics.origin) ?
                t0 + max(0.0, T - t_warmup) :
                manager.window_diagnostics.origin + max(0.0, T - t_warmup)
            _finish_omrf_windows!(manager, stats, final_time, flow)
        end
    end
    initial_refresh_rate = try Float64(PDMPSamplers.refresh_rate(flow)) catch; NaN end
    bridge_profile_after_sample_ns = time_ns()
    result = _pack_result(chains; requested_flow=flow_type,
        adaptive_scheme=adaptive_scheme, initial_refresh_rate=initial_refresh_rate,
        warmup_adaptation_interval=warmup_adaptation_interval, t_warmup=t_warmup)
    subsampling_context_counters = [Dict(
        "requested_full_gradient_backend" => "bridgestan",
        "resolved_full_gradient_backend" => "bridgestan",
        "warmup_gradient_provider_type" => "none",
        "warmup_hvp_provider_type" => "none",
        "native_warmup_gradient_provider_type" => "none",
        "native_warmup_hvp_provider_type" => "none",
        "deterministic_hvp_provider_type" => model.hvp === nothing ? "none" :
            (hasproperty(model.hvp, :f) ? string(typeof(model.hvp.f)) :
             string(typeof(model.hvp))),
        "full_gradient_calls" => ctx.counts.full_gradient,
        "prior_gradient_calls" => ctx.counts.prior_gradient,
        "selected_gradient_calls" => ctx.counts.selected_gradient,
        "persons_evaluated" => ctx.counts.persons_evaluated +
            ctx.counts.analytic_persons,
        "selected_persons_evaluated" => ctx.counts.persons_evaluated,
        "analytic_residual_calls" => ctx.counts.analytic_residual,
        "selected_refinement_calls" => ctx.counts.selected_refinement_calls,
        # This counter is incremented only inside the exact selected-residual
        # branch.  It is deliberately separate from the public flag so a
        # request that is short-circuited by the generic bound remains visible.
        "exact_selected_refinement_branch_calls" =>
            (model.grad.residual_oracle isa OMRFAnalyticSubsamplingOracle &&
             model.grad.residual_oracle.exact_selected_refinement) ?
            ctx.counts.selected_refinement_calls : 0,
        "analytic_persons_evaluated" => ctx.counts.analytic_persons,
        "node_conditionals_evaluated" => ctx.counts.analytic_node_conditionals,
        "analytic_residual_seconds" => ctx.counts.analytic_residual_seconds,
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
        "compact_anchor_likelihood_seconds" =>
            ctx.counts.compact_anchor_likelihood_seconds,
        "direct_factor_bound_scans" =>
            ctx.counts.direct_factor_bound_scans,
        "direct_factors_evaluated" =>
            ctx.counts.direct_factors_evaluated,
        "direct_cell_certificate_seconds" =>
            ctx.counts.direct_cell_certificate_seconds,
        "direct_cached_cell_hits" =>
            ctx.counts.direct_cached_cell_hits,
        "direct_selected_factor_exact_evaluations" =>
            ctx.counts.direct_selected_factor_exact_evaluations,
        "direct_path_preparations" =>
            ctx.counts.direct_path_preparations,
        "direct_path_factors_prepared" =>
            ctx.counts.direct_path_factors_prepared,
        "direct_path_preparation_seconds" =>
            ctx.counts.direct_path_preparation_seconds,
        "direct_adjacent_endpoint_reuses" =>
            ctx.counts.direct_adjacent_endpoint_reuses,
        "direct_left_endpoint_evaluations" =>
            ctx.counts.direct_left_endpoint_evaluations,
        "direct_right_endpoint_evaluations" =>
            ctx.counts.direct_right_endpoint_evaluations,
        "direct_residual_table_constructions" =>
            ctx.counts.direct_residual_table_constructions,
        "direct_internal_cell_uses" =>
            ctx.counts.direct_internal_cell_uses,
        "shared_current_pattern_evaluations" =>
            ctx.counts.shared_current_pattern_evaluations,
        "shared_current_cell_evaluations" =>
            ctx.counts.shared_current_cell_evaluations,
        "shared_current_rebuilt_nodes" =>
            ctx.counts.shared_current_rebuilt_nodes,
        "shared_current_rebuild_seconds" =>
            ctx.counts.shared_current_rebuild_seconds,
        "shared_current_dense_fallback_cells" =>
            ctx.counts.shared_current_dense_fallback_cells,
        "joint_cell_width_bins" =>
            copy(ctx.counts.joint_cell_width_bins),
        "direct_path_anchor_changes" =>
            ctx.counts.direct_path_anchor_changes,
        "direct_path_mask_changes" =>
            ctx.counts.direct_path_mask_changes,
        "direct_path_flow_changes" =>
            ctx.counts.direct_path_flow_changes,
        "direct_path_state_changes" =>
            ctx.counts.direct_path_state_changes,
        "direct_path_forced_rebuilds" =>
            ctx.counts.direct_path_forced_rebuilds,
        "deterministic_hvp_calls" => ctx.counts.deterministic_hvp_calls,
        "deterministic_hvp_seconds" => ctx.counts.deterministic_hvp_seconds,
        "deterministic_source_proposals" =>
            ctx.counts.deterministic_source_proposals,
        "residual_source_proposals" => ctx.counts.residual_source_proposals,
        "deterministic_source_acceptances" =>
            ctx.counts.deterministic_source_acceptances,
        "residual_source_acceptances" =>
            ctx.counts.residual_source_acceptances,
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
        "omrf_factorization" => "person",
        "anchor_cache_gradient_calls" => ctx.counts.anchor_cache_gradient,
        "anchor_cache_persons_evaluated" => ctx.counts.anchor_cache_persons,
        "anchor_likelihood_cache_bytes" => _anchor_cache_bytes(
            model.grad.residual_oracle),
        "model_constructions" => ctx.counts.model_constructions,
        "data_constructions" => ctx.counts.data_constructions,
        "native_full_data_scans" => ctx.counts.native_full_data_scans,
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
        "anchor_recycled_preparations" => manager === nothing ? 0 :
            manager.recycled_preparations,
        "anchor_total_preparation_seconds" => manager === nothing ? 0.0 : manager.preparation_seconds,
        "anchor_activations" => manager === nothing ? 0 : manager.activations,
        "anchor_main_activations" => manager === nothing ? 0 : manager.main_activations,
        "anchor_main_refreshes" => manager === nothing ? 0 : manager.main_refreshes,
        "anchor_main_selections" => manager === nothing ? 0 : manager.main_selections,
        "anchor_selections" => manager === nothing ? 0 : manager.main_selections,
        "anchor_switches" => manager === nothing ? 0 : manager.main_activations,
        "anchor_refreshes" => manager === nothing ? 0 : manager.main_refreshes,
        "anchor_selection_seconds" => manager === nothing ? 0.0 :
            manager.selection_seconds - manager.warmup_finish_selection_seconds,
        "anchor_preparation_seconds" => manager === nothing ? 0.0 :
            manager.preparation_seconds - manager.warmup_finish_preparation_seconds,
        "window_diagnostic_seconds" => manager === nothing ||
            manager.window_diagnostics === nothing ? 0.0 :
            manager.window_diagnostics.diagnostic_seconds,
        "retained_native_full_data_scans" => manager === nothing ?
            ctx.counts.native_full_data_scans - snapshot.native_full_data_scans :
            ctx.counts.native_full_data_scans -
                manager.warmup_finish_native_full_data_scans,
        "anchor_selection_component_evaluations" => manager === nothing ? 0 :
            manager.selection_component_evaluations,
        "anchor_mean_main_distance" => manager === nothing ||
            iszero(manager.main_selections) ? 0.0 :
            manager.main_distance_sum / manager.main_selections,
        "anchor_max_main_distance" => manager === nothing ? 0.0 :
            manager.max_main_distance,
        "anchor_bank_entries" => manager === nothing ? 0 : length(manager.entries),
        "anchor_bank_requested_capacity" => manager === nothing ? 0 :
            manager.capacity,
        "anchor_bank_order" => manager === nothing ? Int[] :
            collect(eachindex(manager.entries)),
        "anchor_selection_policy" => manager === nothing ? "none" :
            String(manager.selection_policy),
        "anchor_selection_boundaries" => manager === nothing ? 0 :
            manager.selection_boundaries,
        "anchor_pending_boundary" => manager === nothing ? 0 :
            manager.pending_boundary_idx,
        "anchor_distance_scales" => manager === nothing ? Float64[] :
            copy(manager.distance_scales),
        "active_anchor" => manager === nothing ? Float64[] : copy(
            manager.entries[manager.active_idx].state.anchor),
        "anchor_positions" => manager === nothing ? zeros(0, 0) :
            reduce(hcat, (entry.state.anchor for entry in manager.entries)),
        "anchor_bank_bytes" => manager === nothing ? 0 :
            Base.summarysize((manager.entries, manager.recycled_state)),
        "analytic_hcv" => model.grad.residual_oracle isa OMRFAnalyticSubsamplingOracle &&
            model.grad.residual_oracle.residual_context.use_hcv,
        "resolved_flow" => flow_type,
        "resolved_factorization" => model.grad.residual_oracle isa
            OMRFAnalyticSubsamplingOracle ? "person" : "unavailable",
        "resolved_geometry" =>
            string(nameof(typeof(model.grad.envelope.scales))),
        "envelope_components" => try
            length(model.grad.envelope.totals)
        catch
            0
        end,
        "exact_selected_refinement" =>
            model.grad.residual_oracle isa OMRFAnalyticSubsamplingOracle &&
            model.grad.residual_oracle.exact_selected_refinement,
        "resolved_subsample_size" => ctx.m,
        "resolved_use_hcv" => model.grad.residual_oracle isa
            OMRFAnalyticSubsamplingOracle &&
            model.grad.residual_oracle.residual_context.use_hcv,
        "analytic_prior" => model.grad.residual_oracle isa
            OMRFAnalyticSubsamplingOracle &&
            model.grad.residual_oracle.analytic_prior !== nothing,
    ) for (model, ctx, snapshot, manager) in
        zip(models, contexts, construction_snapshots, anchor_managers)]
    result["subsampling_context_counters"] = subsampling_context_counters
    # Publish phase-local pipeline fields under the stable research names.
    # These are copied from the sampler statistics object, not reconstructed
    # from provider-level call counters.
    phase_aliases = Dict(
        "cell_roof_candidates" => "main_clock_candidates",
        "grid_bound_evaluations" => "main_grid_bound_evaluations",
        "pointwise_screen_evaluations" => "main_pointwise_screen_evaluations",
        "pointwise_screen_passes" => "main_pointwise_screen_passes",
        "subset_draws" => "main_subset_draws",
        "selected_subset_bound_evaluations" => "main_selected_subset_bound_evaluations",
        "selected_subset_bound_passes" => "main_selected_subset_bound_passes",
        "residual_oracle_evaluations" => "main_residual_oracle_evaluations",
        "final_thinning_acceptances" => "main_final_thinning_acceptances",
        "grid_schedule_builds" => "main_grid_schedule_builds",
        "grid_schedule_searches" => "main_grid_schedule_searches",
        "poisson_time_generations" => "main_poisson_time_generations",
        "candidate_state_copies" => "main_candidate_state_copies",
        "candidate_move_forwards" => "main_candidate_move_forwards",
        "candidate_loop_iterations" => "main_candidate_loop_iterations",
        "statistic_recordings" => "main_statistic_recordings",
        "reflection_events" => "main_reflection_events",
        "sticky_events" => "main_sticky_events",
        "cell_roof_construction_seconds" => "main_grid_bound_seconds",
        "pointwise_screen_seconds" => "main_pointwise_screen_seconds",
        "subset_draw_seconds" => "main_subset_draw_seconds",
        "selected_subset_bound_seconds" => "main_selected_subset_refinement_seconds",
        "residual_oracle_seconds" => "main_residual_oracle_seconds",
        "final_thinning_seconds" => "main_final_thinning_seconds",
        "reflection_update_seconds" => "main_reflection_update_seconds",
        "reflection_gradient_seconds" => "main_reflection_gradient_seconds",
        "reflection_covariance_seconds" => "main_reflection_covariance_seconds",
        "reflection_sticky_schedule_seconds" => "main_reflection_sticky_schedule_seconds",
        "reflection_state_validation_seconds" => "main_reflection_state_validation_seconds",
        "sticky_update_seconds" => "main_sticky_update_seconds",
        "grid_schedule_build_seconds" => "main_grid_schedule_build_seconds",
        "grid_schedule_search_seconds" => "main_grid_schedule_search_seconds",
        "poisson_time_generation_seconds" => "main_poisson_time_generation_seconds",
        "candidate_state_copy_seconds" => "main_candidate_state_copy_seconds",
        "candidate_move_forward_seconds" => "main_candidate_move_forward_seconds",
        "candidate_loop_overhead_seconds" => "main_candidate_loop_overhead_seconds",
        "statistic_recording_seconds" => "main_statistic_recording_seconds",
        "warmup_phase_elapsed_seconds" => "warmup_phase_elapsed_seconds",
        "main_phase_elapsed_seconds" => "main_phase_elapsed_seconds")
    for (public_name, phase_name) in phase_aliases
        result["stats"][public_name] = result["stats"][phase_name]
    end
    result["stats"]["anchor_selections"] = Float64[c["anchor_selections"] for c in subsampling_context_counters]
    result["stats"]["anchor_switches"] = Float64[c["anchor_switches"] for c in subsampling_context_counters]
    result["stats"]["anchor_refreshes"] = Float64[c["anchor_refreshes"] for c in subsampling_context_counters]
    result["stats"]["native_full_data_scans"] = Float64[c["retained_native_full_data_scans"] for c in subsampling_context_counters]
    result["stats"]["anchor_selection_seconds"] = Float64[c["anchor_selection_seconds"] for c in subsampling_context_counters]
    result["stats"]["anchor_preparation_seconds"] = Float64[c["anchor_preparation_seconds"] for c in subsampling_context_counters]
    # The schedule-search and candidate-loop timers are inclusive diagnostic
    # timers.  Only the established pipeline timers enter this residual so
    # the unattributed main-loop time is explicit rather than silently zero.
    exclusive_main = result["stats"]["main_grid_bound_seconds"] .+
        result["stats"]["main_pointwise_screen_seconds"] .+
        result["stats"]["main_subset_draw_seconds"] .+
        result["stats"]["main_selected_subset_refinement_seconds"] .+
        result["stats"]["main_residual_oracle_seconds"] .+
        result["stats"]["main_final_thinning_seconds"] .+
        result["stats"]["main_reflection_update_seconds"] .+
        result["stats"]["main_sticky_update_seconds"]
    result["stats"]["other_main_loop_seconds"] =
        result["stats"]["main_phase_elapsed_seconds"] .- exclusive_main
    window_values = [_omrf_window_matrix(manager) for manager in anchor_managers]
    result["window_diagnostics"] = isone(length(window_values)) ?
        first(window_values) : window_values
    require_no_main_full_data && _enforce_no_main_full_data!(
        result["stats"]["main_full_gradient_calls"],
        subsampling_context_counters)
    result["subsampling"] = true
    if get(ENV, "OMRF_PROFILE_EXCLUSIVE", "") == "1"
        result["bridge_profile"] = Dict(
            "julia_construction_seconds" =>
                (bridge_profile_before_sample_ns - bridge_profile_start_ns) / 1e9,
            "julia_sample_seconds" =>
                (bridge_profile_after_sample_ns - bridge_profile_before_sample_ns) / 1e9,
            "julia_pack_and_counters_seconds" =>
                (time_ns() - bridge_profile_after_sample_ns) / 1e9)
    end
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
        if !iszero(full_gradients) || !iszero(anchor_preparations) ||
                !iszero(native_scans)
            push!(violations, "chain $chain: main_full_gradient_calls=" *
                "$full_gradients, anchor_main_preparations=" *
                "$anchor_preparations, main_native_full_data_scans=" *
                "$native_scans")
        end
    end
    isempty(violations) || throw(ErrorException(
        "retained sampling performed forbidden full-data work (" *
        join(violations, "; ") * ")"))
    return nothing
end

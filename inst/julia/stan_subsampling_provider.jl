# Persistent selected-gradient backend for public custom-Stan subsampling sampling.

mutable struct OMRFSignedBPSCallCounts
    cell_preparations::Int
    endpoint_fills::Int
    person_node_masses_filled::Int
    endpoint_reuses::Int
    interaction_field_preparations::Int
    endpoint_range_fills::Int
    point_screens::Int
    screen_copy_checks::Int
    component_draws::Int
    cell_preparation_seconds::Float64
    endpoint_fill_seconds::Float64
    interaction_field_seconds::Float64
    endpoint_range_seconds::Float64
    point_screen_seconds::Float64
    screen_copy_seconds::Float64
    component_draw_seconds::Float64
    cell_preparation_allocated_bytes::UInt
    endpoint_fill_allocated_bytes::UInt
    interaction_field_allocated_bytes::UInt
    endpoint_range_allocated_bytes::UInt
    point_screen_allocated_bytes::UInt
    screen_copy_allocated_bytes::UInt
    component_draw_allocated_bytes::UInt
    harmonic_cell_mass_constructions::Int
    harmonic_factor_masses_filled::Int
    harmonic_all_factor_point_scans::Int
    harmonic_selected_point_bounds::Int
    harmonic_cell_mass_seconds::Float64
    harmonic_selected_point_seconds::Float64
    harmonic_cell_mass_allocated_bytes::UInt
    harmonic_selected_point_allocated_bytes::UInt
end

OMRFSignedBPSCallCounts() = OMRFSignedBPSCallCounts(
    0, 0, 0, 0, 0, 0, 0, 0, 0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    UInt(0), UInt(0), UInt(0), UInt(0), UInt(0), UInt(0), UInt(0),
    0, 0, 0, 0, 0.0, 0.0, UInt(0), UInt(0))

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
    analytic_factors::Int
    analytic_node_conditionals::Int
    compact_person_bounds::Int
    compact_static_constructions::Int
    compact_aggregate_point_calls::Int
    compact_cell_bound_calls::Int
    compact_term_selection_calls::Int
    compact_person_selection_calls::Int
    compact_observation_bound_calls::Int
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
    compact_static_seconds::Float64
    compact_aggregate_point_seconds::Float64
    compact_cell_bound_seconds::Float64
    compact_term_selection_seconds::Float64
    compact_person_selection_seconds::Float64
    compact_observation_bound_seconds::Float64
    compact_anchor_likelihood_seconds::Float64
    signed_bps::OMRFSignedBPSCallCounts
    harmonic_coefficient_preparations::Int
    harmonic_coefficient_cache_hits::Int
    harmonic_coefficient_seconds::Float64
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
    direct_horizon_table_hits::Int
    direct_internal_cell_uses::Int
    deterministic_hvp_calls::Int
    deterministic_hvp_seconds::Float64
    last_mark_residual::Bool
    deterministic_source_proposals::Int
    residual_source_proposals::Int
    deterministic_source_acceptances::Int
    residual_source_acceptances::Int
    categorical_point_evaluations::Int
    categorical_point_displacement_gt_four::Int
    categorical_point_old_bound_sum::Float64
    categorical_point_new_bound_sum::Float64
    categorical_cell_cache_hits::Int
    categorical_cell_cache_misses::Int
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
    # 4:24 are integer counters.
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    # native/analytic timing, point/cell counters and their timings.
    0.0, 0.0, 0, 0.0, 0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    # proposal count, proposal summaries, lazy counters/timing.
    0, 0.0, 0.0, 0.0, 0, 0, 0, 0.0, 0.0,
    # compact timing counters and signed-BPS diagnostics.
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    OMRFSignedBPSCallCounts(), 0, 0, 0.0,
    0, 0, 0.0, 0, 0,
    0, 0, 0.0, 0, 0, 0, 0, 0, 0,
    0, 0.0, false, 0, 0, 0, 0,
    0, 0, 0.0, 0.0, 0, 0,
    0, 0, 0, 0.0, 0,
    zeros(Int, 10), 0, 0, 0, 0, 0)

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
    hcv_damping::Float64
    displacement::Vector{Float64}
    hcv_product::Vector{Float64}
    hcv_nodewise::Bool
    node_likelihood_hessians::Union{Nothing,Array{Float64,3}}
    node_hcv_product::Matrix{Float64}
    node_hcv_direction_product::Matrix{Float64}
    node_hcv_alphas::Vector{Float64}
    node_supports::Vector{Vector{Int}}
    factorization::Symbol
    touched_positions::Vector{Int}
    touched_mask::BitVector
    counts::StanSubsamplingCallCounts
    active_pattern_compression::Bool
    active_free::BitVector
    active_generation::Int
    zigzag_sign_aware::Bool
    signed_coordinate_bounds::Vector{Float64}
    active_rebuilds::Int
    active_rebuild_seconds::Float64
    active_group_count_sum::Int
    active_group_count_min::Int
    active_group_count_max::Int
    direct_person_node_hcv::Bool
    node_box_first_order::Bool
    hcv_cache::OMRFHCVProductCache
end

# Compatibility constructor for Julia-only tests and small pilots that build
# the pre-nodewise context positionally.
function OMRFResidualContext(X::Matrix{Int}, X_float::Matrix{Float64},
        seen::Vector{Int}, threshold_indices::Vector{Int},
        interaction_indices::Vector{Int}, global_scale_index::Int,
        node_scale_indices::Vector{Int}, threshold_starts::Vector{Int},
        incident_edges::Vector{Vector{Int}},
        incident_neighbours::Vector{Vector{Int}},
        anchor_probabilities::Matrix{Float64},
        anchor_expected_scores::Matrix{Float64},
        expected_scores_workspace::Matrix{Float64},
        interaction_matrix::Matrix{Float64}, linear_fields::Matrix{Float64},
        score_crossproducts::Matrix{Float64},
        observed_threshold_counts::Vector{Float64},
        observed_interaction_sums::Vector{Float64},
        probability_buffer::Vector{Float64}, hcv_buffer::Vector{Float64},
        likelihood_hessian::Union{Nothing,Matrix{Float64}}, use_hcv::Bool,
        hcv_active::Bool, hcv_damping::Float64, displacement::Vector{Float64},
        hcv_product::Vector{Float64}, factorization::Symbol,
        touched_positions::Vector{Int}, touched_mask::BitVector,
        counts::StanSubsamplingCallCounts, active_pattern_compression::Bool,
        active_free::BitVector, active_generation::Int, active_rebuilds::Int,
        active_rebuild_seconds::Float64, active_group_count_sum::Int,
        active_group_count_min::Int, active_group_count_max::Int)
    d = length(displacement)
    P = length(seen)
    context = OMRFResidualContext(X, X_float, seen, threshold_indices,
        interaction_indices, global_scale_index, node_scale_indices,
        threshold_starts, incident_edges, incident_neighbours,
        anchor_probabilities, zeros(P, size(X, 1)),
        anchor_expected_scores, expected_scores_workspace,
        interaction_matrix, linear_fields, score_crossproducts,
        observed_threshold_counts, observed_interaction_sums, probability_buffer,
        hcv_buffer, likelihood_hessian, use_hcv, hcv_active, hcv_damping,
        displacement, hcv_product, false, nothing, zeros(d, P), zeros(d, P),
        zeros(P), Vector{Vector{Int}}(), factorization, touched_positions,
        touched_mask, counts, active_pattern_compression, active_free,
        active_generation, false, zeros(d), active_rebuilds,
        active_rebuild_seconds, active_group_count_sum, active_group_count_min,
        active_group_count_max, false, false, OMRFHCVProductCache(d))
    context.node_supports = _omrf_node_supports(context)
    return context
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
    xlo::Vector{Float64}
    xhi::Vector{Float64}
    vlo::Vector{Float64}
    vhi::Vector{Float64}
    glo::Vector{Float64}
    ghi::Vector{Float64}
end

OMRFFullGradientProvider(context::OMRFResidualContext,
        prior::OMRFAnalyticPrior) = begin
    d = length(context.displacement)
    OMRFFullGradientProvider(context, prior, zeros(d), zeros(d), zeros(d),
        zeros(d), zeros(d), zeros(d))
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

@inline function _omrf_softmax_probability_interval(logit_lo, logit_hi,
        category::Int, q::Int)
    target_lo = category == 0 ? 0.0 : logit_lo[category]
    target_hi = category == 0 ? 0.0 : logit_hi[category]
    lower_den = 1.0
    upper_den = 1.0
    for other in 0:q
        other == category && continue
        other_lo = other == 0 ? 0.0 : logit_lo[other]
        other_hi = other == 0 ? 0.0 : logit_hi[other]
        lo_term = exp(other_lo - target_hi)
        hi_term = exp(other_hi - target_lo)
        # A very wide certified cell can make a logit difference larger than
        # log(floatmax). The probability remains finite; only this direct
        # denominator representation overflows. Use a shifted softmax for
        # that rare case instead of aborting a valid trajectory.
        if !isfinite(lo_term) || !isfinite(hi_term)
            lower_differences = Float64[
                _omrf_up((other_ == 0 ? 0.0 : logit_hi[other_]) - target_lo)
                for other_ in 0:q if other_ != category]
            upper_differences = Float64[
                _omrf_down((other_ == 0 ? 0.0 : logit_lo[other_]) - target_hi)
                for other_ in 0:q if other_ != category]
            lower_max = max(0.0, maximum(lower_differences))
            upper_max = max(0.0, maximum(upper_differences))
            lower_num = _omrf_down(exp(-lower_max))
            upper_num = _omrf_up(exp(-upper_max))
            lower_scaled_den = _omrf_up(exp(-lower_max))
            upper_scaled_den = _omrf_down(exp(-upper_max))
            for difference in lower_differences
                lower_scaled_den = _omrf_up(lower_scaled_den +
                    _omrf_up(exp(_omrf_up(difference - lower_max))))
            end
            for difference in upper_differences
                upper_scaled_den = _omrf_down(upper_scaled_den +
                    _omrf_down(exp(_omrf_down(difference - upper_max))))
            end
            lower = lower_num <= 0.0 ? 0.0 :
                _omrf_down(lower_num / lower_scaled_den)
            upper = upper_num <= 0.0 ? nextfloat(0.0) :
                _omrf_up(upper_num / upper_scaled_den)
            return max(0.0, lower), min(1.0, upper)
        end
        lower_den = _omrf_up(lower_den + hi_term)
        upper_den = _omrf_down(upper_den + lo_term)
    end
    return _omrf_down(inv(lower_den)), _omrf_up(inv(upper_den))
end

function _omrf_full_likelihood_intervals!(glo, ghi,
        context::OMRFResidualContext, xlo, xhi)
    N, P = size(context.X)
    logit_lo, logit_hi = context.probability_buffer, context.hcv_buffer
    @inbounds for node in 1:P
        q = context.seen[node] - 1
        start = context.threshold_starts[node]
        for n in 1:N
            field_lo = 0.0; field_hi = 0.0
            for k in eachindex(context.incident_edges[node])
                edge = context.incident_edges[node][k]
                neighbour = context.incident_neighbours[node][k]
                weight = context.X[n, neighbour]
                index = context.interaction_indices[edge]
                field_lo = _omrf_down(field_lo + weight * xlo[index])
                field_hi = _omrf_up(field_hi + weight * xhi[index])
            end
            for u in 1:q
                ti = context.threshold_indices[start + u - 1]
                logit_lo[u] = _omrf_down(xlo[ti] + u * field_lo)
                logit_hi[u] = _omrf_up(xhi[ti] + u * field_hi)
            end
            expected_lo = 0.0; expected_hi = 0.0
            for u in 1:q
                pl, ph = _omrf_softmax_probability_interval(
                    logit_lo, logit_hi, u, q)
                ti = context.threshold_indices[start + u - 1]
                glo[ti] = _omrf_down(glo[ti] + pl)
                ghi[ti] = _omrf_up(ghi[ti] + ph)
                expected_lo = _omrf_down(expected_lo + u * pl)
                expected_hi = _omrf_up(expected_hi + u * ph)
            end
            for k in eachindex(context.incident_edges[node])
                edge = context.incident_edges[node][k]
                neighbour = context.incident_neighbours[node][k]
                weight = context.X[n, neighbour]
                index = context.interaction_indices[edge]
                glo[index] = _omrf_down(glo[index] + weight * expected_lo)
                ghi[index] = _omrf_up(ghi[index] + weight * expected_hi)
            end
        end
        for u in 1:q
            ti = context.threshold_indices[start + u - 1]
            constant = context.observed_threshold_counts[start + u - 1]
            glo[ti] = _omrf_down(glo[ti] - constant)
            ghi[ti] = _omrf_up(ghi[ti] - constant)
        end
    end
    @inbounds for (edge, index) in pairs(context.interaction_indices)
        constant = context.observed_interaction_sums[edge]
        glo[index] = _omrf_down(glo[index] - constant)
        ghi[index] = _omrf_up(ghi[index] - constant)
    end
    context.counts.event_bound_observation_visits += N * P
    return nothing
end

function _omrf_full_deterministic_cell_bound_general(provider::OMRFFullGradientProvider,
        active, state, flow, left::Real, right::Real)
    _omrf_fill_path_intervals!(provider.xlo, provider.xhi,
        provider.vlo, provider.vhi, state, flow, Float64(left), Float64(right))
    _omrf_prior_intervals!(provider.glo, provider.ghi, provider.prior,
        provider.xlo, provider.xhi, active)
    _omrf_full_likelihood_intervals!(provider.glo, provider.ghi,
        provider.context, provider.xlo, provider.xhi)
    return _omrf_interval_rate_bound(nothing, nothing,
        provider.glo, provider.ghi, provider.xlo, provider.xhi,
        provider.vlo, provider.vhi, state, flow)
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

"""
Certified full-gradient deterministic roof for a linear OMRF BPS cell.

Write the negative log pseudo-likelihood as `L`.  Every conditional
log-sum-exp is convex and the observed sufficient-statistic numerator is
linear, hence along `x(t)=x+t*v`,

    d/dt dot(v, grad(L)(x(t))) = dot(v, Hessian(L)(x(t))*v) >= 0.

The signed likelihood rate is therefore maximized at the right endpoint.  We
evaluate that complete endpoint likelihood gradient (including the exact
numerator) once and add a separate interval upper bound for the prior rate.
The final `eps`-scale padding covers ordinary Float64 summation/evaluation
roundoff; it is scale-aware and is not a substitute for either mathematical
bound.  Nonlinear flows and coordinatewise Zig-Zag rates continue to use the
general coordinate interval construction.
"""
function _omrf_bps_full_deterministic_cell_bound(
        provider::OMRFFullGradientProvider, active, state, flow,
        left::Real, right::Real)
    actual = _omrf_actual_flow(flow)
    actual isa BouncyParticle || throw(ArgumentError(
        "convex endpoint OMRF certificate is restricted to linear BPS"))

    _omrf_fill_path_intervals!(provider.xlo, provider.xhi,
        provider.vlo, provider.vhi, state, flow, Float64(left), Float64(right))
    _omrf_prior_intervals!(provider.glo, provider.ghi, provider.prior,
        provider.xlo, provider.xhi, active)
    prior_upper, prior_scale = _omrf_signed_interval_upper(
        provider.glo, provider.ghi, provider.vlo, provider.vhi)

    # Reuse the provider-owned interval buffers after the prior upper bound has
    # been reduced.  This keeps the prewarmed cell-construction path allocation
    # free without adding another d-vector to every provider.
    endpoint = provider.xlo
    likelihood_gradient = provider.glo
    fill!(likelihood_gradient, 0.0)
    t = Float64(right)
    likelihood_rate = 0.0
    likelihood_scale = 0.0
    @inbounds for i in eachindex(endpoint)
        endpoint[i] = _omrf_is_free(state, i) ?
            muladd(t, state.ξ.θ[i], state.ξ.x[i]) : state.ξ.x[i]
    end
    _omrf_add_likelihood_gradient!(likelihood_gradient, provider.context,
        endpoint)
    @inbounds for i in eachindex(likelihood_gradient)
        term = _omrf_is_free(state, i) ?
            state.ξ.θ[i] * likelihood_gradient[i] : 0.0
        likelihood_rate += term
        likelihood_scale += abs(term)
    end
    provider.context.counts.event_bound_observation_visits +=
        size(provider.context.X, 1) * size(provider.context.X, 2)

    signed_upper = likelihood_rate + prior_upper
    numerical_scale = max(1.0, likelihood_scale + prior_scale,
        abs(signed_upper))
    result = max(0.0, signed_upper + 64eps(Float64) * numerical_scale)
    isfinite(result) || throw(OverflowError(
        "OMRF convex-endpoint BPS rate is not finite; subdivide or reject the state"))
    return nextfloat(result)
end

function _omrf_full_deterministic_cell_bound(provider::OMRFFullGradientProvider,
        active, state, flow, left::Real, right::Real)
    if _omrf_actual_flow(flow) isa BouncyParticle
        return _omrf_bps_full_deterministic_cell_bound(
            provider, active, state, flow, left, right)
    end
    return _omrf_full_deterministic_cell_bound_general(
        provider, active, state, flow, left, right)
end

function PDMPSamplers.deterministic_rate_cell_bound(
        provider::OMRFFullGradientProvider,
        state::PDMPSamplers.AbstractPDMPState,
        flow::PDMPSamplers.ContinuousDynamics, left::Real, right::Real)
    return _omrf_full_deterministic_cell_bound(
        provider, nothing, state, flow, left, right)
end
PDMPSamplers.has_direct_deterministic_rate_cell_bound(
    ::OMRFFullGradientProvider) = true

function PDMPSamplers.deterministic_rate_cell_bound(
        target::PDMPSamplers.DependentSlabTarget{PG,PH,S,O},
        state::PDMPSamplers.AbstractPDMPState,
        flow::PDMPSamplers.ContinuousDynamics, left::Real, right::Real
        ) where {PG<:FullGradient{<:OMRFFullGradientProvider},PH,
                 S<:LogLinearGaussianScaleSlab,O}
    return _omrf_full_deterministic_cell_bound(
        target.posterior_grad!.f, target.free, state, flow, left, right)
end
PDMPSamplers.has_direct_deterministic_rate_cell_bound(
    ::PDMPSamplers.DependentSlabTarget{PG,PH,S,O}
    ) where {PG<:FullGradient{<:OMRFFullGradientProvider},PH,
             S<:LogLinearGaussianScaleSlab,O} = true

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
    nodecv = oracle.node_position_cv
    if context.hcv_active
        context.hcv_nodewise && throw(ArgumentError(
            "direct deterministic HCV certificate requires a global HCV"))
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
    likelihood_baseline = nodecv === nothing ? oracle.full_anchor :
        nodecv.likelihood_baseline
    prior_baseline = nodecv === nothing ? oracle.prior_anchor : nothing
    return _omrf_interval_rate_bound(likelihood_baseline, prior_baseline,
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
    context.hcv_active && !context.hcv_nodewise &&
        isinf(context.hcv_damping) || throw(ArgumentError(
        "joint signed bound requires the installed undamped global HCV"))
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

@inline function _omrf_node_hcv_alpha!(context::OMRFResidualContext,
        x, node::Int)
    support = context.node_supports[node]
    displacement2 = 0.0
    @inbounds for index in support
        delta = context.displacement[index]
        displacement2 += delta * delta
    end
    α = isinf(context.hcv_damping) ? 1.0 :
        context.hcv_damping / (context.hcv_damping + displacement2)
    context.node_hcv_alphas[node] = α
    return α
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
    hcv_damping = use_hcv ? Float64(_rget(envelope_spec, :hcv_damping)) : 1.0
    !isnan(hcv_damping) && hcv_damping > 0 || throw(ArgumentError(
        "OMRF HCV damping must be positive"))
    hcv_scope = use_hcv && _haskey(envelope_spec, :hcv_scope) ?
        Symbol(String(_rget(envelope_spec, :hcv_scope))) : :global
    hcv_scope in (:global, :node) || throw(ArgumentError(
        "unknown OMRF HCV scope: $hcv_scope"))
    hcv_nodewise = use_hcv && hcv_scope === :node
    d = length(anchor)
    factorization = _haskey(envelope_spec, :factorization) ?
        Symbol(String(_rget(envelope_spec, :factorization))) : :person
    factorization in (:person, :person_node) || throw(ArgumentError(
        "unknown OMRF factorization: $factorization"))
    active_pattern_compression = _haskey(
        envelope_spec, :active_set_compression) &&
        Bool(_rget(envelope_spec, :active_set_compression))
    zigzag_sign_aware = _haskey(envelope_spec, :zigzag_sign_aware) &&
        Bool(_rget(envelope_spec, :zigzag_sign_aware))
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
        use_hcv && !hcv_after_warmup, hcv_damping,
        zeros(d), zeros(d), hcv_nodewise, nothing,
        zeros(d, P), zeros(d, P), zeros(P), Vector{Vector{Int}}(),
        factorization, sizehint!(Int[], d), falses(d), counts,
        active_pattern_compression, trues(d), 0, zigzag_sign_aware,
        zeros(d), 0, 0.0, 0, 0, 0,
        _haskey(envelope_spec, :bound_type) &&
            String(_rget(envelope_spec, :bound_type)) in
                ("direct_person_node_hcv", "direct_factor_roof_person_hcv",
                 "direct_person_node_cv"),
        _haskey(envelope_spec, :bound_type) &&
            String(_rget(envelope_spec, :bound_type)) in
                ("node_box_first_order", "cell_marked_node_box_first_order"),
        OMRFHCVProductCache(d))
    context.node_supports = _omrf_node_supports(context)
    _cache_omrf_anchor!(context, anchor)
    if use_hcv
        if hcv_nodewise
            context.likelihood_hessian, context.node_likelihood_hessians =
                _omrf_likelihood_hessians(context, anchor)
        else
            context.likelihood_hessian =
                _omrf_likelihood_hessian(context, anchor)
        end
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

function _omrf_likelihood_hessians(context::OMRFResidualContext, anchor)
    d = length(anchor)
    hessian = zeros(d, d)
    N, P = size(context.X)
    node_hessians = zeros(d, d, P)
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
            node_hessians[support[row], support[col], node] += value
        end
    end
    context.counts.native_full_data_scans += 1
    return hessian, node_hessians
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
    α = if context.hcv_active
        if context.hcv_nodewise
            @inbounds for node in eachindex(context.node_supports)
                _omrf_node_hcv_alpha!(context, x, node)
            end
            0.0
        else
            displacement2 = dot(context.displacement, context.displacement)
            isinf(context.hcv_damping) ? 1.0 :
                context.hcv_damping / (context.hcv_damping + displacement2)
        end
    else
        0.0
    end
    if context.factorization === :person
        # A person touches the support of every node; record it once.
        if !isempty(subset)
            @inbounds for node in 1:P
                _record_omrf_touched_support!(context, node)
            end
        end
        @inbounds for n in subset
            1 <= n <= N || throw(BoundsError(context.X, (n, :)))
            for node in 1:P
                node_α = context.hcv_nodewise && context.hcv_active ?
                    context.node_hcv_alphas[node] : α
                _omrf_node_residual!(context, out, x, n, node, node_α)
            end
        end
        context.counts.analytic_persons += length(subset)
        context.counts.analytic_node_conditionals += P * length(subset)
    else
        @inbounds for factor in subset
            1 <= factor <= N * P || throw(BoundsError(1:(N * P), factor))
            node0, n = divrem(factor - 1, N)
            _record_omrf_touched_support!(context, node0 + 1)
            node = node0 + 1
            node_α = context.hcv_nodewise && context.hcv_active ?
                context.node_hcv_alphas[node] : α
            _omrf_node_residual!(context, out, x, n + 1, node, node_α)
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

"""Frozen ordinary first-order node-position control-variate state."""
mutable struct OMRFNodePositionCV
    bank::OMRFNodeLocalFrozenBank
    supports::Vector{Vector{Int}}
    grams::Vector{Matrix{Float64}}
    linear_terms::Vector{Matrix{Float64}}
    constants::Matrix{Float64}
    sources::Vector{Int}
    requested_sources::Vector{Int}
    pending_sources::Vector{Int}
    likelihood_baseline::Vector{Float64}
    metric_variances::Vector{Float64}
    installation_generation::UInt64
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

function _omrf_node_position_metric_variances(flow, d::Int)
    actual = _omrf_actual_flow(flow)
    actual isa BouncyParticle || throw(ArgumentError(
        "node-position selection is available only for linear BPS flows"))
    if flow isa PDMPSamplers.PreconditionedDynamics
        hasproperty(flow.metric, :scale) || throw(ArgumentError(
            "node-position selection requires a diagonal preconditioner"))
        scale = flow.metric.scale
        length(scale) == d || throw(DimensionMismatch(
            "preconditioner dimension does not match the OMRF parameter order"))
        return abs2.(Float64.(scale))
    end
    return ones(d)
end

"""Build the specification's weighted-logit Gram matrices without dense designs."""
function _prepare_omrf_node_position_cv(envelope_spec, unc_names,
        candidates::AbstractMatrix, flow_or_metric_variances,
        counts::StanSubsamplingCallCounts=StanSubsamplingCallCounts())
    bank = _prepare_omrf_node_local_bank(
        envelope_spec, unc_names, candidates, counts)
    context = bank.contexts[1]
    N, P = size(context.X)
    d, K = size(bank.candidates)
    metric_variances = flow_or_metric_variances isa AbstractVector ?
        collect(Float64, flow_or_metric_variances) :
        _omrf_node_position_metric_variances(flow_or_metric_variances, d)
    length(metric_variances) == d || throw(DimensionMismatch(
        "node-position metric variances have the wrong dimension"))
    supports = Vector{Vector{Int}}(undef, P)
    grams = Vector{Matrix{Float64}}(undef, P)
    linear_terms = Vector{Matrix{Float64}}(undef, P)
    constants = zeros(P, K)
    @inbounds for node in 1:P
        q = context.seen[node] - 1
        threshold_start = context.threshold_starts[node]
        threshold_support = q == 0 ? Int[] : collect(
            context.threshold_indices[threshold_start:(threshold_start + q - 1)])
        edge_support = q == 0 ? Int[] : [context.interaction_indices[edge]
            for edge in context.incident_edges[node]]
        support = q == 0 ? Int[] : vcat(threshold_support, edge_support)
        supports[node] = support
        s = length(support)
        Q = zeros(s, s)
        if q > 0
            L = length(edge_support)
            c2 = sum(abs2, 1:q)
            m0 = 0.0
            m1 = zeros(L)
            M2 = zeros(L, L)
            threshold_variance = sum(metric_variances[index]
                for index in threshold_support)
            for person in 1:N
                weight = threshold_variance
                for local_edge in 1:L
                    y = context.X[person,
                        context.incident_neighbours[node][local_edge]]
                    weight += c2 * metric_variances[edge_support[local_edge]] * y^2
                end
                m0 += weight
                for a in 1:L
                    ya = context.X[person,
                        context.incident_neighbours[node][a]]
                    m1[a] += weight * ya
                    for b in 1:L
                        yb = context.X[person,
                            context.incident_neighbours[node][b]]
                        M2[a, b] += weight * ya * yb
                    end
                end
            end
            m0 /= N
            m1 ./= N
            M2 ./= N
            for u in 1:q
                Q[u, u] = m0
                for edge in 1:L
                    value = u * m1[edge]
                    Q[u, q + edge] = value
                    Q[q + edge, u] = value
                end
            end
            for a in 1:L, b in 1:L
                Q[q + a, q + b] = c2 * M2[a, b]
            end
        end
        grams[node] = Q
        terms = zeros(s, K)
        for source in 1:K
            a = view(bank.candidates, support, source)
            mul!(view(terms, :, source), Q, a)
            constants[node, source] = dot(a, view(terms, :, source))
        end
        linear_terms[node] = terms
    end
    sources = ones(Int, P)
    return OMRFNodePositionCV(bank, supports, grams, linear_terms, constants,
        sources, similar(sources), similar(sources), zeros(d),
        metric_variances, UInt64(1))
end

@inline function _omrf_node_position_score(nodecv::OMRFNodePositionCV,
        node::Int, source::Int, x)
    support = nodecv.supports[node]
    score = nodecv.constants[node, source]
    terms = nodecv.linear_terms[node]
    @inbounds for local_index in eachindex(support)
        score -= 2x[support[local_index]] * terms[local_index, source]
    end
    return score
end

function _omrf_node_position_sources!(out, nodecv::OMRFNodePositionCV, x)
    K = size(nodecv.bank.candidates, 2)
    @inbounds for node in eachindex(out)
        isempty(nodecv.supports[node]) && (out[node] = 1; continue)
        best = 1
        best_score = _omrf_node_position_score(nodecv, node, 1, x)
        for source in 2:K
            score = _omrf_node_position_score(nodecv, node, source, x)
            tolerance = 64eps(Float64) * max(1.0, abs(score), abs(best_score))
            if score < best_score - tolerance
                best = source
                best_score = score
            end
        end
        out[node] = best
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
    node_position_cv::Union{Nothing,OMRFNodePositionCV}
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

# Opt-in research sidecar. It is disabled by default and records only scalar
# proposal quantities, so ordinary sampling has one predictable branch and no
# per-candidate allocation.
const _OMRF_CANDIDATE_RECORD_LIMIT = Ref(0)
const _OMRF_CANDIDATE_RECORDS =
    IdDict{OMRFAnalyticSubsamplingOracle,Vector{NTuple{7,Float64}}}()
const _OMRF_CANDIDATE_RECORD_SECONDS = Ref(0.0)
const _OMRF_CANDIDATE_RECORD_CALLS = Ref(0)

function set_omrf_candidate_record_limit!(limit::Integer)
    limit >= 0 || throw(ArgumentError("candidate record limit must be nonnegative"))
    _OMRF_CANDIDATE_RECORD_LIMIT[] = Int(limit)
    empty!(_OMRF_CANDIDATE_RECORDS)
    _OMRF_CANDIDATE_RECORD_SECONDS[] = 0.0
    _OMRF_CANDIDATE_RECORD_CALLS[] = 0
    return nothing
end

omrf_candidate_record_diagnostics() = (
    seconds=_OMRF_CANDIDATE_RECORD_SECONDS[],
    calls=_OMRF_CANDIDATE_RECORD_CALLS[])

import Profile
const _OMRF_PROFILE_REPORT = Ref("")
function set_omrf_main_profile_enabled!(enabled::Bool)
    if enabled
        start_hook = () -> begin
            Profile.clear()
            Profile.init(n=10^7, delay=0.0005)
            Profile.start_timer()
        end
        stop_hook = () -> begin
            Profile.stop_timer()
            io = IOBuffer()
            Profile.print(io; format=:flat, sortedby=:count, mincount=2)
            _OMRF_PROFILE_REPORT[] = String(take!(io))
        end
        PDMPSamplers.set_main_phase_profile_hooks!(start_hook, stop_hook)
    else
        PDMPSamplers.set_main_phase_profile_hooks!(nothing, nothing)
    end
    return nothing
end
omrf_main_profile_report() = _OMRF_PROFILE_REPORT[]

function PDMPSamplers.record_subsampling_candidate!(
        oracle::OMRFAnalyticSubsamplingOracle,
        cell_roof, aggregate_bound, aggregate_passed, subset_probability,
        selected_subset_bound, exact_selected_intensity, accepted)
    limit = _OMRF_CANDIDATE_RECORD_LIMIT[]
    iszero(limit) && return nothing
    started = time_ns()
    _OMRF_CANDIDATE_RECORD_CALLS[] += 1
    rows = get!(() -> NTuple{7,Float64}[], _OMRF_CANDIDATE_RECORDS, oracle)
    if length(rows) >= limit
        _OMRF_CANDIDATE_RECORD_SECONDS[] += (time_ns() - started) * 1e-9
        return nothing
    end
    push!(rows, (Float64(cell_roof), Float64(aggregate_bound),
        Float64(aggregate_passed), Float64(subset_probability),
        Float64(selected_subset_bound), Float64(exact_selected_intensity),
        Float64(accepted)))
    _OMRF_CANDIDATE_RECORD_SECONDS[] += (time_ns() - started) * 1e-9
    return nothing
end

function _omrf_candidate_record_matrix(oracle::OMRFAnalyticSubsamplingOracle)
    rows = get(_OMRF_CANDIDATE_RECORDS, oracle, NTuple{7,Float64}[])
    isempty(rows) && return zeros(Float64, 0, 7)
    matrix = Matrix{Float64}(undef, length(rows), 7)
    @inbounds for i in eachindex(rows), j in 1:7
        matrix[i, j] = rows[i][j]
    end
    return matrix
end

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
        zeros(d), zeros(d), false, residual_context, analytic_prior, nothing,
        nothing,
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

# Backward-compatible constructor for small Julia-only pilots and tests that
# construct the oracle from already prepared contexts.  The experimental
# exact selected-factor refinement is opt-in in that form.
function OMRFAnalyticSubsamplingOracle(ctx::StanSubsamplingContext,
        anchor::Vector{Float64}, full_anchor::Vector{Float64},
        prior_anchor::Vector{Float64}, prior_x::Vector{Float64},
        prior_x_position::Vector{Float64}, prior_x_valid::Bool,
        residual_context::OMRFResidualContext,
        analytic_prior::Union{Nothing,OMRFAnalyticPrior},
        node_local_hcv::Union{Nothing,OMRFNodeLocalHCVBank},
        lazy_residual_subset_bound::Float64)
    return OMRFAnalyticSubsamplingOracle(ctx, anchor, full_anchor,
        prior_anchor, prior_x, prior_x_position, prior_x_valid,
        residual_context, analytic_prior, node_local_hcv, nothing,
        lazy_residual_subset_bound, NaN, false, zeros(length(anchor)), false,
        trues(length(anchor)), 0,
        zeros(length(anchor)), zeros(length(anchor)), zeros(length(anchor)),
        zeros(length(anchor)), zeros(length(anchor)), zeros(length(anchor)))
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
    if context.active_pattern_compression || context.direct_person_node_hcv
        copyto!(context.active_free, free)
        context.active_generation = oracle.active_generation
    end
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
    if oracle.node_position_cv !== nothing
        nodecv = oracle.node_position_cv
        if oracle.analytic_prior === nothing
            _prior_gradient_at_x!(oracle, x)
            @. out = nodecv.likelihood_baseline + oracle.prior_x
        else
            _omrf_prior_gradient!(out, oracle.analytic_prior, x)
            @. out += nodecv.likelihood_baseline
        end
        return out
    end
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
        if context.hcv_nodewise
            node_hessians = something(context.node_likelihood_hessians)
            @inbounds for node in axes(node_hessians, 3)
                α = _omrf_node_hcv_alpha!(context, x, node)
                support = context.node_supports[node]
                for row in support
                    product = 0.0
                    for col in support
                        product += node_hessians[row, col, node] *
                            context.displacement[col]
                    end
                    context.node_hcv_product[row, node] = product
                    out[row] += α * product
                end
            end
        else
            displacement2 = dot(context.displacement, context.displacement)
            α = isinf(context.hcv_damping) ? 1.0 :
                context.hcv_damping / (context.hcv_damping + displacement2)
            _omrf_haffine!(context.hcv_product, context.hcv_cache,
                something(context.likelihood_hessian), context.displacement)
            @. out += α * context.hcv_product
        end
    end
    return out
end

function (oracle::OMRFAnalyticSubsamplingOracle)(out, x, direction)
    context = oracle.residual_context
    context.hcv_nodewise || throw(ArgumentError(
        "OMRF deterministic HVP is only defined for node-local HCV"))
    prior = oracle.analytic_prior
    prior === nothing && throw(ArgumentError(
        "node-local OMRF HCV requires an analytic prior HVP"))
    _omrf_prior_hvp!(out, prior, x, direction)
    @. context.displacement = x - oracle.anchor
    node_hessians = something(context.node_likelihood_hessians)
    @inbounds for node in axes(node_hessians, 3)
        displacement2 = 0.0
        dot_displacement_direction = 0.0
        for index in context.node_supports[node]
            displacement2 += context.displacement[index] *
                context.displacement[index]
            dot_displacement_direction +=
                context.displacement[index] * direction[index]
        end
        denominator = context.hcv_damping + displacement2
        α = isinf(context.hcv_damping) ? 1.0 :
            context.hcv_damping / denominator
        dα = isinf(context.hcv_damping) ? 0.0 :
            -2 * context.hcv_damping * dot_displacement_direction /
            (denominator * denominator)
        support = context.node_supports[node]
        for row in support
            product = 0.0
            direction_product = 0.0
            for col in support
                hessian = node_hessians[row, col, node]
                product += hessian * context.displacement[col]
                direction_product += hessian * direction[col]
            end
            context.node_hcv_product[row, node] = product
            context.node_hcv_direction_product[row, node] = direction_product
            out[row] += α * direction_product + dα * product
        end
    end
    return out
end

function (oracle::OMRFAnalyticSubsamplingOracle)(out, x, subset, anchor)
    if oracle.node_position_cv !== nothing
        started = time_ns()
        context = oracle.residual_context
        @inbounds for index in context.touched_positions
            context.touched_mask[index] = false
        end
        empty!(context.touched_positions)
        N, P = size(context.X)
        if context.factorization === :person
            isempty(subset) || @inbounds for node in 1:P
                _record_omrf_touched_support!(context, node)
            end
            context.counts.analytic_persons += length(subset)
            context.counts.analytic_node_conditionals += P * length(subset)
        else
            @inbounds for factor in subset
                node0, _ = divrem(factor - 1, N)
                _record_omrf_touched_support!(context, node0 + 1)
            end
            context.counts.analytic_factors += length(subset)
            context.counts.analytic_node_conditionals += length(subset)
        end
        nodecv = oracle.node_position_cv
        _omrf_node_local_residual!(nodecv.bank, out, x, subset,
            nodecv.sources; factorization=context.factorization)
        context.counts.analytic_residual += 1
        context.counts.analytic_residual_seconds += (time_ns() - started) * 1e-9
        return out
    end
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
    context.direct_person_node_hcv &&
        (context.counts.direct_selected_factor_exact_evaluations +=
            length(subset))
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
    node_component_starts::Vector{Int}
    node_component_ends::Vector{Int}
    node_edge_starts::Vector{Int}
    node_edge_ends::Vector{Int}
    node_edge_indices::Vector{Int}
    component_edge_offsets::Vector{Int}
    component_edge_covariates::Vector{Int}
    context::OMRFResidualContext
    active_compression::Bool
    active_component_count::Int
    active_generation::Int
    active_all_free::BitVector
    active_effective_full::Bool
    active_groups::Union{Nothing,Matrix{Int}}
    active_members::Union{Nothing,Vector{Vector{Int}}}
    active_totals::Union{Nothing,Vector{Float64}}
    active_hash_keys::Vector{UInt64}
    active_hash_values::Vector{Int}
    active_hash_stamps::Vector{Int}
    active_hash_epoch::Int
    active_representatives::Vector{Int}
    offsets::Vector{Float64}
    cosines::Vector{Float64}
    sines::Vector{Float64}
    coordinate_offsets::Vector{Float64}
    coordinate_cosines::Vector{Float64}
    coordinate_sines::Vector{Float64}
    cell_left_offsets::Vector{Float64}
    cell_right_offsets::Vector{Float64}
    cell_velocities::Vector{Float64}
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
Research-only first-order person-node certificate over deterministic
node-specific covariate boxes. For each node/box, the callback bounds the
range of every categorical logit displacement and directional velocity,
including the zero reference category, and returns `range(a)*range(b)/4`.
The four coordinate workspaces are owned by the envelope and recycled.
"""
mutable struct OMRFNodeBoxFirstOrderScales
    anchor::Vector{Float64}
    component_nodes::Vector{Int}
    minima::Matrix{Float64}
    maxima::Matrix{Float64}
    context::OMRFResidualContext
    displacement_lower::Vector{Float64}
    displacement_upper::Vector{Float64}
    velocity_lower::Vector{Float64}
    velocity_upper::Vector{Float64}
    node_bank::Union{Nothing,OMRFNodeLocalFrozenBank}
    node_sources::Vector{Int}
    node_installation_generation::UInt64
    categorical_residual_cap::Bool
end

"""
Research-only component-first node-box clock for batch-one person-node
Adaptive Boomerang or linear BPS subsampling. Closed-cell component roofs and their marked
mass distribution are cached during grid construction. At a proposal the
callback retrieves that cell without rebuilding the global point geometry;
only the selected person-node factor is evaluated at the proposal time.
"""
mutable struct OMRFCellMarkedNodeBoxFirstOrderScales
    base::OMRFNodeBoxFirstOrderScales
    component_totals::Vector{Float64}
    cell_lefts::Vector{Float64}
    cell_rights::Vector{Float64}
    cell_scales::Matrix{Float64}
    cell_cumulative_masses::Matrix{Float64}
    cell_total_masses::Vector{Float64}
    n_cached_cells::Int
    cached_time::Float64
    active_cell::Int
    point_x::Vector{Float64}
    point_theta::Vector{Float64}
    point_free::BitVector
    point_is_sticky::Bool
    point_mu::Vector{Float64}
    point_harmonic::Bool
    selected_factor_point_bound::Float64
    component_marks::Int
    factor_point_evaluations::Int
    factor_point_seconds::Float64
    cell_cache_hits::Int
    cell_cache_misses::Int
    point_displacement_gt_four::Int
    point_old_bound_sum::Float64
    point_new_bound_sum::Float64
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
    # Per-envelope scratch: the harmonic closed-cell callback must not create
    # a threshold-coefficient vector for every node and every cell.
    threshold_coefficients::Vector{NTuple{3,Float64}}
    coordinate_displacements::Vector{Float64}
    coordinate_velocities::Vector{Float64}
end

OMRFNormLocalScales(anchor::Vector{Float64}, context::OMRFResidualContext) =
    OMRFNormLocalScales(anchor, context,
        Vector{NTuple{3,Float64}}(undef, maximum(context.seen)),
        zeros(length(anchor)), zeros(length(anchor)))

"""
Research-only compact certificate for the undamped person-node Hessian control
variate. For one categorical conditional, with the reference logit included
as zero, the third derivative of log-sum-exp in directions `(z, delta, delta)`
is their mixed categorical third central moment. Hoeffding's covariance bound
gives `abs(D^3 A[z,delta,delta]) <= R(z)R(delta)^2/4`; the integral Taylor
remainder supplies another factor `1/2`. Hence the directional remainder is
at most `R(z)R(delta)^2/8`.

The four component rows multiply the static person weights `r^0,...,r^3`.
Thus storage is `O(NP + P)` and point/cell callbacks never scan people.
"""
struct OMRFCompactHCVNormLocalScales
    anchor::Vector{Float64}
    context::OMRFResidualContext
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

"""Exact first-order categorical CV directional residual."""
function _omrf_categorical_cv_residual(anchor_logits, delta, z)
    q = length(anchor_logits)
    length(delta) == q == length(z) || throw(DimensionMismatch(
        "categorical CV vectors must have equal lengths"))
    anchor_maximum = 0.0
    displaced_maximum = 0.0
    @inbounds for u in 1:q
        anchor_maximum = max(anchor_maximum, anchor_logits[u])
        displaced_maximum = max(displaced_maximum,
            anchor_logits[u] + delta[u])
    end
    anchor_denominator = exp(-anchor_maximum)
    displaced_denominator = exp(-displaced_maximum)
    anchor_mean_z = 0.0
    displaced_mean_z = 0.0
    @inbounds for u in 1:q
        anchor_weight = exp(anchor_logits[u] - anchor_maximum)
        displaced_weight = exp(
            anchor_logits[u] + delta[u] - displaced_maximum)
        anchor_denominator += anchor_weight
        displaced_denominator += displaced_weight
        anchor_mean_z += anchor_weight * z[u]
        displaced_mean_z += displaced_weight * z[u]
    end
    return displaced_mean_z / displaced_denominator -
        anchor_mean_z / anchor_denominator
end

function _omrf_direct_categorical_cv_cell_bound(
        anchor_logits, delta_left, z, width::Real)
    width >= 0 || throw(ArgumentError("cell width must be nonnegative"))
    left = _omrf_categorical_cv_residual(anchor_logits, delta_left, z)
    right_delta = similar(delta_left)
    @inbounds for index in eachindex(right_delta, delta_left, z)
        right_delta[index] = delta_left[index] + width * z[index]
    end
    right = _omrf_categorical_cv_residual(anchor_logits, right_delta, z)
    return nextfloat(max(abs(left), abs(right)))
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
abstract type OMRFDirectHCVSamplingPolicy end
struct OMRFDirectFactorSampling <: OMRFDirectHCVSamplingPolicy end
struct OMRFDirectPersonSampling <: OMRFDirectHCVSamplingPolicy end

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

mutable struct OMRFDirectPersonNodeHCVEnvelope{S<:OMRFDirectHCVSamplingPolicy} <:
        PDMPSamplers.AbstractResidualEnvelope
    sampling_policy::S
    n_persons::Int
    n_nodes::Int
    anchor::Vector{Float64}
    context::OMRFResidualContext
    hcv_remainder::Bool
    event_horizon_scope::Bool
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
    certified_horizon::Float64
    horizon_valid::Bool
    path_valid::Bool
end

function _new_omrf_direct_hcv_envelope(anchor,
        context::OMRFResidualContext; max_cells::Int=128,
        hcv_remainder::Bool=true, event_horizon_scope::Bool=false,
        person_signed_roof::Bool=false,
        signed_only_person_roof::Bool=false,
        joint_group_size::Int=0,
        shared_current::Bool=false)
    max_cells >= 1 || throw(ArgumentError(
        "direct HCV cache must hold at least one cell"))
    N, P = size(context.X)
    n_factors = N * P
    qmax = maximum(context.seen) - 1
    sampling_policy = context.factorization === :person ?
        OMRFDirectPersonSampling() : OMRFDirectFactorSampling()
    person_signed_roof &&
        !(sampling_policy isa OMRFDirectPersonSampling && hcv_remainder) &&
        throw(ArgumentError(
            "signed person roofs require whole-person undamped HCV"))
    signed_only_person_roof &&
        !(person_signed_roof && !event_horizon_scope) &&
        throw(ArgumentError(
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
        !(shared_current && context.use_hcv && !context.hcv_nodewise &&
          isinf(context.hcv_damping)) && throw(ArgumentError(
            "grouped joint signed bounds require shared-current evaluation and global undamped HCV; got shared=$shared_current use_hcv=$(context.use_hcv) nodewise=$(context.hcv_nodewise) damping=$(context.hcv_damping)"))
    joint_group_size >= 16 && N < 2joint_group_size && throw(ArgumentError(
        "grouped joint signed bounds need at least two full groups (N >= 2m)"))
    shared_current && joint_group_size == 0 && throw(ArgumentError(
        "shared-current categorical evaluation requires joint groups"))
    shared_current && any(!=(4), context.seen) && throw(ArgumentError(
        "shared-current categorical evaluation currently requires four categories at every node"))
    shared_current && !(context.hcv_active && !context.hcv_nodewise &&
        isinf(context.hcv_damping)) && throw(ArgumentError(
        "shared-current JS4 requires active global undamped HCV"))
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
        sampling_policy,
        N, P, anchor, context, hcv_remainder, event_horizon_scope,
        person_signed_roof, signed_only_person_roof,
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
        zeros(n_factors), 0.0, false,
        0.0, false, false)
end

PDMPSamplers.n_observations(
    envelope::OMRFDirectPersonNodeHCVEnvelope{OMRFDirectPersonSampling}) =
    envelope.n_persons
PDMPSamplers.n_observations(
    envelope::OMRFDirectPersonNodeHCVEnvelope{OMRFDirectFactorSampling}) =
    envelope.n_persons * envelope.n_nodes
PDMPSamplers._trajectory_scale_anchor(
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
PDMPSamplers._validate_subsampling_grid_envelope(
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
    if context.hcv_active && !context.hcv_nodewise
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
            if envelope.hcv_remainder
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
            else
                envelope.path_anchor_cov_z_delta0[destination] = 0.0
                envelope.path_anchor_var_z[destination] = 0.0
                envelope.path_curvature_coefficient[destination] = 0.0
                envelope.path_velocity_range[destination] = 0.0
                envelope.path_zigzag_velocity_sum[destination] = 0.0
            end
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
    envelope.horizon_valid = false
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
        PDMPSamplers._inc_gradient_purpose!(
            provider.stats, Val(:deterministic_gradient))
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
    if !context.hcv_active || context.hcv_nodewise || !envelope.path_valid ||
            envelope.path_active_generation != context.active_generation
        return cv.deterministic_gradient!(out, state.ξ.x)
    end
    _omrf_deterministic_gradient!(out, oracle,
        something(oracle.analytic_prior), state.ξ.x)
    _omrf_add_path_terms!(out, envelope, state, flow)
    if slab !== nothing
        x = state.ξ.x
        PDMPSamplers.active_prior_grad!(slab.slab, slab.all_buf, x, slab.all_active)
        PDMPSamplers.active_prior_grad!(slab.slab, slab.active_buf, x, slab.active)
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
    context.hcv_active && !context.hcv_nodewise &&
        oracle.exact_selected_refinement ||
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
    if envelope.hcv_remainder
        residual -= envelope.path_anchor_cov_z_delta0[factor] +
            t * envelope.path_anchor_var_z[factor]
    end
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

function _omrf_direct_hcv_fill_cell_original!(
        envelope::OMRFDirectPersonNodeHCVEnvelope, state, flow,
        left::Float64, right::Float64, cell::Int)
    cumulative = 0.0
    destination = 1
    @inbounds for node in 1:envelope.n_nodes
        for person in 1:envelope.n_persons
            bound = _omrf_direct_hcv_factor_cell_bound(
                envelope, state, flow, person, node, left, right)
            envelope.factor_workspace[destination] = bound
            cumulative += bound
            envelope.cell_cumulative[destination, cell] = cumulative
            destination += 1
        end
    end
    envelope.cell_left[cell] = left
    envelope.cell_right[cell] = right
    envelope.cell_mass[cell] = cumulative
    _omrf_direct_finish_person_cell!(envelope, cell, false)
    counts = envelope.context.counts
    counts.direct_factor_bound_scans += 1
    counts.direct_factors_evaluated += envelope.n_persons * envelope.n_nodes
    counts.direct_residual_table_constructions += 1
    return cumulative
end

@inline _omrf_direct_finish_person_cell!(
    envelope::OMRFDirectPersonNodeHCVEnvelope{OMRFDirectFactorSampling},
    cell::Int, signed_available::Bool) = envelope.cell_mass[cell]

function _omrf_direct_finish_person_cell!(
        envelope::OMRFDirectPersonNodeHCVEnvelope{OMRFDirectPersonSampling},
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
                if envelope.hcv_remainder
                    delta_range = max(
                        _omrf_direct_cached_delta_range(
                            envelope, factor, node, left),
                        _omrf_direct_cached_delta_range(
                            envelope, factor, node, right))
                    velocity_range = envelope.path_velocity_range[factor]
                    global_cap = _omrf_categorical_hcv_bound(
                        velocity_range, delta_range)
                    bound = nextfloat(min(interpolation_bound, global_cap))
                else
                    bound = nextfloat(interpolation_bound)
                end
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
    if envelope.event_horizon_scope
        envelope.horizon_valid || throw(ArgumentError(
            "direct CV horizon table requested before horizon preparation"))
        tolerance = 16eps(max(1.0, abs(envelope.certified_horizon)))
        left >= -tolerance && right <= envelope.certified_horizon + tolerance ||
            throw(ArgumentError(
                "direct CV cell lies outside its certified event horizon"))
        cell = 1
        envelope.context.counts.direct_horizon_table_hits += 1
    else
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

function _omrf_direct_prepare_horizon!(
        envelope::OMRFDirectPersonNodeHCVEnvelope, state,
        flow::BouncyParticle, horizon::Float64, flow_token::UInt)
    envelope.event_horizon_scope || return nothing
    isfinite(horizon) && horizon >= 0 || throw(ArgumentError(
        "direct CV event-horizon certificate requires a finite nonnegative horizon"))
    _omrf_direct_path_matches(envelope, state, flow_token) ||
        _omrf_direct_reset_path!(envelope, state, flow, flow_token)
    if envelope.horizon_valid && envelope.certified_horizon == horizon
        envelope.context.counts.direct_horizon_table_hits += 1
        return nothing
    end
    envelope.last_boundary_valid = false
    started = time_ns()
    _omrf_direct_hcv_fill_cell_cached!(
        envelope, 0.0, horizon, 1, false, false)
    envelope.context.counts.direct_cell_certificate_seconds +=
        (time_ns() - started) * 1e-9
    envelope.cell_count = 1
    envelope.active_cell = 1
    envelope.certified_horizon = horizon
    envelope.horizon_valid = true
    return nothing
end

function PDMPSamplers.prepare_residual_horizon!(
        envelope::OMRFDirectPersonNodeHCVEnvelope, state,
        flow::BouncyParticle, horizon)
    return _omrf_direct_prepare_horizon!(
        envelope, state, flow, Float64(horizon),
        _omrf_direct_flow_token(flow))
end

function PDMPSamplers.prepare_residual_horizon!(
        envelope::OMRFDirectPersonNodeHCVEnvelope, state,
        flow::PDMPSamplers.PreconditionedDynamics, horizon)
    if flow.dynamics isa Union{ZigZag,PDMPSamplers.AnyBoomerang}
        envelope.event_horizon_scope && throw(ArgumentError(
            "direct CV event-horizon cells currently support only linear BPS"))
        return nothing
    end
    flow.dynamics isa BouncyParticle || throw(ArgumentError(
        "direct CV cells require a linear BPS or Zig-Zag flow"))
    return _omrf_direct_prepare_horizon!(envelope, state, flow.dynamics,
        Float64(horizon), _omrf_direct_flow_token(flow))
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

PDMPSamplers._joint_signed_group_enabled(
    envelope::OMRFDirectPersonNodeHCVEnvelope) =
    envelope.joint_group_size > 0

function PDMPSamplers._joint_signed_cell_mass!(cv,
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

function PDMPSamplers._joint_signed_active_mass!(
        envelope::OMRFDirectPersonNodeHCVEnvelope, state, t)
    _omrf_direct_activate_time!(envelope, state, t)
    return @inbounds envelope.joint_cell_mass[envelope.active_cell]
end

function PDMPSamplers._joint_signed_draw_mark!(rng::Random.AbstractRNG,
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

function PDMPSamplers._draw_component(rng::Random.AbstractRNG,
        envelope::OMRFDirectPersonNodeHCVEnvelope{OMRFDirectFactorSampling},
        B::Real)::Int
    cell = envelope.active_cell
    iszero(cell) && throw(ArgumentError(
        "direct HCV distinguished factor requested before cell activation"))
    target = rand(rng) * B
    left = 1
    # The direct table always contains one entry per person-node factor.  A
    # person-gradient route still draws its distinguished *factor* from this
    # table and maps it to the owning person below.
    right = envelope.n_persons * envelope.n_nodes
    @inbounds while left < right
        middle = (left + right) >>> 1
        if envelope.cell_cumulative[middle, cell] < target
            left = middle + 1
        else
            right = middle
        end
    end
    return left
end


function PDMPSamplers._draw_component(rng::Random.AbstractRNG,
        envelope::OMRFDirectPersonNodeHCVEnvelope{OMRFDirectPersonSampling},
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

PDMPSamplers._draw_distinguished(rng::Random.AbstractRNG,
    ::OMRFDirectPersonNodeHCVEnvelope, observation::Int) = observation

@inline function PDMPSamplers.begin_subsampling_mark!(
        oracle::OMRFAnalyticSubsamplingOracle)
    oracle.cached_residual_valid = false
    oracle.cached_deterministic_signed_rate = NaN
    return nothing
end

# Failure-only diagnostics for the direct person-node envelope.  This method
# is reached only when the generic sampler's opt-in immutable failure capture
# is enabled; ordinary proposals retain the allocation-free empty hook.
function PDMPSamplers.subsampling_failure_context(
        oracle::OMRFAnalyticSubsamplingOracle,
        envelope::OMRFDirectPersonNodeHCVEnvelope, subset)
    N = envelope.n_persons
    factors = if envelope.context.factorization === :person
        [((node - 1) * N + Int(person))
            for person in subset for node in 1:envelope.n_nodes]
    else
        collect(Int, subset)
    end
    times = fill(envelope.active_time, length(factors))
    cached_remainders = Vector{Float64}(undef, length(factors))
    point_bounds = similar(cached_remainders)
    cell_bounds = similar(cached_remainders)
    displacement_ranges = similar(cached_remainders)
    velocity_ranges = similar(cached_remainders)
    qmax = maximum(envelope.context.seen) - 1
    anchor_probabilities = fill(NaN, qmax, length(factors))
    path_deltas = fill(NaN, qmax, length(factors))
    path_directional_coefficients = fill(NaN, qmax, length(factors))
    anchor_expected_scores = Vector{Float64}(undef, length(factors))
    cached_anchor_mean_z = Vector{Float64}(undef, length(factors))
    cached_anchor_cov_z_delta = Vector{Float64}(undef, length(factors))
    @inbounds for slot in eachindex(factors)
        factor = factors[slot]
        node = div(factor - 1, N) + 1
        person = mod(factor - 1, N) + 1
        cached_remainders[slot] = _omrf_direct_cached_remainder_at(
            envelope, factor, person, node,
            envelope.active_time)
        point_bounds[slot] = _omrf_direct_factor_point_bound(
            envelope, factor)
        cell_bounds[slot] = envelope.signed_only_person_roof ? NaN :
            PDMPSamplers.observation_residual_bound(envelope, factor)
        displacement_ranges[slot] = _omrf_direct_cached_delta_range(
            envelope, factor, node, envelope.active_time)
        velocity_ranges[slot] = envelope.path_kind == UInt8(1) ?
            envelope.path_zigzag_velocity_sum[factor] :
            envelope.path_velocity_range[factor]
        start = envelope.context.threshold_starts[node]
        q = envelope.context.seen[node] - 1
        for u in 1:q
            anchor_probabilities[u, slot] =
                envelope.context.anchor_probabilities[start + u - 1, person]
            path_deltas[u, slot] = envelope.path_delta0[u, factor] +
                envelope.active_time * envelope.path_z[u, factor]
            path_directional_coefficients[u, slot] =
                envelope.path_z[u, factor]
        end
        anchor_expected_scores[slot] =
            envelope.context.anchor_expected_scores[node, person]
        cached_anchor_mean_z[slot] = envelope.path_anchor_mean_z[factor]
        cached_anchor_cov_z_delta[slot] =
            envelope.path_anchor_cov_z_delta0[factor] +
            envelope.active_time * envelope.path_anchor_var_z[factor]
    end
    return (selected_factors=factors, selected_times=times,
        selected_cached_directional_remainders=cached_remainders,
        selected_point_bounds=point_bounds,
        selected_cell_bounds=cell_bounds,
        selected_displacement_ranges=displacement_ranges,
        selected_velocity_ranges=velocity_ranges,
        selected_anchor_probabilities=anchor_probabilities,
        selected_path_deltas=path_deltas,
        selected_path_directional_coefficients=path_directional_coefficients,
        selected_anchor_expected_scores=anchor_expected_scores,
        selected_cached_anchor_mean_z=cached_anchor_mean_z,
        selected_cached_anchor_cov_z_delta=cached_anchor_cov_z_delta,
        exact_refinement_residual=copy(oracle.cached_residual),
        current_touched_positions=copy(oracle.residual_context.touched_positions),
        oracle_envelope_anchor_max_abs_difference=maximum(abs,
            oracle.anchor .- envelope.anchor),
        envelope_path_anchor_max_abs_difference=maximum(abs,
            envelope.path_anchor .- envelope.anchor),
        oracle_anchor_identity=objectid(oracle.anchor),
        direct_envelope_anchor_identity=objectid(envelope.anchor),
        direct_context_identity=objectid(envelope.context),
        oracle_context_identity=objectid(oracle.residual_context),
        direct_path_kind=Int(envelope.path_kind),
        direct_active_cell=envelope.active_cell,
        direct_active_time=envelope.active_time,
        direct_hcv_active=oracle.residual_context.hcv_active)
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
            E<:OMRFDirectPersonNodeHCVEnvelope{OMRFDirectPersonSampling}}
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
        cv.subset[1] = PDMPSamplers._draw_component(rng, envelope, B)
    else
        PDMPSamplers._draw_base_subset!(rng, cv, cv.subset_design)
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

function PDMPSamplers.draw_subset!(rng::Random.AbstractRNG,
        cv::PDMPSamplers.SubsampledControlVariate{F,O,E},
        D::Real, B::Real) where {F,O,
            E<:OMRFDirectPersonNodeHCVEnvelope{OMRFDirectFactorSampling}}
    envelope = cv.envelope
    oracle = cv.residual_oracle
    oracle isa PDMPSamplers.WithResidualStats && (oracle = oracle.f)
    N = PDMPSamplers.n_observations(envelope)
    m = cv.m
    total = D + B
    residual_source = !(iszero(B) ||
        (!iszero(D) && rand(rng) * total <= D))
    if residual_source
        component = PDMPSamplers._draw_component(rng, envelope, B)
        distinguished = PDMPSamplers._draw_distinguished(
            rng, envelope, component)
        PDMPSamplers._draw_base_subset_conditional!(
            rng, cv, cv.subset_design, distinguished)
    else
        PDMPSamplers._draw_base_subset!(rng, cv, cv.subset_design)
    end
    PDMPSamplers.record_subsampling_mark_source!(
        cv.residual_oracle, residual_source)

    cell_weight = 0.0
    point_weight = 0.0
    @inbounds for factor in cv.subset
        cell_bound = PDMPSamplers.observation_residual_bound(
            envelope, factor)
        point_bound = _omrf_direct_factor_point_bound(envelope, factor)
        cell_weight += cell_bound
        point_weight += min(cell_bound, point_bound)
    end
    scale = N / m
    oracle isa OMRFAnalyticSubsamplingOracle &&
        (oracle.lazy_residual_subset_bound = scale * point_weight)
    return D + scale * cell_weight
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
mutable struct OMRFCovariateLocalScales
    anchor::Vector{Float64}
    context::OMRFResidualContext
    edge_displacements::Vector{Float64}
    edge_velocities::Vector{Float64}
    coordinate_offsets::Vector{Float64}
    coordinate_cosines::Vector{Float64}
    coordinate_sines::Vector{Float64}
    cache_x::Vector{Float64}
    cache_theta::Vector{Float64}
    cache_mu::Vector{Float64}
    cache_anchor::Vector{Float64}
    cache_free::BitVector
    harmonic_cache_flow_hash::UInt
    harmonic_cache_valid::Bool
end

const _OMRFCompactAliasTable =
    Union{Nothing,PDMPSamplers.AliasTables.AliasTable{UInt64,Int}}

"""Immutable data-only state shared by all compact anchor-bank entries."""
struct OMRFCompactStatic
    column_sums::Vector{Float64}
    cross_products::Matrix{Float64}
    covariate_aliases::Vector{_OMRFCompactAliasTable}
    cross_aliases::Matrix{_OMRFCompactAliasTable}
    term_kind::Matrix{UInt8}
    term_k::Matrix{Int}
    term_l::Matrix{Int}
    term_count::Vector{Int}
end

"""Compact person-node covariate-local certificate.

The expanded covariate certificate has one component for every node-local
polynomial term.  This representation keeps the same certified bound, but
aggregates those terms directly into one bound per node and person.  The
node/person masses are reused for distinguished-factor sampling, so no
27,000-by-N weight matrix or per-proposal alias-table construction is needed.
"""
mutable struct OMRFCompactCovariateLocalScales
    base::OMRFCovariateLocalScales
    static::OMRFCompactStatic
    person_bounds::Matrix{Float64}
    person_cumulative::Matrix{Float64}
    term_cumulative::Matrix{Float64}
    prepare_person_bounds::Bool
    node_A::Vector{Float64}
    node_C::Vector{Float64}
    node_edge_displacements::Matrix{Float64}
    node_edge_velocities::Matrix{Float64}
    signed_field_displacement::Matrix{Float64}
    signed_field_velocity::Matrix{Float64}
    signed_field_cosine::Matrix{Float64}
    signed_field_sine::Matrix{Float64}
    signed_threshold_displacement::Vector{Float64}
    signed_threshold_velocity::Vector{Float64}
    signed_threshold_cosine::Vector{Float64}
    signed_threshold_sine::Vector{Float64}
    signed_interaction_displacement::Vector{Float64}
    signed_interaction_velocity::Vector{Float64}
    signed_interaction_cosine::Vector{Float64}
    signed_interaction_sine::Vector{Float64}
    signed_cache_x::Vector{Float64}
    signed_cache_theta::Vector{Float64}
    signed_cache_anchor::Vector{Float64}
    signed_cache_free::BitVector
    signed_cache_flow_hash::UInt
    signed_cache_generation::UInt
    signed_cache_valid::Bool
end

mutable struct OMRFCompactCovariateLocalEnvelope{S,C} <:
        PDMPSamplers.AbstractResidualEnvelope
    n_persons::Int
    n_nodes::Int
    component_scales!::S
    component_cell_scales!::C
    totals::Vector{Float64}
    scales::Vector{Float64}
    cell_scales::Vector{Float64}
    cumulative_masses::Vector{Float64}
    screen_generation::UInt
    screen_pending::Bool
    screen_x::Vector{Float64}
    screen_theta::Vector{Float64}
    cell_term_left::Matrix{Float64}
    cell_term_right::Matrix{Float64}
    cell_node_left::Vector{Float64}
    cell_node_right::Vector{Float64}
    cell_state_id::UInt
    cell_flow_type::DataType
    cell_left::Float64
    cell_right::Float64
    cell_valid::Bool
    cell_residual_left::Float64
    cell_residual_right::Float64
    chord_screen::Bool
    screen_time::Float64
    signed_bps::Bool
    signed_harmonic::Bool
    signed_cached_endpoint_generation::UInt
    signed_cached_endpoint_time::Float64
    signed_cached_endpoint_valid::Bool
    signed_cell_cumulative_left::Vector{Float64}
    signed_cell_cumulative_right::Vector{Float64}
    signed_screen_cumulative::Vector{Float64}
    signed_cell_left_mass::Float64
    signed_cell_right_mass::Float64
    signed_screen_alpha::Float64
    signed_screen_valid::Bool
end

PDMPSamplers.n_observations(envelope::OMRFCompactCovariateLocalEnvelope) =
    envelope.n_persons * envelope.n_nodes

function _omrf_compact_static(context::OMRFResidualContext)
    started = time_ns()
    X = context.X_float
    N, P = size(X)
    column_sums = zeros(Float64, P)
    cross_products = zeros(Float64, P, P)
    @inbounds for person in 1:N, left in 1:P
        x_left = X[person, left]
        column_sums[left] += x_left
        for right in 1:P
            cross_products[left, right] += x_left * X[person, right]
        end
    end
    covariate_aliases = _OMRFCompactAliasTable[nothing for _ in 1:P]
    @inbounds for column in 1:P
        if column_sums[column] > 0
            covariate_aliases[column] = PDMPSamplers.AliasTables.AliasTable(
                view(X, :, column))
        end
    end
    cross_aliases = Matrix{_OMRFCompactAliasTable}(undef, P, P)
    @inbounds for left in 1:P, right in 1:P
        if cross_products[left, right] > 0
            weights = Vector{Float64}(undef, N)
            for person in 1:N
                weights[person] = X[person, left] * X[person, right]
            end
            cross_aliases[left, right] =
                PDMPSamplers.AliasTables.AliasTable(weights)
        else
            cross_aliases[left, right] = nothing
        end
    end
    max_degree = maximum(length, context.incident_neighbours; init=0)
    max_terms = 1 + 2 * max_degree + max_degree^2
    term_kind = zeros(UInt8, P, max_terms)
    term_k = zeros(Int, P, max_terms)
    term_l = zeros(Int, P, max_terms)
    term_count = zeros(Int, P)
    @inbounds for node in 1:P
        neighbours = context.incident_neighbours[node]
        term = 1
        term_kind[node, term] = 0x01
        term_count[node] = term
        for neighbour in neighbours
            term += 1
            term_kind[node, term] = 0x02
            term_k[node, term] = neighbour
            term_count[node] = term
        end
        for neighbour in neighbours
            term += 1
            term_kind[node, term] = 0x03
            term_l[node, term] = neighbour
            term_count[node] = term
        end
        for left in neighbours, right in neighbours
            term += 1
            term_kind[node, term] = 0x04
            term_k[node, term] = left
            term_l[node, term] = right
            term_count[node] = term
        end
    end
    context.counts.compact_static_constructions += 1
    context.counts.compact_static_seconds += (time_ns() - started) * 1e-9
    return OMRFCompactStatic(column_sums, cross_products,
        covariate_aliases, cross_aliases, term_kind, term_k, term_l,
        term_count)
end

function _new_omrf_compact_envelope(anchor,
        residual_context::OMRFResidualContext,
        static::Union{Nothing,OMRFCompactStatic}=nothing;
        signed_bps::Bool=false, signed_harmonic::Bool=false)
    N, P = size(residual_context.X)
    max_degree = maximum(length, residual_context.incident_edges; init=0)
    base = OMRFCovariateLocalScales(anchor, residual_context,
        zeros(max_degree), zeros(max_degree), zeros(length(anchor)),
        zeros(length(anchor)), zeros(length(anchor)),
        zeros(length(anchor)), zeros(length(anchor)), zeros(length(anchor)),
        zeros(length(anchor)), falses(length(anchor)), UInt(0), false)
    static_state = static === nothing ? _omrf_compact_static(residual_context) : static
    max_terms = size(static_state.term_kind, 2)
    scales = OMRFCompactCovariateLocalScales(base, static_state,
        zeros(Float64, P, N), zeros(Float64, P, N),
        zeros(Float64, P, max_terms), false,
        zeros(Float64, P), zeros(Float64, P),
        zeros(Float64, P, max_degree), zeros(Float64, P, max_degree),
        zeros(Float64, N, P), zeros(Float64, N, P),
        zeros(Float64, N, P), zeros(Float64, N, P),
        zeros(Float64, length(residual_context.threshold_indices)),
        zeros(Float64, length(residual_context.threshold_indices)),
        zeros(Float64, length(residual_context.threshold_indices)),
        zeros(Float64, length(residual_context.threshold_indices)),
        zeros(Float64, length(residual_context.interaction_indices)),
        zeros(Float64, length(residual_context.interaction_indices)),
        zeros(Float64, length(residual_context.interaction_indices)),
        zeros(Float64, length(residual_context.interaction_indices)),
        zeros(Float64, length(anchor)), zeros(Float64, length(anchor)),
        zeros(Float64, length(anchor)), falses(length(anchor)),
        UInt(0), UInt(0), false)
    return OMRFCompactCovariateLocalEnvelope(
        N, P, scales, scales, ones(Float64, P), zeros(Float64, P),
        zeros(Float64, P), zeros(Float64, P), UInt(0), false,
        zeros(Float64, length(anchor)), zeros(Float64, length(anchor)),
        zeros(Float64, P, max_terms), zeros(Float64, P, max_terms),
        zeros(Float64, P), zeros(Float64, P), UInt(0), typeof(nothing),
        0.0, 0.0, false, 0.0, 0.0, false, 0.0, signed_bps,
        signed_harmonic,
        UInt(0), 0.0, false,
        zeros(Float64, N * P), zeros(Float64, N * P),
        zeros(Float64, N * P), 0.0, 0.0, 0.0, false)
end

@inline function _omrf_compact_fill_node!(out, scales,
        node::Int, q::Int, A::Float64, C::Float64,
        edge_displacements, edge_velocities)
    context = scales.base.context
    X = context.X
    neighbours = context.incident_neighbours[node]
    static = scales.static
    scales.node_A[node] = A
    scales.node_C[node] = C
    @inbounds for k in eachindex(neighbours)
        scales.node_edge_displacements[node, k] = edge_displacements[k]
        scales.node_edge_velocities[node, k] = edge_velocities[k]
    end
    if !scales.prepare_person_bounds
        cumulative = 0.0
        term = 1
        cumulative += 0.25 * size(X, 1) * A * C
        scales.term_cumulative[node, term] = cumulative
        @inbounds for k in eachindex(neighbours)
            term += 1
            cumulative += 0.25 * q * C * edge_displacements[k] *
                static.column_sums[neighbours[k]]
            scales.term_cumulative[node, term] = cumulative
        end
        @inbounds for l in eachindex(neighbours)
            term += 1
            cumulative += 0.25 * q * A * edge_velocities[l] *
                static.column_sums[neighbours[l]]
            scales.term_cumulative[node, term] = cumulative
        end
        @inbounds for k in eachindex(neighbours), l in eachindex(neighbours)
            term += 1
            cumulative += 0.25 * q^2 * edge_displacements[k] *
                edge_velocities[l] * static.cross_products[
                    neighbours[k], neighbours[l]]
            scales.term_cumulative[node, term] = cumulative
        end
        out[node] = cumulative
        return nothing
    end
    total = 0.0
    @inbounds for person in axes(X, 1)
        displacement_sum = 0.0
        velocity_sum = 0.0
        for k in eachindex(neighbours)
            x = X[person, neighbours[k]]
            displacement_sum += x * edge_displacements[k]
            velocity_sum += x * edge_velocities[k]
        end
        bound = 0.25 * (A + q * displacement_sum) *
            (C + q * velocity_sum)
        scales.person_bounds[node, person] = bound
        total += bound
        scales.person_cumulative[node, person] = total
    end
    out[node] = total
    return nothing
end

@inline function _omrf_compact_prepare_person_node!(scales,
        node::Int)
    context = scales.base.context
    X = context.X
    neighbours = context.incident_neighbours[node]
    q = context.seen[node] - 1
    context.counts.compact_person_bounds += size(X, 1)
    A = scales.node_A[node]
    C = scales.node_C[node]
    total = 0.0
    @inbounds for person in axes(X, 1)
        displacement_sum = 0.0
        velocity_sum = 0.0
        for k in eachindex(neighbours)
            x = X[person, neighbours[k]]
            displacement_sum += x * scales.node_edge_displacements[node, k]
            velocity_sum += x * scales.node_edge_velocities[node, k]
        end
        bound = 0.25 * (A + q * displacement_sum) *
            (C + q * velocity_sum)
        scales.person_bounds[node, person] = bound
        total += bound
        scales.person_cumulative[node, person] = total
    end
    return nothing
end

@inline function _omrf_compact_person_bound(scales, node::Int,
        person::Int)
    context = scales.base.context
    neighbours = context.incident_neighbours[node]
    q = context.seen[node] - 1
    displacement_sum = 0.0
    velocity_sum = 0.0
    @inbounds for k in eachindex(neighbours)
        x = context.X[person, neighbours[k]]
        displacement_sum += x * scales.node_edge_displacements[node, k]
        velocity_sum += x * scales.node_edge_velocities[node, k]
    end
    return 0.25 * (scales.node_A[node] + q * displacement_sum) *
        (scales.node_C[node] + q * velocity_sum)
end

function PDMPSamplers.prepare_residual_sampling!(
        envelope::OMRFCompactCovariateLocalEnvelope, state, flow, t)
    scales = envelope.component_scales!
    envelope.screen_pending = false
    envelope.chord_screen = false
    scales.prepare_person_bounds = false
    actual_flow = flow isa PDMPSamplers.PreconditionedDynamics ?
        flow.dynamics : flow
    scales(envelope.scales, state, actual_flow, t)
    cumulative = 0.0
    @inbounds for node in axes(envelope.scales, 1)
        cumulative += envelope.scales[node] * envelope.totals[node]
        envelope.cumulative_masses[node] = cumulative
    end
    return nothing
end

@inline function _omrf_signed_factor_bound_coordinates(
        context::OMRFResidualContext, anchor, x, theta,
        person::Int, node::Int)
    q = context.seen[node] - 1
    iszero(q) && return 0.0
    field_displacement = 0.0
    field_velocity = 0.0
    weighted_edge_velocity = 0.0
    @inbounds for k in eachindex(context.incident_edges[node])
        covariate = context.X[person, context.incident_neighbours[node][k]]
        iszero(covariate) && continue
        edge_index = context.interaction_indices[
            context.incident_edges[node][k]]
        displacement = x[edge_index] - anchor[edge_index]
        velocity = theta[edge_index]
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
        displacement = x[index] - anchor[index]
        velocity = theta[index]
        logit = displacement + u * field_displacement
        minimum_logit = min(minimum_logit, logit)
        maximum_logit = max(maximum_logit, logit)
        logit_velocity = velocity + u * field_velocity
        minimum_logit_velocity = min(minimum_logit_velocity, logit_velocity)
        maximum_logit_velocity = max(maximum_logit_velocity, logit_velocity)
        threshold_velocity += abs(velocity)
    end
    return 0.25 * (maximum_logit - minimum_logit) *
        (maximum_logit_velocity - minimum_logit_velocity)
end

@inline function _omrf_signed_fill_interaction_matrix!(context, values)
    matrix = context.interaction_matrix
    fill!(matrix, 0.0)
    edge = 1
    P = size(matrix, 1)
    @inbounds for left in 1:(P - 1), right in (left + 1):P
        value = values[edge]
        matrix[left, right] = value
        matrix[right, left] = value
        edge += 1
    end
    return matrix
end

@inline function _omrf_signed_cache_matches(scales,
        state, flow, anchor, flow_hash::UInt)
    scales.signed_cache_valid || return false
    scales.signed_cache_flow_hash == flow_hash || return false
    @inbounds for index in eachindex(anchor, scales.signed_cache_anchor)
        anchor[index] == scales.signed_cache_anchor[index] || return false
    end
    @inbounds for index in eachindex(state.ξ.x, scales.signed_cache_x)
        state.ξ.x[index] == scales.signed_cache_x[index] || return false
        state.ξ.θ[index] == scales.signed_cache_theta[index] || return false
    end
    if state isa StickyPDMPState
        @inbounds for index in eachindex(state.free, scales.signed_cache_free)
            state.free[index] == scales.signed_cache_free[index] || return false
        end
    end
    return true
end

function _omrf_signed_prepare_flight_cache!(scales, state, flow)
    context = scales.base.context
    anchor = scales.base.anchor
    flow_generation = flow isa PDMPSamplers.PreconditionedDynamics ?
        (hasproperty(flow.metric, :generation) ? flow.metric.generation : UInt(0)) :
        UInt(0)
    flow_hash = hash((typeof(flow), flow_generation, hash(flow)))
    _omrf_signed_cache_matches(scales, state, flow, anchor, flow_hash) &&
        return scales.signed_cache_generation
    started = time_ns()
    allocated = _omrf_allocated_bytes()
    counts = context.counts.signed_bps
    interaction_started = time_ns()
    interaction_allocated = _omrf_allocated_bytes()
    harmonic_flow = flow isa PDMPSamplers.AnyBoomerang
    edge = 1
    @inbounds for left in 1:(size(context.X, 2) - 1), right in (left + 1):size(context.X, 2)
        index = context.interaction_indices[edge]
        if harmonic_flow
            scales.signed_interaction_displacement[edge] =
                flow.μ[index] - anchor[index]
            scales.signed_interaction_cosine[edge] =
                _omrf_coordinate_free(state, index) ? state.ξ.x[index] - flow.μ[index] :
                state.ξ.x[index] - anchor[index]
            scales.signed_interaction_sine[edge] =
                _omrf_coordinate_free(state, index) ? state.ξ.θ[index] : 0.0
        else
            scales.signed_interaction_displacement[edge] =
                state.ξ.x[index] - anchor[index]
            scales.signed_interaction_velocity[edge] =
                _omrf_coordinate_free(state, index) ? state.ξ.θ[index] : 0.0
        end
        edge += 1
    end
    _omrf_signed_fill_interaction_matrix!(
        context, scales.signed_interaction_displacement)
    mul!(scales.signed_field_displacement,
        context.X_float, context.interaction_matrix)
    if harmonic_flow
        _omrf_signed_fill_interaction_matrix!(
            context, scales.signed_interaction_cosine)
        mul!(scales.signed_field_cosine,
            context.X_float, context.interaction_matrix)
        _omrf_signed_fill_interaction_matrix!(
            context, scales.signed_interaction_sine)
        mul!(scales.signed_field_sine,
            context.X_float, context.interaction_matrix)
    else
        _omrf_signed_fill_interaction_matrix!(
            context, scales.signed_interaction_velocity)
        mul!(scales.signed_field_velocity,
            context.X_float, context.interaction_matrix)
    end
    interaction_elapsed = (time_ns() - interaction_started) * 1e-9
    counts.interaction_field_preparations += 1
    counts.interaction_field_seconds += interaction_elapsed
    counts.interaction_field_allocated_bytes +=
        _omrf_allocated_bytes() - interaction_allocated

    @inbounds for local_index in eachindex(context.threshold_indices)
        index = context.threshold_indices[local_index]
        if harmonic_flow
            scales.signed_threshold_displacement[local_index] =
                flow.μ[index] - anchor[index]
            scales.signed_threshold_cosine[local_index] =
                _omrf_coordinate_free(state, index) ? state.ξ.x[index] - flow.μ[index] :
                state.ξ.x[index] - anchor[index]
            scales.signed_threshold_sine[local_index] =
                _omrf_coordinate_free(state, index) ? state.ξ.θ[index] : 0.0
        else
            scales.signed_threshold_displacement[local_index] =
                state.ξ.x[index] - anchor[index]
            scales.signed_threshold_velocity[local_index] =
                _omrf_coordinate_free(state, index) ? state.ξ.θ[index] : 0.0
        end
    end
    copyto!(scales.signed_cache_x, state.ξ.x)
    copyto!(scales.signed_cache_theta, state.ξ.θ)
    copyto!(scales.signed_cache_anchor, anchor)
    if state isa StickyPDMPState
        copyto!(scales.signed_cache_free, state.free)
    else
        fill!(scales.signed_cache_free, true)
    end
    scales.signed_cache_flow_hash = flow_hash
    scales.signed_cache_generation += UInt(1)
    scales.signed_cache_valid = true
    counts.endpoint_range_fills += 1
    counts.endpoint_range_seconds += (time_ns() - started) * 1e-9 -
        interaction_elapsed
    counts.endpoint_range_allocated_bytes +=
        _omrf_allocated_bytes() - allocated
    return scales.signed_cache_generation
end

function _omrf_signed_fill_endpoint!(envelope, state, flow, t, cumulative)
    scales = envelope.component_scales!
    context = scales.base.context
    people, nodes = size(context.X)
    started = time_ns()
    allocated = _omrf_allocated_bytes()
    _omrf_signed_prepare_flight_cache!(scales, state, flow)
    total = 0.0
    @inbounds for node in 1:nodes, person in 1:people
        q = context.seen[node] - 1
        minimum_logit = 0.0
        maximum_logit = 0.0
        minimum_velocity = 0.0
        maximum_velocity = 0.0
        field_displacement = scales.signed_field_displacement[person, node] +
            t * scales.signed_field_velocity[person, node]
        field_velocity = scales.signed_field_velocity[person, node]
        threshold_start = context.threshold_starts[node]
        for u in 1:q
            local_index = threshold_start + u - 1
            displacement = scales.signed_threshold_displacement[local_index] +
                t * scales.signed_threshold_velocity[local_index]
            velocity = scales.signed_threshold_velocity[local_index]
            logit = displacement + u * field_displacement
            logit_velocity = velocity + u * field_velocity
            minimum_logit = min(minimum_logit, logit)
            maximum_logit = max(maximum_logit, logit)
            minimum_velocity = min(minimum_velocity, logit_velocity)
            maximum_velocity = max(maximum_velocity, logit_velocity)
        end
        bound = 0.25 * (maximum_logit - minimum_logit) *
            (maximum_velocity - minimum_velocity)
        total += bound
        cumulative[(node - 1) * people + person] = total
    end
    counts = context.counts.signed_bps
    counts.endpoint_fills += 1
    counts.person_node_masses_filled += people * nodes
    elapsed = (time_ns() - started) * 1e-9
    counts.endpoint_fill_seconds += elapsed
    counts.endpoint_fill_allocated_bytes += _omrf_allocated_bytes() - allocated
    return total
end

"""Fill a signed factor mass vector for a harmonic Boomerang endpoint."""
function _omrf_signed_harmonic_fill_endpoint!(envelope, state, flow, t, cumulative)
    scales = envelope.component_scales!
    context = scales.base.context
    people, nodes = size(context.X)
    counts = context.counts.signed_bps
    started = time_ns()
    allocated = _omrf_allocated_bytes()
    _omrf_signed_prepare_flight_cache!(scales, state, flow)
    sn, cs = sincos(t)
    total = 0.0
    @inbounds for node in 1:nodes, person in 1:people
        q = context.seen[node] - 1
        min_eta = 0.0
        max_eta = 0.0
        min_vel = 0.0
        max_vel = 0.0
        for u in 1:q
            index = context.threshold_indices[
                context.threshold_starts[node] + u - 1]
            local_index = context.threshold_starts[node] + u - 1
            offset = scales.signed_threshold_displacement[local_index] +
                u * scales.signed_field_displacement[person, node]
            cosine = scales.signed_threshold_cosine[local_index] +
                u * scales.signed_field_cosine[person, node]
            sine = scales.signed_threshold_sine[local_index] +
                u * scales.signed_field_sine[person, node]
            eta = offset + cosine * cs + sine * sn
            velocity = -cosine * sn + sine * cs
            min_eta = min(min_eta, eta)
            max_eta = max(max_eta, eta)
            min_vel = min(min_vel, velocity)
            max_vel = max(max_vel, velocity)
        end
        total += 0.25 * (max_eta - min_eta) * (max_vel - min_vel)
        cumulative[(node - 1) * people + person] = total
    end
    counts.harmonic_all_factor_point_scans += people * nodes
    counts.endpoint_fills += 1
    counts.endpoint_fill_seconds += (time_ns() - started) * 1e-9
    counts.endpoint_fill_allocated_bytes += _omrf_allocated_bytes() - allocated
    return total
end

"""Certified harmonic cell mass for the signed Boomerang research route."""
function _omrf_signed_harmonic_fill_cell!(envelope, state, flow,
        left, right, cumulative)
    scales = envelope.component_scales!
    context = scales.base.context
    people, nodes = size(context.X)
    total = 0.0
    midpoint_sine, midpoint_cosine = sincos((left + right) / 2)
    half_width = (right - left) / 2
    @inbounds for node in 1:nodes, person in 1:people
        q = context.seen[node] - 1
        displacement_range = 0.0
        velocity_range = 0.0
        # Include the baseline category u=0, then all threshold categories.
        for u in 1:q, v in 0:(u - 1)
            iu = context.threshold_indices[context.threshold_starts[node] + u - 1]
            offset_u = scales.base.coordinate_offsets[iu]
            cosine_u = scales.base.coordinate_cosines[iu]
            sine_u = scales.base.coordinate_sines[iu]
            cu = (scales.signed_threshold_displacement[context.threshold_starts[node] + u - 1] +
                u * scales.signed_field_displacement[person, node],
                scales.signed_threshold_cosine[context.threshold_starts[node] + u - 1] +
                u * scales.signed_field_cosine[person, node],
                scales.signed_threshold_sine[context.threshold_starts[node] + u - 1] +
                u * scales.signed_field_sine[person, node])
            cv = (0.0, 0.0, 0.0)
            if v > 0
                iv = context.threshold_indices[context.threshold_starts[node] + v - 1]
                cv = (scales.signed_threshold_displacement[context.threshold_starts[node] + v - 1] +
                    v * scales.signed_field_displacement[person, node],
                    scales.signed_threshold_cosine[context.threshold_starts[node] + v - 1] +
                    v * scales.signed_field_cosine[person, node],
                    scales.signed_threshold_sine[context.threshold_starts[node] + v - 1] +
                    v * scales.signed_field_sine[person, node])
            end
            displacement_range = max(displacement_range,
                _omrf_harmonic_abs_cell_bound(cu[1]-cv[1], cu[2]-cv[2], cu[3]-cv[3],
                    midpoint_sine, midpoint_cosine, half_width))
            velocity_range = max(velocity_range,
                _omrf_harmonic_abs_cell_bound(0.0, cu[3]-cv[3], -(cu[2]-cv[2]),
                    midpoint_sine, midpoint_cosine, half_width))
        end
        total += 0.25 * displacement_range * velocity_range
        cumulative[(node - 1) * people + person] = total
    end
    return total
end

function _omrf_signed_prepare_cell!(out, envelope, state, flow, left, right)
    scales = envelope.component_scales!
    counts = scales.base.context.counts.signed_bps
    cell_started = time_ns()
    cell_allocated = _omrf_allocated_bytes()
    generation = _omrf_signed_prepare_flight_cache!(scales, state, flow)
    reuse_left = envelope.signed_cached_endpoint_valid &&
        envelope.signed_cached_endpoint_generation == generation &&
        envelope.signed_cached_endpoint_time == Float64(left)
    harmonic_route = envelope.signed_harmonic &&
        flow isa PDMPSamplers.AnyBoomerang
    left_total = if harmonic_route
        harmonic_started = time_ns()
        harmonic_allocated = _omrf_allocated_bytes()
        total = _omrf_signed_harmonic_fill_cell!(
            envelope, state, flow, left, right,
            envelope.signed_cell_cumulative_left)
        # The harmonic certificate is constant over the closed cell.  Keep
        # one mass distribution and expose it through both endpoint slots for
        # the generic endpoint-mixture interface.
        copyto!(envelope.signed_cell_cumulative_right,
            envelope.signed_cell_cumulative_left)
        counts.harmonic_cell_mass_constructions += 1
        counts.harmonic_factor_masses_filled += envelope.n_persons *
            envelope.n_nodes
        counts.harmonic_cell_mass_seconds +=
            (time_ns() - harmonic_started) * 1e-9
        counts.harmonic_cell_mass_allocated_bytes +=
            _omrf_allocated_bytes() - harmonic_allocated
        total
    elseif reuse_left
        copyto!(envelope.signed_cell_cumulative_left,
            envelope.signed_cell_cumulative_right)
        counts.endpoint_reuses += 1
        envelope.signed_cell_right_mass
    else
        _omrf_signed_fill_endpoint!(
            envelope, state, flow, left, envelope.signed_cell_cumulative_left)
    end
    right_total = harmonic_route ?
        left_total :
        _omrf_signed_fill_endpoint!(
            envelope, state, flow, right, envelope.signed_cell_cumulative_right)
    envelope.signed_cell_left_mass = left_total
    envelope.signed_cell_right_mass = right_total
    envelope.cell_state_id = objectid(state)
    envelope.cell_flow_type = typeof(flow)
    envelope.cell_left = Float64(left)
    envelope.cell_right = Float64(right)
    envelope.cell_valid = true
    envelope.cell_residual_left = left_total
    envelope.cell_residual_right = right_total
    envelope.signed_cached_endpoint_generation = generation
    envelope.signed_cached_endpoint_time = Float64(right)
    envelope.signed_cached_endpoint_valid = true
    counts.cell_preparations += 1
    counts.cell_preparation_seconds += (time_ns() - cell_started) * 1e-9
    counts.cell_preparation_allocated_bytes +=
        _omrf_allocated_bytes() - cell_allocated
    # Keep the generic component vector valid for diagnostics and for code
    # paths that inspect it, although the signed BPS hook samples factors
    # directly from the endpoint cumulative vectors.
    people = envelope.n_persons
    @inbounds for node in 1:envelope.n_nodes
        first = (node - 1) * people + 1
        last = node * people
        out[node] = envelope.signed_cell_cumulative_right[last] -
            (first == 1 ? 0.0 : envelope.signed_cell_cumulative_right[first - 1])
    end
    return PDMPSamplers._validate_component_scales(
        out, envelope.totals, "component cell")
end

function PDMPSamplers.component_cell_scales!(out,
        envelope::OMRFCompactCovariateLocalEnvelope, state,
        flow::Union{ZigZag,BouncyParticle}, left, right)
    if envelope.signed_bps && flow isa BouncyParticle
        return _omrf_signed_prepare_cell!(out, envelope, state, flow, left, right)
    end
    started = time_ns()
    scales = envelope.component_scales!
    scales(out, state, flow, left)
    copyto!(envelope.cell_node_left, out)
    copyto!(envelope.cell_term_left, scales.term_cumulative)
    scales(out, state, flow, right)
    copyto!(envelope.cell_node_right, out)
    copyto!(envelope.cell_term_right, scales.term_cumulative)
    # Retain the established closed-cell calculation for validation and for
    # fallback diagnostics, while the residual-affine hook below exposes the
    # tighter endpoint chord to the actual proposal clock.
    scales(out, state, flow, left, right)
    envelope.cell_state_id = objectid(state)
    envelope.cell_flow_type = typeof(flow)
    envelope.cell_left = Float64(left)
    envelope.cell_right = Float64(right)
    envelope.cell_valid = true
    envelope.cell_residual_left = dot(envelope.cell_node_left,
        envelope.totals)
    envelope.cell_residual_right = dot(envelope.cell_node_right,
        envelope.totals)
    context = scales.base.context
    context.counts.compact_cell_bound_calls += 1
    context.counts.compact_cell_bound_seconds +=
        (time_ns() - started) * 1e-9
    context.counts.envelope_cell_calls += 1
    context.counts.envelope_cell_seconds +=
        (time_ns() - started) * 1e-9
    return PDMPSamplers._validate_component_scales(
        out, envelope.totals, "component cell")
end

function PDMPSamplers.component_cell_scales!(out,
        envelope::OMRFCompactCovariateLocalEnvelope, state,
        flow::PDMPSamplers.AnyBoomerang, left, right)
    if envelope.signed_harmonic
        return _omrf_signed_prepare_cell!(out, envelope, state, flow, left, right)
    end
    scales = envelope.component_scales!
    # Adaptive Boomerang needs the certified harmonic closed-cell callback;
    # the right endpoint alone is only a point bound and is not a cell roof.
    scales(out, state, flow, left, right)
    return PDMPSamplers._validate_component_scales(
        out, envelope.totals, "component cell")
end

function PDMPSamplers.component_cell_scales!(out,
        envelope::OMRFCompactCovariateLocalEnvelope,
        state, flow::PDMPSamplers.PreconditionedDynamics, left, right)
    if envelope.signed_bps && flow.dynamics isa BouncyParticle
        return _omrf_signed_prepare_cell!(out, envelope, state,
            flow, left, right)
    end
    throw(ArgumentError("compact preconditioned OMRF cell bounds require an explicit non-signed route"))
end

function PDMPSamplers.residual_affine_cell_bound(
        envelope::OMRFCompactCovariateLocalEnvelope, state,
        flow::Union{ZigZag,BouncyParticle}, left, right)
    envelope.cell_valid || return nothing
    envelope.cell_state_id == objectid(state) || return nothing
    typeof(flow) === envelope.cell_flow_type || return nothing
    envelope.cell_left == Float64(left) || return nothing
    envelope.cell_right == Float64(right) || return nothing
    width = envelope.cell_right - envelope.cell_left
    width > 0 || return nothing
    if envelope.signed_bps && flow isa BouncyParticle
        return PDMPSamplers.ResidualAffineCell(
            envelope.signed_cell_left_mass,
            (envelope.signed_cell_right_mass -
             envelope.signed_cell_left_mass) / width,
            envelope.signed_cell_cumulative_left,
            envelope.signed_cell_cumulative_right,
            nothing, nothing)
    end
    return PDMPSamplers.ResidualAffineCell(
        envelope.cell_residual_left,
        (envelope.cell_residual_right - envelope.cell_residual_left) / width,
        envelope.cell_node_left, envelope.cell_node_right,
        envelope.cell_term_left, envelope.cell_term_right)
end

function PDMPSamplers.residual_affine_cell_bound(
        envelope::OMRFCompactCovariateLocalEnvelope, state,
        flow::PDMPSamplers.PreconditionedDynamics, left, right)
    envelope.signed_bps && flow.dynamics isa BouncyParticle || return nothing
    envelope.cell_valid || return nothing
    envelope.cell_state_id == objectid(state) || return nothing
    envelope.cell_left == Float64(left) || return nothing
    envelope.cell_right == Float64(right) || return nothing
    width = envelope.cell_right - envelope.cell_left
    width > 0 || return nothing
    return PDMPSamplers.ResidualAffineCell(
        envelope.signed_cell_left_mass,
        (envelope.signed_cell_right_mass - envelope.signed_cell_left_mass) / width,
        envelope.signed_cell_cumulative_left,
        envelope.signed_cell_cumulative_right, nothing, nothing)
end

function PDMPSamplers.residual_affine_cell_bound(
        envelope::OMRFCompactCovariateLocalEnvelope, state,
        flow::PDMPSamplers.AnyBoomerang, left, right)
    envelope.signed_harmonic && envelope.cell_valid || return nothing
    envelope.cell_state_id == objectid(state) || return nothing
    envelope.cell_left == Float64(left) || return nothing
    envelope.cell_right == Float64(right) || return nothing
    width = envelope.cell_right - envelope.cell_left
    width > 0 || return nothing
    return PDMPSamplers.ResidualAffineCell(
        envelope.signed_cell_left_mass,
        (envelope.signed_cell_right_mass -
         envelope.signed_cell_left_mass) / width,
        envelope.signed_cell_cumulative_left,
        envelope.signed_cell_cumulative_right, nothing, nothing)
end

function PDMPSamplers.screening_residual_bound(
        envelope::OMRFCompactCovariateLocalEnvelope, state, flow, t)
    scales = envelope.component_scales!
    scales.prepare_person_bounds = false
    actual_flow = flow isa PDMPSamplers.PreconditionedDynamics ?
        flow.dynamics : flow
    signed_route = envelope.signed_bps && actual_flow isa BouncyParticle
    harmonic_route = envelope.signed_harmonic &&
        actual_flow isa PDMPSamplers.AnyBoomerang
    if harmonic_route
        envelope.cell_valid || throw(ArgumentError(
            "signed harmonic screening requires a prepared closed cell"))
        objectid(state) == envelope.cell_state_id || throw(ArgumentError(
            "signed harmonic screening cell state is stale"))
        t >= envelope.cell_left && t <= envelope.cell_right ||
            throw(ArgumentError("signed harmonic screening time is outside the prepared cell"))
        total = envelope.signed_cell_right_mass
        @inbounds for j in eachindex(envelope.screen_x, state.ξ.x)
            displacement, velocity = _omrf_coordinate_at(
                state, actual_flow, envelope.component_scales!.base.anchor, j, t)
            envelope.screen_x[j] = envelope.component_scales!.base.anchor[j] + displacement
            envelope.screen_theta[j] = velocity
        end
        envelope.signed_screen_valid = true
        envelope.signed_screen_alpha = 0.0
        envelope.screen_generation += UInt(1)
        envelope.screen_pending = true
        envelope.screen_time = Float64(t)
        scales.base.context.counts.signed_bps.point_screens += 1
        return total
    end
    if signed_route
        signed_screen_started = time_ns()
        signed_screen_allocated = _omrf_allocated_bytes()
        same_cell = envelope.cell_valid && objectid(state) == envelope.cell_state_id &&
            typeof(flow) === envelope.cell_flow_type &&
            t >= envelope.cell_left && t <= envelope.cell_right &&
            envelope.cell_right > envelope.cell_left
        if same_cell
            α = (t - envelope.cell_left) /
                (envelope.cell_right - envelope.cell_left)
            envelope.signed_screen_alpha = α
            envelope.signed_cell_left_mass >= 0 ||
                throw(ArgumentError("signed BPS left mass is negative"))
            envelope.signed_cell_right_mass >= 0 ||
                throw(ArgumentError("signed BPS right mass is negative"))
            total = (1.0 - α) * envelope.signed_cell_left_mass +
                α * envelope.signed_cell_right_mass
            envelope.signed_screen_valid = true
        else
            total = _omrf_signed_fill_endpoint!(
                envelope, state, flow, t, envelope.signed_screen_cumulative)
            envelope.signed_cell_left_mass = total
            envelope.signed_cell_right_mass = total
            envelope.signed_screen_alpha = 0.0
            envelope.signed_screen_valid = true
            copyto!(envelope.signed_cell_cumulative_left,
                envelope.signed_screen_cumulative)
            copyto!(envelope.signed_cell_cumulative_right,
                envelope.signed_screen_cumulative)
        end
        copy_started = time_ns()
        copy_allocated = _omrf_allocated_bytes()
        @inbounds for j in eachindex(envelope.screen_x, state.ξ.x)
            displacement, velocity = _omrf_coordinate_at(
                state, actual_flow, scales.base.anchor, j, t)
            envelope.screen_x[j] = scales.base.anchor[j] + displacement
            envelope.screen_theta[j] = velocity
        end
        screen_counts = scales.base.context.counts.signed_bps
        screen_counts.screen_copy_checks += length(envelope.screen_x)
        screen_counts.screen_copy_seconds +=
            (time_ns() - copy_started) * 1e-9
        screen_counts.screen_copy_allocated_bytes +=
            _omrf_allocated_bytes() - copy_allocated
        envelope.chord_screen = false
        envelope.screen_generation += UInt(1)
        envelope.screen_pending = true
        envelope.screen_time = Float64(t)
        screen_counts.point_screens += 1
        screen_counts.point_screen_seconds +=
            (time_ns() - signed_screen_started) * 1e-9
        screen_counts.point_screen_allocated_bytes +=
            _omrf_allocated_bytes() - signed_screen_allocated
        return total
    end
    envelope.signed_screen_valid = false
    chord = actual_flow isa Union{ZigZag,BouncyParticle} &&
        envelope.cell_valid && objectid(state) == envelope.cell_state_id &&
        typeof(actual_flow) === envelope.cell_flow_type &&
        t >= envelope.cell_left && t <= envelope.cell_right &&
        envelope.cell_right > envelope.cell_left
    if chord
        α = (t - envelope.cell_left) /
            (envelope.cell_right - envelope.cell_left)
        @inbounds for node in axes(envelope.scales, 1)
            envelope.scales[node] = (1.0 - α) *
                envelope.cell_node_left[node] + α *
                envelope.cell_node_right[node]
        end
        envelope.chord_screen = true
    else
        envelope.chord_screen = false
        scales(envelope.scales, state, actual_flow, t)
    end
    @inbounds for j in eachindex(envelope.screen_x, state.ξ.x)
        displacement, velocity = _omrf_coordinate_at(
            state, actual_flow, scales.base.anchor, j, t)
        envelope.screen_x[j] = scales.base.anchor[j] + displacement
        envelope.screen_theta[j] = velocity
    end
    cumulative = 0.0
    @inbounds for node in axes(envelope.scales, 1)
        cumulative += envelope.scales[node] * envelope.totals[node]
        envelope.cumulative_masses[node] = cumulative
    end
    envelope.screen_generation += UInt(1)
    envelope.screen_pending = true
    envelope.screen_time = Float64(t)
    return cumulative
end

PDMPSamplers.residual_sampling_generation(
        envelope::OMRFCompactCovariateLocalEnvelope) =
    envelope.screen_generation

function PDMPSamplers.reuse_screened_residual_sampling!(
        envelope::OMRFCompactCovariateLocalEnvelope, state, flow, t, token)
    signed_check = envelope.signed_bps
    check_started = signed_check ? time_ns() : 0
    check_allocated = signed_check ? _omrf_allocated_bytes() : UInt(0)
    envelope.screen_pending || return false
    token == envelope.screen_generation || return false
    @inbounds for j in eachindex(envelope.screen_x, state.ξ.x)
        x = state.ξ.x[j]
        θ = state.ξ.θ[j]
        target_x = envelope.screen_x[j]
        target_θ = envelope.screen_theta[j]
        # The candidate is produced by the same flow integrator.  The
        # tolerance only absorbs harmless reassociation/rounding differences;
        # a different state or flow trajectory falls back to recomputation.
        tol_x = 1e-12 * (1.0 + abs(target_x))
        tol_θ = 1e-12 * (1.0 + abs(target_θ))
        if abs(x - target_x) > tol_x || abs(θ - target_θ) > tol_θ
            envelope.screen_pending = false
            if signed_check
                counts = envelope.component_scales!.base.context.counts.signed_bps
                counts.screen_copy_checks += length(envelope.screen_x)
                counts.screen_copy_seconds +=
                    (time_ns() - check_started) * 1e-9
                counts.screen_copy_allocated_bytes +=
                    _omrf_allocated_bytes() - check_allocated
            end
            return false
        end
    end
    envelope.screen_pending = false
    if signed_check
        counts = envelope.component_scales!.base.context.counts.signed_bps
        counts.screen_copy_checks += length(envelope.screen_x)
        counts.screen_copy_seconds += (time_ns() - check_started) * 1e-9
        counts.screen_copy_allocated_bytes +=
            _omrf_allocated_bytes() - check_allocated
    end
    return true
end

function (scales::OMRFCompactCovariateLocalScales)(out, state, flow, t)
    started = time_ns()
    base = scales.base
    context = base.context
    _cache_omrf_pattern_point_coordinates!(base, state, flow, t)
    offsets = base.coordinate_offsets
    velocities = base.coordinate_cosines
    @inbounds for node in eachindex(context.seen)
        q = context.seen[node] - 1
        minimum_displacement = 0.0
        maximum_displacement = 0.0
        minimum_velocity = 0.0
        maximum_velocity = 0.0
        threshold_start = context.threshold_starts[node]
        for u in 1:q
            index = context.threshold_indices[threshold_start + u - 1]
            displacement = offsets[index]
            velocity = velocities[index]
            minimum_displacement = min(minimum_displacement, displacement)
            maximum_displacement = max(maximum_displacement, displacement)
            minimum_velocity = min(minimum_velocity, velocity)
            maximum_velocity = max(maximum_velocity, velocity)
        end
        edge_displacements = base.edge_displacements
        edge_velocities = base.edge_velocities
        edges = context.incident_edges[node]
        @inbounds for k in eachindex(edges)
            index = context.interaction_indices[edges[k]]
            edge_displacements[k] = abs(offsets[index])
            edge_velocities[k] = abs(velocities[index])
        end
        _omrf_compact_fill_node!(out, scales, node, q,
            maximum_displacement - minimum_displacement,
            maximum_velocity - minimum_velocity,
            edge_displacements, edge_velocities)
    end
    context.counts.compact_aggregate_point_calls += 1
    context.counts.compact_aggregate_point_seconds +=
        (time_ns() - started) * 1e-9
    context.counts.envelope_point_calls += 1
    context.counts.envelope_point_seconds +=
        (time_ns() - started) * 1e-9
    return out
end

function (scales::OMRFCompactCovariateLocalScales)(out, state,
        flow::Union{ZigZag,BouncyParticle}, left, right)
    base = scales.base
    context = base.context
    @inbounds for node in eachindex(context.seen)
        q = context.seen[node] - 1
        threshold_start = context.threshold_starts[node]
        minimum_velocity = 0.0
        maximum_velocity = 0.0
        A = 0.0
        for u in 1:q, v in 0:(u - 1)
            index = context.threshold_indices[threshold_start + u - 1]
            displacement_left, velocity = _omrf_coordinate_at(
                state, flow, base.anchor, index, left)
            displacement_right, _ = _omrf_coordinate_at(
                state, flow, base.anchor, index, right)
            minimum_velocity = min(minimum_velocity, velocity)
            maximum_velocity = max(maximum_velocity, velocity)
            if v > 0
                index = context.threshold_indices[threshold_start + v - 1]
                other_left, _ = _omrf_coordinate_at(
                    state, flow, base.anchor, index, left)
                other_right, _ = _omrf_coordinate_at(
                    state, flow, base.anchor, index, right)
                displacement_left -= other_left
                displacement_right -= other_right
            end
            A = max(A, abs(displacement_left), abs(displacement_right))
        end
        edge_displacements = base.edge_displacements
        edge_velocities = base.edge_velocities
        edges = context.incident_edges[node]
        @inbounds for k in eachindex(edges)
            index = context.interaction_indices[edges[k]]
            left_displacement, velocity = _omrf_coordinate_at(
                state, flow, base.anchor, index, left)
            right_displacement, _ = _omrf_coordinate_at(
                state, flow, base.anchor, index, right)
            edge_displacements[k] = max(abs(left_displacement),
                abs(right_displacement))
            edge_velocities[k] = abs(velocity)
        end
        _omrf_compact_fill_node!(out, scales, node, q, A,
            maximum_velocity - minimum_velocity,
            edge_displacements, edge_velocities)
    end
    return out
end

function (scales::OMRFCompactCovariateLocalScales)(out, state,
        flow::PDMPSamplers.AnyBoomerang, left, right)
    started = time_ns()
    base = scales.base
    context = base.context
    _cache_omrf_pattern_harmonic_coefficients!(base, state, flow)
    midpoint = (left + right) / 2
    midpoint_sine, midpoint_cosine = sincos(midpoint)
    half_width = (right - left) / 2
    @inbounds for node in eachindex(context.seen)
        q = context.seen[node] - 1
        threshold_start = context.threshold_starts[node]
        A = 0.0
        C = 0.0
        for u in 1:q, v in 0:(u - 1)
            index = context.threshold_indices[threshold_start + u - 1]
            offset = base.coordinate_offsets[index]
            cosine = base.coordinate_cosines[index]
            sine = base.coordinate_sines[index]
            if v > 0
                index = context.threshold_indices[threshold_start + v - 1]
                offset -= base.coordinate_offsets[index]
                cosine -= base.coordinate_cosines[index]
                sine -= base.coordinate_sines[index]
            end
            A = max(A, _omrf_harmonic_abs_cell_bound(offset, cosine, sine,
                midpoint_sine, midpoint_cosine, half_width))
            C = max(C, _omrf_harmonic_abs_cell_bound(0.0, sine, -cosine,
                midpoint_sine, midpoint_cosine, half_width))
        end
        edge_displacements = base.edge_displacements
        edge_velocities = base.edge_velocities
        edges = context.incident_edges[node]
        @inbounds for k in eachindex(edges)
            index = context.interaction_indices[edges[k]]
            offset = base.coordinate_offsets[index]
            cosine = base.coordinate_cosines[index]
            sine = base.coordinate_sines[index]
            edge_displacements[k] = _omrf_harmonic_abs_cell_bound(
                offset, cosine, sine, midpoint_sine, midpoint_cosine,
                half_width)
            edge_velocities[k] = _omrf_harmonic_abs_cell_bound(
                0.0, sine, -cosine, midpoint_sine, midpoint_cosine,
                half_width)
        end
        _omrf_compact_fill_node!(out, scales, node, q, A, C,
            edge_displacements, edge_velocities)
    end
    context.counts.compact_cell_bound_calls += 1
    context.counts.compact_cell_bound_seconds +=
        (time_ns() - started) * 1e-9
    context.counts.envelope_cell_calls += 1
    context.counts.envelope_cell_seconds +=
        (time_ns() - started) * 1e-9
    return out
end

PDMPSamplers.n_observations(envelope::OMRFCompactCovariateLocalEnvelope) =
    envelope.n_persons * envelope.n_nodes

function PDMPSamplers.observation_residual_bound(
        envelope::OMRFCompactCovariateLocalEnvelope, observation::Integer)
    started = time_ns()
    selected_allocated = envelope.signed_harmonic ?
        _omrf_allocated_bytes() : UInt(0)
    block0, person0 = divrem(observation - 1, envelope.n_persons)
    node = block0 + 1
    person = person0 + 1
    1 <= node <= envelope.n_nodes || throw(BoundsError(
        1:PDMPSamplers.n_observations(envelope), observation))
    result = if envelope.signed_bps || envelope.signed_harmonic
        envelope.signed_screen_valid || throw(ArgumentError(
            "signed BPS factor bound requested before point screening"))
        bound = _omrf_signed_factor_bound_coordinates(
            envelope.component_scales!.base.context,
            envelope.component_scales!.base.anchor,
            envelope.screen_x, envelope.screen_theta, person, node)
        if envelope.signed_harmonic
            counts = envelope.component_scales!.base.context.counts.signed_bps
            counts.harmonic_selected_point_bounds += 1
            counts.harmonic_selected_point_seconds +=
                (time_ns() - started) * 1e-9
            counts.harmonic_selected_point_allocated_bytes +=
                _omrf_allocated_bytes() - selected_allocated
        end
        bound
    else
        _omrf_compact_person_bound(envelope.component_scales!, node, person)
    end
    context = envelope.component_scales!.base.context
    context.counts.compact_observation_bound_calls += 1
    context.counts.compact_observation_bound_seconds +=
        (time_ns() - started) * 1e-9
    return result
end

function PDMPSamplers._draw_component(rng::Random.AbstractRNG,
        envelope::OMRFCompactCovariateLocalEnvelope, B::Real)
    if envelope.signed_bps || envelope.signed_harmonic
        draw_started = time_ns()
        draw_allocated = _omrf_allocated_bytes()
        envelope.signed_screen_valid || throw(ArgumentError(
            "signed BPS component draw requested before point screening"))
        left_mass = (1.0 - envelope.signed_screen_alpha) *
            envelope.signed_cell_left_mass
        right_mass = envelope.signed_screen_alpha *
            envelope.signed_cell_right_mass
        total = left_mass + right_mass
        if total <= 0
            result = 1
        elseif rand(rng) * total <= left_mass
            u = rand(rng) * envelope.signed_cell_left_mass
            result = min(searchsortedlast(
                envelope.signed_cell_cumulative_left, u) + 1,
                envelope.n_persons * envelope.n_nodes)
        else
            u = rand(rng) * envelope.signed_cell_right_mass
            result = min(searchsortedlast(
                envelope.signed_cell_cumulative_right, u) + 1,
                envelope.n_persons * envelope.n_nodes)
        end
        counts = envelope.component_scales!.base.context.counts.signed_bps
        counts.component_draws += 1
        counts.component_draw_seconds += (time_ns() - draw_started) * 1e-9
        counts.component_draw_allocated_bytes +=
            _omrf_allocated_bytes() - draw_allocated
        return result
    end
    u = rand(rng) * B
    return min(searchsortedlast(envelope.cumulative_masses, u) + 1,
        envelope.n_nodes)
end

@inline function _omrf_compact_search_term(cumulative::Matrix{Float64},
        node::Int, n_terms::Int, u::Float64)
    left = 1
    right = n_terms
    while left < right
        middle = (left + right) >>> 1
        if cumulative[node, middle] < u
            left = middle + 1
        else
            right = middle
        end
    end
    return left
end

function PDMPSamplers._draw_distinguished(rng::Random.AbstractRNG,
        envelope::OMRFCompactCovariateLocalEnvelope, node::Int)
    if envelope.signed_bps || envelope.signed_harmonic
        envelope.signed_screen_valid || throw(ArgumentError(
            "signed BPS distinguished draw requested before point screening"))
        return node
    end
    scales = envelope.component_scales!
    context = scales.base.context
    static = scales.static
    if envelope.chord_screen
        # The point screen only needs node totals.  Build the selected node's
        # term cumulative masses lazily, after component selection, so the
        # roof loop does not repeat the O(P^3) term traversal.
        α = clamp(envelope.screen_time - envelope.cell_left,
            0.0, envelope.cell_right - envelope.cell_left) /
            (envelope.cell_right - envelope.cell_left)
        n_terms = static.term_count[node]
        @inbounds for term in 1:n_terms
            scales.term_cumulative[node, term] =
                (1.0 - α) * envelope.cell_term_left[node, term] +
                α * envelope.cell_term_right[node, term]
        end
        envelope.chord_screen = false
    end
    selection_started = time_ns()
    total = envelope.scales[node]
    u = rand(rng) * total
    term = _omrf_compact_search_term(scales.term_cumulative, node,
        static.term_count[node], u)
    context.counts.compact_term_selection_calls += 1
    context.counts.compact_term_selection_seconds +=
        (time_ns() - selection_started) * 1e-9
    person_started = time_ns()
    kind = static.term_kind[node, term]
    person = if kind == 0x01
        rand(rng, 1:envelope.n_persons)
    elseif kind == 0x02
        table = static.covariate_aliases[static.term_k[node, term]]
        table === nothing ? rand(rng, 1:envelope.n_persons) : rand(rng, table)
    elseif kind == 0x03
        table = static.covariate_aliases[static.term_l[node, term]]
        table === nothing ? rand(rng, 1:envelope.n_persons) : rand(rng, table)
    else
        table = static.cross_aliases[static.term_k[node, term],
            static.term_l[node, term]]
        table === nothing ? rand(rng, 1:envelope.n_persons) : rand(rng, table)
    end
    context.counts.compact_person_selection_calls += 1
    context.counts.compact_person_selection_seconds +=
        (time_ns() - person_started) * 1e-9
    return (node - 1) * envelope.n_persons + person
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
        flow::PDMPSamplers.AnyBoomerang, t)
    # Every node uses the same harmonic phase.  Compute sin/cos and each
    # coordinate's displacement/velocity once, then reuse interaction
    # coordinates at both incident nodes.  This is algebraically identical to
    # `_omrf_norm_node_terms_at` and changes only repeated arithmetic.
    s, c = sincos(t)
    displacement = scales.coordinate_displacements
    velocity = scales.coordinate_velocities
    @inbounds for index in eachindex(scales.anchor)
        if _omrf_coordinate_free(state, index)
            centered = state.ξ.x[index] - flow.μ[index]
            displacement[index] = flow.μ[index] - scales.anchor[index] +
                centered * c + state.ξ.θ[index] * s
            velocity[index] = -centered * s + state.ξ.θ[index] * c
        else
            displacement[index] = state.ξ.x[index] - scales.anchor[index]
            velocity[index] = 0.0
        end
    end
    context = scales.context
    P = length(context.seen)
    @inbounds for node in 1:P
        q = context.seen[node] - 1
        if iszero(q)
            _fill_omrf_norm_scales!(out, node, P, 0.0, 0.0, 0.0, 0.0)
            continue
        end
        minimum_threshold = 0.0
        maximum_threshold = 0.0
        threshold_velocity = 0.0
        threshold_start = context.threshold_starts[node]
        for u in 1:q
            index = context.threshold_indices[threshold_start + u - 1]
            value = displacement[index]
            minimum_threshold = min(minimum_threshold, value)
            maximum_threshold = max(maximum_threshold, value)
            threshold_velocity += abs(velocity[index])
        end
        edge_displacement2 = 0.0
        edge_velocity2 = 0.0
        for edge in context.incident_edges[node]
            index = context.interaction_indices[edge]
            edge_displacement2 += displacement[index]^2
            edge_velocity2 += velocity[index]^2
        end
        _fill_omrf_norm_scales!(out, node, P,
            maximum_threshold - minimum_threshold,
            q * sqrt(edge_displacement2), threshold_velocity,
            q * sqrt(edge_velocity2))
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

function _omrf_compact_hcv_node_terms_at(
        scales::OMRFCompactHCVNormLocalScales,
        state, flow, node::Int, t)
    context = scales.context
    q = context.seen[node] - 1
    iszero(q) && return 0.0, 0.0, 0.0, 0.0
    threshold_start = context.threshold_starts[node]
    minimum_displacement = 0.0
    maximum_displacement = 0.0
    minimum_velocity = 0.0
    maximum_velocity = 0.0
    @inbounds for u in 1:q
        index = context.threshold_indices[threshold_start + u - 1]
        displacement, velocity = _omrf_coordinate_at(
            state, flow, scales.anchor, index, t)
        minimum_displacement = min(minimum_displacement, displacement)
        maximum_displacement = max(maximum_displacement, displacement)
        minimum_velocity = min(minimum_velocity, velocity)
        maximum_velocity = max(maximum_velocity, velocity)
    end
    edge_displacement2 = 0.0
    edge_velocity2 = 0.0
    @inbounds for edge in context.incident_edges[node]
        index = context.interaction_indices[edge]
        displacement, velocity = _omrf_coordinate_at(
            state, flow, scales.anchor, index, t)
        edge_displacement2 += displacement * displacement
        edge_velocity2 += velocity * velocity
    end
    A = maximum_displacement - minimum_displacement
    C = maximum_velocity - minimum_velocity
    return A, C, sqrt(edge_displacement2), sqrt(edge_velocity2)
end

@inline function _fill_omrf_compact_hcv_scales!(
        out, node, P, q, A, C, D, V)
    qD = q * D
    qV = q * V
    out[node] = C * A * A / 8
    out[P + node] = (2 * C * A * qD + qV * A * A) / 8
    out[2P + node] = (C * qD * qD + 2 * qV * A * qD) / 8
    out[3P + node] = qV * qD * qD / 8
    return nothing
end

function (scales::OMRFCompactHCVNormLocalScales)(out, state,
        flow::Union{ZigZag,BouncyParticle}, t)
    P = length(scales.context.seen)
    @inbounds for node in 1:P
        A, C, D, V = _omrf_compact_hcv_node_terms_at(
            scales, state, flow, node, t)
        _fill_omrf_compact_hcv_scales!(out, node, P,
            scales.context.seen[node] - 1, A, C, D, V)
    end
    return out
end

function (scales::OMRFCompactHCVNormLocalScales)(out, state,
        flow::Union{ZigZag,BouncyParticle}, left, right)
    P = length(scales.context.seen)
    @inbounds for node in 1:P
        Al, C, Dl, V = _omrf_compact_hcv_node_terms_at(
            scales, state, flow, node, left)
        Ar, _, Dr, _ = _omrf_compact_hcv_node_terms_at(
            scales, state, flow, node, right)
        _fill_omrf_compact_hcv_scales!(out, node, P,
            scales.context.seen[node] - 1,
            max(Al, Ar), C, max(Dl, Dr), V)
    end
    return out
end

@inline (scales::OMRFCompactHCVNormLocalScales)(out, state,
    flow::PDMPSamplers.PreconditionedDynamics, args...) =
    scales(out, state, flow.dynamics, args...)

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
        threshold_coefficients = scales.threshold_coefficients
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
    cell_started = time_ns()
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
    context.counts.envelope_cell_calls += 1
    context.counts.envelope_cell_seconds += (time_ns() - cell_started) * 1e-9
    return out
end

(scales::OMRFCovariateLocalScales)(out, state,
    flow::PDMPSamplers.PreconditionedDynamics, args...) =
    scales(out, state, flow.dynamics, args...)

(scales::OMRFNormLocalScales)(out, state,
    flow::PDMPSamplers.PreconditionedDynamics, args...) =
    scales(out, state, flow.dynamics, args...)

function _build_omrf_pattern_packed_metadata(
        component_nodes::Vector{Int}, component_covariates::Matrix{Int},
        context::OMRFResidualContext)
    n_nodes = length(context.seen)
    n_components = length(component_nodes)
    node_component_starts = Vector{Int}(undef, n_nodes)
    node_component_ends = Vector{Int}(undef, n_nodes)
    node_edge_starts = Vector{Int}(undef, n_nodes)
    node_edge_ends = Vector{Int}(undef, n_nodes)
    n_edges = sum(length, context.incident_edges)
    node_edge_indices = Vector{Int}(undef, n_edges)
    component_edge_offsets = Vector{Int}(undef, n_components + 1)
    component_edge_count = sum(
        length(context.incident_edges[component_nodes[component]])
        for component in eachindex(component_nodes))
    component_edge_covariates = Vector{Int}(undef, component_edge_count)

    edge_destination = 1
    component_edge_destination = 1
    @inbounds for node in 1:n_nodes
        first_component = findfirst(==(node), component_nodes)
        last_component = findlast(==(node), component_nodes)
        if first_component === nothing
            node_component_starts[node] = 1
            node_component_ends[node] = 0
        else
            node_component_starts[node] = first_component
            node_component_ends[node] = last_component
            all(==(node), view(component_nodes, first_component:last_component)) ||
                throw(ArgumentError(
                    "OMRF pattern components must be contiguous by node"))
        end

        edges = context.incident_edges[node]
        node_edge_starts[node] = edge_destination
        for edge in edges
            node_edge_indices[edge_destination] =
                context.interaction_indices[edge]
            edge_destination += 1
        end
        node_edge_ends[node] = edge_destination - 1

        if first_component !== nothing
            for component in first_component:last_component
                component_edge_offsets[component] = component_edge_destination
                for k in eachindex(edges)
                    component_edge_covariates[component_edge_destination] =
                        component_covariates[
                            component, context.incident_neighbours[node][k]]
                    component_edge_destination += 1
                end
            end
        end
    end
    component_edge_offsets[n_components + 1] = component_edge_destination
    return node_component_starts, node_component_ends,
        node_edge_starts, node_edge_ends, node_edge_indices,
        component_edge_offsets, component_edge_covariates
end

@inline function _omrf_active_pattern_hash(
        context::OMRFResidualContext, node::Int, person::Int,
        free::BitVector)
    hash = UInt64(1469598103934665603)
    @inbounds for k in eachindex(context.incident_edges[node])
        edge = context.incident_edges[node][k]
        index = context.interaction_indices[edge]
        free[index] || continue
        hash ⊻= UInt64(context.X[person,
            context.incident_neighbours[node][k]] + 1)
        hash *= UInt64(1099511628211)
    end
    return hash
end

@inline function _omrf_active_pattern_equal(
        context::OMRFResidualContext, node::Int, person::Int,
        representative::Int, free::BitVector)
    @inbounds for k in eachindex(context.incident_edges[node])
        edge = context.incident_edges[node][k]
        index = context.interaction_indices[edge]
        free[index] || continue
        neighbour = context.incident_neighbours[node][k]
        context.X[person, neighbour] == context.X[representative, neighbour] ||
            return false
    end
    return true
end

function _rebuild_omrf_active_pattern!(scales::OMRFPatternLocalScales)
    _rebuild_omrf_active_pattern!(scales, false)
end

function _rebuild_omrf_active_pattern!(scales::OMRFPatternLocalScales,
        use_all_free::Bool)
    started = time_ns()
    context = scales.context
    free = use_all_free ? scales.active_all_free : context.active_free
    groups = scales.active_groups::Matrix{Int}
    members = scales.active_members::Vector{Vector{Int}}
    totals = scales.active_totals::Vector{Float64}
    X = context.X
    N, P = size(X)
    fill!(totals, 0.0)
    @inbounds for member in members
        empty!(member)
    end

    epoch = scales.active_hash_epoch
    component_base = 0
    component_edge_destination = 1
    factorized = context.factorization === :person_node

    @inbounds for node in 1:P
        epoch = epoch == typemax(Int) ? 1 : epoch + 1
        epoch == 1 && fill!(scales.active_hash_stamps, 0)
        first_component = component_base + 1
        group_count = 0
        for person in 1:N
            key = _omrf_active_pattern_hash(context, node, person, free)
            slot = Int(key & UInt64(length(scales.active_hash_keys) - 1)) + 1
            local_group = 0
            while scales.active_hash_stamps[slot] == epoch
                if scales.active_hash_keys[slot] == key
                    candidate = scales.active_hash_values[slot]
                    representative = scales.active_representatives[candidate]
                    if _omrf_active_pattern_equal(
                            context, node, person, representative, free)
                        local_group = candidate
                        break
                    end
                end
                slot = slot == length(scales.active_hash_keys) ? 1 : slot + 1
            end
            if iszero(local_group)
                group_count += 1
                local_group = group_count
                scales.active_hash_stamps[slot] = epoch
                scales.active_hash_keys[slot] = key
                scales.active_hash_values[slot] = local_group
                scales.active_representatives[local_group] = person
                component = component_base + local_group
                scales.component_nodes[component] = node
                for column in axes(scales.component_covariates, 2)
                    scales.component_covariates[component, column] = 0
                end
                for k in eachindex(context.incident_edges[node])
                    edge = context.incident_edges[node][k]
                    index = context.interaction_indices[edge]
                    free[index] || continue
                    neighbour = context.incident_neighbours[node][k]
                    scales.component_covariates[component, neighbour] =
                        X[person, neighbour]
                end
            end
            component = component_base + local_group
            observation = factorized ? (node - 1) * N + person : person
            partition = factorized ? 1 : node
            groups[partition, observation] = component
            push!(members[component], observation)
            totals[component] += 1.0
        end

        node_component_starts = first_component
        node_component_ends = component_base + group_count
        scales.node_component_starts[node] = node_component_starts
        scales.node_component_ends[node] = node_component_ends
        edges = context.incident_edges[node]
        for component in node_component_starts:node_component_ends
            scales.component_edge_offsets[component] =
                component_edge_destination
            for k in eachindex(edges)
                edge = edges[k]
                index = context.interaction_indices[edge]
                scales.component_edge_covariates[component_edge_destination] =
                    free[index] ? scales.component_covariates[
                        component, context.incident_neighbours[node][k]] : 0
                component_edge_destination += 1
            end
        end
        component_base = node_component_ends
    end
    scales.component_edge_offsets[component_base + 1] =
        component_edge_destination
    scales.active_component_count = component_base
    scales.active_hash_epoch = epoch
    scales.active_generation = context.active_generation
    scales.active_effective_full = use_all_free
    context.active_rebuilds += 1
    context.active_rebuild_seconds += (time_ns() - started) * 1e-9
    context.active_group_count_sum += component_base
    if context.active_rebuilds == 1
        context.active_group_count_min = component_base
        context.active_group_count_max = component_base
    else
        context.active_group_count_min = min(
            context.active_group_count_min, component_base)
        context.active_group_count_max = max(
            context.active_group_count_max, component_base)
    end
    return nothing
end

@inline function _ensure_omrf_active_pattern!(
        scales::OMRFPatternLocalScales)
    scales.active_compression || return nothing
    context = scales.context
    scales.active_generation == context.active_generation ||
        _rebuild_omrf_active_pattern!(scales)
    return nothing
end

@inline function _omrf_active_pattern_requires_full(
        scales::OMRFPatternLocalScales, state)
    state isa StickyPDMPState || return false
    context = scales.context
    @inbounds for index in context.interaction_indices
        !state.free[index] && state.ξ.x[index] != scales.anchor[index] &&
            return true
    end
    return false
end

@inline function _ensure_omrf_active_pattern!(
        scales::OMRFPatternLocalScales, state)
    scales.active_compression || return nothing
    context = scales.context
    use_all_free = _omrf_active_pattern_requires_full(scales, state)
    if scales.active_generation != context.active_generation ||
            scales.active_effective_full != use_all_free
        _rebuild_omrf_active_pattern!(scales, use_all_free)
    end
    return nothing
end

function _attach_omrf_active_pattern!(
        envelope::PDMPSamplers.GroupedResidualEnvelope,
        scales::OMRFPatternLocalScales)
    scales.active_groups = envelope.groups
    scales.active_members = envelope.members
    scales.active_totals = envelope.totals
    N = size(scales.context.X, 1)
    @inbounds for member in envelope.members
        sizehint!(member, N)
    end
    return envelope
end

function PDMPSamplers.set_active_set!(
        scales::OMRFPatternLocalScales, free::BitVector)
    scales.active_compression || return nothing
    context = scales.context
    length(free) == length(context.active_free) || throw(DimensionMismatch(
        "active OMRF set has the wrong dimension"))
    changed = false
    @inbounds for index in eachindex(free, context.active_free)
        changed |= free[index] != context.active_free[index]
    end
    changed || return nothing
    copyto!(context.active_free, free)
    context.active_generation += 1
    return nothing
end

function PDMPSamplers.set_active_set!(
        envelope::PDMPSamplers.GroupedResidualEnvelope, free::BitVector)
    PDMPSamplers.set_active_set!(envelope.component_scales!, free)
    return nothing
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

@inline function _omrf_pattern_zigzag_scale(
        scales::OMRFPatternLocalScales, component::Int,
        q::Int, threshold_start::Int, node_edge_start::Int,
        node_edge_end::Int, displacements, velocities)
    context = scales.context
    field_displacement = 0.0
    field_velocity = 0.0
    weighted_edge_velocity = 0.0
    component_edge_start = scales.component_edge_offsets[component]
    @inbounds for edge_position in node_edge_start:node_edge_end
        k = edge_position - node_edge_start + 1
        covariate = scales.component_edge_covariates[
            component_edge_start + k - 1]
        iszero(covariate) && continue
        index = scales.node_edge_indices[edge_position]
        field_displacement += covariate * displacements[index]
        field_velocity += covariate * velocities[index]
        weighted_edge_velocity += covariate * abs(velocities[index])
    end
    minimum_logit = 0.0
    maximum_logit = 0.0
    minimum_logit_velocity = 0.0
    maximum_logit_velocity = 0.0
    threshold_velocity = 0.0
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
    return 0.25 * (maximum_logit - minimum_logit) *
        (threshold_velocity + q * weighted_edge_velocity)
end

@inline function _omrf_pattern_bouncy_scale(
        scales::OMRFPatternLocalScales, component::Int,
        q::Int, threshold_start::Int, node_edge_start::Int,
        node_edge_end::Int, displacements, velocities)
    context = scales.context
    field_displacement = 0.0
    field_velocity = 0.0
    component_edge_start = scales.component_edge_offsets[component]
    @inbounds for edge_position in node_edge_start:node_edge_end
        k = edge_position - node_edge_start + 1
        covariate = scales.component_edge_covariates[
            component_edge_start + k - 1]
        iszero(covariate) && continue
        index = scales.node_edge_indices[edge_position]
        field_displacement += covariate * displacements[index]
        field_velocity += covariate * velocities[index]
    end
    minimum_logit = 0.0
    maximum_logit = 0.0
    minimum_logit_velocity = 0.0
    maximum_logit_velocity = 0.0
    @inbounds for u in 1:q
        index = context.threshold_indices[threshold_start + u - 1]
        logit = displacements[index] + u * field_displacement
        minimum_logit = min(minimum_logit, logit)
        maximum_logit = max(maximum_logit, logit)
        logit_velocity = velocities[index] + u * field_velocity
        minimum_logit_velocity = min(minimum_logit_velocity, logit_velocity)
        maximum_logit_velocity = max(maximum_logit_velocity, logit_velocity)
    end
    return 0.25 * (maximum_logit - minimum_logit) *
        (maximum_logit_velocity - minimum_logit_velocity)
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

@inline function _omrf_add_person_node_zigzag_bounds!(
        context::OMRFResidualContext, out, state, flow,
        anchor, person::Int, node::Int, t)
    q = context.seen[node] - 1
    iszero(q) && return nothing
    field_displacement = 0.0
    @inbounds for k in eachindex(context.incident_edges[node])
        covariate = context.X[person,
            context.incident_neighbours[node][k]]
        edge_index = context.interaction_indices[
            context.incident_edges[node][k]]
        displacement, _ = _omrf_coordinate_at(
            state, flow, anchor, edge_index, t)
        field_displacement += covariate * displacement
    end
    minimum_logit = 0.0
    maximum_logit = 0.0
    threshold_start = context.threshold_starts[node]
    @inbounds for u in 1:q
        index = context.threshold_indices[threshold_start + u - 1]
        displacement, _ = _omrf_coordinate_at(
            state, flow, anchor, index, t)
        logit = displacement + u * field_displacement
        minimum_logit = min(minimum_logit, logit)
        maximum_logit = max(maximum_logit, logit)
    end
    coordinate_scale = 0.25 * (maximum_logit - minimum_logit)
    @inbounds for u in 1:q
        index = context.threshold_indices[threshold_start + u - 1]
        out[index] += coordinate_scale
    end
    @inbounds for k in eachindex(context.incident_edges[node])
        covariate = context.X[person,
            context.incident_neighbours[node][k]]
        edge_index = context.interaction_indices[
            context.incident_edges[node][k]]
        out[edge_index] += coordinate_scale * q * abs(covariate)
    end
    return nothing
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
    oracle.exact_selected_refinement || return PDMPSamplers._signed_subset_bound(
        state, gradient, flow, D, M)
    supported_flow = flow isa Union{ZigZag,BouncyParticle,
        PDMPSamplers.AnyBoomerang,PDMPSamplers.PreconditionedDynamics}
    supported_flow || return PDMPSamplers._signed_subset_bound(
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
    if isfinite(oracle.lazy_residual_subset_bound)
        bound = D + oracle.lazy_residual_subset_bound
        oracle.lazy_residual_subset_bound = NaN
        return min(M, bound)
    end
    context = oracle.residual_context
    if context.node_box_first_order
        _omrf_prepare_exact_selected_residual!(oracle, state, subset)
        residual_rate_bound = scale *
            _omrf_exact_selected_residual_abs_rate(
                state, oracle.cached_residual, flow,
                context.touched_positions)
        return min(M, D + residual_rate_bound)
    end
    zigzag = flow isa ZigZag ||
        (flow isa PDMPSamplers.PreconditionedDynamics &&
         flow.dynamics isa ZigZag)
    if context.zigzag_sign_aware && zigzag && !context.use_hcv &&
            gradient !== nothing
        fill!(context.signed_coordinate_bounds, 0.0)
        N = size(context.X, 1)
        if context.factorization === :person
            @inbounds for person in subset, node in eachindex(context.seen)
                _omrf_add_person_node_zigzag_bounds!(context,
                    context.signed_coordinate_bounds, state, flow,
                    oracle.anchor, person, node, 0.0)
            end
        else
            @inbounds for factor in subset
                node0, person0 = divrem(factor - 1, N)
                _omrf_add_person_node_zigzag_bounds!(context,
                    context.signed_coordinate_bounds, state, flow,
                    oracle.anchor, person0 + 1, node0 + 1, 0.0)
            end
        end
        bound = D
        @inbounds for index in eachindex(context.signed_coordinate_bounds)
            residual_bound = scale *
                context.signed_coordinate_bounds[index]
            signed_deterministic = state.ξ.θ[index] * gradient[index]
            bound += max(0.0, signed_deterministic + residual_bound) -
                max(0.0, signed_deterministic)
        end
        return min(M, bound)
    end
    factor_pattern_dynamics = flow isa Union{
        ZigZag,BouncyParticle,PDMPSamplers.AnyBoomerang} ||
        (flow isa PDMPSamplers.PreconditionedDynamics &&
         flow.dynamics isa Union{ZigZag,BouncyParticle,
             PDMPSamplers.AnyBoomerang})
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

PDMPSamplers.needs_subsampling_residual_gradient(
    oracle::OMRFAnalyticSubsamplingOracle) =
    oracle.residual_context.zigzag_sign_aware &&
    !oracle.residual_context.use_hcv

function PDMPSamplers.subsampling_residual_subset_bound(
        oracle::OMRFAnalyticSubsamplingOracle, state, flow,
        D, M, subset, scale)
    return PDMPSamplers.subsampling_residual_subset_bound(
        oracle, state, nothing, flow, D, M, subset, scale)
end

function _cache_omrf_pattern_linear_cell_coordinates!(
        scales::OMRFPatternLocalScales, state,
        flow::Union{ZigZag,BouncyParticle}, left, right)
    left_offsets = scales.cell_left_offsets
    right_offsets = scales.cell_right_offsets
    velocities = scales.cell_velocities
    context = scales.context
    @inbounds for index in context.threshold_indices
        left_offsets[index], velocity = _omrf_coordinate_at(
            state, flow, scales.anchor, index, left)
        right_offsets[index], _ = _omrf_coordinate_at(
            state, flow, scales.anchor, index, right)
        velocities[index] = velocity
    end
    @inbounds for index in context.interaction_indices
        left_offsets[index], velocity = _omrf_coordinate_at(
            state, flow, scales.anchor, index, left)
        right_offsets[index], _ = _omrf_coordinate_at(
            state, flow, scales.anchor, index, right)
        velocities[index] = velocity
    end
    return nothing
end

function _cache_omrf_pattern_point_coordinates!(
        scales::OMRFPatternLocalScales, state,
        flow::Union{ZigZag,BouncyParticle}, t)
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

function _omrf_pattern_point_scales!(out, scales::OMRFPatternLocalScales,
        ::ZigZag)
    context = scales.context
    displacements = scales.coordinate_offsets
    velocities = scales.coordinate_cosines
    @inbounds for node in eachindex(scales.node_component_starts)
        q = context.seen[node] - 1
        first_component = scales.node_component_starts[node]
        last_component = scales.node_component_ends[node]
        if q == 0
            for component in first_component:last_component
                out[component] = 0.0
            end
            continue
        end
        threshold_start = context.threshold_starts[node]
        edge_start = scales.node_edge_starts[node]
        edge_end = scales.node_edge_ends[node]
        for component in first_component:last_component
            out[component] = _omrf_pattern_zigzag_scale(
                scales, component, q, threshold_start,
                edge_start, edge_end, displacements, velocities)
        end
    end
    if scales.active_compression
        @inbounds for component in (scales.active_component_count + 1):length(out)
            out[component] = 0.0
        end
    end
    return out
end

function _omrf_pattern_point_scales!(out, scales::OMRFPatternLocalScales,
        ::BouncyParticle)
    context = scales.context
    displacements = scales.coordinate_offsets
    velocities = scales.coordinate_cosines
    @inbounds for node in eachindex(scales.node_component_starts)
        q = context.seen[node] - 1
        first_component = scales.node_component_starts[node]
        last_component = scales.node_component_ends[node]
        if q == 0
            for component in first_component:last_component
                out[component] = 0.0
            end
            continue
        end
        threshold_start = context.threshold_starts[node]
        edge_start = scales.node_edge_starts[node]
        edge_end = scales.node_edge_ends[node]
        for component in first_component:last_component
            out[component] = _omrf_pattern_bouncy_scale(
                scales, component, q, threshold_start,
                edge_start, edge_end, displacements, velocities)
        end
    end
    if scales.active_compression
        @inbounds for component in (scales.active_component_count + 1):length(out)
            out[component] = 0.0
        end
    end
    return out
end

@inline _omrf_pattern_point_scales_from_nodes!(
        out, scales::OMRFPatternLocalScales, flow::ZigZag) =
    _omrf_pattern_point_scales!(out, scales, flow)

@inline _omrf_pattern_point_scales_from_nodes!(
        out, scales::OMRFPatternLocalScales, flow::BouncyParticle) =
    _omrf_pattern_point_scales!(out, scales, flow)

@inline _omrf_pattern_point_scales_from_nodes!(
        out, scales::OMRFPatternLocalScales,
        flow::PDMPSamplers.AnyBoomerang) =
    _omrf_pattern_point_scales!(out, scales, flow)

function (scales::OMRFPatternLocalScales)(out, state,
        flow::ZigZag, t)
    _ensure_omrf_active_pattern!(scales, state)
    _cache_omrf_pattern_point_coordinates!(scales, state, flow, t)
    _omrf_pattern_point_scales!(out, scales, flow)
    return out
end

function (scales::OMRFPatternLocalScales)(out, state,
        flow::BouncyParticle, t)
    _ensure_omrf_active_pattern!(scales, state)
    _cache_omrf_pattern_point_coordinates!(scales, state, flow, t)
    _omrf_pattern_point_scales!(out, scales, flow)
    return out
end

function _omrf_pattern_cell_scales!(out, scales::OMRFPatternLocalScales,
        ::ZigZag)
    context = scales.context
    left_offsets = scales.cell_left_offsets
    right_offsets = scales.cell_right_offsets
    velocities = scales.cell_velocities
    @inbounds for node in eachindex(scales.node_component_starts)
        q = context.seen[node] - 1
        first_component = scales.node_component_starts[node]
        last_component = scales.node_component_ends[node]
        if q == 0
            for component in first_component:last_component
                out[component] = 0.0
            end
            continue
        end
        threshold_start = context.threshold_starts[node]
        edge_start = scales.node_edge_starts[node]
        edge_end = scales.node_edge_ends[node]
        for component in first_component:last_component
            left_scale = _omrf_pattern_zigzag_scale(
                scales, component, q, threshold_start,
                edge_start, edge_end, left_offsets, velocities)
            right_scale = _omrf_pattern_zigzag_scale(
                scales, component, q, threshold_start,
                edge_start, edge_end, right_offsets, velocities)
            out[component] = max(left_scale, right_scale)
        end
    end
    if scales.active_compression
        @inbounds for component in (scales.active_component_count + 1):length(out)
            out[component] = 0.0
        end
    end
    return out
end

function _omrf_pattern_cell_scales!(out, scales::OMRFPatternLocalScales,
        ::BouncyParticle)
    context = scales.context
    left_offsets = scales.cell_left_offsets
    right_offsets = scales.cell_right_offsets
    velocities = scales.cell_velocities
    @inbounds for node in eachindex(scales.node_component_starts)
        q = context.seen[node] - 1
        first_component = scales.node_component_starts[node]
        last_component = scales.node_component_ends[node]
        if q == 0
            for component in first_component:last_component
                out[component] = 0.0
            end
            continue
        end
        threshold_start = context.threshold_starts[node]
        edge_start = scales.node_edge_starts[node]
        edge_end = scales.node_edge_ends[node]
        for component in first_component:last_component
            left_scale = _omrf_pattern_bouncy_scale(
                scales, component, q, threshold_start,
                edge_start, edge_end, left_offsets, velocities)
            right_scale = _omrf_pattern_bouncy_scale(
                scales, component, q, threshold_start,
                edge_start, edge_end, right_offsets, velocities)
            out[component] = max(left_scale, right_scale)
        end
    end
    if scales.active_compression
        @inbounds for component in (scales.active_component_count + 1):length(out)
            out[component] = 0.0
        end
    end
    return out
end

function (scales::OMRFPatternLocalScales)(out, state,
        flow::ZigZag, left, right)
    _ensure_omrf_active_pattern!(scales, state)
    _cache_omrf_pattern_linear_cell_coordinates!(
        scales, state, flow, left, right)
    _omrf_pattern_cell_scales!(out, scales, flow)
    return out
end

function (scales::OMRFPatternLocalScales)(out, state,
        flow::BouncyParticle, left, right)
    _ensure_omrf_active_pattern!(scales, state)
    _cache_omrf_pattern_linear_cell_coordinates!(
        scales, state, flow, left, right)
    _omrf_pattern_cell_scales!(out, scales, flow)
    return out
end

function (scales::OMRFPatternLocalScales)(out, state,
        flow::PDMPSamplers.AnyBoomerang, t)
    _ensure_omrf_active_pattern!(scales, state)
    _cache_omrf_pattern_point_coordinates!(scales, state, flow, t)
    _omrf_pattern_point_scales!(out, scales, flow)
    return out
end

function _omrf_pattern_point_scales!(out, scales::OMRFPatternLocalScales,
        ::PDMPSamplers.AnyBoomerang)
    context = scales.context
    displacements = scales.coordinate_offsets
    velocities = scales.coordinate_cosines
    @inbounds for node in eachindex(scales.node_component_starts)
        q = context.seen[node] - 1
        first_component = scales.node_component_starts[node]
        last_component = scales.node_component_ends[node]
        if q == 0
            for component in first_component:last_component
                out[component] = 0.0
            end
            continue
        end
        threshold_start = context.threshold_starts[node]
        edge_start = scales.node_edge_starts[node]
        edge_end = scales.node_edge_ends[node]
        for component in first_component:last_component
            out[component] = _omrf_pattern_bouncy_scale(
                scales, component, q, threshold_start,
                edge_start, edge_end, displacements, velocities)
        end
    end
    if scales.active_compression
        @inbounds for component in (scales.active_component_count + 1):length(out)
            out[component] = 0.0
        end
    end
    return out
end

# The generic GroupedResidualEnvelope screening path materializes every
# component scale through a separate callback and constructs cumulative
# masses.  The pattern certificate still visits its structural components,
# but the pointwise aggregate can be accumulated in the same pass while
# retaining the scales needed if the aggregate screen accepts.
function _omrf_pattern_screening_bound!(
        envelope::PDMPSamplers.GroupedResidualEnvelope{F,C}, state,
        flow::ZigZag, t) where {F<:OMRFPatternLocalScales,C}
    scales = envelope.component_scales!
    _ensure_omrf_active_pattern!(scales, state)
    _cache_omrf_pattern_point_coordinates!(scales, state, flow, t)
    component_scales = envelope.scales
    cumulative_masses = envelope.cumulative_masses
    result = 0.0
    _omrf_pattern_point_scales_from_nodes!(
        component_scales, scales, flow)
    @inbounds for component in eachindex(scales.component_nodes, envelope.totals)
        result += component_scales[component] * envelope.totals[component]
        cumulative_masses[component] = result
    end
    return nextfloat(result)
end

function _omrf_pattern_screening_bound!(
        envelope::PDMPSamplers.GroupedResidualEnvelope{F,C}, state,
        flow::BouncyParticle, t) where {F<:OMRFPatternLocalScales,C}
    scales = envelope.component_scales!
    _ensure_omrf_active_pattern!(scales, state)
    _cache_omrf_pattern_point_coordinates!(scales, state, flow, t)
    component_scales = envelope.scales
    cumulative_masses = envelope.cumulative_masses
    result = 0.0
    _omrf_pattern_point_scales!(component_scales, scales, flow)
    @inbounds for component in eachindex(scales.component_nodes, envelope.totals)
        result += component_scales[component] * envelope.totals[component]
        cumulative_masses[component] = result
    end
    return nextfloat(result)
end

PDMPSamplers.screening_residual_bound(
        envelope::PDMPSamplers.GroupedResidualEnvelope{F,C}, state,
        flow::ZigZag, t) where {F<:OMRFPatternLocalScales,C} =
    _omrf_pattern_screening_bound!(envelope, state, flow, t)

PDMPSamplers.screening_residual_bound(
        envelope::PDMPSamplers.GroupedResidualEnvelope{F,C}, state,
        flow::BouncyParticle, t) where {F<:OMRFPatternLocalScales,C} =
    _omrf_pattern_screening_bound!(envelope, state, flow, t)

function PDMPSamplers.prepare_residual_sampling!(
        envelope::PDMPSamplers.GroupedResidualEnvelope{F,C}, state,
        flow::ZigZag, t) where
        {F<:OMRFPatternLocalScales,C}
    # `screening_residual_bound` has already populated both buffers in the
    # same pass.  Keeping this hook explicit documents that no second scale
    # evaluation is needed before drawing a subset.
    return nothing
end

PDMPSamplers.prepare_residual_sampling!(
        envelope::PDMPSamplers.GroupedResidualEnvelope{F,C}, state,
        flow::BouncyParticle, t) where {F<:OMRFPatternLocalScales,C} = nothing

PDMPSamplers.screening_residual_bound(
        envelope::PDMPSamplers.GroupedResidualEnvelope{F,C}, state,
        flow::PDMPSamplers.PreconditionedDynamics, t) where
        {F<:OMRFPatternLocalScales,C} =
    PDMPSamplers.screening_residual_bound(
        envelope, state, flow.dynamics, t)

PDMPSamplers.prepare_residual_sampling!(
        envelope::PDMPSamplers.GroupedResidualEnvelope{F,C}, state,
        flow::PDMPSamplers.PreconditionedDynamics, t) where
        {F<:OMRFPatternLocalScales,C} =
    PDMPSamplers.prepare_residual_sampling!(
        envelope, state, flow.dynamics, t)

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

@inline function _omrf_interval_product_lower(xl, xu, yl, yu)
    return min(xl * yl, xl * yu, xu * yl, xu * yu)
end

@inline function _omrf_interval_product_upper(xl, xu, yl, yu)
    return max(xl * yl, xl * yu, xu * yl, xu * yu)
end

function _cache_omrf_node_box_point!(scales::OMRFNodeBoxFirstOrderScales,
        state, flow, t)
    @inbounds for index in eachindex(scales.anchor)
        value, velocity = _omrf_coordinate_at_anchor_value(
            state, flow, 0.0, index, t)
        scales.displacement_lower[index] = value
        scales.displacement_upper[index] = value
        scales.velocity_lower[index] = velocity
        scales.velocity_upper[index] = velocity
    end
    return scales
end

function _cache_omrf_node_box_cell!(scales::OMRFNodeBoxFirstOrderScales,
        state, flow::PDMPSamplers.AnyBoomerang, left, right)
    @inbounds for index in eachindex(scales.anchor)
        offset, cosine, sine = _omrf_harmonic_coefficients_anchor_value(
            state, flow, 0.0, index)
        scales.displacement_lower[index], scales.displacement_upper[index] =
            _omrf_harmonic_cell_interval(
                offset, cosine, sine, left, right)
        scales.velocity_lower[index], scales.velocity_upper[index] =
            _omrf_harmonic_cell_interval(
                0.0, sine, -cosine, left, right)
    end
    return scales
end

function _cache_omrf_node_box_cell!(scales::OMRFNodeBoxFirstOrderScales,
        state, flow::BouncyParticle, left, right)
    @inbounds for index in eachindex(scales.anchor)
        left_displacement, velocity = _omrf_coordinate_at_anchor_value(
            state, flow, 0.0, index, left)
        right_displacement, _ = _omrf_coordinate_at_anchor_value(
            state, flow, 0.0, index, right)
        scales.displacement_lower[index] = min(
            left_displacement, right_displacement)
        scales.displacement_upper[index] = max(
            left_displacement, right_displacement)
        scales.velocity_lower[index] = velocity
        scales.velocity_upper[index] = velocity
    end
    return scales
end

@inline function _omrf_node_reference(scales::OMRFNodeBoxFirstOrderScales,
        node::Int, index::Int)
    bank = scales.node_bank
    return bank === nothing ? scales.anchor[index] :
        bank.candidates[index, scales.node_sources[node]]
end

"""Certified categorical residual bound, optionally using total variation."""
@inline function _omrf_categorical_residual_bound(
        displacement_range::Float64, velocity_range::Float64,
        capped::Bool)
    multiplier = capped ? min(0.25 * displacement_range, 1.0) :
        0.25 * displacement_range
    value = velocity_range * multiplier
    return iszero(value) ? 0.0 : nextfloat(value)
end

@inline function _omrf_node_box_first_order_bound(
        scales::OMRFNodeBoxFirstOrderScales, component::Int, node::Int)
    context = scales.context
    q = context.seen[node] - 1
    iszero(q) && return 0.0
    field_displacement_lower = 0.0
    field_displacement_upper = 0.0
    field_velocity_lower = 0.0
    field_velocity_upper = 0.0
    @inbounds for k in eachindex(context.incident_edges[node])
        edge_index = context.interaction_indices[
            context.incident_edges[node][k]]
        reference = _omrf_node_reference(scales, node, edge_index)
        neighbour = context.incident_neighbours[node][k]
        covariate_lower = scales.minima[component, neighbour]
        covariate_upper = scales.maxima[component, neighbour]
        field_displacement_lower += _omrf_interval_product_lower(
            covariate_lower, covariate_upper,
            scales.displacement_lower[edge_index] - reference,
            scales.displacement_upper[edge_index] - reference)
        field_displacement_upper += _omrf_interval_product_upper(
            covariate_lower, covariate_upper,
            scales.displacement_lower[edge_index] - reference,
            scales.displacement_upper[edge_index] - reference)
        field_velocity_lower += _omrf_interval_product_lower(
            covariate_lower, covariate_upper,
            scales.velocity_lower[edge_index],
            scales.velocity_upper[edge_index])
        field_velocity_upper += _omrf_interval_product_upper(
            covariate_lower, covariate_upper,
            scales.velocity_lower[edge_index],
            scales.velocity_upper[edge_index])
    end
    displacement_minimum = 0.0
    displacement_maximum = 0.0
    velocity_minimum = 0.0
    velocity_maximum = 0.0
    threshold_start = context.threshold_starts[node]
    @inbounds for category in 1:q
        threshold_index = context.threshold_indices[
            threshold_start + category - 1]
        threshold_reference = _omrf_node_reference(
            scales, node, threshold_index)
        category_displacement_lower =
            scales.displacement_lower[threshold_index] - threshold_reference +
            category * field_displacement_lower
        category_displacement_upper =
            scales.displacement_upper[threshold_index] - threshold_reference +
            category * field_displacement_upper
        category_velocity_lower = scales.velocity_lower[threshold_index] +
            category * field_velocity_lower
        category_velocity_upper = scales.velocity_upper[threshold_index] +
            category * field_velocity_upper
        displacement_minimum = min(
            displacement_minimum, category_displacement_lower)
        displacement_maximum = max(
            displacement_maximum, category_displacement_upper)
        velocity_minimum = min(velocity_minimum, category_velocity_lower)
        velocity_maximum = max(velocity_maximum, category_velocity_upper)
    end
    displacement_range = displacement_maximum - displacement_minimum
    velocity_range = velocity_maximum - velocity_minimum
    # For categorical p, Hoeffding's covariance inequality yields
    # |b'*(p(x)-p(anchor))| <= range(a)*range(b)/4.  The zero reference
    # category is included by the zero initial extrema above.
    return _omrf_categorical_residual_bound(displacement_range,
        velocity_range, scales.categorical_residual_cap)
end

@inline function _omrf_node_box_two_bounds!(out,
        scales::OMRFNodeBoxFirstOrderScales, node::Int,
        first_component::Int, second_component::Int)
    context = scales.context
    q = context.seen[node] - 1
    if iszero(q)
        out[first_component] = 0.0
        out[second_component] = 0.0
        return nothing
    end
    first_dl = 0.0; first_du = 0.0
    first_vl = 0.0; first_vu = 0.0
    second_dl = 0.0; second_du = 0.0
    second_vl = 0.0; second_vu = 0.0
    # G=2 specialization: every incident edge is loaded once and contributes
    # to both contiguous node-box rows. This is algebraically identical to
    # two calls of `_omrf_node_box_first_order_bound`.
    @inbounds for k in eachindex(context.incident_edges[node])
        edge_index = context.interaction_indices[
            context.incident_edges[node][k]]
        reference = _omrf_node_reference(scales, node, edge_index)
        neighbour = context.incident_neighbours[node][k]
        displacement_lower = scales.displacement_lower[edge_index] - reference
        displacement_upper = scales.displacement_upper[edge_index] - reference
        velocity_lower = scales.velocity_lower[edge_index]
        velocity_upper = scales.velocity_upper[edge_index]

        lower = scales.minima[first_component, neighbour]
        upper = scales.maxima[first_component, neighbour]
        first_dl += _omrf_interval_product_lower(
            lower, upper, displacement_lower, displacement_upper)
        first_du += _omrf_interval_product_upper(
            lower, upper, displacement_lower, displacement_upper)
        first_vl += _omrf_interval_product_lower(
            lower, upper, velocity_lower, velocity_upper)
        first_vu += _omrf_interval_product_upper(
            lower, upper, velocity_lower, velocity_upper)

        lower = scales.minima[second_component, neighbour]
        upper = scales.maxima[second_component, neighbour]
        second_dl += _omrf_interval_product_lower(
            lower, upper, displacement_lower, displacement_upper)
        second_du += _omrf_interval_product_upper(
            lower, upper, displacement_lower, displacement_upper)
        second_vl += _omrf_interval_product_lower(
            lower, upper, velocity_lower, velocity_upper)
        second_vu += _omrf_interval_product_upper(
            lower, upper, velocity_lower, velocity_upper)
    end
    first_dmin = 0.0; first_dmax = 0.0
    first_vmin = 0.0; first_vmax = 0.0
    second_dmin = 0.0; second_dmax = 0.0
    second_vmin = 0.0; second_vmax = 0.0
    threshold_start = context.threshold_starts[node]
    @inbounds for category in 1:q
        threshold_index = context.threshold_indices[
            threshold_start + category - 1]
        threshold_reference = _omrf_node_reference(
            scales, node, threshold_index)
        threshold_dl = scales.displacement_lower[threshold_index] -
            threshold_reference
        threshold_du = scales.displacement_upper[threshold_index] -
            threshold_reference
        threshold_vl = scales.velocity_lower[threshold_index]
        threshold_vu = scales.velocity_upper[threshold_index]
        first_dmin = min(first_dmin, threshold_dl + category * first_dl)
        first_dmax = max(first_dmax, threshold_du + category * first_du)
        first_vmin = min(first_vmin, threshold_vl + category * first_vl)
        first_vmax = max(first_vmax, threshold_vu + category * first_vu)
        second_dmin = min(second_dmin, threshold_dl + category * second_dl)
        second_dmax = max(second_dmax, threshold_du + category * second_du)
        second_vmin = min(second_vmin, threshold_vl + category * second_vl)
        second_vmax = max(second_vmax, threshold_vu + category * second_vu)
    end
    out[first_component] = _omrf_categorical_residual_bound(
        first_dmax - first_dmin, first_vmax - first_vmin,
        scales.categorical_residual_cap)
    out[second_component] = _omrf_categorical_residual_bound(
        second_dmax - second_dmin, second_vmax - second_vmin,
        scales.categorical_residual_cap)
    return nothing
end

@inline function _omrf_node_box_fill_bounds!(out,
        scales::OMRFNodeBoxFirstOrderScales)
    node_count = length(scales.context.seen)
    if length(out) == 2node_count
        @inbounds for node in 1:node_count
            first_component = 2node - 1
            second_component = first_component + 1
            if scales.component_nodes[first_component] != node ||
                    scales.component_nodes[second_component] != node
                # The compact construction normally stores exactly two
                # consecutive rows per node. Retain the generic fallback for
                # malformed/legacy layouts without changing their behavior.
                for component in eachindex(out)
                    out[component] = _omrf_node_box_first_order_bound(
                        scales, component,
                        scales.component_nodes[component])
                end
                return out
            end
            _omrf_node_box_two_bounds!(out, scales, node,
                first_component, second_component)
        end
        return out
    end
    @inbounds for component in eachindex(out)
        out[component] = _omrf_node_box_first_order_bound(
            scales, component, scales.component_nodes[component])
    end
    return out
end

function (scales::OMRFNodeBoxFirstOrderScales)(out, state,
        flow::PDMPSamplers.AnyBoomerang, t)
    started = time_ns()
    _cache_omrf_node_box_point!(scales, state, flow, t)
    _omrf_node_box_fill_bounds!(out, scales)
    scales.context.counts.envelope_point_calls += 1
    scales.context.counts.envelope_point_seconds +=
        (time_ns() - started) * 1e-9
    return out
end

function (scales::OMRFNodeBoxFirstOrderScales)(out, state,
        flow::BouncyParticle, t)
    started = time_ns()
    _cache_omrf_node_box_point!(scales, state, flow, t)
    _omrf_node_box_fill_bounds!(out, scales)
    scales.context.counts.envelope_point_calls += 1
    scales.context.counts.envelope_point_seconds +=
        (time_ns() - started) * 1e-9
    return out
end

function (scales::OMRFNodeBoxFirstOrderScales)(out, state,
        flow::PDMPSamplers.AnyBoomerang, left, right)
    started = time_ns()
    _cache_omrf_node_box_cell!(scales, state, flow, left, right)
    _omrf_node_box_fill_bounds!(out, scales)
    scales.context.counts.envelope_cell_calls += 1
    scales.context.counts.envelope_cell_seconds +=
        (time_ns() - started) * 1e-9
    return out
end

function (scales::OMRFNodeBoxFirstOrderScales)(out, state,
        flow::BouncyParticle, left, right)
    started = time_ns()
    _cache_omrf_node_box_cell!(scales, state, flow, left, right)
    _omrf_node_box_fill_bounds!(out, scales)
    scales.context.counts.envelope_cell_calls += 1
    scales.context.counts.envelope_cell_seconds +=
        (time_ns() - started) * 1e-9
    return out
end

(scales::OMRFNodeBoxFirstOrderScales)(out, state,
    flow::PDMPSamplers.PreconditionedDynamics, args...) =
    scales(out, state, flow.dynamics, args...)

PDMPSamplers._trajectory_scale_anchor(scales::OMRFNodeBoxFirstOrderScales) =
    scales.anchor

function _omrf_cell_marked_cache_slot!(
        scales::OMRFCellMarkedNodeBoxFirstOrderScales,
        left::Float64, right::Float64)
    slot = scales.n_cached_cells + 1
    if slot > length(scales.cell_lefts)
        old_capacity = length(scales.cell_lefts)
        new_capacity = max(slot, max(8, 2old_capacity))
        resize!(scales.cell_lefts, new_capacity)
        resize!(scales.cell_rights, new_capacity)
        resize!(scales.cell_total_masses, new_capacity)
        n_components = length(scales.component_totals)
        old_scales = scales.cell_scales
        old_cumulative = scales.cell_cumulative_masses
        scales.cell_scales = zeros(n_components, new_capacity)
        scales.cell_cumulative_masses = zeros(n_components, new_capacity)
        if old_capacity > 0
            copyto!(view(scales.cell_scales, :, 1:old_capacity), old_scales)
            copyto!(view(scales.cell_cumulative_masses, :, 1:old_capacity),
                old_cumulative)
        end
    end
    scales.cell_lefts[slot] = left
    scales.cell_rights[slot] = right
    scales.n_cached_cells = slot
    return slot
end

@inline function _omrf_cell_marked_lookup(
        scales::OMRFCellMarkedNodeBoxFirstOrderScales, time::Float64)
    # Grid cells are appended in increasing order. Find the last cell whose
    # left endpoint does not exceed the proposal time in logarithmic time;
    # the cache is reset for every new event-search schedule.
    lo = 1
    hi = scales.n_cached_cells
    slot = 0
    @inbounds while lo <= hi
        mid = (lo + hi) >>> 1
        if scales.cell_lefts[mid] <= time
            slot = mid
            lo = mid + 1
        else
            hi = mid - 1
        end
    end
    return slot > 0 && time <= scales.cell_rights[slot] ? slot : 0
end

function (scales::OMRFCellMarkedNodeBoxFirstOrderScales)(out, state,
        flow::PDMPSamplers.AnyBoomerang, left, right)
    iszero(left) && (scales.n_cached_cells = 0)
    scales.base(out, state, flow, left, right)
    slot = _omrf_cell_marked_cache_slot!(
        scales, Float64(left), Float64(right))
    cumulative = 0.0
    @inbounds for component in eachindex(out, scales.component_totals)
        value = out[component]
        cumulative += value * scales.component_totals[component]
        scales.cell_scales[component, slot] = value
        scales.cell_cumulative_masses[component, slot] = cumulative
    end
    scales.cell_total_masses[slot] = cumulative
    return out
end

function (scales::OMRFCellMarkedNodeBoxFirstOrderScales)(out, state,
        flow::BouncyParticle, left, right)
    iszero(left) && (scales.n_cached_cells = 0)
    scales.base(out, state, flow, left, right)
    slot = _omrf_cell_marked_cache_slot!(
        scales, Float64(left), Float64(right))
    cumulative = 0.0
    @inbounds for component in eachindex(out, scales.component_totals)
        value = out[component]
        cumulative += value * scales.component_totals[component]
        scales.cell_scales[component, slot] = value
        scales.cell_cumulative_masses[component, slot] = cumulative
    end
    scales.cell_total_masses[slot] = cumulative
    return out
end

@inline function _omrf_cell_marked_set_point!(scales, state, flow, time,
        harmonic::Bool)
    scales.point_x = state.ξ.x
    scales.point_theta = state.ξ.θ
    state isa StickyPDMPState && (scales.point_free = state.free)
    scales.point_is_sticky = state isa StickyPDMPState
    scales.point_mu = flow.μ
    scales.point_harmonic = harmonic
    scales.cached_time = time
    return nothing
end

function (scales::OMRFCellMarkedNodeBoxFirstOrderScales)(out, state,
        flow::PDMPSamplers.AnyBoomerang, t)
    time = Float64(t)
    slot = _omrf_cell_marked_lookup(scales, time)
    if iszero(slot)
        scales.cell_cache_misses += 1
        scales.base.context.counts.categorical_cell_cache_misses += 1
        throw(ErrorException(
            "cell-marked node-box proposal lies outside the active cell cache"))
    end
    @inbounds for component in eachindex(out)
        out[component] = scales.cell_scales[component, slot]
    end
    _omrf_cell_marked_set_point!(scales, state, flow, time, true)
    scales.active_cell = slot
    scales.cell_cache_hits += 1
    scales.base.context.counts.categorical_cell_cache_hits += 1
    return out
end


function (scales::OMRFCellMarkedNodeBoxFirstOrderScales)(out, state,
        flow::BouncyParticle, t)
    time = Float64(t)
    slot = _omrf_cell_marked_lookup(scales, time)
    if iszero(slot)
        scales.cell_cache_misses += 1
        scales.base.context.counts.categorical_cell_cache_misses += 1
        throw(ErrorException(
            "cell-marked node-box proposal lies outside the active cell cache"))
    end
    @inbounds for component in eachindex(out)
        out[component] = scales.cell_scales[component, slot]
    end
    _omrf_cell_marked_set_point!(scales, state, flow, time, false)
    scales.active_cell = slot
    scales.cell_cache_hits += 1
    scales.base.context.counts.categorical_cell_cache_hits += 1
    return out
end

(scales::OMRFCellMarkedNodeBoxFirstOrderScales)(out, state,
    flow::PDMPSamplers.PreconditionedDynamics, args...) =
    scales(out, state, flow.dynamics, args...)

PDMPSamplers._trajectory_scale_anchor(
    scales::OMRFCellMarkedNodeBoxFirstOrderScales) = scales.base.anchor

function PDMPSamplers.subsampling_failure_context(
        oracle::OMRFAnalyticSubsamplingOracle,
        envelope::PDMPSamplers.GroupedResidualEnvelope{
            <:OMRFCellMarkedNodeBoxFirstOrderScales}, subset)
    scales = envelope.component_scales!
    factor = isempty(subset) ? 0 : subset[1]
    component = iszero(factor) ? 0 : envelope.groups[1, factor]
    context = oracle.residual_context
    count_names = fieldnames(typeof(context.counts))
    count_values = ntuple(i -> getfield(context.counts, count_names[i]),
        length(count_names))
    count_snapshot = NamedTuple{count_names}(count_values)
    N = size(context.X, 1)
    node0, person0 = iszero(factor) ? (-1, -1) : divrem(factor - 1, N)
    slot = scales.active_cell
    return (
        selected_factor=factor,
        selected_factors=copy(subset),
        selected_person=person0 + 1,
        selected_node=node0 + 1,
        selected_component=component,
        component_population=iszero(component) ? 0.0 :
            scales.component_totals[component],
        component_cell_certificate=iszero(component) ? NaN :
            envelope.scales[component],
        selected_point_certificate=scales.selected_factor_point_bound,
        active_cell=slot,
        cell_left=iszero(slot) ? NaN : scales.cell_lefts[slot],
        cell_right=iszero(slot) ? NaN : scales.cell_rights[slot],
        cache_generation=scales.n_cached_cells,
        cached_point_time=scales.cached_time,
        cached_point_position=copy(scales.point_x),
        cached_point_velocity=copy(scales.point_theta),
        cached_point_free=copy(scales.point_free),
        cached_point_mean=copy(scales.point_mu),
        cached_point_harmonic=scales.point_harmonic,
        oracle_anchor=copy(oracle.anchor),
        oracle_anchor_identity=objectid(oracle.anchor),
        envelope_anchor=copy(scales.base.anchor),
        envelope_anchor_identity=objectid(scales.base.anchor),
        envelope_weights=copy(envelope.weights),
        envelope_groups=copy(envelope.groups),
        envelope_totals=copy(envelope.totals),
        live_component_scales=copy(envelope.scales),
        live_component_cumulative_masses=copy(envelope.cumulative_masses),
        component_membership_counts=copy(scales.component_totals),
        cached_cell_component_scales=iszero(slot) ? Float64[] :
            copy(view(scales.cell_scales, :, slot)),
        envelope_wrapper_type=string(typeof(envelope)),
        component_provider_type=string(typeof(scales)),
        anchor_probabilities=copy(context.anchor_probabilities),
        touched_coordinates=copy(context.touched_positions),
        active_generation=context.active_generation,
        oracle_event_history=count_snapshot)
end

function PDMPSamplers.screening_residual_bound(
        envelope::PDMPSamplers.GroupedResidualEnvelope{
            <:OMRFCellMarkedNodeBoxFirstOrderScales}, state, flow, t)
    scales = envelope.component_scales!
    scales(envelope.scales, state, flow, t)
    slot = scales.active_cell
    @inbounds for component in eachindex(envelope.cumulative_masses)
        envelope.cumulative_masses[component] =
            scales.cell_cumulative_masses[component, slot]
    end
    return scales.cell_total_masses[slot]
end

@inline function _omrf_cell_marked_coordinate(scales,
        node::Int, index::Int, sine::Float64, cosine::Float64)
    reference = _omrf_node_reference(scales.base, node, index)
    if scales.point_is_sticky && !scales.point_free[index]
        return scales.point_x[index] - reference, 0.0
    end
    if !scales.point_harmonic
        displacement = scales.point_x[index] - reference +
            scales.cached_time * scales.point_theta[index]
        return displacement, scales.point_theta[index]
    end
    centered = scales.point_x[index] - scales.point_mu[index]
    displacement = scales.point_mu[index] - reference +
        centered * cosine + scales.point_theta[index] * sine
    velocity = -centered * sine + scales.point_theta[index] * cosine
    return displacement, velocity
end

@inline function _omrf_cell_marked_factor_point_bound(
        scales::OMRFCellMarkedNodeBoxFirstOrderScales,
        factor::Int)
    context = scales.base.context
    N = size(context.X, 1)
    node0, person0 = divrem(factor - 1, N)
    node = node0 + 1
    person = person0 + 1
    sine, cosine = scales.point_harmonic ? sincos(scales.cached_time) :
        (0.0, 1.0)
    field_displacement = 0.0
    field_velocity = 0.0
    @inbounds for k in eachindex(context.incident_edges[node])
        covariate = context.X[person, context.incident_neighbours[node][k]]
        iszero(covariate) && continue
        index = context.interaction_indices[context.incident_edges[node][k]]
        displacement, velocity = _omrf_cell_marked_coordinate(
            scales, node, index, sine, cosine)
        field_displacement += covariate * displacement
        field_velocity += covariate * velocity
    end
    displacement_minimum = 0.0
    displacement_maximum = 0.0
    velocity_minimum = 0.0
    velocity_maximum = 0.0
    threshold_start = context.threshold_starts[node]
    q = context.seen[node] - 1
    @inbounds for category in 1:q
        index = context.threshold_indices[threshold_start + category - 1]
        displacement, velocity = _omrf_cell_marked_coordinate(
            scales, node, index, sine, cosine)
        logit_displacement = displacement + category * field_displacement
        logit_velocity = velocity + category * field_velocity
        displacement_minimum = min(displacement_minimum, logit_displacement)
        displacement_maximum = max(displacement_maximum, logit_displacement)
        velocity_minimum = min(velocity_minimum, logit_velocity)
        velocity_maximum = max(velocity_maximum, logit_velocity)
    end
    displacement_range = displacement_maximum - displacement_minimum
    velocity_range = velocity_maximum - velocity_minimum
    old_bound = _omrf_categorical_residual_bound(
        displacement_range, velocity_range, false)
    new_bound = _omrf_categorical_residual_bound(
        displacement_range, velocity_range, true)
    scales.point_displacement_gt_four += displacement_range > 4.0
    scales.point_old_bound_sum += old_bound
    scales.point_new_bound_sum += new_bound
    counts = scales.base.context.counts
    counts.categorical_point_evaluations += 1
    counts.categorical_point_displacement_gt_four += displacement_range > 4.0
    counts.categorical_point_old_bound_sum += old_bound
    counts.categorical_point_new_bound_sum += new_bound
    return scales.base.categorical_residual_cap ? new_bound : old_bound
end

function PDMPSamplers.draw_subset!(rng::Random.AbstractRNG,
        cv::PDMPSamplers.SubsampledControlVariate{<:Any,<:Any,
            <:PDMPSamplers.GroupedResidualEnvelope{
                <:OMRFCellMarkedNodeBoxFirstOrderScales}},
        D::Real, B::Real)
    envelope = cv.envelope
    scales = envelope.component_scales!
    N = PDMPSamplers.n_observations(envelope)
    m = cv.m
    total = D + B
    residual_source = !iszero(B) &&
        (iszero(D) || rand(rng) * total > D)
    if residual_source
        component = PDMPSamplers._draw_component(rng, envelope, B)
        distinguished = PDMPSamplers._draw_distinguished(
            rng, envelope, component)
        if isone(m)
            # Preserve the established batch-one draw stream exactly.
            cv.subset[1] = distinguished
        else
            PDMPSamplers._draw_base_subset_conditional!(
                rng, cv, cv.subset_design, distinguished)
        end
        scales.component_marks += 1
    else
        PDMPSamplers._draw_base_subset!(rng, cv, cv.subset_design)
    end
    PDMPSamplers.record_subsampling_mark_source!(
        cv.residual_oracle, residual_source)
    started = time_ns()
    cell_bound = 0.0
    point_bound = 0.0
    @inbounds for factor in cv.subset
        component = envelope.groups[1, factor]
        factor_cell_bound = envelope.scales[component]
        factor_point_bound = _omrf_cell_marked_factor_point_bound(
            scales, factor)
        factor_point_bound <= factor_cell_bound * (1 + 64eps(Float64)) ||
            throw(ErrorException(
                "selected person-node point bound exceeds its closed-cell certificate"))
        cell_bound += factor_cell_bound
        point_bound += factor_point_bound
    end
    scales.factor_point_evaluations += m
    scales.factor_point_seconds += (time_ns() - started) * 1e-9
    scale = N / m
    scales.selected_factor_point_bound = scale * point_bound
    oracle = cv.residual_oracle
    oracle isa PDMPSamplers.WithResidualStats && (oracle = oracle.f)
    oracle isa OMRFAnalyticSubsamplingOracle &&
        (oracle.lazy_residual_subset_bound = scale * point_bound)
    return D + scale * cell_bound
end

function _cache_omrf_pattern_harmonic_coefficients!(
        scales::Union{OMRFPatternLocalScales,OMRFPatternLocalHCVScales,
            OMRFCovariateLocalScales},
        state, flow)
    if scales isa OMRFCovariateLocalScales
        flow_hash = hash((typeof(flow), hash(flow.Γ)))
        cache_matches = scales.harmonic_cache_valid &&
            scales.harmonic_cache_flow_hash == flow_hash &&
            length(scales.cache_x) == length(state.ξ.x) &&
            length(scales.cache_theta) == length(state.ξ.θ) &&
            all(@inbounds(scales.cache_x[i] == state.ξ.x[i])
                for i in eachindex(state.ξ.x)) &&
            all(@inbounds(scales.cache_theta[i] == state.ξ.θ[i])
                for i in eachindex(state.ξ.θ)) &&
            all(@inbounds(scales.cache_mu[i] == flow.μ[i])
                for i in eachindex(flow.μ)) &&
            all(@inbounds(scales.cache_anchor[i] == scales.anchor[i])
                for i in eachindex(scales.anchor)) &&
            all(@inbounds(scales.cache_free[i] ==
                _omrf_coordinate_free(state, i)) for i in eachindex(scales.anchor))
        if cache_matches
            scales.context.counts.harmonic_coefficient_cache_hits += 1
            return nothing
        end
        scales.context.counts.harmonic_coefficient_preparations += 1
        started = time_ns()
        scales.harmonic_cache_flow_hash = flow_hash
        copyto!(scales.cache_x, state.ξ.x)
        copyto!(scales.cache_theta, state.ξ.θ)
        copyto!(scales.cache_mu, flow.μ)
        copyto!(scales.cache_anchor, scales.anchor)
        @inbounds for i in eachindex(scales.cache_free)
            scales.cache_free[i] = _omrf_coordinate_free(state, i)
        end
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
        scales.harmonic_cache_valid = true
        context.counts.harmonic_coefficient_seconds +=
            (time_ns() - started) * 1e-9
        return nothing
    end
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
    _ensure_omrf_active_pattern!(scales, state)
    context = scales.context
    _cache_omrf_pattern_harmonic_coefficients!(scales, state, flow)
    midpoint_sine, midpoint_cosine = sincos((left + right) / 2)
    half_width = (right - left) / 2
    coordinate_offsets = scales.coordinate_offsets
    coordinate_cosines = scales.coordinate_cosines
    coordinate_sines = scales.coordinate_sines
    component_limit = scales.active_compression ?
        scales.active_component_count : length(scales.component_nodes)
    @inbounds for component in 1:component_limit
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
        _omrf_categorical_hcv_bound(velocity_range, displacement_range) :
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
        _omrf_categorical_hcv_bound(velocity_range, displacement_range) :
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
        _omrf_categorical_hcv_bound(velocity_range, displacement_range) :
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
        _omrf_categorical_hcv_bound(velocity_range, displacement_range) :
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
            _omrf_categorical_hcv_bound(
                velocity_range, displacement_range) :
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
        _omrf_categorical_hcv_bound(velocity_range, displacement_range) :
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
        _omrf_categorical_hcv_bound(velocity_range, displacement_range) :
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
        state, flow::Union{PDMPSamplers.AnyBoomerang,
            ZigZag,BouncyParticle}, t::Real)
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
        state, flow::Union{PDMPSamplers.AnyBoomerang,
            ZigZag,BouncyParticle}, node::Int, anchor)
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

# Anchor-bank selection is called through the installed preconditioned flow,
# whereas the node geometry is expressed in physical BPS/Zig-Zag/Boomerang
# coordinates.  Mirror the existing whole-envelope forwarding methods so the
# research whole-person HCV bank uses the same physical flow as its event
# certificate.
_cache_omrf_node_grouped_point_node!(
    scales::OMRFNodeGroupedPersonHCVBoxScales, state,
    flow::PDMPSamplers.PreconditionedDynamics, node::Int, anchor) =
    _cache_omrf_node_grouped_point_node!(
        scales, state, flow.dynamics, node, anchor)

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
    return _omrf_categorical_hcv_bound(
        velocity_range, displacement_range)
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
    return _omrf_categorical_hcv_bound(
        velocity_range, displacement_range)
end

@inline function _omrf_node_grouped_zigzag_point_bound(
        scales::OMRFNodeGroupedPersonHCVBoxScales,
        component::Int, node::Int, state, flow, anchor)
    context = scales.context
    q = context.seen[node] - 1
    q == 0 && return 0.0
    edges = context.incident_edges[node]
    neighbours = context.incident_neighbours[node]
    maximum_logit_abs = 0.0
    start = context.threshold_starts[node]
    @inbounds for u in 1:q
        index = context.threshold_indices[start + u - 1]
        displacement, _ = _omrf_coordinate_at(
            state, flow, anchor, index, 0.0)
        field_center = 0.0
        field_uncertainty = 0.0
        for k in eachindex(edges)
            edge_index = context.interaction_indices[edges[k]]
            lower = scales.minima[component, neighbours[k]]
            upper = scales.maxima[component, neighbours[k]]
            midpoint = (lower + upper) / 2
            radius = (upper - lower) / 2
            edge_displacement, _ = _omrf_coordinate_at(
                state, flow, anchor, edge_index, 0.0)
            field_center += midpoint * edge_displacement
            field_uncertainty += radius * abs(edge_displacement)
        end
        logit_abs = abs(displacement + u * field_center) +
            u * field_uncertainty
        maximum_logit_abs = max(maximum_logit_abs, logit_abs)
    end
    weighted_velocity = 0.0
    @inbounds for u in 1:q
        index = context.threshold_indices[start + u - 1]
        _, velocity = _omrf_coordinate_at(
            state, flow, anchor, index, 0.0)
        weighted_velocity += abs(velocity)
    end
    @inbounds for k in eachindex(edges)
        edge_index = context.interaction_indices[edges[k]]
        _, velocity = _omrf_coordinate_at(
            state, flow, anchor, edge_index, 0.0)
        lower = scales.minima[component, neighbours[k]]
        upper = scales.maxima[component, neighbours[k]]
        weighted_velocity += q * abs(velocity) * max(abs(lower), abs(upper))
    end
    # Every represented logit lies in [-maximum_logit_abs,
    # maximum_logit_abs], so twice that value encloses the categorical range.
    return _omrf_categorical_hcv_bound(
        weighted_velocity, 2maximum_logit_abs)
end

@inline function _omrf_node_grouped_person_zigzag_point_bound(
        scales::OMRFNodeGroupedPersonHCVBoxScales,
        person::Int, node::Int, state, flow, anchor)
    context = scales.context
    q = context.seen[node] - 1
    q == 0 && return 0.0
    edges = context.incident_edges[node]
    neighbours = context.incident_neighbours[node]
    minimum_logit = 0.0
    maximum_logit = 0.0
    start = context.threshold_starts[node]
    @inbounds for u in 1:q
        index = context.threshold_indices[start + u - 1]
        displacement, _ = _omrf_coordinate_at(
            state, flow, anchor, index, 0.0)
        field = 0.0
        for k in eachindex(edges)
            edge_index = context.interaction_indices[edges[k]]
            edge_displacement, _ = _omrf_coordinate_at(
                state, flow, anchor, edge_index, 0.0)
            field += context.X[person, neighbours[k]] * edge_displacement
        end
        logit = displacement + u * field
        minimum_logit = min(minimum_logit, logit)
        maximum_logit = max(maximum_logit, logit)
    end
    weighted_velocity = 0.0
    @inbounds for u in 1:q
        index = context.threshold_indices[start + u - 1]
        _, velocity = _omrf_coordinate_at(
            state, flow, anchor, index, 0.0)
        weighted_velocity += abs(velocity)
    end
    @inbounds for k in eachindex(edges)
        edge_index = context.interaction_indices[edges[k]]
        _, velocity = _omrf_coordinate_at(
            state, flow, anchor, edge_index, 0.0)
        weighted_velocity += q * abs(velocity) *
            abs(context.X[person, neighbours[k]])
    end
    return _omrf_categorical_hcv_bound(
        weighted_velocity, maximum_logit - minimum_logit)
end

@inline function _omrf_node_grouped_zigzag_cell_bound(
        scales::OMRFNodeGroupedPersonHCVBoxScales,
        component::Int, node::Int, state, flow, anchor,
        left::Real, right::Real)
    context = scales.context
    q = context.seen[node] - 1
    q == 0 && return 0.0
    edges = context.incident_edges[node]
    neighbours = context.incident_neighbours[node]
    maximum_logit_abs = 0.0
    start = context.threshold_starts[node]
    @inbounds for u in 1:q
        index = context.threshold_indices[start + u - 1]
        left_displacement, _ = _omrf_coordinate_at(
            state, flow, anchor, index, left)
        right_displacement, _ = _omrf_coordinate_at(
            state, flow, anchor, index, right)
        left_field = 0.0
        right_field = 0.0
        field_uncertainty = 0.0
        for k in eachindex(edges)
            edge_index = context.interaction_indices[edges[k]]
            left_edge, _ = _omrf_coordinate_at(
                state, flow, anchor, edge_index, left)
            right_edge, _ = _omrf_coordinate_at(
                state, flow, anchor, edge_index, right)
            lower = scales.minima[component, neighbours[k]]
            upper = scales.maxima[component, neighbours[k]]
            midpoint = (lower + upper) / 2
            radius = (upper - lower) / 2
            left_field += midpoint * left_edge
            right_field += midpoint * right_edge
            field_uncertainty += radius * max(abs(left_edge), abs(right_edge))
        end
        left_abs = abs(left_displacement + u * left_field) +
            u * field_uncertainty
        right_abs = abs(right_displacement + u * right_field) +
            u * field_uncertainty
        maximum_logit_abs = max(maximum_logit_abs, left_abs, right_abs)
    end
    weighted_velocity = 0.0
    @inbounds for u in 1:q
        index = context.threshold_indices[start + u - 1]
        _, velocity = _omrf_coordinate_at(
            state, flow, anchor, index, left)
        weighted_velocity += abs(velocity)
    end
    @inbounds for k in eachindex(edges)
        edge_index = context.interaction_indices[edges[k]]
        _, velocity = _omrf_coordinate_at(
            state, flow, anchor, edge_index, left)
        lower = scales.minima[component, neighbours[k]]
        upper = scales.maxima[component, neighbours[k]]
        weighted_velocity += q * abs(velocity) * max(abs(lower), abs(upper))
    end
    return _omrf_categorical_hcv_bound(
        weighted_velocity, 2maximum_logit_abs)
end

@inline function _omrf_node_grouped_person_point_bound(
        scales::OMRFNodeGroupedPersonHCVBoxScales,
        person::Int, node::Int, state, flow, anchor)
    if flow isa ZigZag
        node_anchor = scales.node_anchors === nothing ? scales.anchor :
            view(scales.node_anchors, :, scales.node_sources[node])
        return _omrf_node_grouped_person_zigzag_point_bound(
            scales, person, node, state, flow, node_anchor)
    end
    return _omrf_node_grouped_person_point_bound(scales, person, node)
end

function (scales::OMRFNodeGroupedPersonHCVBoxScales)(out, state,
        flow::Union{PDMPSamplers.AnyBoomerang,BouncyParticle}, t)
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

function (scales::OMRFNodeGroupedPersonHCVBoxScales)(out, state,
        flow::ZigZag, t)
    started = time_ns()
    @inbounds for component in eachindex(out)
        node = scales.component_nodes[component]
        anchor = scales.node_anchors === nothing ? scales.anchor :
            view(scales.node_anchors, :, scales.node_sources[node])
        out[component] = _omrf_node_grouped_zigzag_point_bound(
            scales, component, node, state, flow, anchor)
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
    return _omrf_categorical_hcv_bound(
        velocity_range, displacement_range)
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

function _omrf_node_grouped_linear_cell_bound(
        scales::OMRFNodeGroupedPersonHCVBoxScales,
        component::Int, node::Int, state, flow, left::Real, right::Real)
    context = scales.context
    q = context.seen[node] - 1
    q == 0 && return 0.0
    anchor = scales.node_anchors === nothing ? scales.anchor :
        view(scales.node_anchors, :, scales.node_sources[node])
    if flow isa ZigZag
        return _omrf_node_grouped_zigzag_cell_bound(
            scales, component, node, state, flow, anchor, left, right)
    end
    edges = context.incident_edges[node]
    neighbours = context.incident_neighbours[node]
    displacement_range = 0.0
    velocity_range = 0.0
    start = context.threshold_starts[node]
    @inbounds for u in 1:q, v in 0:(u - 1)
        u_index = context.threshold_indices[start + u - 1]
        left_displacement, left_velocity = _omrf_coordinate_at(
            state, flow, anchor, u_index, left)
        right_displacement, _ = _omrf_coordinate_at(
            state, flow, anchor, u_index, right)
        if iszero(v)
            left_pair = left_displacement
            right_pair = right_displacement
            pair_velocity = left_velocity
        else
            v_index = context.threshold_indices[start + v - 1]
            v_left, v_velocity = _omrf_coordinate_at(
                state, flow, anchor, v_index, left)
            v_right, _ = _omrf_coordinate_at(
                state, flow, anchor, v_index, right)
            left_pair = left_displacement - v_left
            right_pair = right_displacement - v_right
            pair_velocity = left_velocity - v_velocity
        end
        displacement_uncertainty = 0.0
        velocity_uncertainty = 0.0
        multiplier = u - v
        for k in eachindex(edges)
            index = context.interaction_indices[edges[k]]
            left_edge, edge_velocity = _omrf_coordinate_at(
                state, flow, anchor, index, left)
            right_edge, _ = _omrf_coordinate_at(
                state, flow, anchor, index, right)
            lower = scales.minima[component, neighbours[k]]
            upper = scales.maxima[component, neighbours[k]]
            midpoint = (lower + upper) / 2
            radius = (upper - lower) / 2
            left_pair += multiplier * midpoint * left_edge
            right_pair += multiplier * midpoint * right_edge
            pair_velocity += multiplier * midpoint * edge_velocity
            displacement_uncertainty += multiplier * radius *
                max(abs(left_edge), abs(right_edge))
            velocity_uncertainty += multiplier * radius * abs(edge_velocity)
        end
        displacement_range = max(displacement_range,
            abs(left_pair) + displacement_uncertainty,
            abs(right_pair) + displacement_uncertainty)
        velocity_range = max(velocity_range,
            abs(pair_velocity) + velocity_uncertainty)
    end
    return _omrf_categorical_hcv_bound(
        velocity_range, displacement_range)
end

function (scales::OMRFNodeGroupedPersonHCVBoxScales)(out, state,
        flow::BouncyParticle, left, right)
    started = time_ns()
    @inbounds for component in eachindex(out)
        node = scales.component_nodes[component]
        out[component] = _omrf_node_grouped_linear_cell_bound(
            scales, component, node, state, flow, left, right)
    end
    scales.context.counts.envelope_cell_calls += 1
    scales.context.counts.envelope_cell_seconds +=
        (time_ns() - started) * 1e-9
    return out
end

function (scales::OMRFNodeGroupedPersonHCVBoxScales)(out, state,
        flow::ZigZag, left, right)
    started = time_ns()
    @inbounds for component in eachindex(out)
        node = scales.component_nodes[component]
        anchor = scales.node_anchors === nothing ? scales.anchor :
            view(scales.node_anchors, :, scales.node_sources[node])
        out[component] = _omrf_node_grouped_zigzag_cell_bound(
            scales, component, node, state, flow, anchor, left, right)
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
        flow::Union{PDMPSamplers.AnyBoomerang,BouncyParticle,ZigZag}, left, right)
    iszero(left) && (scales.n_cached_cells = 0)
    scales.base(out, state, flow, left, right)
    slot = _omrf_cell_cache_slot!(
        scales, Float64(left), Float64(right), length(out))
    copyto!(view(scales.cell_scales, :, slot), out)
    return out
end

function (scales::OMRFCellCachedNodeGroupedPersonHCVBoxScales)(out, state,
        flow::Union{PDMPSamplers.AnyBoomerang,BouncyParticle,ZigZag}, t)
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
                scales.base, observation, partition, scales.cached_state,
                base_flow, scales.base.anchor)
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

@inline function _omrf_lazy_component_point_bound(
        scales::OMRFLazyNodeGroupedPersonHCVBoxScales,
        component::Int, flow)
    node = scales.base.component_nodes[component]
    anchor = scales.base.node_anchors === nothing ? scales.base.anchor :
        view(scales.base.node_anchors, :, scales.base.node_sources[node])
    if flow isa ZigZag
        return _omrf_node_grouped_zigzag_point_bound(
            scales.base, component, node, scales.cached_state,
            flow, anchor)
    end
    return _omrf_node_grouped_point_bound(scales.base, component, node)
end

function (scales::OMRFLazyNodeGroupedPersonHCVBoxScales)(out, state,
        flow::Union{PDMPSamplers.AnyBoomerang,BouncyParticle,ZigZag}, t)
    started = time_ns()
    _cache_omrf_lazy_point!(scales, state, flow, t)
    @inbounds for root_index in eachindex(scales.roots, out)
        out[root_index] = _omrf_lazy_component_point_bound(
            scales, scales.roots[root_index], _omrf_lazy_base_flow(flow))
    end
    scales.base.context.counts.envelope_point_calls += 1
    scales.base.context.counts.envelope_point_seconds +=
        (time_ns() - started) * 1e-9
    return out
end

function (scales::OMRFLazyNodeGroupedPersonHCVBoxScales)(out, state,
        flow::Union{PDMPSamplers.AnyBoomerang,BouncyParticle,ZigZag}, left, right)
    started = time_ns()
    midpoint = (left + right) / 2
    midpoint_sine, midpoint_cosine = sincos(midpoint)
    half_width = (right - left) / 2
    if flow isa Union{BouncyParticle,ZigZag}
        scales.base(scales.leaf_scales, state, flow, left, right)
        @inbounds for root_index in eachindex(scales.roots, out)
            out[root_index] = scales.leaf_scales[scales.roots[root_index]]
        end
    else
        _cache_omrf_node_grouped_cell!(scales.base, state, flow,
            midpoint_sine, midpoint_cosine, half_width)
        @inbounds for root_index in eachindex(scales.roots, out)
            component = scales.roots[root_index]
            node = scales.base.component_nodes[component]
            out[root_index] = _omrf_node_grouped_cell_bound(
                scales.base, component, node, midpoint_sine, midpoint_cosine,
                half_width)
        end
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
    base_flow = _omrf_lazy_base_flow(scales.cached_flow)
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
            left_scale = _omrf_lazy_component_point_bound(
                scales, left, base_flow)
            right_scale = _omrf_lazy_component_point_bound(
                scales, right, base_flow)
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
            value = _omrf_lazy_component_point_bound(
                scales, component, base_flow)
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

function _cache_omrf_batched_linear_fields!(
        scales::OMRFBatchedFactorHCVScales, state, flow::ZigZag)
    context = scales.context
    fill!(scales.interaction_offsets, 0.0)
    fill!(scales.interaction_cosines, 0.0)
    P = length(context.seen)
    edge = 1
    @inbounds for left in 1:(P - 1), right in (left + 1):P
        index = context.interaction_indices[edge]
        left_anchor = _omrf_factor_anchor(scales, index, left)
        right_anchor = _omrf_factor_anchor(scales, index, right)
        left_delta, velocity = _omrf_coordinate_at_anchor_value(
            state, flow, left_anchor, index, 0.0)
        right_delta, _ = _omrf_coordinate_at_anchor_value(
            state, flow, right_anchor, index, 0.0)
        scales.interaction_offsets[right, left] = left_delta
        scales.interaction_offsets[left, right] = right_delta
        scales.interaction_cosines[left, right] = velocity
        scales.interaction_cosines[right, left] = velocity
        # Reuse the otherwise harmonic-only sine buffers to cache the
        # coordinatewise Zig-Zag absolute edge-velocity contribution.  This
        # removes an O(N P^2) incident-edge traversal from every point/cell
        # certificate without changing the bound or summation over factors.
        scales.interaction_sines[left, right] = abs(velocity)
        scales.interaction_sines[right, left] = abs(velocity)
        edge += 1
    end
    mul!(scales.field_offsets, context.X_float,
        scales.interaction_offsets)
    mul!(scales.field_cosines, context.X_float,
        scales.interaction_cosines)
    mul!(scales.field_sines, context.X_float,
        scales.interaction_sines)
    return nothing
end

@inline function _omrf_batched_zigzag_hcv_scale(
        scales::OMRFBatchedFactorHCVScales, state, flow::ZigZag,
        person::Int, node::Int, t::Float64)
    context = scales.context
    q = context.seen[node] - 1
    iszero(q) && return 0.0
    field_delta = scales.field_offsets[person, node] +
        t * scales.field_cosines[person, node]
    minimum_delta = 0.0
    maximum_delta = 0.0
    coordinate_velocity_range_sum = 0.0
    start = context.threshold_starts[node]
    @inbounds for u in 1:q
        index = context.threshold_indices[start + u - 1]
        anchor = _omrf_factor_anchor(scales, index, node)
        delta, velocity = _omrf_coordinate_at_anchor_value(
            state, flow, anchor, index, t)
        logit_delta = delta + u * field_delta
        minimum_delta = min(minimum_delta, logit_delta)
        maximum_delta = max(maximum_delta, logit_delta)
        # Coordinatewise Zig-Zag rates cannot use cancellation in B*v.
        # Sum the range of each coordinate's induced logit direction.
        coordinate_velocity_range_sum += abs(velocity)
    end
    coordinate_velocity_range_sum += q * scales.field_sines[person, node]
    return context.hcv_active ? _omrf_categorical_hcv_bound(
        coordinate_velocity_range_sum, maximum_delta - minimum_delta) :
        coordinate_velocity_range_sum * (maximum_delta - minimum_delta) / 4
end

function (scales::OMRFBatchedFactorHCVScales)(out, state,
        flow::ZigZag, t)
    started = time_ns()
    _cache_omrf_batched_linear_fields!(scales, state, flow)
    destination = 1
    N, P = size(scales.context.X)
    @inbounds for node in 1:P, person in 1:N
        out[destination] = _omrf_batched_zigzag_hcv_scale(
            scales, state, flow, person, node, Float64(t))
        destination += 1
    end
    scales.context.counts.event_bound_observation_visits += N * P
    scales.context.counts.envelope_point_calls += 1
    scales.context.counts.envelope_point_seconds +=
        (time_ns() - started) * 1e-9
    return out
end

function (scales::OMRFBatchedFactorHCVScales)(out, state,
        flow::ZigZag, left, right)
    started = time_ns()
    _cache_omrf_batched_linear_fields!(scales, state, flow)
    destination = 1
    N, P = size(scales.context.X)
    @inbounds for node in 1:P, person in 1:N
        left_scale = _omrf_batched_zigzag_hcv_scale(
            scales, state, flow, person, node, Float64(left))
        right_scale = _omrf_batched_zigzag_hcv_scale(
            scales, state, flow, person, node, Float64(right))
        # Along a linear stratum every pairwise logit displacement is affine;
        # its absolute maximum, and hence Delta, occurs at a cell endpoint.
        # V_ZZ is constant on the stratum.  The capped bound is monotone in
        # Delta, so the larger endpoint certificate encloses the whole cell.
        out[destination] = max(left_scale, right_scale)
        destination += 1
    end
    scales.context.counts.event_bound_observation_visits += 2N * P
    scales.context.counts.envelope_cell_calls += 1
    scales.context.counts.envelope_cell_seconds +=
        (time_ns() - started) * 1e-9
    return out
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
                _omrf_categorical_hcv_bound(
                    velocity_range, displacement_range) :
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
                _omrf_categorical_hcv_bound(
                    velocity_range, displacement_range) :
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
    if scales.context.hcv_nodewise
        @inbounds for node in 1:P
            displacement, velocity = _omrf_node_geometry(
                state, flow, scales.anchor, scales.supports[node], t)
            displacement2 = displacement * displacement
            α = isinf(scales.context.hcv_damping) ? 1.0 :
                scales.context.hcv_damping /
                (scales.context.hcv_damping + displacement2)
            _fill_omrf_damped_hcv_scales!(out, node, P,
                displacement, velocity, α)
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
    if scales.context.hcv_nodewise
        damping = scales.context.hcv_damping
        @inbounds for node in 1:P
            displacement, velocity = _omrf_node_cell_geometry(
                state, flow, scales.anchor, scales.supports[node], left, right)
            displacement2 = displacement * displacement
            one_minus_α_bound = isinf(damping) ? 0.0 :
                displacement2 / (damping + displacement2)
            out[node] = velocity * displacement * one_minus_α_bound
            out[P + node] = velocity * displacement2
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
    # Compact norm-local specifications intentionally omit the generic dense
    # P-by-(NP) matrices. Dispatch on their specialized payload first; legacy
    # specifications carrying those matrices reach exactly the same branch.
    if bound_type in ("node_box_first_order",
            "cell_marked_node_box_first_order")
        residual_context.factorization === :person_node ||
            throw(ArgumentError(
                "node-box first-order geometry requires person-node factorization"))
        !residual_context.use_hcv && !residual_context.hcv_active ||
            throw(ArgumentError(
                "node-box first-order geometry requires an ordinary first-order CV"))
        groups_by_node = Matrix{Int}(
            _rget(envelope_spec, :node_box_membership))
        N, P = size(residual_context.X)
        size(groups_by_node) == (P, N) || throw(DimensionMismatch(
            "node-box membership must be a P-by-N matrix"))
        component_nodes = _as_int_vector(
            _rget(envelope_spec, :node_box_nodes))
        minima = _as_float_matrix(_rget(envelope_spec, :node_box_minima))
        maxima = _as_float_matrix(_rget(envelope_spec, :node_box_maxima))
        n_components = length(component_nodes)
        size(minima) == size(maxima) == (n_components, P) ||
            throw(DimensionMismatch(
                "node-box extrema must contain one P-vector per component"))
        all(component -> 1 <= component <= n_components, groups_by_node) ||
            throw(ArgumentError("node-box membership is out of range"))
        all(minima .<= maxima) || throw(ArgumentError(
            "node-box minima must not exceed maxima"))
        if bound_type == "cell_marked_node_box_first_order"
            n_components == 2P || throw(ArgumentError(
                "cell-marked node-box geometry currently requires G=2"))
            exact_refinement = _haskey(envelope_spec,
                :exact_selected_refinement) && Bool(_rget(
                    envelope_spec, :exact_selected_refinement))
            exact_refinement || throw(ArgumentError(
                "cell-marked node-box geometry requires exact selected-factor refinement"))
        end
        @inbounds for node in 1:P, person in 1:N
            component_nodes[groups_by_node[node, person]] == node ||
                throw(ArgumentError(
                    "node-box membership refers to a component for another node"))
        end
        # Person-node factor indices are node-major: all N people for node 1,
        # then all N people for node 2. One grouped partition therefore
        # represents exactly NP factors without a dense P-by-NP matrix.
        groups = reshape(vec(permutedims(groups_by_node)), 1, N * P)
        d = length(anchor)
        categorical_residual_cap = _haskey(envelope_spec,
            :categorical_residual_cap) && Bool(_rget(
                envelope_spec, :categorical_residual_cap))
        scales = OMRFNodeBoxFirstOrderScales(copy(anchor), component_nodes,
            minima, maxima, residual_context, zeros(d), zeros(d), zeros(d),
            zeros(d), nothing, Int[], UInt64(0), categorical_residual_cap)
        if bound_type == "cell_marked_node_box_first_order"
            component_totals = zeros(n_components)
            @inbounds for factor in axes(groups, 2)
                component_totals[groups[1, factor]] += 1.0
            end
            # A schedule normally uses only a handful of cells. Preallocating
            # 128 keeps ordinary cell construction allocation-free; unusual
            # larger schedules grow geometrically during construction.
            capacity = 128
            cell_marked = OMRFCellMarkedNodeBoxFirstOrderScales(
                scales, component_totals, zeros(capacity), zeros(capacity),
                zeros(n_components, capacity),
                zeros(n_components, capacity), zeros(capacity), 0,
                0.0, 0, zeros(d), zeros(d), trues(d), false, zeros(d), false, NaN,
                0, 0, 0.0, 0, 0, 0, 0.0, 0.0)
            return PDMPSamplers.GroupedResidualEnvelope(
                groups, n_components, cell_marked;
                component_cell_scales! = cell_marked)
        end
        return PDMPSamplers.GroupedResidualEnvelope(
            groups, n_components, scales; component_cell_scales! = scales)
    end
    if bound_type == "norm_local_hcv_compact"
        residual_context.factorization === :person_node ||
            throw(ArgumentError(
                "compact HCV norm-local geometry requires person-node factorization"))
        residual_context.use_hcv && residual_context.hcv_active &&
            isinf(residual_context.hcv_damping) || throw(ArgumentError(
                "compact HCV norm-local geometry requires active undamped HCV"))
        residual_context.hcv_nodewise && throw(ArgumentError(
            "compact HCV norm-local geometry requires a global HCV"))
        weights = _as_float_matrix(_rget(envelope_spec, :norm_weights))
        N, P = size(residual_context.X)
        size(weights) == (4P, N) || throw(DimensionMismatch(
            "compact HCV norm-local weights must be a 4P-by-N matrix"))
        scales = OMRFCompactHCVNormLocalScales(anchor, residual_context)
        component_blocks = repeat(collect(1:P), 4)
        return PDMPSamplers.BlockSeparableResidualEnvelope(
            weights, component_blocks, P, scales;
            component_cell_scales! = scales)
    end
    if bound_type == "node_local_norm_low_rank" &&
            _haskey(envelope_spec, :norm_weights) &&
            !residual_context.use_hcv
        weights = _as_float_matrix(_rget(envelope_spec, :norm_weights))
        N, P = size(residual_context.X)
        size(weights) == (3P, N) || throw(DimensionMismatch(
            "OMRF norm-local weights must be a 3P-by-N matrix"))
        scales = OMRFNormLocalScales(anchor, residual_context)
        if residual_context.factorization === :person_node
            component_blocks = repeat(collect(1:P), 3)
            return PDMPSamplers.BlockSeparableResidualEnvelope(
                weights, component_blocks, P, scales;
                component_cell_scales! = scales)
        end
        return SeparableResidualEnvelope(weights, scales;
            component_cell_scales! = scales)
    end
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
            packed_metadata = _build_omrf_pattern_packed_metadata(
                component_nodes, component_covariates, residual_context)
            workspace_length = maximum(residual_context.seen)
            active_compression = residual_context.active_pattern_compression
            hash_capacity = 1
            while hash_capacity < 2 * size(residual_context.X, 1)
                hash_capacity <<= 1
            end
            scales = OMRFPatternLocalScales(anchor, component_nodes,
                component_covariates, packed_metadata..., residual_context,
                active_compression, n_components,
                active_compression ? -1 : 0,
                active_compression ? trues(length(anchor)) : BitVector(),
                false, nothing, nothing, nothing,
                active_compression ? zeros(UInt64, hash_capacity) :
                    zeros(UInt64, 1),
                active_compression ? zeros(Int, hash_capacity) : zeros(Int, 1),
                active_compression ? zeros(Int, hash_capacity) : zeros(Int, 1),
                0,
                active_compression ? zeros(Int, size(residual_context.X, 1)) :
                    Int[],
                zeros(workspace_length), zeros(workspace_length),
                zeros(workspace_length), zeros(length(anchor)),
                zeros(length(anchor)), zeros(length(anchor)),
                zeros(length(anchor)), zeros(length(anchor)),
                zeros(length(anchor)))
            envelope = PDMPSamplers.GroupedResidualEnvelope(
                groups, n_components, scales;
                component_cell_scales! = scales)
            if active_compression
                _attach_omrf_active_pattern!(envelope, scales)
                # Build the initial active partition once at construction so
                # the envelope is valid before the first sampler callback.
                # Subsequent sticky transitions use the generation check in
                # `_ensure_omrf_active_pattern!` and rebuild lazily at the
                # next point/cell query.
                _rebuild_omrf_active_pattern!(scales)
            end
            return envelope
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
                zeros(length(anchor)), zeros(length(anchor)),
                zeros(length(anchor)), zeros(length(anchor)), zeros(length(anchor)),
                zeros(length(anchor)), falses(length(anchor)), UInt(0), false)
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
        if bound_type == "covariate_local_compact" &&
                residual_context.factorization === :person_node &&
                !residual_context.use_hcv
            N = size(residual_context.X, 1)
            P = size(residual_context.X, 2)
            signed_bps = _haskey(envelope_spec, :signed_bps) &&
                Bool(_rget(envelope_spec, :signed_bps))
            signed_harmonic = _haskey(envelope_spec, :signed_harmonic) &&
                Bool(_rget(envelope_spec, :signed_harmonic))
            return _new_omrf_compact_envelope(anchor, residual_context;
                signed_bps=signed_bps, signed_harmonic=signed_harmonic)
        end
        if bound_type in ("direct_person_node_hcv",
                "direct_factor_roof_person_hcv")
            expected_factorization = bound_type == "direct_person_node_hcv" ?
                :person_node : :person
            residual_context.factorization === expected_factorization ||
                throw(ArgumentError("$bound_type requires " *
                    "$(String(expected_factorization)) factorization"))
            residual_context.use_hcv && residual_context.hcv_active &&
                isinf(residual_context.hcv_damping) || throw(ArgumentError(
                    "direct HCV geometry requires active undamped HCV"))
            residual_context.hcv_nodewise && throw(ArgumentError(
                "direct HCV geometry requires a global HCV"))
            person_signed_roof = _haskey(envelope_spec,
                :person_signed_roof) &&
                Bool(_rget(envelope_spec, :person_signed_roof))
            signed_only_person_roof = _haskey(envelope_spec,
                :signed_only_person_roof) &&
                Bool(_rget(envelope_spec, :signed_only_person_roof))
            person_signed_roof && expected_factorization !== :person &&
                throw(ArgumentError(
                    "signed person roofs require person factorization"))
            cell_scope = _haskey(envelope_spec, :direct_cv_cell_scope) ?
                String(_rget(envelope_spec, :direct_cv_cell_scope)) : "cell"
            cell_scope in ("cell", "event_horizon") || throw(ArgumentError(
                "direct HCV cell scope must be cell or event_horizon"))
            cell_scope == "event_horizon" &&
                bound_type != "direct_factor_roof_person_hcv" &&
                throw(ArgumentError(
                    "HCV event-horizon cells require whole-person sampling"))
            return _new_omrf_direct_hcv_envelope(anchor, residual_context;
                person_signed_roof=person_signed_roof,
                signed_only_person_roof=signed_only_person_roof,
                joint_group_size=_haskey(envelope_spec,
                    :joint_signed_group_size) ?
                    Int(_rget(envelope_spec, :joint_signed_group_size)) : 0,
                shared_current=_haskey(envelope_spec,
                    :shared_current_categorical) &&
                    Bool(_rget(envelope_spec, :shared_current_categorical)),
                event_horizon_scope=cell_scope == "event_horizon")
        end
        if bound_type == "direct_person_node_cv"
            residual_context.factorization === :person_node ||
                throw(ArgumentError(
                    "direct first-order CV geometry requires person-node factorization"))
            !residual_context.use_hcv && !residual_context.hcv_active ||
                throw(ArgumentError(
                    "direct first-order CV geometry does not use a Hessian CV"))
            cell_scope = _haskey(envelope_spec, :direct_cv_cell_scope) ?
                String(_rget(envelope_spec, :direct_cv_cell_scope)) : "cell"
            cell_scope in ("cell", "event_horizon") || throw(ArgumentError(
                "direct CV cell scope must be cell or event_horizon"))
            return _new_omrf_direct_hcv_envelope(anchor, residual_context;
                hcv_remainder=false,
                event_horizon_scope=cell_scope == "event_horizon")
        end
        if bound_type == "norm_local_hcv_compact"
            residual_context.factorization === :person_node ||
                throw(ArgumentError(
                    "compact HCV norm-local geometry requires person-node factorization"))
            residual_context.use_hcv && residual_context.hcv_active &&
                isinf(residual_context.hcv_damping) || throw(ArgumentError(
                    "compact HCV norm-local geometry requires active undamped HCV"))
            residual_context.hcv_nodewise && throw(ArgumentError(
                "compact HCV norm-local geometry requires a global HCV"))
            weights = _as_float_matrix(_rget(envelope_spec, :norm_weights))
            P = size(residual_context.X, 2)
            size(weights) == (4P, size(residual_context.X, 1)) ||
                throw(DimensionMismatch(
                    "compact HCV norm-local weights must contain four rows per node"))
            scales = OMRFCompactHCVNormLocalScales(anchor, residual_context)
            component_blocks = repeat(collect(1:P), 4)
            return PDMPSamplers.BlockSeparableResidualEnvelope(
                weights, component_blocks, P, scales;
                component_cell_scales! = scales)
        end
        if _haskey(envelope_spec, :norm_weights) &&
                !residual_context.use_hcv
            weights = _as_float_matrix(_rget(envelope_spec, :norm_weights))
            size(weights, 1) == 3 * size(residual_context.X, 2) ||
                throw(DimensionMismatch(
                    "OMRF norm-local weights must have three rows per node"))
            scales = OMRFNormLocalScales(anchor, residual_context)
            if residual_context.factorization === :person_node
                P = size(residual_context.X, 2)
                component_blocks = vcat(collect(1:P), collect(1:P),
                    collect(1:P))
                return PDMPSamplers.BlockSeparableResidualEnvelope(
                    weights, component_blocks, P, scales;
                    component_cell_scales! = scales)
            end
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
    warmup_finish_event_bound_observation_visits::Int
    selection_seconds::Float64
    selection_component_evaluations::Int
    node_source_check_distance::Float64
    node_source_last_position::Vector{Float64}
    node_source_checks::Int
    node_source_skips::Int
    node_source_changes::Int
    window_diagnostics::Any
    selection_policy::Symbol
    pending_boundary_idx::Int
    pending_boundary_time::Float64
    boundary_generation::UInt64
    pending_boundary_generation::UInt64
    selection_boundaries::Int
    node_position_policy_version::Int
    warmup_bank_learning_policy::Symbol
    warmup_bank_candidate_count::Int
    warmup_bank_training_state_count::Int
    warmup_bank_objective_single::Float64
    warmup_bank_objective_selected::Float64
    warmup_bank_learning_seconds::Float64
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
    spec = manager.envelope_spec
    bound_type = _haskey(spec, :bound_type) ?
        String(_rget(spec, :bound_type)) : ""
    if bound_type in ("direct_person_node_hcv",
            "direct_factor_roof_person_hcv", "direct_person_node_cv") &&
            (bound_type != "direct_person_node_cv") == residual_context.use_hcv
        event_horizon_scope =
            bound_type in ("direct_person_node_cv",
                           "direct_factor_roof_person_hcv") &&
            _haskey(spec, :direct_cv_cell_scope) &&
            String(_rget(spec, :direct_cv_cell_scope)) == "event_horizon"
        return _new_omrf_direct_hcv_envelope(requested, residual_context;
            hcv_remainder=bound_type != "direct_person_node_cv",
            event_horizon_scope=event_horizon_scope,
            person_signed_roof=_haskey(spec, :person_signed_roof) &&
                Bool(_rget(spec, :person_signed_roof)),
            signed_only_person_roof=
                _haskey(spec, :signed_only_person_roof) &&
                Bool(_rget(spec, :signed_only_person_roof)),
            joint_group_size=_haskey(spec, :joint_signed_group_size) ?
                Int(_rget(spec, :joint_signed_group_size)) : 0,
            shared_current=_haskey(spec, :shared_current_categorical) &&
                Bool(_rget(spec, :shared_current_categorical)))
    end
    if bound_type == "norm_local_hcv_compact" &&
            residual_context.use_hcv && !isempty(manager.entries)
        template = manager.entries[1].state.envelope
        if template isa PDMPSamplers.BlockSeparableResidualEnvelope
            P = size(residual_context.X, 2)
            scales = OMRFCompactHCVNormLocalScales(
                requested, residual_context)
            return PDMPSamplers.BlockSeparableResidualEnvelope(
                template.weights, template.component_blocks,
                template.block_components, template.n_blocks, scales, scales,
                template.totals, template.alias_tables,
                zeros(length(template.scales)),
                zeros(length(template.cell_scales)),
                zeros(length(template.cumulative_masses)))
        end
    end
    if bound_type == "covariate_local_expansion" &&
            !residual_context.use_hcv && !isempty(manager.entries)
        template = manager.entries[1].state.envelope
        if template isa Union{SeparableResidualEnvelope,
                PDMPSamplers.BlockSeparableResidualEnvelope}
            max_degree = maximum(length,
                residual_context.incident_edges; init=0)
            scales = OMRFCovariateLocalScales(requested, residual_context,
                zeros(max_degree), zeros(max_degree), zeros(length(requested)),
                zeros(length(requested)), zeros(length(requested)),
                zeros(length(requested)), zeros(length(requested)),
                zeros(length(requested)), zeros(length(requested)),
                falses(length(requested)), UInt(0), false)
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
    if bound_type == "covariate_local_compact" &&
            !residual_context.use_hcv && !isempty(manager.entries)
        template = manager.entries[1].state.envelope
        if template isa OMRFCompactCovariateLocalEnvelope
            static = template.component_scales!.static
            signed_bps = _haskey(spec, :signed_bps) &&
                Bool(_rget(spec, :signed_bps))
            signed_harmonic = _haskey(spec, :signed_harmonic) &&
                Bool(_rget(spec, :signed_harmonic))
            return _new_omrf_compact_envelope(
                requested, residual_context, static;
                signed_bps=signed_bps, signed_harmonic=signed_harmonic)
        end
    end
    return _build_omrf_residual_envelope(spec, requested, residual_context)
end

function _prepare_omrf_anchor(manager::OMRFAnchorManager,
        anchor::AbstractVector)
    started = time_ns()
    oracle = manager.oracle
    recycled = manager.recycled_state
    compact_hcv_recycle = recycled !== nothing &&
        manager.hcv_active && recycled.residual_context.use_hcv &&
        ((recycled.envelope isa PDMPSamplers.BlockSeparableResidualEnvelope &&
          recycled.envelope.component_scales! isa
            OMRFCompactHCVNormLocalScales) ||
         recycled.envelope isa OMRFDirectPersonNodeHCVEnvelope)
    if recycled !== nothing && oracle.analytic_prior !== nothing &&
            ((!manager.hcv_active && !recycled.residual_context.use_hcv) ||
             compact_hcv_recycle) &&
            recycled.envelope isa Union{SeparableResidualEnvelope,
                PDMPSamplers.GroupedResidualEnvelope,
                PDMPSamplers.BlockSeparableResidualEnvelope,
                OMRFCompactCovariateLocalEnvelope,
                OMRFDirectPersonNodeHCVEnvelope} &&
            (recycled.envelope isa OMRFDirectPersonNodeHCVEnvelope ||
             recycled.envelope.component_scales! isa Union{
                OMRFNodeLocalScales,
                OMRFNormLocalScales,
                OMRFNodeBoxFirstOrderScales,
                OMRFCellMarkedNodeBoxFirstOrderScales,
                OMRFCovariateLocalScales,
                OMRFCompactCovariateLocalScales,
                OMRFCompactHCVNormLocalScales})
        # The recycled state is not installed in the bank, so it can be
        # mutated without exposing a partially prepared anchor. Its envelope
        # callbacks retain the same anchor/context array identities.
        manager.recycled_state = nothing
        copyto!(recycled.anchor, anchor)
        scales = recycled.envelope.component_scales!
        anchor_storage = recycled.envelope isa
                OMRFDirectPersonNodeHCVEnvelope ?
            recycled.envelope.anchor :
            (scales isa OMRFCompactCovariateLocalScales ?
                scales.base.anchor :
             scales isa OMRFCellMarkedNodeBoxFirstOrderScales ?
                scales.base.anchor : scales.anchor)
        anchor_storage === recycled.anchor || copyto!(anchor_storage, anchor)
        if scales isa OMRFCellMarkedNodeBoxFirstOrderScales
            scales.n_cached_cells = 0
            scales.active_cell = 0
            scales.selected_factor_point_bound = NaN
        end
        recycled.residual_context.hcv_active = manager.hcv_active &&
            recycled.residual_context.use_hcv
        likelihood_started = time_ns()
        _omrf_prior_gradient!(
            recycled.prior_anchor, oracle.analytic_prior, recycled.anchor)
        _cache_omrf_anchor_and_full_gradient!(recycled.residual_context,
            recycled.full_anchor, oracle.analytic_prior, recycled.anchor)
        if compact_hcv_recycle
            hessian = recycled.residual_context.likelihood_hessian::Matrix{Float64}
            _omrf_likelihood_hessian!(hessian, recycled.residual_context)
            if recycled.envelope isa OMRFDirectPersonNodeHCVEnvelope
                recycled.envelope.cell_count = 0
                recycled.envelope.active_cell = 0
                recycled.envelope.path_valid = false
            end
        end
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
    if state.residual_context.active_pattern_compression ||
            state.residual_context.direct_person_node_hcv
        copyto!(state.residual_context.active_free, oracle.active_free)
        state.residual_context.active_generation = oracle.active_generation
    end
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
    bank_learning = _haskey(envelope_spec, :hcv_bank_learning) ?
        Symbol(String(_rget(envelope_spec, :hcv_bank_learning))) : :geometric
    bank_learning in (:geometric, :exact_variance,
            :selector_coupled_variance) || throw(ArgumentError(
        "unknown OMRF HCV warmup-bank learning policy: $bank_learning"))
    manager = OMRFAnchorManager(oracle, envelope_spec, String.(unc_names),
        OMRFAnchorEntry[], 0, Int(capacity), zeros(length(oracle.full_anchor)),
        ones(length(oracle.full_anchor)),
        0, 0.0, 0, 0, 0, 0, 0.0, 0.0, Inf,
        oracle.residual_context.hcv_active, nothing, 0, 0, 0.0, 0.0,
        oracle.ctx.counts.native_full_data_scans,
        oracle.ctx.counts.event_bound_observation_visits, 0.0, 0,
        0.0, copy(oracle.anchor), 0, 0, 0, nothing,
        :adaptive, 0, Inf, UInt64(0), UInt64(0), 0, 0,
        bank_learning, 0, 0, NaN, NaN, 0.0,
        Float64[], zeros(0, 0), Float64[], :scaled, false)
    initial = OMRFPreparedAnchorState(oracle.anchor, oracle.full_anchor,
        oracle.prior_anchor, oracle.residual_context, initial_envelope)
    manager.active_idx = _insert_omrf_anchor!(manager, initial)
    return manager
end

function _install_omrf_node_position_sources!(manager::OMRFAnchorManager,
        cv::SubsampledControlVariate, requested::AbstractVector{<:Integer})
    nodecv = manager.oracle.node_position_cv::OMRFNodePositionCV
    scales = cv.envelope.component_scales!
    base = scales isa OMRFCellMarkedNodeBoxFirstOrderScales ? scales.base :
        (scales isa OMRFNodeBoxFirstOrderScales ? scales : nothing)
    base === nothing && throw(ArgumentError(
        "node-position requires the first-order node-box certificate"))
    sources_changed = false
    @inbounds for node in eachindex(nodecv.sources, requested)
        sources_changed |= nodecv.sources[node] != requested[node]
    end
    installed_sources_match = length(base.node_sources) == length(requested)
    if installed_sources_match
        @inbounds for node in eachindex(base.node_sources, requested)
            installed_sources_match &= base.node_sources[node] == requested[node]
        end
    end
    authoritative_installation = base.node_bank === nodecv.bank &&
        base.node_installation_generation == nodecv.installation_generation &&
        installed_sources_match
    if !sources_changed && authoritative_installation
        return false
    end
    sources_changed && (nodecv.installation_generation += UInt64(1))
    copyto!(nodecv.sources, requested)
    _omrf_node_local_baseline!(nodecv.likelihood_baseline,
        nodecv.bank, nodecv.sources)
    # Keep the legacy full/prior fields coherent for the direct deterministic
    # certificate, but never interpret this sum as one global anchor.
    @. manager.oracle.full_anchor = nodecv.likelihood_baseline +
        manager.oracle.prior_anchor
    manager.oracle.prior_x_valid = false
    manager.oracle.cached_residual_valid = false
    base.node_bank === nodecv.bank || (base.node_bank = nodecv.bank)
    resize!(base.node_sources, length(nodecv.sources))
    copyto!(base.node_sources, nodecv.sources)
    base.node_installation_generation = nodecv.installation_generation
    base.context = manager.oracle.residual_context
    if scales isa OMRFCellMarkedNodeBoxFirstOrderScales
        scales.n_cached_cells = 0
        scales.active_cell = 0
        scales.selected_factor_point_bound = NaN
    end
    # `true` means downstream proposal caches must be invalidated.  This is
    # also required for same-index installation into a replacement wrapper.
    return true
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

@inline function _omrf_hcv_velocity_variance(flow, coordinate::Int)
    sigma = PDMPSamplers._physical_velocity_sd(flow, coordinate)
    isfinite(sigma) && sigma > 0 || throw(DomainError(sigma,
        "HCV bank learning requires finite positive velocity scales"))
    return sigma * sigma
end

"""Exact full-stratum HCV residual variance under the installed velocity law.

The score is the trace of the finite-population residual covariance against
the diagonal physical-velocity covariance.  Constants depending only on batch
size and population size are omitted because they cannot change the selected
bank.  This routine is warmup-only and intentionally enumerates the complete
factor population.
"""
function _omrf_hcv_expected_directional_variance!(
        context::OMRFResidualContext, anchor, x, velocity_variances,
        work, residual_mean, singleton)
    context.use_hcv && context.hcv_active && isinf(context.hcv_damping) ||
        throw(ArgumentError(
            "exact-variance bank learning requires active undamped HCV"))
    N, P = size(context.X)
    population = context.factorization === :person ? N : N * P
    @. context.displacement = x - anchor
    fill!(residual_mean, 0.0)
    sum_metric_norm2 = 0.0
    @inbounds for factor in 1:population
        singleton[1] = factor
        omrf_residual!(context, work, x, singleton)
        for coordinate in eachindex(work, velocity_variances)
            value = work[coordinate]
            residual_mean[coordinate] += value
            sum_metric_norm2 += velocity_variances[coordinate] * value * value
        end
    end
    inverse_population = inv(Float64(population))
    centered = sum_metric_norm2 * inverse_population
    @inbounds for coordinate in eachindex(residual_mean, velocity_variances)
        mean_value = residual_mean[coordinate] * inverse_population
        centered -= velocity_variances[coordinate] * mean_value * mean_value
    end
    tolerance = 256eps(Float64) * max(sum_metric_norm2 * inverse_population, 1.0)
    centered >= -tolerance || throw(ErrorException(
        "negative HCV finite-population variance: $centered"))
    return max(0.0, centered)
end

"""Choose a frozen HCV bank from this chain's own warmup trajectory.

Candidate positions and training positions are the same bounded time-uniform
warmup grid plus the authoritative endpoint.  Every candidate is scored with
the exact factor residual covariance under the final adapted velocity metric.
The endpoint is installed first; subsequent anchors greedily minimize the
time-uniform total of the best available score.  No retained state is used.
"""
function _omrf_selector_assigned_objective(costs, candidates, scales,
        selected)
    total = 0.0
    @inbounds for state_index in axes(candidates, 2)
        best_slot = 1
        best_distance = Inf
        for slot in eachindex(selected)
            candidate = selected[slot]
            distance = 0.0
            for coordinate in axes(candidates, 1)
                delta = (candidates[coordinate, state_index] -
                    candidates[coordinate, candidate]) / scales[coordinate]
                distance += delta * delta
            end
            if distance < best_distance
                best_slot = slot
                best_distance = distance
            end
        end
        total += costs[selected[best_slot], state_index]
    end
    return total / size(candidates, 2)
end

"""Exact subset search for a bank evaluated through the retained selector.

Candidate one is the authoritative endpoint and is always retained. Candidate
order is chronological, so the ordinary lowest-index tie convention is
preserved. This is warmup-only work over at most twice the requested bank size.
"""
function _omrf_selector_coupled_indices(costs, candidates, scales,
        capacity::Integer)
    n = size(candidates, 2)
    k = min(Int(capacity), n)
    k == 1 && return [1], _omrf_selector_assigned_objective(
        costs, candidates, scales, [1])
    k == n && return collect(1:n), _omrf_selector_assigned_objective(
        costs, candidates, scales, collect(1:n))
    tail = collect(2:k)
    selected = Vector{Int}(undef, k)
    selected[1] = 1
    best = zeros(Int, k)
    best_objective = Inf
    while true
        copyto!(selected, 2, tail, 1, length(tail))
        objective = _omrf_selector_assigned_objective(
            costs, candidates, scales, selected)
        if objective < best_objective
            best_objective = objective
            copyto!(best, selected)
        end
        position = length(tail)
        while position >= 1 && tail[position] == n - length(tail) + position
            position -= 1
        end
        position == 0 && break
        tail[position] += 1
        if position < length(tail)
            for j in (position + 1):length(tail)
                tail[j] = tail[j - 1] + 1
            end
        end
    end
    return best, best_objective
end

function _omrf_exact_variance_anchors!(manager::OMRFAnchorManager,
        trace::PDMPSamplers.AbstractPDMPTrace, position::AbstractVector,
        capacity::Integer, flow; selector_coupled::Bool=false)
    started = time_ns()
    # Twice the requested capacity is enough to expose dispersed warmup
    # regions without turning finalization into a second simulation phase.
    candidate_capacity = max(1, 2 * Int(capacity))
    candidates = _omrf_warmup_candidates(
        trace, position, candidate_capacity; retained_fraction=0.5,
        n_grid=candidate_capacity)
    d, n = size(candidates)
    velocity_variances = [_omrf_hcv_velocity_variance(flow, coordinate)
        for coordinate in 1:d]
    costs = Matrix{Float64}(undef, n, n)
    work = zeros(d)
    residual_mean = zeros(d)
    singleton = [1]

    # Reuse one uninstalled prepared state while constructing the score
    # matrix. Selected anchors are prepared afresh below and become the only
    # durable bank entries.
    manager.recycled_state = nothing
    temporary = nothing
    @inbounds for anchor_index in 1:n
        temporary === nothing || (manager.recycled_state = temporary)
        temporary = _prepare_omrf_anchor(
            manager, view(candidates, :, anchor_index))
        context = temporary.residual_context
        context.hcv_active = true
        for state_index in 1:n
            costs[anchor_index, state_index] =
                _omrf_hcv_expected_directional_variance!(context,
                    temporary.anchor, view(candidates, :, state_index),
                    velocity_variances,
                    work, residual_mean, singleton)
        end
    end
    manager.recycled_state = temporary

    # Preserve the established positive distance scales for the retained
    # position selector. Only bank learning changes in this policy.
    means = vec(sum(candidates; dims=2)) ./ n
    scales = zeros(d)
    if n > 1
        @inbounds for column in axes(candidates, 2), coordinate in 1:d
            delta = candidates[coordinate, column] - means[coordinate]
            scales[coordinate] += delta * delta
        end
        scales .= sqrt.(scales ./ (n - 1))
    end
    positive = filter(scale -> scale > sqrt(eps(Float64)), scales)
    scale_floor = isempty(positive) ? 1.0 :
        0.1 * Statistics.median(positive)
    @. scales = max(scales, scale_floor)
    objective_single = sum(view(costs, 1, :)) / n
    k = min(Int(capacity), n)
    selected_indices = Vector{Int}(undef, k)
    selected_objective = NaN
    if selector_coupled
        selected_indices, selected_objective =
            _omrf_selector_coupled_indices(
                costs, candidates, scales, capacity)
    else
        # Historical research comparator: oracle assignment by minimum exact
        # variance. This is deliberately distinct from retained selection.
        selected_indices[1] = 1
        minimum_cost = copy(view(costs, 1, :))
        available = trues(n)
        available[1] = false
        @inbounds for selection in 2:k
            best_candidate = 0
            best_objective = Inf
            for candidate in 1:n
                available[candidate] || continue
                objective = 0.0
                for state_index in 1:n
                    objective += min(minimum_cost[state_index],
                        costs[candidate, state_index])
                end
                if objective < best_objective
                    best_candidate = candidate
                    best_objective = objective
                end
            end
            best_candidate > 0 || throw(ErrorException(
                "exact-variance HCV bank selection exhausted its candidates"))
            selected_indices[selection] = best_candidate
            available[best_candidate] = false
            for state_index in 1:n
                minimum_cost[state_index] = min(minimum_cost[state_index],
                    costs[best_candidate, state_index])
            end
        end
        selected_objective = sum(minimum_cost) / n
    end
    selected = candidates[:, selected_indices]

    manager.warmup_bank_candidate_count = n
    manager.warmup_bank_training_state_count = n
    manager.warmup_bank_objective_single = objective_single
    manager.warmup_bank_objective_selected = selected_objective
    manager.warmup_bank_learning_seconds += (time_ns() - started) * 1e-9
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
    N, P = size(bank.contexts[1].X)
    K = size(bank.candidates, 2)
    0 < retained_fraction <= 1 || throw(ArgumentError(
        "warmup scoring retained fraction must lie in (0, 1]"))
    n_grid = min(128, max(8K, 32))
    pdmp_trace = trace isa PDMPSamplers.StreamingPDMPTrace ? nothing :
        (trace isa PDMPSamplers.PDMPTrace ? trace :
            PDMPSamplers.PDMPTrace(trace))
    pdmp_trace !== nothing && length(pdmp_trace.times) < 2 &&
        return _select_omrf_node_local_sources!(manager, cv, state, flow)
    streaming_positions, streaming_velocities =
        trace isa PDMPSamplers.StreamingPDMPTrace ?
        PDMPSamplers.streaming_time_uniform_states(
            trace, retained_fraction, n_grid) : (nothing, nothing)
    final_time = trace isa PDMPSamplers.StreamingPDMPTrace ?
        trace.terminal_state.t : last(pdmp_trace.times)
    initial_time = trace isa PDMPSamplers.StreamingPDMPTrace ?
        trace.initial_state.t : first(pdmp_trace.times)
    grid_start = final_time - retained_fraction * (final_time - initial_time)
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
        if trace isa PDMPSamplers.StreamingPDMPTrace
            copyto!(position, view(streaming_positions, :, grid_index))
            copyto!(velocity, view(streaming_velocities, :, grid_index))
        else
            n_events = length(pdmp_trace.times)
            idx = clamp(searchsortedlast(pdmp_trace.times, time), 1, n_events)
            copyto!(position, view(pdmp_trace.positions, :, idx))
            copyto!(velocity, view(pdmp_trace.velocities, :, idx))
            if idx < n_events && time > pdmp_trace.times[idx]
                PDMPSamplers.move_forward_time!(
                    scoring_state.ξ, time - pdmp_trace.times[idx], flow)
            end
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
                            node_scales, person, node, scoring_state, flow,
                            anchor)
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
                        node_scales, persons[pattern], node, state, flow,
                        anchor)
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
    # The position policy owns a fixed bank explicitly.  The legacy `Inf`
    # refresh sentinel remains meaningful only for the adaptive comparator.
    frozen_bank = manager.selection_policy in (:position, :node_position) ||
        (isinf(manager.main_refresh_distance) && manager.capacity > 1)
    trace = trace_manager === nothing ? nothing :
        PDMPSamplers.get_warmup_trace(trace_manager)
    if frozen_bank && trace !== nothing &&
            PDMPSamplers._has_integrable_segment(trace)
        representatives, scales = if manager.warmup_bank_learning_policy in
                (:exact_variance, :selector_coupled_variance)
            manager.hcv_active || throw(ArgumentError(
                "exact-variance HCV bank learning requires HCV during warmup"))
            flow === nothing && throw(ArgumentError(
                "exact-variance HCV bank learning requires the final warmup flow"))
            _omrf_exact_variance_anchors!(
                manager, trace, position, manager.capacity, flow;
                selector_coupled=manager.warmup_bank_learning_policy ===
                    :selector_coupled_variance)
        else
            _omrf_representative_anchors(
                trace, position, manager.capacity)
        end
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
        if manager.selection_policy === :node_position
            manager.active_idx = 1
        else
            _select_omrf_anchor!(manager, cv, position; phase=:warmup)
        end
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
    if manager.selection_policy === :node_position
        manager.hcv_active && throw(ArgumentError(
            "node-position is an ordinary first-order CV policy; HCV is disabled"))
        candidates = reduce(hcat,
            (entry.state.anchor for entry in manager.entries))
        manager.oracle.node_position_cv = _prepare_omrf_node_position_cv(
            manager.envelope_spec, manager.unc_names, candidates, flow,
            manager.oracle.ctx.counts)
        manager.node_position_policy_version = 1
        nodecv = manager.oracle.node_position_cv
        _omrf_node_position_sources!(nodecv.requested_sources, nodecv, position)
        _install_omrf_node_position_sources!(
            manager, cv, nodecv.requested_sources)
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
    envelope_anchor = PDMPSamplers._trajectory_scale_anchor(cv.envelope)
    needs_sync = manager.selection_policy === :node_position ? false :
        cv.anchor != active_state.anchor ||
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
    manager.warmup_finish_event_bound_observation_visits =
        manager.oracle.ctx.counts.event_bound_observation_visits
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
    policy in ("position", "node_position") || throw(ArgumentError(
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
    if policy == "node_position"
        Int(_rget(spec, :policy_version)) == 1 || throw(ArgumentError(
            "unsupported node-position mathematical policy version"))
        metric_variances = _as_float_vector(
            _rget(spec, :physical_velocity_variances))
        manager.oracle.node_position_cv = _prepare_omrf_node_position_cv(
            manager.envelope_spec, manager.unc_names, anchors,
            metric_variances, manager.oracle.ctx.counts)
        manager.node_position_policy_version = Int(_rget(spec, :policy_version))
        nodecv = manager.oracle.node_position_cv
        saved_supports = [Int.(support)
            for support in _rget(spec, :node_support_order)]
        saved_supports == nodecv.supports || throw(ArgumentError(
            "restored node-position support ordering does not match the provider"))
        _omrf_node_position_sources!(
            nodecv.requested_sources, nodecv, position)
        _install_omrf_node_position_sources!(
            manager, cv, nodecv.requested_sources)
        manager.active_idx = 1
    else
        manager.active_idx = _omrf_position_anchor_index(manager, position)
        PDMPSamplers.refresh_anchor_owned!(cv,
            manager.entries[manager.active_idx].state.anchor)
    end
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

@inline function _omrf_node_position_score_slope(nodecv::OMRFNodePositionCV,
        node::Int, source::Int, velocity)
    support = nodecv.supports[node]
    terms = nodecv.linear_terms[node]
    slope = 0.0
    @inbounds for local_index in eachindex(support)
        slope -= 2velocity[support[local_index]] * terms[local_index, source]
    end
    return slope
end

function _omrf_node_position_outgoing_sources!(out,
        nodecv::OMRFNodePositionCV, x, velocity)
    K = size(nodecv.bank.candidates, 2)
    @inbounds for node in eachindex(out)
        isempty(nodecv.supports[node]) && (out[node] = 1; continue)
        minimum_score = Inf
        for source in 1:K
            minimum_score = min(minimum_score,
                _omrf_node_position_score(nodecv, node, source, x))
        end
        best = 1
        best_slope = Inf
        for source in 1:K
            score = _omrf_node_position_score(nodecv, node, source, x)
            tolerance = 64eps(Float64) * max(1.0, abs(score),
                abs(minimum_score))
            abs(score - minimum_score) <= tolerance || continue
            slope = _omrf_node_position_score_slope(
                nodecv, node, source, velocity)
            if slope < best_slope || (slope == best_slope && source < best)
                best = source
                best_slope = slope
            end
        end
        out[node] = best
    end
    return out
end

function _omrf_node_position_horizon_impl!(manager::OMRFAnchorManager,
        state)
    nodecv = manager.oracle.node_position_cv::OMRFNodePositionCV
    manager.boundary_generation += UInt64(1)
    manager.pending_boundary_idx = 0
    manager.pending_boundary_time = Inf
    manager.pending_boundary_generation = manager.boundary_generation
    x, velocity = state.ξ.x, state.ξ.θ
    _omrf_node_position_outgoing_sources!(
        nodecv.pending_sources, nodecv, x, velocity)
    if nodecv.pending_sources != nodecv.sources
        manager.pending_boundary_idx = 1
        manager.pending_boundary_time = state.t[]
        return 0.0
    end
    earliest = Inf
    K = size(nodecv.bank.candidates, 2)
    @inbounds for node in eachindex(nodecv.sources)
        current = nodecv.sources[node]
        isempty(nodecv.supports[node]) && continue
        current_score = _omrf_node_position_score(nodecv, node, current, x)
        current_slope = _omrf_node_position_score_slope(
            nodecv, node, current, velocity)
        for competitor in 1:K
            competitor == current && continue
            alpha = _omrf_node_position_score(
                nodecv, node, competitor, x) - current_score
            beta = _omrf_node_position_score_slope(
                nodecv, node, competitor, velocity) - current_slope
            tolerance = 64eps(Float64) * max(1.0, abs(alpha))
            if alpha >= -tolerance && beta < 0.0
                crossing = max(0.0, alpha / -beta)
                earliest = min(earliest, crossing)
            end
        end
    end
    manager.pending_boundary_idx = isfinite(earliest) ? 1 : 0
    manager.pending_boundary_time = isfinite(earliest) ?
        state.t[] + earliest : Inf
    return earliest
end

function _omrf_node_position_anchor_horizon!(manager::OMRFAnchorManager,
        state, flow, phase::Symbol)
    manager.selection_policy === :node_position || return Inf
    phase === :main || return Inf
    _omrf_actual_flow(flow) isa BouncyParticle || throw(ArgumentError(
        "node-position anchor banks require linear BPS dynamics"))
    started = time_ns()
    horizon = _omrf_node_position_horizon_impl!(manager, state)
    manager.selection_seconds += (time_ns() - started) * 1e-9
    nodecv = manager.oracle.node_position_cv::OMRFNodePositionCV
    manager.selection_component_evaluations +=
        length(nodecv.sources) * size(nodecv.bank.candidates, 2)
    return horizon
end

function _select_omrf_node_position_core!(manager::OMRFAnchorManager,
        cv::SubsampledControlVariate, x, phase::Symbol,
        state, boundary_won::Bool)
    nodecv = manager.oracle.node_position_cv::OMRFNodePositionCV
    owns_boundary = boundary_won && phase === :main &&
        manager.pending_boundary_idx > 0 && state !== nothing &&
        manager.pending_boundary_generation == manager.boundary_generation &&
        state.t[] == manager.pending_boundary_time
    started = time_ns()
    if owns_boundary
        _omrf_node_position_outgoing_sources!(
            nodecv.requested_sources, nodecv, x, state.ξ.θ)
        manager.selection_boundaries += 1
    else
        _omrf_node_position_sources!(nodecv.requested_sources, nodecv, x)
    end
    changes = 0
    @inbounds for node in eachindex(nodecv.sources)
        changes += nodecv.sources[node] != nodecv.requested_sources[node]
    end
    changed = _install_omrf_node_position_sources!(
        manager, cv, nodecv.requested_sources)
    manager.selection_seconds += (time_ns() - started) * 1e-9
    manager.selection_component_evaluations +=
        length(nodecv.sources) * size(nodecv.bank.candidates, 2)
    manager.pending_boundary_idx = 0
    manager.pending_boundary_time = Inf
    manager.pending_boundary_generation = UInt64(0)
    if phase === :main
        manager.main_selections += 1
        changed && (manager.main_activations += 1)
        manager.node_source_checks += 1
        changed && (manager.node_source_changes += changes)
    end
    return changed
end

function _select_omrf_node_position!(manager::OMRFAnchorManager,
        cv::SubsampledControlVariate, x; phase::Symbol=:unknown,
        state=nothing, boundary_won::Bool=false)
    _attach_omrf_live_envelope_context!(manager, cv)
    return _select_omrf_node_position_core!(
        manager, cv, x, phase, state, boundary_won)
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
    if cv.envelope isa OMRFDirectPersonNodeHCVEnvelope
        if cv.envelope.context !== authoritative
            cv.envelope.context = authoritative
            return true
        end
        return false
    end
    scales = cv.envelope.component_scales!
    base = scales isa OMRFCellMarkedNodeBoxFirstOrderScales ? scales.base :
        (scales isa OMRFNodeBoxFirstOrderScales ? scales : nothing)
    if base !== nothing && base.context !== authoritative
        base.context = authoritative
        return true
    end
    if scales isa Union{OMRFBatchedFactorHCVScales,
            OMRFPatternLocalHCVScales} &&
            scales.context !== authoritative
        scales.context = authoritative
        return true
    end
    return false
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
    if manager.selection_policy === :node_position
        return _select_omrf_node_position_core!(
            manager, cv, x, phase, state, boundary_won)
    end
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
    PDMPSamplers.active_prior_grad!(target.slab, target.all_buf, x,
        target.all_active)
    PDMPSamplers.active_prior_grad!(target.slab, target.active_buf, x,
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
    if oracle.residual_context.hcv_nodewise
        oracle(out, x, direction)
    else
        _omrf_prior_hvp!(out, something(oracle.analytic_prior), x, direction)
    end
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

function _build_stan_subsampling_model(ctx::StanSubsamplingContext,
        unc_names, subsampling, anchor, slab_prior, can_stick;
        anchor_capacity::Integer=0, manage_during_warmup::Bool=true)
    N = Int(_rget(subsampling, :n_observations))
    m = Int(_rget(subsampling, :subsample_size))
    envelope_spec = _rget(subsampling, :residual_envelope)
    is_omrf = _haskey(envelope_spec, :type) &&
        String(_rget(envelope_spec, :type)) == "omrf"
    backend = is_omrf && _haskey(envelope_spec, :backend) ?
        String(_rget(envelope_spec, :backend)) : "stan"
    backend in ("analytic", "stan") || throw(ArgumentError(
        "unknown OMRF residual backend: $backend"))
    # Generic and Stan-backed envelopes retain the matrix-column contract.
    # Analytic OMRF validates its contribution population from X and the
    # factorization below, without materializing an unused generic matrix.
    if !is_omrf || backend != "analytic"
        weights = _as_float_matrix(_rget(envelope_spec, :weights))
        size(weights, 2) == N || throw(DimensionMismatch(
            "residual envelope must have N observation columns"))
    end
    d = ctx.full.d
    length(anchor) == d || throw(DimensionMismatch("subsampling anchor has the wrong dimension"))
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
    deterministic_hvp = if oracle isa OMRFAnalyticSubsamplingOracle &&
            oracle.analytic_prior !== nothing
        if slab_prior === nothing
            oracle.residual_context.hcv_nodewise ? oracle : nothing
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
        (live_grad, state, flow, phase, stats, event_type,
                anchor_boundary_won) -> begin
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
                phase, state, flow, boundary_won=anchor_boundary_won,
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
            _finish_omrf_warmup!(manager, active_cv, state.ξ.x, trace_manager;
                state, flow)
        end,
        (live_grad, state, flow, phase, stats) -> begin
            active_cv = live_grad isa PDMPSamplers.SubsampledControlVariate ?
                live_grad : cv
            # The horizon callback is invoked before event search, including
            # the first retained proposal after `with_stats` or restoration.
            _attach_omrf_live_envelope_context!(manager, active_cv)
            min(_omrf_position_anchor_horizon!(manager, state, flow, phase),
                _omrf_node_position_anchor_horizon!(
                    manager, state, flow, phase))
        end,
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

"""Enumerate the batch-one accepted marked BPS generator at one complete state."""
function r_omrf_bps_matched_generator(prepared, subsampling, position,
        velocity, free, flow_cov, anchor, model_prior, slab_prior, can_stick)
    x = _as_float_vector(position)
    v = _as_float_vector(velocity)
    free_vec = BitVector(_as_bool_vector(free))
    can_stick_vec = _as_bool_vector(can_stick)
    covariance = _as_float_matrix(flow_cov)
    names = prepared.unc_names
    d = length(x)
    length(v) == d == length(free_vec) == length(names) ||
        throw(DimensionMismatch("matched BPS state dimensions disagree"))
    size(covariance) == (d, d) || throw(DimensionMismatch(
        "matched BPS covariance has the wrong dimension"))
    m = Int(_rget(subsampling, :subsample_size))
    m == 1 || throw(ArgumentError("matched BPS enumeration requires batch one"))
    N = Int(_rget(subsampling, :n_observations))
    scale_idx = findall(name -> name == "global_scale" ||
        startswith(name, "node_scale."), names)
    structural_idx = setdiff(collect(1:d), scale_idx)

    function reflected_increment(g)
        gm = copy(g)
        @inbounds for i in eachindex(gm, free_vec)
            free_vec[i] || (gm[i] = 0.0)
        end
        rate = max(0.0, dot(v, gm))
        Σg = covariance * gm
        denominator = dot(gm, Σg)
        dv = zeros(d)
        if denominator > 0
            dv .= (-2 * dot(v, gm) / denominator) .* Σg
            @inbounds for i in eachindex(dv, free_vec)
                free_vec[i] || (dv[i] = 0.0)
            end
        end
        num_scale = dot(view(v, scale_idx), view(gm, scale_idx))
        num_struct = dot(view(v, structural_idx), view(gm, structural_idx))
        den_scale = dot(view(gm, scale_idx), view(Σg, scale_idx))
        den_struct = dot(view(gm, structural_idx), view(Σg, structural_idx))
        return rate, dv, num_scale, num_struct, den_scale, den_struct
    end

    spec = _rget(subsampling, :residual_envelope)
    full_model = _build_sticky_target(
        build_omrf_full_model(spec, names, d), model_prior, slab_prior,
        can_stick_vec, names)
    PDMPSamplers.set_active_set!(full_model, free_vec)
    g_full = zeros(d)
    PDMPSamplers.compute_gradient_for_reflection!(full_model.grad, x, g_full)
    full_rate, full_dv, full_ns, full_nt, full_ds, full_dt =
        reflected_increment(g_full)

    function enumerate_anchor(anchor_value)
        built = _build_stan_subsampling_model(first(prepared.contexts), names,
            subsampling, _as_float_vector(anchor_value), slab_prior,
            can_stick_vec; anchor_capacity=0)
        model = first(built)
        PDMPSamplers.set_active_set!(model, free_vec)
        cv = model.grad
        deterministic = zeros(d)
        PDMPSamplers.deterministic_gradient!(deterministic, cv, x)
        residual = zeros(d)
        g = similar(deterministic)
        b = zeros(length(scale_idx))
        Q = zeros(length(scale_idx))
        intensity = 0.0
        weighted_ns = weighted_nt = weighted_ds = weighted_dt = 0.0
        max_scale_residual = 0.0
        for mark in 1:N
            fill!(residual, 0.0)
            cv.residual_oracle(residual, x, [mark], cv.anchor)
            max_scale_residual = max(max_scale_residual,
                maximum(abs, view(residual, scale_idx); init=0.0))
            @. g = deterministic + N * residual
            rate, dv, ns, nt, ds, dt = reflected_increment(g)
            a = rate / N
            intensity += a
            @inbounds for (j, index) in pairs(scale_idx)
                b[j] += a * dv[index]
                Q[j] += a * dv[index]^2
            end
            weighted_ns += a * ns
            weighted_nt += a * nt
            weighted_ds += a * ds
            weighted_dt += a * dt
        end
        return (intensity=intensity, scale_drift=b, scale_squared_jump=Q,
            weighted_scale_numerator=weighted_ns,
            weighted_structural_numerator=weighted_nt,
            weighted_scale_denominator=weighted_ds,
            weighted_structural_denominator=weighted_dt,
            max_selected_scale_residual=max_scale_residual)
    end

    current = enumerate_anchor(anchor)
    anchor_current = enumerate_anchor(x)
    return (parameter_names=names, scale_indices=scale_idx,
        full=(intensity=full_rate,
            scale_drift=full_rate .* full_dv[scale_idx],
            scale_squared_jump=full_rate .* abs2.(full_dv[scale_idx]),
            weighted_scale_numerator=full_rate * full_ns,
            weighted_structural_numerator=full_rate * full_nt,
            weighted_scale_denominator=full_rate * full_ds,
            weighted_structural_denominator=full_rate * full_dt),
        subsampled=current, anchor_current=anchor_current)
end

"""Uniform-base-subset fixed-state intensity audit for BPS/Boomerang.

This samples subsets from the base law, not from the size-biased proposal.
The returned rate and canonical jump-intensity Monte Carlo errors therefore
refer directly to the marked event law in the minibatch specification.
"""
function r_omrf_minibatch_fixed_state(prepared, subsampling, position,
        velocity, free, flow_cov, flow_mean, flow_type::String, anchor,
        model_prior, slab_prior, can_stick, batches_raw::AbstractVector,
        draws::Integer, seed::Integer)
    x = _as_float_vector(position); v = _as_float_vector(velocity)
    free_vec = BitVector(_as_bool_vector(free))
    covariance = _as_float_matrix(flow_cov)
    mean = _as_float_vector(flow_mean)
    can_stick_vec = _as_bool_vector(can_stick)
    names = prepared.unc_names; d = length(x)
    length(v) == d == length(free_vec) == length(mean) == length(names) ||
        throw(DimensionMismatch("fixed-state dimensions disagree"))
    size(covariance) == (d, d) || throw(DimensionMismatch(
        "fixed-state covariance has the wrong dimension"))
    spec = _rget(subsampling, :residual_envelope)
    full_model = _build_sticky_target(build_omrf_full_model(spec, names, d),
        model_prior, slab_prior, can_stick_vec, names)
    PDMPSamplers.set_active_set!(full_model, free_vec)
    g_full = zeros(d)
    PDMPSamplers.compute_gradient_for_reflection!(full_model.grad, x, g_full)
    harmonic = occursin("Boomerang", flow_type)
    precision = harmonic ? inv(Symmetric(covariance)) : zeros(0, 0)
    correction = harmonic ? precision * (x - mean) : zeros(d)
    harmonic && (g_full .-= correction)

    function rate_and_jump(g)
        @inbounds for i in eachindex(g)
            free_vec[i] || (g[i] = 0.0)
        end
        signed = dot(v, g); rate = max(0.0, signed)
        iszero(rate) && return 0.0, 0.0
        Cg = covariance * g
        denominator = dot(g, Cg)
        denominator > 0 || return rate, 0.0
        coefficient = -2signed / denominator
        jump = 0.0
        @inbounds for i in eachindex(g)
            if free_vec[i]
                delta = coefficient * Cg[i]
                jump += delta^2 / max(covariance[i, i], eps(Float64))
            end
        end
        return rate, jump
    end
    full_gradient = copy(g_full)
    full_rate, full_jump = rate_and_jump(full_gradient)

    built = _build_stan_subsampling_model(first(prepared.contexts), names,
        subsampling, _as_float_vector(anchor), slab_prior, can_stick_vec;
        anchor_capacity=0)
    model = first(built)
    PDMPSamplers.set_active_set!(model, free_vec)
    cv = model.grad
    deterministic = zeros(d)
    PDMPSamplers.deterministic_gradient!(deterministic, cv, x)
    harmonic && (deterministic .-= correction)
    M = PDMPSamplers.n_observations(cv.envelope)
    rng = Xoshiro(Int(seed)); scratch = Dict{Int,Int}()
    residual = zeros(d); gradient = zeros(d); gradient_mean = zeros(d)
    batches = Int.(batches_raw)
    rates = zeros(length(batches)); rate_se = zeros(length(batches))
    jump_intensities = zeros(length(batches)); jump_se = zeros(length(batches))
    gradient_errors = zeros(length(batches))
    for (bi, m) in pairs(batches)
        1 <= m <= M || throw(ArgumentError("batch size lies outside 1:M"))
        subset = Vector{Int}(undef, m)
        n_draws = m == M ? 1 : Int(draws)
        rate_mean = rate_m2 = jump_mean = jump_m2 = 0.0
        fill!(gradient_mean, 0.0)
        for draw in 1:n_draws
            if m == M
                @inbounds for j in 1:M; subset[j] = j; end
            else
                PDMPSamplers._draw_uniform_without_replacement!(rng, subset,
                    M, scratch)
            end
            fill!(residual, 0.0)
            cv.residual_oracle(residual, x, subset, cv.anchor)
            @. gradient = deterministic + (M / m) * residual
            @. gradient_mean += (gradient - gradient_mean) / draw
            rate, jump = rate_and_jump(gradient)
            jump_intensity = rate * jump
            delta = rate - rate_mean; rate_mean += delta / draw
            rate_m2 += delta * (rate - rate_mean)
            delta_jump = jump_intensity - jump_mean
            jump_mean += delta_jump / draw
            jump_m2 += delta_jump * (jump_intensity - jump_mean)
        end
        rates[bi] = rate_mean; jump_intensities[bi] = jump_mean
        if n_draws > 1
            rate_se[bi] = sqrt(rate_m2 / (n_draws - 1) / n_draws)
            jump_se[bi] = sqrt(jump_m2 / (n_draws - 1) / n_draws)
        end
        gradient_errors[bi] = maximum(abs.(gradient_mean .- g_full))
    end
    return (batch=batches, draws=fill(Int(draws), length(batches)),
        factor_population=fill(M, length(batches)), intensity=rates,
        intensity_mcse=rate_se, canonical_jump_intensity=jump_intensities,
        canonical_jump_intensity_mcse=jump_se,
        averaged_gradient_max_abs_error=gradient_errors,
        full_intensity=fill(full_rate, length(batches)),
        full_canonical_jump_intensity=fill(full_rate * full_jump,
            length(batches)))
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
        free0=nothing,
        stored_velocity0=nothing,
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
        boomerang_variance_floor_fraction::Float64=0.0,
        sticky::Bool=false, can_stick=nothing, model_prior=nothing,
        parameter_prior=nothing, slab_prior=nothing,
        show_progress::Bool=true, n_chains::Int=1, threaded::Bool=false,
        seed::Union{Integer,Nothing}=nothing,
        adaptive_scheme::String="diagonal",
        refresh_rate::Float64=1.0,
        warmup_mode::String="subsampled",
        n_anchor_updates::Int=0,
        use_anchor_bank::Bool=false,
        anchor_bank_capacity::Int=8,
        anchor_selection_policy::String="adaptive",
        position_anchor_bank_spec=nothing,
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
        support_boundary_min_safe_time::Float64=1e-12,
        diagnostic_window_width::Float64=0.0,
        streaming_trace_directory::Union{Nothing,String}=nothing,
        streaming_trace_buffer_events::Int=4096)
    bridge_profile_start_ns = time_ns()
    algorithm_type in ("GridThinningStrategy", "ThinningStrategy") ||
        throw(ArgumentError("subsampling custom-Stan sampling requires GridThinningStrategy or ThinningStrategy"))
    warmup_mode in ("subsampled", "full") || throw(ArgumentError(
        "warmup_mode must be \"subsampled\" or \"full\""))
    warmup_mode == "full" && t_warmup <= 0 && throw(ArgumentError(
        "full-gradient warmup requires positive t_warmup"))
    x0_vec = _as_float_vector(x0)
    envelope_specification = _rget(subsampling, :residual_envelope)
    hcv_bank_learning = _haskey(envelope_specification,
            :hcv_bank_learning) ?
        Symbol(String(_rget(envelope_specification, :hcv_bank_learning))) :
        :geometric
    hcv_bank_learning in (:exact_variance, :selector_coupled_variance) &&
        t_warmup > 0 &&
        warmup_mode != "full" &&
        throw(ArgumentError(
            "exact-variance HCV bank learning requires full-gradient warmup"))
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
    anchor_selection_policy in ("adaptive", "position", "node_position") ||
        throw(ArgumentError(
            "anchor_selection_policy must be \"adaptive\", \"position\", or \"node_position\""))
    anchor_selection_policy in ("position", "node_position") && !use_anchor_bank &&
        throw(ArgumentError(
            "position-based anchor selection requires use_anchor_bank=true"))
    anchor_selection_policy == "position" &&
        !(flow_type in ("BouncyParticle", "PreconditionedBPS",
            "ZigZag", "PreconditionedZigZag", "Boomerang",
            "AdaptiveBoomerang")) &&
        throw(ArgumentError(
            "position anchor selection requires linear Zig-Zag/BPS or harmonic Boomerang dynamics"))
    anchor_selection_policy == "node_position" &&
        !(flow_type in ("BouncyParticle", "PreconditionedBPS")) &&
        throw(ArgumentError(
            "node-position anchor selection is available only for linear BPS flows"))
    !isnothing(position_anchor_bank_spec) &&
        !(anchor_selection_policy in ("position", "node_position")) && throw(ArgumentError(
            "a position anchor bank specification requires position selection"))
    !isnothing(position_anchor_bank_spec) && t_warmup > 0 &&
        throw(ArgumentError(
            "restored position anchor banks require zero warmup"))
    if anchor_selection_policy == "node_position"
        envelope_spec = _rget(subsampling, :residual_envelope)
        String(_rget(envelope_spec, :bound_type)) ==
            "cell_marked_node_box_first_order" || throw(ArgumentError(
            "node_position requires the cell-marked first-order node-box certificate"))
    end
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
    isfinite(boomerang_variance_floor_fraction) &&
        boomerang_variance_floor_fraction >= 0 || throw(ArgumentError(
            "boomerang_variance_floor_fraction must be finite and nonnegative"))
    for ctx in contexts
        built = _build_stan_subsampling_model(ctx, unc_names,
            subsampling, anchor, slab_prior, can_stick_vec; anchor_capacity,
            manage_during_warmup=warmup_mode != "full")
        if length(built) == 5
            model, _, _, adapter, manager = built
            adapter.update_dt = warmup_mode == "full" ? Inf : n_anchor_updates > 0 ?
                t_warmup / n_anchor_updates : Inf
            adapter.last_update = t0
            manager.main_refresh_distance = main_anchor_refresh_distance
            manager.node_source_check_distance = node_source_check_distance
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
        event_bound_observation_visits=
            ctx.counts.event_bound_observation_visits,
    ) for ctx in contexts]
    d = first(models).d
    flow = build_flow(flow_type, _to_precision(_as_flow_cov(flow_cov, d), d),
        _as_flow_mean(flow_mean, d); adaptive_scheme, refresh_rate)
    boom_options = PDMPSamplers.BoomerangAdaptationOptions(
        initial_variance_floor_fraction=boomerang_variance_floor_fraction)
    boom_floor_mask = _boomerang_logscale_floor_mask(unc_names, d)
    alg0 = build_algorithm(algorithm_type; c0, d, grid_n, grid_t_max,
        use_fd_hvp, curvature_backend,
        post_warmup_simplify, grid_bound, grid_curvature_bound,
        linear_area_threshold, linear_min_area_gain,
        lazy_low_tightness_threshold, lazy_max_low_tightness_rejections,
        lazy_max_rejections)
    parameter_prior_vec = isnothing(parameter_prior) ? nothing :
        _as_float_vector(parameter_prior)
    alg = _build_sticky_event_algorithm(alg0, sticky, model_prior,
        parameter_prior_vec, slab_prior, can_stick_vec, flow, unc_names)
    # The retained model is an active-slab target.  The full-data warmup
    # model must represent the same target, otherwise adaptation follows the
    # all-active analytic prior and then switches target at the warmup/main
    # boundary.  This is especially important for node-shared slabs: freezing
    # one interaction also changes its node-scale contribution.
    native_warmup_gradient_types = String[]
    native_warmup_hvp_types = String[]
    warmup_models = warmup_mode == "full" ?
        [begin
            native = build_omrf_full_model(
                _rget(subsampling, :residual_envelope), unc_names, d)
            push!(native_warmup_gradient_types, string(typeof(native.grad.f)))
            push!(native_warmup_hvp_types, string(typeof(native.hvp)))
            _build_sticky_target(native, model_prior, slab_prior,
                can_stick_vec, unc_names)
        end for _ in 1:n_chains] : nothing
    # This invariant must be checked while the native warmup providers still
    # exist.  After the dependent-slab wrapper is applied, the public model
    # type is intentionally different and is not evidence about the native
    # provider used for warmup.
    if warmup_mode == "full"
        all(occursin("OMRFFullGradientProvider", provider_type) for
            provider_type in native_warmup_gradient_types) ||
            throw(ArgumentError(
                "full-gradient warmup did not construct native " *
                "OMRFFullGradientProvider providers"))
        all(occursin("OMRFFullHVPProvider", provider_type) for
            provider_type in native_warmup_hvp_types) ||
            throw(ArgumentError(
                "full-gradient warmup did not construct native " *
                "OMRFFullHVPProvider providers"))
    end
    warmup_alg0 = warmup_mode == "full" ? build_algorithm(algorithm_type;
        c0, d, grid_n, grid_t_max, use_fd_hvp=false,
        curvature_backend="exact", post_warmup_simplify,
        grid_bound, grid_curvature_bound, linear_area_threshold,
        linear_min_area_gain, lazy_low_tightness_threshold,
        lazy_max_low_tightness_rejections, lazy_max_rejections) : nothing
    warmup_alg = if warmup_mode == "full"
        _build_sticky_event_algorithm(warmup_alg0, sticky, model_prior,
            parameter_prior_vec, slab_prior, can_stick_vec, flow, unc_names)
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
    # Adaptation must see the same gradient provider as the warmup event
    # process.  Keep this selection identical for both the anchor-manager and
    # no-anchor-manager paths; the retained subsampling provider is not a
    # substitute for the full-gradient warmup target.
    warmup_adaptation_gradient(i) = isnothing(warmup_models) ? models[i].grad :
        warmup_models[i].grad
    bridge_profile_before_sample_ns = time_ns()
    chains = if all(isnothing, anchor_managers)
        adaptation_grad = warmup_adaptation_gradient(1)
        adapter = PDMPSamplers.default_warmup_adapter(
            flow, adaptation_grad, t_warmup, t0;
            warmup_adaptation_interval,
            can_stick=sticky ? can_stick_vec : nothing,
            boomerang_options=boom_options,
            boomerang_variance_floor_mask=boom_floor_mask)
        pdmp_sample(initial, flow, models, alg, t0, T, t_warmup;
            progress=show_progress, threaded, seed,
            adapter,
            warmup_models, warmup_algorithm=warmup_alg,
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
            adaptation_grad = warmup_adaptation_gradient(i)
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
                warmup_models=isnothing(warmup_models) ? nothing :
                    [warmup_models[i]], warmup_algorithm=warmup_alg,
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
        "requested_full_gradient_backend" =>
            warmup_mode == "full" ? "native" : "bridgestan",
        "resolved_full_gradient_backend" => if warmup_mode == "full" &&
                !isempty(native_warmup_gradient_types)
            provider_type = first(native_warmup_gradient_types)
            occursin("OMRFFullGradientProvider", provider_type) ?
                "native" : "unresolved"
        else
            "bridgestan"
        end,
        "warmup_gradient_provider_type" => if warmup_mode == "full" &&
                !isempty(native_warmup_gradient_types)
            string(typeof(first(warmup_models).grad.f))
        else
            "none"
        end,
        "warmup_hvp_provider_type" => if warmup_mode == "full" &&
                !isempty(native_warmup_hvp_types)
            string(typeof(first(warmup_models).hvp))
        else
            "none"
        end,
        "native_warmup_gradient_provider_type" =>
            warmup_mode == "full" && !isempty(native_warmup_gradient_types) ?
            first(native_warmup_gradient_types) : "none",
        "native_warmup_hvp_provider_type" =>
            warmup_mode == "full" && !isempty(native_warmup_hvp_types) ?
            first(native_warmup_hvp_types) : "none",
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
        "analytic_factors_evaluated" => ctx.counts.analytic_factors,
        "node_conditionals_evaluated" => ctx.counts.analytic_node_conditionals,
        "compact_person_bounds_constructed" =>
            ctx.counts.compact_person_bounds,
        "compact_static_constructions" =>
            ctx.counts.compact_static_constructions,
        "compact_aggregate_point_calls" =>
            ctx.counts.compact_aggregate_point_calls,
        "compact_cell_bound_calls" => ctx.counts.compact_cell_bound_calls,
        "compact_term_selection_calls" =>
            ctx.counts.compact_term_selection_calls,
        "compact_person_selection_calls" =>
            ctx.counts.compact_person_selection_calls,
        "compact_observation_bound_calls" =>
            ctx.counts.compact_observation_bound_calls,
        "native_full_gradient_seconds" =>
            ctx.counts.native_full_gradient_seconds,
        "analytic_residual_seconds" => ctx.counts.analytic_residual_seconds,
        "envelope_point_calls" => ctx.counts.envelope_point_calls,
        "envelope_point_seconds" => ctx.counts.envelope_point_seconds,
        "envelope_cell_calls" => ctx.counts.envelope_cell_calls,
        "envelope_cell_seconds" => ctx.counts.envelope_cell_seconds,
        "harmonic_coefficient_preparations" =>
            ctx.counts.harmonic_coefficient_preparations,
        "harmonic_coefficient_cache_hits" =>
            ctx.counts.harmonic_coefficient_cache_hits,
        "harmonic_coefficient_seconds" =>
            ctx.counts.harmonic_coefficient_seconds,
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
        "cell_marked_component_draws" =>
            (model.grad.envelope isa PDMPSamplers.GroupedResidualEnvelope &&
             model.grad.envelope.component_scales! isa
                OMRFCellMarkedNodeBoxFirstOrderScales) ?
            model.grad.envelope.component_scales!.component_marks : 0,
        "cell_marked_factor_point_evaluations" =>
            (model.grad.envelope isa PDMPSamplers.GroupedResidualEnvelope &&
             model.grad.envelope.component_scales! isa
                OMRFCellMarkedNodeBoxFirstOrderScales) ?
            model.grad.envelope.component_scales!.factor_point_evaluations : 0,
        "cell_marked_factor_point_seconds" =>
            (model.grad.envelope isa PDMPSamplers.GroupedResidualEnvelope &&
             model.grad.envelope.component_scales! isa
                OMRFCellMarkedNodeBoxFirstOrderScales) ?
            model.grad.envelope.component_scales!.factor_point_seconds : 0.0,
        "categorical_residual_cap" =>
            (model.grad.envelope isa PDMPSamplers.GroupedResidualEnvelope &&
             model.grad.envelope.component_scales! isa
                OMRFCellMarkedNodeBoxFirstOrderScales) ?
            model.grad.envelope.component_scales!.base.
                categorical_residual_cap : false,
        "cell_marked_cell_cache_hits" =>
            (model.grad.envelope isa PDMPSamplers.GroupedResidualEnvelope &&
             model.grad.envelope.component_scales! isa
                OMRFCellMarkedNodeBoxFirstOrderScales) ?
            model.grad.envelope.component_scales!.cell_cache_hits : 0,
        "cell_marked_cell_cache_misses" =>
            (model.grad.envelope isa PDMPSamplers.GroupedResidualEnvelope &&
             model.grad.envelope.component_scales! isa
                OMRFCellMarkedNodeBoxFirstOrderScales) ?
            model.grad.envelope.component_scales!.cell_cache_misses : 0,
        "cell_marked_point_displacement_gt_four" =>
            (model.grad.envelope isa PDMPSamplers.GroupedResidualEnvelope &&
             model.grad.envelope.component_scales! isa
                OMRFCellMarkedNodeBoxFirstOrderScales) ?
            model.grad.envelope.component_scales!.
                point_displacement_gt_four : 0,
        "cell_marked_point_old_bound_sum" =>
            (model.grad.envelope isa PDMPSamplers.GroupedResidualEnvelope &&
             model.grad.envelope.component_scales! isa
                OMRFCellMarkedNodeBoxFirstOrderScales) ?
            model.grad.envelope.component_scales!.point_old_bound_sum : 0.0,
        "cell_marked_point_new_bound_sum" =>
            (model.grad.envelope isa PDMPSamplers.GroupedResidualEnvelope &&
             model.grad.envelope.component_scales! isa
                OMRFCellMarkedNodeBoxFirstOrderScales) ?
            model.grad.envelope.component_scales!.point_new_bound_sum : 0.0,
        "categorical_point_evaluations" =>
            ctx.counts.categorical_point_evaluations,
        "categorical_point_displacement_gt_four" =>
            ctx.counts.categorical_point_displacement_gt_four,
        "categorical_point_old_bound_sum" =>
            ctx.counts.categorical_point_old_bound_sum,
        "categorical_point_new_bound_sum" =>
            ctx.counts.categorical_point_new_bound_sum,
        "categorical_cell_cache_hits" =>
            ctx.counts.categorical_cell_cache_hits,
        "categorical_cell_cache_misses" =>
            ctx.counts.categorical_cell_cache_misses,
        "compact_static_seconds" => ctx.counts.compact_static_seconds,
        "compact_aggregate_point_seconds" =>
            ctx.counts.compact_aggregate_point_seconds,
        "compact_cell_bound_seconds" =>
            ctx.counts.compact_cell_bound_seconds,
        "compact_term_selection_seconds" =>
            ctx.counts.compact_term_selection_seconds,
        "compact_person_selection_seconds" =>
            ctx.counts.compact_person_selection_seconds,
        "compact_observation_bound_seconds" =>
            ctx.counts.compact_observation_bound_seconds,
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
        "direct_horizon_table_hits" =>
            ctx.counts.direct_horizon_table_hits,
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
        "signed_bps" => model.grad.envelope isa
            OMRFCompactCovariateLocalEnvelope && model.grad.envelope.signed_bps,
        "signed_harmonic" => model.grad.envelope isa
            OMRFCompactCovariateLocalEnvelope && model.grad.envelope.signed_harmonic,
        "signed_cell_preparations" =>
            ctx.counts.signed_bps.cell_preparations,
        "signed_endpoint_fills" => ctx.counts.signed_bps.endpoint_fills,
        "signed_person_node_masses_filled" =>
            ctx.counts.signed_bps.person_node_masses_filled,
        "signed_endpoint_reuses" => ctx.counts.signed_bps.endpoint_reuses,
        "signed_interaction_field_preparations" =>
            ctx.counts.signed_bps.interaction_field_preparations,
        "signed_endpoint_range_fills" =>
            ctx.counts.signed_bps.endpoint_range_fills,
        "signed_point_screens" => ctx.counts.signed_bps.point_screens,
        "signed_screen_copy_checks" =>
            ctx.counts.signed_bps.screen_copy_checks,
        "signed_component_draws" => ctx.counts.signed_bps.component_draws,
        "signed_cell_preparation_seconds" =>
            ctx.counts.signed_bps.cell_preparation_seconds,
        "signed_endpoint_fill_seconds" =>
            ctx.counts.signed_bps.endpoint_fill_seconds,
        "signed_interaction_field_seconds" =>
            ctx.counts.signed_bps.interaction_field_seconds,
        "signed_endpoint_range_seconds" =>
            ctx.counts.signed_bps.endpoint_range_seconds,
        "signed_point_screen_seconds" =>
            ctx.counts.signed_bps.point_screen_seconds,
        "signed_screen_copy_seconds" =>
            ctx.counts.signed_bps.screen_copy_seconds,
        "signed_component_draw_seconds" =>
            ctx.counts.signed_bps.component_draw_seconds,
        "signed_cell_preparation_allocated_bytes" =>
            ctx.counts.signed_bps.cell_preparation_allocated_bytes,
        "signed_endpoint_fill_allocated_bytes" =>
            ctx.counts.signed_bps.endpoint_fill_allocated_bytes,
        "signed_interaction_field_allocated_bytes" =>
            ctx.counts.signed_bps.interaction_field_allocated_bytes,
        "signed_endpoint_range_allocated_bytes" =>
            ctx.counts.signed_bps.endpoint_range_allocated_bytes,
        "signed_point_screen_allocated_bytes" =>
            ctx.counts.signed_bps.point_screen_allocated_bytes,
        "signed_screen_copy_allocated_bytes" =>
            ctx.counts.signed_bps.screen_copy_allocated_bytes,
        "signed_component_draw_allocated_bytes" =>
            ctx.counts.signed_bps.component_draw_allocated_bytes,
        "harmonic_cell_mass_constructions" =>
            ctx.counts.signed_bps.harmonic_cell_mass_constructions,
        "harmonic_factor_masses_filled" =>
            ctx.counts.signed_bps.harmonic_factor_masses_filled,
        "harmonic_all_factor_point_scans" =>
            ctx.counts.signed_bps.harmonic_all_factor_point_scans,
        "harmonic_selected_point_bounds" =>
            ctx.counts.signed_bps.harmonic_selected_point_bounds,
        "harmonic_cell_mass_seconds" =>
            ctx.counts.signed_bps.harmonic_cell_mass_seconds,
        "harmonic_selected_point_seconds" =>
            ctx.counts.signed_bps.harmonic_selected_point_seconds,
        "harmonic_cell_mass_allocated_bytes" =>
            ctx.counts.signed_bps.harmonic_cell_mass_allocated_bytes,
        "harmonic_selected_point_allocated_bytes" =>
            ctx.counts.signed_bps.harmonic_selected_point_allocated_bytes,
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
        "anchor_bank_requested_capacity" => manager === nothing ? 0 :
            manager.capacity,
        "anchor_bank_order" => manager === nothing ? Int[] :
            collect(eachindex(manager.entries)),
        "anchor_selection_policy" => manager === nothing ? "none" :
            String(manager.selection_policy),
        "anchor_selection_boundaries" => manager === nothing ? 0 :
            manager.selection_boundaries,
        "node_position_policy_version" => manager === nothing ? 0 :
            manager.node_position_policy_version,
        "warmup_bank_learning_policy" => manager === nothing ? "none" :
            String(manager.warmup_bank_learning_policy),
        "warmup_bank_candidate_count" => manager === nothing ? 0 :
            manager.warmup_bank_candidate_count,
        "warmup_bank_training_state_count" => manager === nothing ? 0 :
            manager.warmup_bank_training_state_count,
        "warmup_bank_objective_single" => manager === nothing ? NaN :
            manager.warmup_bank_objective_single,
        "warmup_bank_objective_selected" => manager === nothing ? NaN :
            manager.warmup_bank_objective_selected,
        "warmup_bank_objective_reduction" => manager === nothing ||
                !isfinite(manager.warmup_bank_objective_single) ||
                iszero(manager.warmup_bank_objective_selected) ? NaN :
            manager.warmup_bank_objective_single /
                manager.warmup_bank_objective_selected,
        "warmup_bank_learning_seconds" => manager === nothing ? 0.0 :
            manager.warmup_bank_learning_seconds,
        "node_position_sources" => manager === nothing ||
            manager.oracle.node_position_cv === nothing ? Int[] :
            copy(manager.oracle.node_position_cv.sources),
        "node_position_metric_variances" => manager === nothing ||
            manager.oracle.node_position_cv === nothing ? Float64[] :
            copy(manager.oracle.node_position_cv.metric_variances),
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
        "resolved_signed_bps" => model.grad.envelope isa
            OMRFCompactCovariateLocalEnvelope &&
            model.grad.envelope.signed_bps,
        "resolved_signed_harmonic" => model.grad.envelope isa
            OMRFCompactCovariateLocalEnvelope &&
            model.grad.envelope.signed_harmonic,
        "resolved_flow" => flow_type,
        "resolved_factorization" => model.grad.residual_oracle isa
            OMRFAnalyticSubsamplingOracle ? String(
                model.grad.residual_oracle.residual_context.factorization) :
            "unavailable",
        "resolved_geometry" => model.grad.envelope isa
            OMRFCompactCovariateLocalEnvelope ? "covariate_local_compact" :
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
        "hcv_damping" => model.grad.residual_oracle isa OMRFAnalyticSubsamplingOracle &&
            model.grad.residual_oracle.residual_context.use_hcv ?
            model.grad.residual_oracle.residual_context.hcv_damping : 0.0,
    ) for (model, ctx, snapshot, manager) in
        zip(models, contexts, construction_snapshots, anchor_managers)]
    result["subsampling_context_counters"] = subsampling_context_counters
    candidate_values = [model.grad.residual_oracle isa
        OMRFAnalyticSubsamplingOracle ?
        _omrf_candidate_record_matrix(model.grad.residual_oracle) :
        zeros(Float64, 0, 7) for model in models]
    result["candidate_records"] = isone(length(candidate_values)) ?
        first(candidate_values) : candidate_values
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

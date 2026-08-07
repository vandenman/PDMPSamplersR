# Analytic marked-family BridgeStan provider. Included inside
# PDMPSamplersRBridge after the shared BridgeStan and sampler helpers.

mutable struct MarkedBridgeCallCounts
    analytic_residual::Int
    prior_gradient::Int
    prior_hvp::Int
    full_gradient::Int
    anchor_preparations::Int
    anchor_activations::Int
end
MarkedBridgeCallCounts() = MarkedBridgeCallCounts(0, 0, 0, 0, 0, 0)

mutable struct MarkedBridgeStanContext
    sm_full::BridgeStan.StanModel
    sm_prior::BridgeStan.StanModel
    family::Symbol
    predictor_designs::Vector{Matrix{Float64}}
    predictor_indices::Vector{Vector{Int}}
    predictor_operator_norms::Vector{Float64}
    d::Int
    offsets::Matrix{Float64}
    response::Matrix{Float64}
    known_se::Vector{Float64}
    observation_multipliers::Vector{Float64}
    N::Int
    m::Int
    anchor::Vector{Float64}
    anchor_predictors::Matrix{Float64}
    anchor_kernel_values::Matrix{Float64}
    anchor_envelope::Any
    likelihood_hessian::Union{Nothing,Matrix{Float64}}
    use_hcv::Bool
    hcv_damping::Float64
    predictor_buffer::Vector{Float64}
    hcv_displacement::Vector{Float64}
    hcv_product::Vector{Float64}
    hcv_direction_product::Vector{Float64}
    full_anchor::Vector{Float64}
    prior_anchor::Vector{Float64}
    prior_x::Vector{Float64}
    calls::MarkedBridgeCallCounts
end

struct MarkedBridgeAnchorState
    anchor::Vector{Float64}
    anchor_predictors::Matrix{Float64}
    anchor_kernel_values::Matrix{Float64}
    full_anchor::Vector{Float64}
    prior_anchor::Vector{Float64}
    envelope::SeparableResidualEnvelope
    likelihood_hessian::Union{Nothing,Matrix{Float64}}
end

mutable struct MarkedBridgeAnchorEntry
    state::MarkedBridgeAnchorState
    age::Int
end

mutable struct MarkedBridgeAnchorManager
    ctx::MarkedBridgeStanContext
    entries::Vector{MarkedBridgeAnchorEntry}
    active_idx::Int
    capacity::Int
    proposal::Vector{Float64}
    main_activations::Int
end

function _marked_predictor_operator_norms(
        designs::Vector{Matrix{Float64}}, indices::Vector{Vector{Int}},
        N::Int, d::Int)
    K = length(designs)
    operator_norms = Vector{Float64}(undef, N)
    if K == 1
        # The operator norm of a d×1 predictor is its Euclidean row norm.
        design = only(designs)
        @inbounds for i in 1:N
            operator_norms[i] = sqrt(sum(abs2, @view design[i, :]))
        end
        return operator_norms
    end
    workspace = zeros(d, K)
    @inbounds for i in 1:N
        fill!(workspace, 0.0)
        for k in 1:K
            design = designs[k]
            active = indices[k]
            for column in eachindex(active)
                workspace[active[column], k] = design[i, column]
            end
        end
        operator_norms[i] = opnorm(workspace)
    end
    return operator_norms
end

function _fill_marked_anchor_cache!(ctx::MarkedBridgeStanContext, anchor,
        predictors::AbstractMatrix, kernels::AbstractMatrix)
    N, K = ctx.N, length(ctx.predictor_designs)
    @inbounds for i in 1:N, k in 1:K
        predictors[i, k] = _marked_eta(ctx, i, k, anchor)
    end
    if ctx.family in (:bernoulli, :binomial, :independent_bernoulli)
        @inbounds for i in 1:N, k in 1:K
            kernels[i, k] = _marked_logistic(predictors[i, k])
        end
    elseif ctx.family === :poisson
        @inbounds for i in 1:N
            kernels[i, 1] = exp(predictors[i, 1])
        end
    elseif ctx.family in (:categorical, :multinomial)
        @inbounds for i in 1:N
            maximum_eta = max(0.0, maximum(@view predictors[i, :]))
            denominator = exp(-maximum_eta)
            for k in 1:K
                denominator += exp(predictors[i, k] - maximum_eta)
            end
            for k in 1:K
                kernels[i, k] = exp(predictors[i, k] - maximum_eta) / denominator
            end
        end
    end
    return nothing
end

function _new_marked_bss_context(lib_path_std::String,
        data_full_file::String, data_prior_file::String,
        family::Symbol, predictor_designs::Vector{Matrix{Float64}},
        predictor_indices::Vector{Vector{Int}}, d::Int,
        offsets::Matrix{Float64}, response::Matrix{Float64}, known_se::Vector{Float64},
        observation_multipliers::Vector{Float64},
        N::Int, m::Int, anchor::Vector{Float64};
        use_hcv::Bool=false, hcv_damping::Real=d)
    K = length(predictor_designs)
    K == length(predictor_indices) || throw(DimensionMismatch(
        "predictor designs and active-coordinate indices must have equal length"))
    all(size(design, 1) == N for design in predictor_designs) ||
        throw(DimensionMismatch("each compact predictor design must have N rows"))
    all(size(predictor_designs[k], 2) == length(predictor_indices[k]) for k in 1:K) ||
        throw(DimensionMismatch("compact predictor columns must match active indices"))
    all(all(i -> 1 <= i <= d, active) && issorted(active) && allunique(active)
        for active in predictor_indices) || throw(ArgumentError(
        "predictor indices must be sorted unique coordinates in 1:d"))
    size(offsets) == (N, K) || throw(DimensionMismatch(
        "offset matrix must have size (N, K)"))
    size(response, 1) == N || throw(DimensionMismatch("response must have N rows"))
    length(observation_multipliers) == N || throw(DimensionMismatch(
        "observation multiplier vector must have length N"))
    any(x -> !isfinite(x) || x < 0, observation_multipliers) &&
        throw(ArgumentError("observation multipliers must be finite and nonnegative"))
    if use_hcv && !(family in (:bernoulli, :binomial) && K == 1)
        throw(ArgumentError(
            "analytic marked HCV requires a one-predictor Bernoulli or binomial likelihood"))
    end
    isfinite(hcv_damping) && hcv_damping > 0 || throw(ArgumentError(
        "analytic marked HCV damping must be finite and positive"))
    operator_norms = _marked_predictor_operator_norms(
        predictor_designs, predictor_indices, N, d)
    ctx = MarkedBridgeStanContext(
        BridgeStan.StanModel(lib_path_std, data_full_file; warn=false),
        BridgeStan.StanModel(lib_path_std, data_prior_file; warn=false),
        family, predictor_designs, predictor_indices, operator_norms, d,
        offsets, response, known_se,
        observation_multipliers, N, m, copy(anchor),
        zeros(N, K), zeros(N, K), nothing, nothing, use_hcv,
        Float64(hcv_damping), zeros(K), zeros(d), zeros(d), zeros(d),
        zeros(d), zeros(d), zeros(d), MarkedBridgeCallCounts())
    initial = _prepare_marked_anchor(ctx, anchor; count_preparation=false)
    _install_marked_anchor_state!(ctx, initial)
    return ctx
end

@inline function _marked_logistic(x::Float64)
    if x >= 0
        return inv(1 + exp(-x))
    end
    ex = exp(x)
    return ex / (1 + ex)
end

function _marked_eta(ctx, i, k, x)
    value = ctx.offsets[i, k]
    design = ctx.predictor_designs[k]
    indices = ctx.predictor_indices[k]
    @inbounds for column in eachindex(indices)
        value += design[i, column] * x[indices[column]]
    end
    return value
end


function _add_marked_predictor!(out, ctx, i, k, value)
    design = ctx.predictor_designs[k]
    indices = ctx.predictor_indices[k]
    @inbounds for column in eachindex(indices)
        out[indices[column]] += design[i, column] * value
    end
    return out
end

function _marked_growth_envelope(weights_by_observation, growth, anchor)
    N = length(growth)
    positive = findall(i -> growth[i] > 0 && weights_by_observation[i] > 0, 1:N)
    if isempty(positive)
        return TrajectoryResidualEnvelope(reshape(weights_by_observation, 1, N), anchor)
    end
    n_bins = min(16, max(1, ceil(Int, sqrt(length(positive)))))
    order = sort(positive; by=i -> growth[i])
    cumulative_weight = cumsum(weights_by_observation[order])
    total_weight = cumulative_weight[end]
    edges = Float64[]
    for target in range(total_weight / n_bins, total_weight; length=n_bins)
        index = searchsortedfirst(cumulative_weight, target)
        push!(edges, growth[order[min(index, length(order))]])
    end
    edges[end] = growth[order[end]]
    unique!(edges)
    weights = zeros(length(edges), N)
    for i in 1:N
        bin = searchsortedfirst(edges, growth[i])
        bin = min(max(bin, 1), length(edges))
        weights[bin, i] = weights_by_observation[i]
    end
    return TrajectoryResidualEnvelope(weights, anchor; growth_rates=edges)
end

function _marked_family_envelope(ctx::MarkedBridgeStanContext,
        anchor=ctx.anchor, anchor_predictors=ctx.anchor_predictors)
    designs = ctx.predictor_designs
    N, K = ctx.N, length(designs)
    weights = zeros(N)
    growth = zeros(N)
    hcv_remainder_weights = ctx.use_hcv ? zeros(N) : nothing
    for i in 1:N
        Bnorm = ctx.predictor_operator_norms[i]
        multiplier = ctx.observation_multipliers[i]
        if ctx.family in (:bernoulli, :binomial, :independent_bernoulli)
            weights[i] = multiplier * Bnorm^2 / 4
            ctx.use_hcv && (hcv_remainder_weights[i] =
                multiplier * Bnorm^3 / (12sqrt(3.0)))
        elseif ctx.family in (:categorical, :multinomial)
            weights[i] = multiplier * Bnorm^2 / 2
        elseif ctx.family === :poisson
            ηa = anchor_predictors[i, 1]
            weights[i] = multiplier * exp(ηa) * Bnorm^2
            growth[i] = Bnorm
        elseif ctx.family === :gaussian && K == 1
            isempty(ctx.known_se) && throw(ArgumentError(
                "fixed-scale Gaussian marked sampling requires known standard errors"))
            weights[i] = multiplier * Bnorm^2 / ctx.known_se[i]^2
        elseif ctx.family === :gaussian && K == 2
            μa = anchor_predictors[i, 1]
            sa = anchor_predictors[i, 2]
            q0 = abs(ctx.response[i, 1] - μa) + 1
            coefficient = 0.5 + exp(-2sa) *
                (1 + 2sqrt(2.0) * q0 + 2.5q0^2)
            weights[i] = multiplier * Bnorm^2 * coefficient
            bμ = norm(@view designs[1][i, :])
            bs = norm(@view designs[2][i, :])
            growth[i] = 2 * (bμ + bs)
        else
            throw(ArgumentError("no marked family envelope for $(ctx.family) with $K predictors"))
        end
    end
    ctx.use_hcv && return PDMPSamplers.DampedHCVResidualEnvelope(
        weights, hcv_remainder_weights, anchor; damping=ctx.hcv_damping)
    return _marked_growth_envelope(weights, growth, anchor)
end

function _marked_likelihood_hessian(ctx::MarkedBridgeStanContext,
        anchor_kernels::AbstractMatrix)
    ctx.use_hcv || return nothing
    hessian = zeros(ctx.d, ctx.d)
    design = only(ctx.predictor_designs)
    indices = only(ctx.predictor_indices)
    @inbounds for i in 1:ctx.N
        p = anchor_kernels[i, 1]
        coefficient = ctx.observation_multipliers[i] * p * (1 - p)
        for column in eachindex(indices), row in eachindex(indices)
            hessian[indices[row], indices[column]] += coefficient *
                design[i, row] * design[i, column]
        end
    end
    return hessian
end

function _prepare_marked_anchor(ctx::MarkedBridgeStanContext,
        anchor::AbstractVector; count_preparation::Bool=true)
    length(anchor) == ctx.d || throw(DimensionMismatch(
        "marked anchor has the wrong dimension"))
    requested = similar(ctx.anchor)
    copyto!(requested, anchor)
    K = length(ctx.predictor_designs)
    predictors = zeros(ctx.N, K)
    kernels = zeros(ctx.N, K)
    _fill_marked_anchor_cache!(ctx, requested, predictors, kernels)
    likelihood_hessian = _marked_likelihood_hessian(ctx, kernels)
    full_anchor = zeros(ctx.d)
    prior_anchor = zeros(ctx.d)
    BridgeStan.log_density_gradient!(ctx.sm_full, requested, full_anchor)
    ctx.calls.full_gradient += 1
    full_anchor .*= -1
    BridgeStan.log_density_gradient!(ctx.sm_prior, requested, prior_anchor)
    ctx.calls.prior_gradient += 1
    envelope = _marked_family_envelope(ctx, requested, predictors)
    count_preparation && (ctx.calls.anchor_preparations += 1)
    return MarkedBridgeAnchorState(requested, predictors, kernels,
        full_anchor, prior_anchor, envelope, likelihood_hessian)
end

function _anchor_state_from_context(ctx::MarkedBridgeStanContext)
    return MarkedBridgeAnchorState(ctx.anchor, ctx.anchor_predictors,
        ctx.anchor_kernel_values, ctx.full_anchor, ctx.prior_anchor,
        ctx.anchor_envelope, ctx.likelihood_hessian)
end

function _install_marked_anchor_state!(ctx::MarkedBridgeStanContext,
        state::MarkedBridgeAnchorState)
    ctx.anchor = state.anchor
    ctx.anchor_predictors = state.anchor_predictors
    ctx.anchor_kernel_values = state.anchor_kernel_values
    ctx.anchor_envelope = state.envelope
    ctx.likelihood_hessian = state.likelihood_hessian
    ctx.full_anchor = state.full_anchor
    ctx.prior_anchor = state.prior_anchor
    return ctx
end

@inline function _marked_anchor_distance2(x::AbstractVector,
        state::MarkedBridgeAnchorState)
    distance = zero(eltype(x))
    @inbounds for j in eachindex(x, state.anchor)
        delta = x[j] - state.anchor[j]
        distance += delta * delta
    end
    return distance
end

function _active_marked_anchor(manager::MarkedBridgeAnchorManager)
    1 <= manager.active_idx <= length(manager.entries) || throw(ArgumentError(
        "marked anchor bank has no active entry"))
    return manager.entries[manager.active_idx]
end

function _insert_marked_anchor!(manager::MarkedBridgeAnchorManager,
        state::MarkedBridgeAnchorState)
    idx = if length(manager.entries) < manager.capacity
        push!(manager.entries, MarkedBridgeAnchorEntry(state, 0))
        length(manager.entries)
    else
        ages = map(entry -> entry.age, manager.entries)
        replace_idx = argmax(ages)
        manager.entries[replace_idx] = MarkedBridgeAnchorEntry(state, 0)
        replace_idx
    end
    iszero(manager.active_idx) && (manager.active_idx = idx)
    return idx
end

function _activate_marked_anchor!(manager::MarkedBridgeAnchorManager,
        requested::AbstractVector)
    state = _active_marked_anchor(manager).state
    state.anchor == requested || throw(ArgumentError(
        "active marked bank entry does not match the requested anchor"))
    ctx = manager.ctx
    # Every operation below is an assignment of already validated, prepared
    # state. No full gradient, observation scan, or envelope construction is
    # performed when an existing bank entry is selected.
    _install_marked_anchor_state!(ctx, state)
    ctx.calls.anchor_activations += 1
    return state.envelope
end

function _new_marked_anchor_manager(ctx::MarkedBridgeStanContext,
        capacity::Integer)
    capacity >= 1 || throw(ArgumentError("marked anchor-bank capacity must be positive"))
    manager = MarkedBridgeAnchorManager(ctx, MarkedBridgeAnchorEntry[], 0,
        Int(capacity), zeros(ctx.d), 0)
    _insert_marked_anchor!(manager, _anchor_state_from_context(ctx))
    return manager
end

function _select_marked_anchor!(manager::MarkedBridgeAnchorManager,
        cv::MarkedControlVariate, x::AbstractVector; phase::Symbol=:unknown)
    previous = manager.active_idx
    best_idx = 1
    best_distance = Inf
    @inbounds for idx in eachindex(manager.entries)
        entry = manager.entries[idx]
        entry.age += 1
        distance = _marked_anchor_distance2(x, entry.state)
        if distance < best_distance
            best_distance = distance
            best_idx = idx
        end
    end
    manager.entries[best_idx].age = 0
    manager.active_idx = best_idx
    if best_idx != previous
        PDMPSamplers.refresh_anchor!(cv,
            manager.entries[best_idx].state.anchor)
        phase === :main && (manager.main_activations += 1)
    end
    return nothing
end

function _add_marked_anchor!(manager::MarkedBridgeAnchorManager,
        cv::MarkedControlVariate, trace)
    Statistics.mean!(manager.proposal, trace)
    prepared = _prepare_marked_anchor(manager.ctx, manager.proposal)
    return _store_prepared_marked_anchor!(manager, cv, prepared)
end

function _store_prepared_marked_anchor!(manager::MarkedBridgeAnchorManager,
        cv::MarkedControlVariate, prepared::MarkedBridgeAnchorState)
    previous = manager.active_idx
    idx = _insert_marked_anchor!(manager, prepared)
    # If LRU replacement overwrote the active slot, immediately switch all
    # marked state to the replacement. Otherwise selection occurs at the next
    # event boundary.
    if idx == previous
        manager.active_idx = idx
        PDMPSamplers.refresh_anchor!(cv, prepared.anchor)
    end
    return idx
end

function _build_marked_bss_model(ctx::MarkedBridgeStanContext;
        use_fd_hvp::Bool=false, anchor_capacity::Integer=0)
    designs = ctx.predictor_designs
    N, d, K = ctx.N, ctx.d, length(designs)
    d == length(ctx.anchor) || throw(DimensionMismatch("likelihood design column count must equal the unconstrained dimension"))
    envelope = ctx.anchor_envelope
    envelope isa SeparableResidualEnvelope || throw(ArgumentError(
        "marked context has no prepared initial envelope"))
    manager = anchor_capacity > 0 ?
        _new_marked_anchor_manager(ctx, anchor_capacity) : nothing

    function deterministic_gradient!(out, x)
        BridgeStan.log_density_gradient!(ctx.sm_prior, x, ctx.prior_x)
        ctx.calls.prior_gradient += 1
        @. out = ctx.full_anchor - ctx.prior_x + ctx.prior_anchor
        if ctx.use_hcv
            likelihood_hessian = something(ctx.likelihood_hessian)
            displacement = ctx.hcv_displacement
            @. displacement = x - ctx.anchor
            displacement2 = dot(displacement, displacement)
            α = ctx.hcv_damping / (ctx.hcv_damping + displacement2)
            mul!(ctx.hcv_product, likelihood_hessian, displacement)
            @. out += α * ctx.hcv_product
        end
        return out
    end

    function deterministic_hvp!(out, x, v)
        BridgeStan.log_density_hessian_vector_product!(ctx.sm_prior, x, v, out)
        ctx.calls.prior_hvp += 1
        out .*= -1
        if ctx.use_hcv
            likelihood_hessian = something(ctx.likelihood_hessian)
            displacement = ctx.hcv_displacement
            @. displacement = x - ctx.anchor
            denominator = ctx.hcv_damping + dot(displacement, displacement)
            α = ctx.hcv_damping / denominator
            derivative_α = -2 * ctx.hcv_damping * dot(displacement, v) /
                denominator^2
            mul!(ctx.hcv_product, likelihood_hessian, displacement)
            mul!(ctx.hcv_direction_product, likelihood_hessian, v)
            @. out += α * ctx.hcv_direction_product + derivative_α * ctx.hcv_product
        end
        return out
    end

    function residual_oracle!(out, x, subset, anchor)
        fill!(out, 0.0)
        α = 0.0
        if ctx.use_hcv
            displacement = ctx.hcv_displacement
            @. displacement = x - ctx.anchor
            displacement2 = dot(displacement, displacement)
            α = ctx.hcv_damping / (ctx.hcv_damping + displacement2)
        end
        @inbounds for i in subset
            multiplier = ctx.observation_multipliers[i]
            if ctx.family in (:bernoulli, :binomial, :independent_bernoulli)
                predictors = ctx.family === :independent_bernoulli ? (1:K) : (1:1)
                for k in predictors
                    Δη = multiplier * (_marked_logistic(_marked_eta(ctx, i, k, x)) -
                        ctx.anchor_kernel_values[i, k])
                    if ctx.use_hcv
                        predictor_delta = _marked_eta(ctx, i, k, x) -
                            ctx.anchor_predictors[i, k]
                        p_anchor = ctx.anchor_kernel_values[i, k]
                        Δη -= α * multiplier * p_anchor * (1 - p_anchor) *
                            predictor_delta
                    end
                    _add_marked_predictor!(out, ctx, i, k, Δη)
                end
            elseif ctx.family === :poisson
                Δη = multiplier * (exp(_marked_eta(ctx, i, 1, x)) -
                    ctx.anchor_kernel_values[i, 1])
                _add_marked_predictor!(out, ctx, i, 1, Δη)
            elseif ctx.family === :gaussian && K == 1
                Δη = multiplier * (_marked_eta(ctx, i, 1, x) -
                    ctx.anchor_predictors[i, 1]) / ctx.known_se[i]^2
                _add_marked_predictor!(out, ctx, i, 1, Δη)
            elseif ctx.family === :gaussian && K == 2
                μx = _marked_eta(ctx, i, 1, x); sx = _marked_eta(ctx, i, 2, x)
                μa = ctx.anchor_predictors[i, 1]; sa = ctx.anchor_predictors[i, 2]
                qx = ctx.response[i, 1] - μx; qa = ctx.response[i, 1] - μa
                se2 = isempty(ctx.known_se) ? 0.0 : ctx.known_se[i]^2
                ex = exp(2sx); ea = exp(2sa)
                vx = ex + se2; va = ea + se2
                Δμ = multiplier * ((μx - ctx.response[i, 1]) / vx -
                    (μa - ctx.response[i, 1]) / va)
                Δs = multiplier * ((ex / vx) * (1 - qx^2 / vx) -
                    (ea / va) * (1 - qa^2 / va))
                _add_marked_predictor!(out, ctx, i, 1, Δμ)
                _add_marked_predictor!(out, ctx, i, 2, Δs)
            elseif ctx.family in (:categorical, :multinomial)
                predictor_buffer = ctx.predictor_buffer
                max_x = 0.0
                for k in 1:K
                    predictor_buffer[k] = _marked_eta(ctx, i, k, x)
                    max_x = max(max_x, predictor_buffer[k])
                end
                # Include the reference category's zero logit and normalize
                # after subtracting the maximum. This is the ordinary stable
                # softmax and remains finite for extreme affine predictors.
                denom_x = exp(-max_x)
                for k in 1:K
                    denom_x += exp(predictor_buffer[k] - max_x)
                end
                for k in 1:K
                    Δp = multiplier * (
                        exp(predictor_buffer[k] - max_x) / denom_x -
                        ctx.anchor_kernel_values[i, k])
                    _add_marked_predictor!(out, ctx, i, k, Δp)
                end
            end
        end
        ctx.calls.analytic_residual += 1
        return out
    end

    refresh_callback = manager === nothing ? nothing :
        (anchor -> _activate_marked_anchor!(manager, anchor))
    cv = MarkedControlVariate(deterministic_gradient!, residual_oracle!, envelope,
        ctx.anchor, ctx.m;
        deterministic_hvp! = use_fd_hvp ? nothing : deterministic_hvp!,
        refresh_anchor! = refresh_callback)
    model = PDMPModel(d, cv)
    if manager === nothing
        return model
    end
    adapter = PDMPSamplers.MarkedAnchorBankAdapter(
        (grad, x, phase) -> _select_marked_anchor!(manager, grad, x; phase),
        (grad, trace) -> _add_marked_anchor!(manager, grad, trace),
        Inf, 0.0)
    return model, adapter, manager
end

function _marked_call_counts(ctx::MarkedBridgeStanContext,
        manager::Union{Nothing,MarkedBridgeAnchorManager}=nothing)
    return Dict(
        "analytic_residual_calls" => ctx.calls.analytic_residual,
        "prior_gradient_calls" => ctx.calls.prior_gradient,
        "prior_hvp_calls" => ctx.calls.prior_hvp,
        "full_gradient_calls" => ctx.calls.full_gradient,
        "anchor_preparations" => ctx.calls.anchor_preparations,
        "anchor_activations" => ctx.calls.anchor_activations,
        "anchor_bank_entries" => manager === nothing ? 0 : length(manager.entries),
        "anchor_bank_bytes" => manager === nothing ? 0 : Base.summarysize(manager.entries),
        "anchor_main_activations" => manager === nothing ? 0 : manager.main_activations,
        "analytic_hcv" => ctx.use_hcv,
        "hcv_damping" => ctx.hcv_damping,
    )
end

function r_marked_family_closure_diagnostics(stan_file::String,
        data_full_file::String, data_prior_file::String, family_name,
        predictor_designs, predictor_indices, predictor_dimension::Integer,
        offsets, response, known_se, observation_multipliers,
        anchor, theta, use_hcv::Bool=false)
    designs = Matrix{Float64}[Matrix{Float64}(design)
        for design in _marked_collection_values(predictor_designs)]
    indices = Vector{Int}[_as_marked_predictor_indices(active)
        for active in _marked_collection_values(predictor_indices)]
    N = size(first(designs), 1)
    d = Int(predictor_dimension)
    anchor_vec = _marked_float_vector(anchor)
    theta_vec = _marked_float_vector(theta)
    ctx = _new_marked_bss_context(_compile_model(stan_file),
        data_full_file, data_prior_file, Symbol(family_name), designs, indices, d,
        Matrix{Float64}(offsets), Matrix{Float64}(response),
        _marked_float_vector(known_se), _marked_float_vector(observation_multipliers),
        N, 1, anchor_vec; use_hcv)
    model = _build_marked_bss_model(ctx)

    deterministic = zeros(d)
    residual = zeros(d)
    estimate = zeros(d)
    model.grad.deterministic_gradient!(deterministic, theta_vec)
    copyto!(estimate, deterministic)

    displacement = theta_vec - anchor_vec
    state = PDMPState(0.0,
        SkeletonPoint(copy(anchor_vec), copy(displacement)))
    flow = BouncyParticle(d, 0.0)
    scales = PDMPSamplers.component_scales!(model.grad.envelope.scales,
        model.grad.envelope, state, flow, 1.0)
    max_bound_ratio = 0.0
    for i in 1:N
        model.grad.residual_oracle(residual, theta_vec, [i], anchor_vec)
        estimate .+= residual
        bound = dot(scales, @view model.grad.envelope.weights[:, i])
        actual = norm(displacement) * norm(residual)
        ratio = iszero(bound) ? (iszero(actual) ? 0.0 : Inf) : actual / bound
        max_bound_ratio = max(max_bound_ratio, ratio)
    end

    full = zeros(d)
    BridgeStan.log_density_gradient!(ctx.sm_full, theta_vec, full)
    full .*= -1
    return Dict(
        "closure_error" => maximum(abs, estimate - full),
        "max_bound_ratio" => max_bound_ratio,
        "analytic_residual_calls" => ctx.calls.analytic_residual,
    )
end

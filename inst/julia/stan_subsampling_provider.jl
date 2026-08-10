# Persistent selected-gradient backend for public custom-Stan subsampling sampling.

mutable struct StanSubsamplingCallCounts
    full_gradient::Int
    prior_gradient::Int
    selected_gradient::Int
    persons_evaluated::Int
    anchor_cache_gradient::Int
    anchor_cache_persons::Int
    model_constructions::Int
    data_constructions::Int
end
StanSubsamplingCallCounts() = StanSubsamplingCallCounts(0, 0, 0, 0, 0, 0, 0, 0)

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
        unc_names, subsampling, anchor, slab_prior, can_stick)
    N = Int(_rget(subsampling, :n_observations))
    m = Int(_rget(subsampling, :subsample_size))
    envelope_spec = _rget(subsampling, :residual_envelope)
    weights = _as_float_matrix(_rget(envelope_spec, :weights))
    size(weights, 2) == N || throw(DimensionMismatch(
        "residual envelope must have N observation columns"))
    growth = _as_float_vector(_rget(envelope_spec, :growth_rates))
    d = ctx.full.d
    length(anchor) == d || throw(DimensionMismatch("subsampling anchor has the wrong dimension"))
    envelope = TrajectoryResidualEnvelope(weights, anchor; growth_rates=growth)
    oracle = StanSubsamplingOracle(ctx, anchor, N)
    deterministic = if isnothing(slab_prior)
        oracle
    else
        provider = build_slab_provider(slab_prior, unc_names, can_stick, d)
        StanSubsamplingSlabDeterministic(oracle, provider, d)
    end
    cv = SubsampledControlVariate(deterministic, oracle, envelope, anchor, m)
    return PDMPModel(d, cv), ctx, unc_names
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
        flow_type::String, flow_mean, flow_cov)
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
    oracle = StanSubsamplingOracle(ctx, anchor, N)
    deterministic = zeros(d)
    residual = zeros(d)
    full = zeros(d)
    oracle(deterministic, x)
    oracle(residual, x, subset_vec, anchor)
    prior_gradient = copy(oracle.prior_x)
    selected_gradient = copy(oracle.selected_x)
    selected_likelihood = oracle.selected_x - oracle.prior_x
    selected_anchor_likelihood = oracle.selected_anchor - oracle.prior_anchor
    _clear_gradient!(ctx, :full, full, x)
    subsampled_gradient = deterministic + (N / m) * residual
    anchor_deterministic = zeros(d)
    anchor_residual = zeros(d)
    oracle(anchor_deterministic, anchor)
    oracle(anchor_residual, anchor, subset_vec, anchor)
    anchor_subsampling = anchor_deterministic + (N / m) * anchor_residual
    anchor_closure_error = maximum(abs, anchor_subsampling - oracle.full_anchor)

    envelope_spec = _rget(subsampling, :residual_envelope)
    weights = _as_float_matrix(_rget(envelope_spec, :weights))
    growth = _as_float_vector(_rget(envelope_spec, :growth_rates))
    envelope = TrajectoryResidualEnvelope(weights, anchor; growth_rates=growth)
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
    envelope_rate = (N / m) * sum(i -> dot(scales, @view(weights[:, i])), subset_vec)
    scaled_residual = (N / m) * residual
    residual_rate = PDMPSamplers.λ(state, scaled_residual, flow) +
        PDMPSamplers.λ(state, -scaled_residual, flow)

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
        sticky::Bool=false, can_stick=nothing, model_prior=nothing,
        parameter_prior=nothing, slab_prior=nothing,
        show_progress::Bool=true, n_chains::Int=1, threaded::Bool=false,
        seed::Union{Integer,Nothing}=nothing,
        adaptive_scheme::String="diagonal",
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
    x0_vec = _as_float_vector(x0)
    subsampling_anchor = _rget(subsampling, :anchor)
    anchor = isnothing(subsampling_anchor) ? copy(x0_vec) :
        _as_float_vector(subsampling_anchor)
    can_stick_vec = _as_bool_vector(can_stick)
    models = PDMPModel[]
    contexts = prepared.contexts
    length(contexts) == n_chains || throw(DimensionMismatch(
        "prepared subsampling context count must equal n_chains"))
    unc_names = prepared.unc_names
    for ctx in contexts
        model, _, _ = _build_stan_subsampling_model(ctx, unc_names,
            subsampling, anchor, slab_prior, can_stick_vec)
        push!(models, model)
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
    chains = pdmp_sample(initial, flow, models, alg, t0, T, t_warmup;
        progress=show_progress, threaded, seed,
        warmup_stop=build_warmup_stop(t_warmup),
        support_boundary_options=sbopts,
        statistic_counter=PDMPSamplers.DevelStatisticCounter)
    result = _pack_result(chains)
    result["subsampling_context_counters"] = [Dict(
        "full_gradient_calls" => ctx.counts.full_gradient,
        "prior_gradient_calls" => ctx.counts.prior_gradient,
        "selected_gradient_calls" => ctx.counts.selected_gradient,
        "persons_evaluated" => ctx.counts.persons_evaluated,
        "anchor_cache_gradient_calls" => ctx.counts.anchor_cache_gradient,
        "anchor_cache_persons_evaluated" => ctx.counts.anchor_cache_persons,
        "anchor_likelihood_cache_bytes" => Base.summarysize(
            model.grad.residual_oracle.anchor_likelihood_cache),
        "model_constructions" => ctx.counts.model_constructions,
        "data_constructions" => ctx.counts.data_constructions,
        "initialization_model_constructions" => snapshot.model,
        "initialization_data_constructions" => snapshot.data,
        "initialization_full_gradient_calls" => snapshot.full_gradient,
        "sampling_model_constructions" => ctx.counts.model_constructions - snapshot.model,
        "sampling_data_constructions" => ctx.counts.data_constructions - snapshot.data,
        "sampling_full_gradient_calls" => ctx.counts.full_gradient - snapshot.full_gradient,
    ) for (model, ctx, snapshot) in zip(models, contexts, construction_snapshots)]
    result["subsampling"] = true
    return result
end

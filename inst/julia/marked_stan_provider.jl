# Persistent selected-gradient backend for public custom-Stan marked sampling.

const _STAN_MARKED_BRIDGE_LOCK = ReentrantLock()

mutable struct StanMarkedCallCounts
    full_gradient::Int
    prior_gradient::Int
    selected_gradient::Int
    persons_evaluated::Int
    model_constructions::Int
end
StanMarkedCallCounts() = StanMarkedCallCounts(0, 0, 0, 0, 3)

mutable struct StanMarkedContext
    full::PDMPModel
    prior::PDMPModel
    selected::PDMPModel
    set_fn::Ptr{Nothing}
    clear_fn::Ptr{Nothing}
    subset0::Vector{Int32}
    m::Int
    counts::StanMarkedCallCounts
end

function _resolve_subset_hooks(sm::BridgeStan.StanModel)
    set_fn = Libdl.dlsym(sm.lib, :pdmp_set_subsample_indices; throw_error=false)
    clear_fn = Libdl.dlsym(sm.lib, :pdmp_clear_subsample_indices; throw_error=false)
    (set_fn == C_NULL || clear_fn == C_NULL) && throw(ArgumentError(
        "the compiled Stan model does not export the PDMPSamplersR subset hook; compile it with compile_pdmp_stan_model()"))
    return set_fn, clear_fn
end

@inline function _clear_stan_subset!(ctx::StanMarkedContext)
    @ccall $(ctx.clear_fn)()::Cvoid
    return nothing
end

@inline function _install_stan_subset!(ctx::StanMarkedContext, subset)
    length(subset) == ctx.m || throw(DimensionMismatch("selected subset has the wrong size"))
    @inbounds for j in eachindex(subset)
        ctx.subset0[j] = Int32(subset[j] - 1)
    end
    subset0 = ctx.subset0
    set_fn = ctx.set_fn
    m = Cint(ctx.m)
    GC.@preserve subset0 begin
        @ccall $set_fn(pointer(subset0)::Ptr{Int32}, m::Cint)::Cvoid
    end
    return nothing
end

function _clear_gradient!(ctx::StanMarkedContext, which::Symbol, out, x)
    lock(_STAN_MARKED_BRIDGE_LOCK) do
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

function _selected_gradient!(ctx::StanMarkedContext, out, x, subset)
    lock(_STAN_MARKED_BRIDGE_LOCK) do
        _install_stan_subset!(ctx, subset)
        try
            compute_gradient!(ctx.selected.grad, x, out)
            ctx.counts.selected_gradient += 1
            ctx.counts.persons_evaluated += ctx.m
        finally
            _clear_stan_subset!(ctx)
        end
    end
    return out
end

mutable struct StanMarkedOracle
    ctx::StanMarkedContext
    full_anchor::Vector{Float64}
    prior_anchor::Vector{Float64}
    prior_x::Vector{Float64}
    selected_x::Vector{Float64}
    selected_anchor::Vector{Float64}
end

function StanMarkedOracle(ctx::StanMarkedContext, anchor::Vector{Float64})
    d = length(anchor)
    oracle = StanMarkedOracle(ctx, zeros(d), zeros(d), zeros(d), zeros(d), zeros(d))
    _clear_gradient!(ctx, :full, oracle.full_anchor, anchor)
    _clear_gradient!(ctx, :prior, oracle.prior_anchor, anchor)
    return oracle
end

function (oracle::StanMarkedOracle)(out, x)
    _clear_gradient!(oracle.ctx, :prior, oracle.prior_x, x)
    @. out = oracle.full_anchor + oracle.prior_x - oracle.prior_anchor
    return out
end

function (oracle::StanMarkedOracle)(out, x, subset, anchor)
    _selected_gradient!(oracle.ctx, oracle.selected_x, x, subset)
    _selected_gradient!(oracle.ctx, oracle.selected_anchor, anchor, subset)
    _clear_gradient!(oracle.ctx, :prior, oracle.prior_x, x)
    @. out = (oracle.selected_x - oracle.prior_x) -
        (oracle.selected_anchor - oracle.prior_anchor)
    return out
end

mutable struct StanMarkedSlabDeterministic{O,S}
    oracle::O
    slab::S
    active::BitVector
    all_active::BitVector
    all_buf::Vector{Float64}
    active_buf::Vector{Float64}
end

function StanMarkedSlabDeterministic(oracle, slab, d)
    p = length(PDMPSamplers.beta_indices(slab))
    return StanMarkedSlabDeterministic(oracle, slab, trues(p), trues(p),
        zeros(d), zeros(d))
end

function (target::StanMarkedSlabDeterministic)(out, x)
    target.oracle(out, x)
    PDMPSamplers.active_prior_grad!(target.slab, target.all_buf, x,
        target.all_active)
    PDMPSamplers.active_prior_grad!(target.slab, target.active_buf, x,
        target.active)
    @. out = out - target.all_buf + target.active_buf
    return out
end

function PDMPSamplers.set_active_set!(target::StanMarkedSlabDeterministic,
        free::BitVector)
    indices = PDMPSamplers.beta_indices(target.slab)
    @inbounds for j in eachindex(indices)
        target.active[j] = free[indices[j]]
    end
    return nothing
end

function _new_stan_marked_context(lib_path::String, full_data::String,
        prior_data::String, m::Int)
    sm_full = BridgeStan.StanModel(lib_path, full_data; warn=false)
    sm_prior = BridgeStan.StanModel(lib_path, prior_data; warn=false)
    sm_selected = BridgeStan.StanModel(lib_path, full_data; warn=false)
    names = BridgeStan.param_unc_names(sm_full)
    BridgeStan.param_unc_names(sm_prior) == names || throw(ArgumentError(
        "full and prior-only data must produce identical unconstrained parameter names"))
    BridgeStan.param_unc_names(sm_selected) == names || throw(ArgumentError(
        "selected and full contexts must have identical unconstrained parameter names"))
    set_fn, clear_fn = _resolve_subset_hooks(sm_selected)
    ctx = StanMarkedContext(PDMPModel(sm_full), PDMPModel(sm_prior),
        PDMPModel(sm_selected), set_fn, clear_fn, Vector{Int32}(undef, m), m,
        StanMarkedCallCounts())
    _clear_stan_subset!(ctx)
    return ctx, String.(names)
end

function _build_stan_marked_model(lib_path, full_data, prior_data,
        marked, anchor, slab_prior, model_prior, can_stick)
    N = Int(_rget(marked, :n_observations))
    m = Int(_rget(marked, :subsample_size))
    envelope_spec = _rget(marked, :residual_envelope)
    weights = _as_float_matrix(_rget(envelope_spec, :weights))
    size(weights, 2) == N || throw(DimensionMismatch(
        "residual envelope must have N observation columns"))
    growth = _as_float_vector(_rget(envelope_spec, :growth_rates))
    ctx, unc_names = _new_stan_marked_context(lib_path, full_data, prior_data, m)
    d = ctx.full.d
    length(anchor) == d || throw(DimensionMismatch("marked anchor has the wrong dimension"))
    envelope = TrajectoryResidualEnvelope(weights, anchor; growth_rates=growth)
    oracle = StanMarkedOracle(ctx, anchor)
    deterministic = if isnothing(slab_prior)
        oracle
    else
        provider = build_slab_provider(slab_prior, unc_names, can_stick, d)
        StanMarkedSlabDeterministic(oracle, provider, d)
    end
    cv = MarkedControlVariate(deterministic, oracle, envelope, anchor, m)
    return PDMPModel(d, cv), ctx, unc_names
end

function r_pdmp_stan_marked(lib_path::String, full_data::String,
        prior_data::String, marked, x0, flow_type::String,
        algorithm_type::String, flow_mean, flow_cov;
        c0::Float64=1e-2, grid_n::Int=30, grid_t_max::Float64=2.0,
        grid_bound::String="constant",
        grid_curvature_bound::Union{Nothing,Float64}=nothing,
        linear_area_threshold::Float64=0.95,
        linear_min_area_gain::Float64=0.0,
        t0::Float64=0.0, T::Float64=10000.0, t_warmup::Float64=0.0,
        sticky::Bool=false, can_stick=nothing, model_prior=nothing,
        parameter_prior=nothing, slab_prior=nothing,
        show_progress::Bool=true, n_chains::Int=1, threaded::Bool=false,
        seed::Union{Integer,Nothing}=nothing,
        adaptive_scheme::String="diagonal")
    algorithm_type in ("GridThinningStrategy", "ThinningStrategy") ||
        throw(ArgumentError("marked custom-Stan sampling requires GridThinningStrategy or ThinningStrategy"))
    x0_vec = _as_float_vector(x0)
    marked_anchor = _rget(marked, :anchor)
    anchor = isnothing(marked_anchor) ? copy(x0_vec) : _as_float_vector(marked_anchor)
    can_stick_vec = _as_bool_vector(can_stick)
    models = PDMPModel[]
    contexts = StanMarkedContext[]
    unc_names = String[]
    for _ in 1:n_chains
        model, ctx, names = _build_stan_marked_model(lib_path, full_data,
            prior_data, marked, anchor, slab_prior, model_prior, can_stick_vec)
        push!(models, model)
        push!(contexts, ctx)
        unc_names = names
    end
    d = first(models).d
    flow = build_flow(flow_type, _to_precision(_as_flow_cov(flow_cov, d), d),
        _as_flow_mean(flow_mean, d); adaptive_scheme)
    alg0 = build_algorithm(algorithm_type; c0, d, grid_n, grid_t_max,
        grid_bound, grid_curvature_bound, linear_area_threshold,
        linear_min_area_gain)
    alg = isnothing(slab_prior) ?
        wrap_sticky(alg0, sticky, model_prior,
            isnothing(parameter_prior) ? nothing : _as_float_vector(parameter_prior),
            can_stick_vec) :
        wrap_dependent_sticky(alg0, sticky, model_prior, slab_prior,
            can_stick_vec, flow_type, unc_names)
    chains = pdmp_sample(x0_vec, flow, models, alg, t0, T, t_warmup;
        progress=show_progress, threaded, seed,
        statistic_counter=PDMPSamplers.DevelStatisticCounter)
    result = _pack_result(chains)
    result["marked_context_counters"] = [Dict(
        "full_gradient_calls" => ctx.counts.full_gradient,
        "prior_gradient_calls" => ctx.counts.prior_gradient,
        "selected_gradient_calls" => ctx.counts.selected_gradient,
        "persons_evaluated" => ctx.counts.persons_evaluated,
        "model_constructions" => ctx.counts.model_constructions) for ctx in contexts]
    result["marked_subsampling"] = true
    return result
end

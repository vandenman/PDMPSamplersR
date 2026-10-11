module PDMPSamplersRBridge

using PDMPSamplers, LinearAlgebra, BridgeStan, Random, SparseArrays, Statistics, Libdl
using FileWatching
const LogExpFunctions = PDMPSamplers.LogExpFunctions

# Research-only timing/allocation diagnostics are deliberately opt-in.  The
# ordinary production path still retains event/provenance counters, but does
# not query Julia's GC allocator counters at every instrumented operation.
const _OMRF_RESEARCH_DIAGNOSTICS = lowercase(get(ENV,
    "OMRF_RESEARCH_DIAGNOSTICS", "false")) in ("1", "true", "yes", "on")

# BridgeStan/Stan Math autodiff state is serialized consistently for every
# Stan-backed provider used by this bridge.
const _BRIDGESTAN_CALL_LOCK = ReentrantLock()

export build_flow, build_algorithm, wrap_sticky
export build_model_prior_odds, build_slab_provider, wrap_dependent_sticky
export r_discretize, r_mean, r_var, r_std, r_cov, r_cor, r_quantile, r_median, r_cdf, r_ess, r_summary_all
export r_inclusion_probs, extract_stats
export r_chain_times, r_chain_positions, r_chain_velocities, r_chain_free_masks, r_chain_is_boomerang, r_chain_is_mutable_boomerang, r_chain_mu, r_chain_flow_cov
export r_chain_terminal_time, r_chain_terminal_position,
    r_chain_terminal_physical_velocity, r_chain_terminal_stored_velocity,
    r_chain_terminal_free_mask, r_chain_terminal_flow_spec
export r_chain_initial_time, r_chain_initial_position,
    r_chain_initial_physical_velocity, r_chain_initial_stored_velocity,
    r_chain_initial_free_mask, r_chain_initial_flow_spec
export r_chain_retained_initial_time, r_chain_retained_initial_position,
    r_chain_retained_initial_physical_velocity,
    r_chain_retained_initial_stored_velocity,
    r_chain_retained_initial_free_mask, r_chain_retained_initial_flow_spec
export r_resolve_flow_spec, r_chain_flow_spec
export r_from_skeleton
export r_chain_is_factorized
export r_chain_is_streaming, r_chain_streaming_chunk_paths,
    r_chain_streaming_chunk_events, r_chain_streaming_chunk_bytes,
    r_chain_streaming_directory, r_chain_streaming_chunk_count,
    r_chain_streaming_output_bytes, r_chain_streaming_flow_cov,
    r_chain_streaming_physical_events,
    r_chain_streaming_computational_boundaries,
    r_chain_streaming_buffer_maximum,
    r_chain_streaming_reflection_jump_summary
export r_chain_sparse_initial_time, r_chain_sparse_initial_position, r_chain_sparse_initial_velocity
export r_chain_sparse_event_indices, r_chain_sparse_event_times, r_chain_sparse_event_positions, r_chain_sparse_event_velocities
export r_chain_sparse_initial_free_mask, r_chain_sparse_event_free_values
export r_from_sparse_skeleton
export r_from_streaming_skeleton
export r_streaming_scientific_reduction, r_streaming_replay_validation,
    r_streaming_dense_test_reference
export r_chain_has_live_state_provenance
export r_pdmp_stan, r_pdmp_custom
export r_stan_subsampling_diagnostics, prepare_stan_subsampling,
    r_pdmp_stan_subsampling
export write_cmdstan_csv, r_constrain_and_write_csv
export r_pdmp_brms_subsampling, r_pdmp_stan_for_brms
export r_get_param_unc_names
export r_threading_available
export r_cleanup_worker_state

# Condition-specific objects are assigned in Main by JuliaCall for historical
# API compatibility.  A Mirai daemon is reused across conditions, so those
# roots must be cleared explicitly after the checkpoint is written.  Keeping
# the list here gives the R scheduler a stable cleanup API without duplicating
# Julia variable knowledge in simulation code.
const _WORKER_STATE_ROOTS = (
    :_stan_model, :_pdmp_model, :_prepared_stan_subsampling, :_subsampling,
    :_path_to_stan_model, :_path_to_stan_data,
    :_path_to_prior_stan_data_full, :_n_chains, :_curvature_backend,
    :flow, :algorithm, :x0, :theta0, :free0, :stored_velocity0,
    :flow_cov, :flow_mean, :c0, :grid_n, :grid_t_max,
    :grid_bound, :grid_curvature_bound, :t0, :T, :t_warmup, :f, :hessian_f, :hessian, :sticky, :can_stick, :model_prior,
    :parameter_prior, :slab_prior, :show_progress, :seed, :adaptive_scheme,
    :refresh_rate, :unc_names, :_unc_names, :subsampling_anchor_updates,
    :subsampling_use_anchor_bank,
    :subsampling_anchor_bank_capacity, :subsampling_main_anchor_refresh_distance,
    :subsampling_require_no_main_full_data, :subsampling_diagnostic_window,
    :warmup_adaptation_interval, :subsampling_diagnostic_window,
    :boomerang_variance_floor_fraction,
    :support_boundary_mode, :support_boundary_max_bisection_steps,
    :support_boundary_time_rtol, :support_boundary_time_atol,
    :support_boundary_clip_fraction, :support_boundary_max_refresh_attempts,
    :support_boundary_refresh_probe_time, :support_boundary_min_safe_time,
    :_diagnostic_model, :_diagnostic_full_data, :_diagnostic_prior_data,
    :_diagnostic_subsampling, :_diagnostic_position, :_diagnostic_subset,
    :_diagnostic_velocity, :_diagnostic_flow, :_diagnostic_flow_mean,
    :_diagnostic_flow_cov, :_diagnostic_phase)

function r_cleanup_worker_state()
    for name in _WORKER_STATE_ROOTS
        if isdefined(Main, name)
            Core.eval(Main, Expr(:(=), name, nothing))
        end
    end
    GC.gc(true)
    return nothing
end

function _boomerang_logscale_floor_mask(unc_names, d::Integer)
    names = String.(unc_names)
    length(names) == d || throw(DimensionMismatch(
        "parameter-name count $(length(names)) does not match dimension $d"))
    # The research variance safeguard applies to the hierarchical log-scale
    # parameters of the OMRF Stan models.
    return BitVector(name == "global_scale" || startswith(name, "node_scale.")
        for name in names)
end

_subsampling_collection_values(values::Union{AbstractVector,Tuple}) = values
_subsampling_collection_values(values::NamedTuple) = Base.values(values)
_subsampling_collection_values(values::AbstractDict) = Base.values(values)
_as_subsampling_predictor_indices(index::Number) = [Int(index)]
_as_subsampling_predictor_indices(indices) = collect(Int, vec(indices))
_subsampling_float_vector(value::Number) = [Float64(value)]
_subsampling_float_vector(values) = collect(Float64, vec(values))

function build_flow(flow_type::String, prec::AbstractMatrix{Float64}, flow_mean::AbstractVector{Float64};
                    adaptive_scheme::String="diagonal", refresh_rate::Float64=1.0)
    if flow_type == "ZigZag"
        return ZigZag(prec, flow_mean)
    elseif flow_type == "BouncyParticle"
        return BouncyParticle(prec, flow_mean, refresh_rate)
    elseif flow_type == "Boomerang"
        return Boomerang(prec, flow_mean)
    elseif flow_type == "AdaptiveBoomerang"
        λref = parse(Float64, get(ENV, "PDMP_ADAPTIVE_BOOMERANG_LAMBDA_REF", "0.1"))
        # `flow_mean` and `prec` are an initial adaptive specification, not a
        # promise that warmup leaves the flow unchanged.  Starting from the
        # supplied values is nevertheless essential for honest provenance and
        # for reproducible adaptive initialization.
        return AdaptiveBoomerang(prec, flow_mean; λref=λref)
    elseif flow_type == "PreconditionedZigZag"
        scale = prec isa Diagonal ? inv.(sqrt.(prec.diag)) :
            ones(Float64, length(flow_mean))
        return PreconditionedZigZag(prec, flow_mean; scale)
    elseif flow_type == "PreconditionedBPS"
        scale = prec isa Diagonal ? inv.(sqrt.(prec.diag)) :
            ones(Float64, length(flow_mean))
        return PreconditionedBPS(prec, flow_mean; refresh_rate, scale)
    elseif flow_type == "DensePreconditionedZigZag"
        flow = DensePreconditionedZigZag(prec, flow_mean)
        set_dense_preconditioner!(flow.metric,
            cholesky(inv(Symmetric(prec))).L)
        return flow
    elseif flow_type == "DensePreconditionedBPS"
        flow = DensePreconditionedBPS(prec, flow_mean)
        set_dense_preconditioner!(flow.metric,
            cholesky(inv(Symmetric(prec))).L)
        return flow
    else
        throw(ArgumentError("Unknown flow type: $flow_type"))
    end
end

# Full-gradient warmup and ordinary full sampling must resolve the same target
# and sticky event wrapper.  Keep these constructors independent of the
# retained subsampling representation so the switching-phase path cannot
# silently acquire a different prior or boundary clock.
function _build_sticky_target(base_model::PDMPModel, model_prior, slab_prior,
        can_stick_vec, unc_names_vec)
    return isnothing(slab_prior) ? base_model :
        _build_dependent_slab_model(base_model, model_prior, slab_prior,
            can_stick_vec, unc_names_vec)
end

function _build_sticky_event_algorithm(alg0, sticky::Bool, model_prior,
        parameter_prior_vec, slab_prior, can_stick_vec, flow, unc_names_vec)
    return isnothing(slab_prior) ?
        wrap_sticky(alg0, sticky, model_prior, parameter_prior_vec,
            can_stick_vec) :
        wrap_dependent_sticky(alg0, sticky, model_prior, slab_prior,
            can_stick_vec, flow, unc_names_vec)
end

function _resolved_flow_mean(flow)
    base = PDMPSamplers.underlying_flow(flow)
    return hasproperty(base, :μ) ? Vector{Float64}(getproperty(base, :μ)) : Float64[]
end


function _resolved_flow_covariance(flow)
    d = length(_resolved_flow_mean(flow))
    if flow isa PDMPSamplers.PreconditionedDynamics
        metric = flow.metric
        if metric isa PDMPSamplers.DiagonalPreconditioner
            return Matrix{Float64}(Diagonal(abs2.(metric.scale)))
        elseif metric isa PDMPSamplers.DensePreconditioner
            return Matrix{Float64}(metric.L * transpose(metric.L))
        elseif metric isa PDMPSamplers.IdentityPreconditioner
            return Matrix{Float64}(I, d, d)
        end
    end
    base = PDMPSamplers.underlying_flow(flow)
    if hasproperty(base, :ΣL)
        factor = getproperty(base, :ΣL)
        return Matrix{Float64}(factor * transpose(factor))
    elseif hasproperty(base, :Γ)
        return Matrix{Float64}(inv(Symmetric(Matrix(getproperty(base, :Γ)))))
    end
    return Matrix{Float64}(I, d, d)
end


function _flow_spec(flow)
    Dict("flow_type" => string(typeof(flow)),
         "inner_flow_type" => string(typeof(PDMPSamplers.underlying_flow(flow))),
         "mean" => _resolved_flow_mean(flow),
         "covariance" => _resolved_flow_covariance(flow))
end


function r_resolve_flow_spec(flow_type::String, flow_mean, flow_cov;
        adaptive_scheme::String="diagonal", refresh_rate::Float64=1.0)
    d = length(flow_mean)
    covariance = _as_flow_cov(flow_cov, d)
    mean = _as_flow_mean(flow_mean, d)
    flow = build_flow(flow_type, _to_precision(covariance, d), mean;
                      adaptive_scheme, refresh_rate)
    return _flow_spec(flow)
end

function build_algorithm(algorithm_type::String; c0::Float64, d::Integer, grid_n::Int, grid_t_max::Float64,
        curvature_backend::String="auto",
        grid_bound::String="value_quadratic",
        grid_curvature_bound::Union{Nothing,Float64}=nothing,
        lazy_low_tightness_threshold::Float64=0.1,
        lazy_max_low_tightness_rejections::Int=3,
        lazy_max_rejections::Int=0)
    if algorithm_type == "ThinningStrategy"
        return ThinningStrategy(GlobalBounds(c0, d))
    elseif algorithm_type == "GridThinningStrategy"
        grid_n_min = parse(Int, get(ENV, "PDMP_GRID_N_MIN", string(min(grid_n, 5))))
        bound_violation_value = get(ENV, "PDMP_GRID_BOUND_VIOLATION", "")
        bound_violation = isempty(bound_violation_value) ? nothing :
            Symbol(bound_violation_value)
        return GridThinningStrategy(; N = grid_n, t_max = grid_t_max,
            N_min = grid_n_min,
            curvature_backend = Symbol(curvature_backend),
            bound = Symbol(grid_bound),
            bound_violation = bound_violation,
            curvature_bound = grid_curvature_bound,
            lazy_low_tightness_threshold = lazy_low_tightness_threshold,
            lazy_max_low_tightness_rejections = lazy_max_low_tightness_rejections,
            lazy_max_rejections = lazy_max_rejections)
    elseif algorithm_type == "RootsPoissonStrategy"
        return RootsPoissonTimeStrategy()
    else
        throw(ArgumentError("Unknown algorithm type: $algorithm_type"))
    end
end

function wrap_sticky(alg::PDMPSamplers.PoissonTimeStrategy, sticky::Bool, model_prior, parameter_prior, can_stick)
    !sticky && return alg

    if haskey(model_prior, :prob)
        w = model_prior[:prob] ./ (1 .- model_prior[:prob])
        κ = w .* parameter_prior
    else
        κ = BetaBernoulliKappa(model_prior[:a]::Float64,
            model_prior[:b]::Float64, parameter_prior, can_stick)
    end

    return Sticky(alg, κ, BitVector(can_stick))
end

function build_warmup_stop(t_warmup::Float64)
    enabled = parse(Bool, get(ENV, "PDMP_ADAPTIVE_WARMUP_STOP", "false"))
    enabled || return nothing
    min_time = parse(Float64, get(ENV, "PDMP_ADAPTIVE_WARMUP_MIN_TIME", string(0.25 * t_warmup)))
    max_time = parse(Float64, get(ENV, "PDMP_ADAPTIVE_WARMUP_MAX_TIME", string(t_warmup)))
    stable_time = parse(Float64, get(ENV, "PDMP_ADAPTIVE_WARMUP_STABLE_TIME", string(0.25 * t_warmup)))
    min_events = parse(Int, get(ENV, "PDMP_ADAPTIVE_WARMUP_MIN_EVENTS", "100"))
    check_every = parse(Int, get(ENV, "PDMP_ADAPTIVE_WARMUP_CHECK_EVERY", "25"))
    return PDMPSamplers.AdaptiveWarmupCriterion(;
        min_time, max_time, stable_time, min_events, check_every)
end

_haskey(x, key::Symbol) = haskey(x, key) || haskey(x, String(key))
_rget(x, key::Symbol) = haskey(x, key) ? x[key] : x[String(key)]
_as_float_vector(x::Number) = [Float64(x)]
_as_float_vector(x) = Vector{Float64}(x)
_as_int_vector(x::Number) = [Int(x)]
_as_int_vector(x) = Vector{Int}(x)
_as_bool_vector(x::Bool) = Bool[x]
_as_bool_vector(x::Nothing) = nothing
_as_bool_vector(x) = Bool.(x)
_as_string_vector(x::AbstractString) = String[String(x)]
_as_string_vector(x::Nothing) = String[]
_as_string_vector(x) = String.(x)
_as_float_matrix(x::Number) = reshape([Float64(x)], 1, 1)
_as_float_matrix(x) = Matrix{Float64}(x)

function _as_float_sparse_rows(spec)
    dims = _as_int_vector(_rget(spec, :dims))
    length(dims) == 2 || throw(DimensionMismatch("sparse design dims must have length two"))
    rows = _as_int_vector(_rget(spec, :i))
    columns = _as_int_vector(_rget(spec, :j))
    values = _as_float_vector(_rget(spec, :x))
    length(rows) == length(columns) == length(values) || throw(DimensionMismatch(
        "sparse design i, j, and x vectors must have equal length"))
    rowptr = zeros(Int, dims[1] + 1)
    for row in rows
        1 <= row <= dims[1] || throw(ArgumentError("sparse design row index is out of bounds"))
        rowptr[row + 1] += 1
    end
    cumsum!(rowptr, rowptr)
    rowptr .+= 1
    colidx = Vector{Int}(undef, length(columns))
    nzval = Vector{Float64}(undef, length(values))
    nextptr = copy(rowptr)
    for k in eachindex(rows, columns, values)
        1 <= columns[k] <= dims[2] || throw(ArgumentError(
            "sparse design column index is out of bounds"))
        destination = nextptr[rows[k]]
        colidx[destination] = columns[k]
        nzval[destination] = values[k]
        nextptr[rows[k]] += 1
    end
    return dims, rowptr, colidx, nzval
end

function _slab_type(slab_prior)
    t = _rget(slab_prior, :type)
    return String(t)
end


function _resolve_unc_spec(value, unc_names::AbstractVector{<:AbstractString}, label)
    values = value isa AbstractString ? [String(value)] : String.(value)
    names = String.(unc_names)
    idx = Int[]
    for requested in values
        exact = findall(==(requested), names)
        matches = isempty(exact) ? findall(name ->
            startswith(name, requested * ".") || startswith(name, requested * "["), names) : exact
        isempty(matches) && throw(ArgumentError(
            "$label name or block prefix $requested was not found among unconstrained parameter names"))
        append!(idx, matches)
    end
    allunique(idx) || throw(ArgumentError("resolved $label coordinates contain duplicates"))
    return idx
end

function _coef_indices(slab_prior, unc_names::AbstractVector{<:AbstractString},
        can_stick::AbstractVector{Bool}, d::Integer)
    coef = _rget(slab_prior, :coef)
    if isnothing(coef)
        idx = findall(Bool.(can_stick))
        isempty(idx) && throw(ArgumentError("slab_prior requires at least one stickable coefficient"))
        return Int.(idx)
    end
    if coef isa AbstractString
        isempty(unc_names) && throw(ArgumentError("character slab coef requires unconstrained parameter names"))
        idx = _resolve_unc_spec(coef, unc_names, "slab coef")
        any(i -> !can_stick[i], idx) && throw(ArgumentError("explicit slab coef must be a subset of can_stick coordinates"))
        return idx
    end
    if coef isa AbstractVector{<:AbstractString}
        isempty(unc_names) && throw(ArgumentError("character slab coef requires unconstrained parameter names"))
        idx = _resolve_unc_spec(coef, unc_names, "slab coef")
        any(i -> !can_stick[i], idx) &&
            throw(ArgumentError("explicit slab coef must be a subset of can_stick coordinates"))
        return idx
    end
    idx = coef isa Integer ? [Int(coef)] : Int.(coef)
    any(i -> i < 1 || i > d, idx) && throw(ArgumentError("integer slab coef indices must lie in 1:d"))
    any(i -> !can_stick[i], idx) &&
        throw(ArgumentError("explicit slab coef must be a subset of can_stick coordinates"))
    return idx
end

function _state_indices(value, unc_names::AbstractVector{<:AbstractString}, d::Integer, label::AbstractString)
    if value isa AbstractString
        isempty(unc_names) && throw(ArgumentError("character $label requires unconstrained parameter names"))
        return _resolve_unc_spec(value, unc_names, label)
    elseif value isa AbstractVector{<:AbstractString}
        isempty(unc_names) && throw(ArgumentError("character $label requires unconstrained parameter names"))
        return _resolve_unc_spec(value, unc_names, label)
    end
    idx = value isa Integer ? [Int(value)] : Int.(value)
    any(i -> i < 1 || i > d, idx) && throw(ArgumentError("$label indices must lie in 1:d"))
    return idx
end

function build_model_prior_odds(model_prior, beta_indices::AbstractVector{Int}, d::Integer)
    m = length(beta_indices)
    if _haskey(model_prior, :prob)
        prob = _as_float_vector(_rget(model_prior, :prob))
        if length(prob) == 1
            return PDMPSamplers.BernoulliModelPrior(fill(prob[1], m))
        elseif length(prob) == m
            return PDMPSamplers.BernoulliModelPrior(prob)
        elseif length(prob) == d
            return PDMPSamplers.BernoulliModelPrior(prob[beta_indices])
        else
            throw(DimensionMismatch("Bernoulli model prior length must be 1, beta dimension, or full dimension"))
        end
    elseif _haskey(model_prior, :a) && _haskey(model_prior, :b)
        return PDMPSamplers.BetaBernoulliModelPrior(
            m, Float64(_rget(model_prior, :a)), Float64(_rget(model_prior, :b)))
    else
        throw(ArgumentError("unknown model_prior specification"))
    end
end

function build_slab_provider(slab_prior, unc_names::AbstractVector{<:AbstractString},
        can_stick::AbstractVector{Bool}, d::Integer)
    beta_idx = _coef_indices(slab_prior, unc_names, can_stick, d)
    t = _slab_type(slab_prior)
    if t == "dense_gaussian"
        return DenseGaussianSlab(_as_float_vector(_rget(slab_prior, :mean)),
            _as_float_matrix(_rget(slab_prior, :cov)), beta_idx)
    elseif t == "exchangeable_gaussian"
        if Bool(_rget(slab_prior, :zero_mean))
            return ZeroMeanExchangeableGaussianSlab(beta_idx,
                Float64(_rget(slab_prior, :u)), Float64(_rget(slab_prior, :v)))
        else
            return ExchangeableGaussianSlab(beta_idx, Float64(_rget(slab_prior, :mean)),
                Float64(_rget(slab_prior, :u)), Float64(_rget(slab_prior, :v)))
        end
    elseif t == "loglinear_gaussian_scale"
        log_base_scales = _as_float_vector(_rget(slab_prior, :log_base_scales))
        length(log_base_scales) == 1 &&
            (log_base_scales = fill(log_base_scales[1], length(beta_idx)))
        logscale_idx = _state_indices(_rget(slab_prior, :logscale), unc_names, d, "logscale")
        design = _rget(slab_prior, :logscale_design)
        storage = String(_rget(design, :storage))
        dims = _as_int_vector(_rget(design, :dims))
        dims == [length(beta_idx), length(logscale_idx)] ||
            throw(DimensionMismatch("logscale_design must have size (beta dimension, logscale dimension)"))
        isempty(intersect(beta_idx, logscale_idx)) || throw(ArgumentError(
            "loglinear_gaussian_scale_slab logscale coordinates must be disjoint from slab coef coordinates"))
        any(i -> can_stick[i], logscale_idx) && throw(ArgumentError(
            "loglinear_gaussian_scale_slab logscale coordinates must be non-stickable"))
        if storage == "dense"
            values = _as_float_vector(_rget(design, :x))
            length(values) == prod(dims) || throw(DimensionMismatch(
                "dense logscale_design value count does not match its dimensions"))
            return LogLinearGaussianScaleSlab(beta_idx, logscale_idx,
                log_base_scales, reshape(values, dims...))
        elseif storage == "sparse_rows"
            _, rowptr, colidx, nzval = _as_float_sparse_rows(design)
            return LogLinearGaussianScaleSlab(beta_idx, logscale_idx,
                log_base_scales, rowptr, colidx, nzval)
        end
        throw(ArgumentError("unknown logscale_design storage $storage"))
    else
        throw(ArgumentError("unknown slab_prior type: $t"))
    end
end

"""
Wrap a dependent-slab sampler using the concrete constructed dynamics.

Clock selection is flow-aware: in particular, Boomerang-family dynamics with
log-linear Gaussian scale slabs use the harmonic clock. Passing an R flow name
here would discard preconditioning and adaptive-dynamics information.
"""
function wrap_dependent_sticky(alg::PDMPSamplers.PoissonTimeStrategy, sticky::Bool,
        model_prior, slab_prior, can_stick,
        flow::PDMPSamplers.ContinuousDynamics,
        unc_names::AbstractVector{<:AbstractString}=String[])
    sticky || throw(ArgumentError("wrap_dependent_sticky requires sticky=true"))
    isnothing(slab_prior) && throw(ArgumentError("wrap_dependent_sticky requires slab_prior"))
    d = length(can_stick)
    provider = build_slab_provider(slab_prior, unc_names, can_stick, d)
    odds = build_model_prior_odds(model_prior, PDMPSamplers.beta_indices(provider), d)
    clock = default_aggregate_unstick_clock(provider, odds, flow)
    return AggregateSticky(alg, clock, BitVector(can_stick))
end

function _build_dependent_slab_model(
        posterior_model::PDMPModel,
        model_prior,
        slab_prior,
        can_stick,
        unc_names::AbstractVector{<:AbstractString}=String[])
    d = posterior_model.d
    provider = build_slab_provider(slab_prior, unc_names, can_stick, d)
    odds = build_model_prior_odds(model_prior, PDMPSamplers.beta_indices(provider), d)
    target = DependentSlabTarget(d, posterior_model.grad, provider, odds;
        posterior_hvp=posterior_model.hvp)
    return PDMPModel(target; hvp=true)
end

function _to_precision(flow_cov::AbstractMatrix{Float64}, d::Int)
    if isdiag(flow_cov)
        diagonal = Vector{Float64}(undef, d)
        @inbounds for i in 1:d
            diagonal[i] = inv(flow_cov[i, i])
        end
        return Diagonal(diagonal)
    end
    return inv(Symmetric(flow_cov))
end

function _as_flow_mean(flow_mean, d::Int)
    if flow_mean isa Number
        return [Float64(flow_mean)]
    end
    isempty(flow_mean) && return zeros(d)
    return Vector{Float64}(flow_mean)
end

function _as_flow_cov(flow_cov, d::Int)
    if flow_cov isa Number
        return reshape([Float64(flow_cov)], 1, 1)
    end
    isempty(flow_cov) && return Matrix{Float64}(I, d, d)
    return Matrix{Float64}(flow_cov)
end

_compile_env_enabled(name::String) = lowercase(get(ENV, name, "false")) in ("1", "true", "yes", "y")

function _compile_cache_checksum(bytes::Vector{UInt8})
    value = UInt64(0xcbf29ce484222325)
    @inbounds for byte in bytes
        value = (value ⊻ UInt64(byte)) * UInt64(0x100000001b3)
    end
    return value
end

# The compiled BridgeStan library is an ABI artifact, not just a function of
# the Stan source.  Include an explicit schema marker so caches produced by an
# older BridgeStan/PDMPSamplersR build cannot be loaded after the toolchain or
# generated-model ABI changes.
const _BRIDGESTAN_CACHE_SCHEMA = "pdmpsamplersr-bridgestan-cache-v2"

function _cache_build_dir(path_to_stan_model::String,
        stanc_args::Vector{String}, make_args::Vector{String})
    cache_root = get(ENV, "PDMPSAMPLERSR_BRIDGESTAN_CACHE_DIR",
        joinpath(first(DEPOT_PATH), "pdmpsamplersr", "bridgestan"))
    source_bytes = read(path_to_stan_model)
    option_bytes = Vector{UInt8}(codeunits(join([stanc_args; make_args], '\0')))
    schema_bytes = Vector{UInt8}(codeunits(_BRIDGESTAN_CACHE_SCHEMA))
    checksum = _compile_cache_checksum([schema_bytes; source_bytes; option_bytes])
    return joinpath(cache_root, string(checksum; base=16, pad=16))
end

function _with_compile_lock(f::Function, lock_path::String)
    mkpath(dirname(lock_path))
    # Compilation can take minutes on a cold cache.  Refresh the pidfile while
    # waiting/compiling, and recover locks left by a killed worker after one
    # hour.  This lock is shared by independent R/mirai worker processes.
    FileWatching.Pidfile.mkpidlock(lock_path;
        stale_age=3600.0, refresh=600.0) do
        f()
    end
end

function _cached_stan_source(path_to_stan_model::String,
        stanc_args::Vector{String}, make_args::Vector{String},
        build_dir::String=_cache_build_dir(path_to_stan_model, stanc_args, make_args))
    mkpath(build_dir)
    cached_source = joinpath(build_dir, basename(path_to_stan_model))
    isfile(cached_source) || cp(path_to_stan_model, cached_source)
    return cached_source
end

function _compile_model(path_to_stan_model::String)
    !endswith(path_to_stan_model, ".stan") && return path_to_stan_model
    optimize = _compile_env_enabled("PDMPSAMPLERSR_BRIDGESTAN_STANC_O1")
    stanc_flags = optimize ? "--O1" : ""
    stanc_args = isempty(stanc_flags) ? String[] : [stanc_flags]
    make_args = String[]
    _compile_env_enabled("PDMPSAMPLERSR_BRIDGESTAN_NATIVE") &&
        push!(make_args, "CXXFLAGS+=-march=native")
    if _compile_env_enabled("PDMPSAMPLERSR_BRIDGESTAN_CACHE")
        build_dir = _cache_build_dir(path_to_stan_model, stanc_args, make_args)
        return _with_compile_lock(joinpath(build_dir, ".compile.lock")) do
            source = _cached_stan_source(path_to_stan_model, stanc_args, make_args, build_dir)
            library = replace(source, r"\.stan$" => "_model.$(Libdl.dlext)")
            isfile(library) && return library
            BridgeStan.compile_model(source; stanc_args, make_args)
        end
    end
    return _with_compile_lock(string(path_to_stan_model, ".compile.lock")) do
        library = replace(path_to_stan_model, r"\.stan$" => "_model.$(Libdl.dlext)")
        isfile(library) && return library
        BridgeStan.compile_model(path_to_stan_model; stanc_args, make_args)
    end
end

function _compile_model_with_header(path_to_stan_model::String, hpp_path::String)
    !endswith(path_to_stan_model, ".stan") && return path_to_stan_model
    header = replace(abspath(normpath(hpp_path)), "\\" => "/")
    stanc_args = ["--allow-undefined"]
    _compile_env_enabled("PDMPSAMPLERSR_BRIDGESTAN_STANC_O1") &&
        push!(stanc_args, "--O1")
    make_args = ["USER_HEADER=$(header)"]
    _compile_env_enabled("PDMPSAMPLERSR_BRIDGESTAN_NATIVE") &&
        push!(make_args, "CXXFLAGS+=-march=native")
    if _compile_env_enabled("PDMPSAMPLERSR_BRIDGESTAN_CACHE")
        build_dir = _cache_build_dir(path_to_stan_model, stanc_args, make_args)
        return _with_compile_lock(joinpath(build_dir, ".compile.lock")) do
            source = _cached_stan_source(path_to_stan_model, stanc_args, make_args, build_dir)
            library = replace(source, r"\.stan$" => "_model.$(Libdl.dlext)")
            isfile(library) && mtime(library) >=
                max(mtime(source), mtime(hpp_path)) && return library
            BridgeStan.compile_model(source; stanc_args, make_args)
        end
    end
    return _with_compile_lock(string(path_to_stan_model, ".compile.lock")) do
        library = replace(path_to_stan_model, r"\.stan$" => "_model.$(Libdl.dlext)")
        isfile(library) && mtime(library) >=
            max(mtime(path_to_stan_model), mtime(hpp_path)) && return library
        BridgeStan.compile_model(path_to_stan_model; stanc_args, make_args)
    end
end

function _flow_provenance(chains::PDMPChains, requested_flow::String,
        adaptive_scheme::String, initial_refresh_rate::Float64,
        warmup_adaptation_interval::Union{Nothing,Float64}=nothing,
        t_warmup::Float64=0.0)
    trace = chains.traces[1]
    flow = trace.flow
    underlying = PDMPSamplers.underlying_flow(flow)
    preconditioner = flow isa PDMPSamplers.PreconditionedDynamics ?
        string(typeof(flow.metric)) : "None"
    adaptation_class = if flow isa PDMPSamplers.PreconditionedDynamics &&
            flow.metric isa PDMPSamplers.DiagonalPreconditioner
        "diagonal_preconditioner_warmup"
    elseif underlying isa PDMPSamplers.MutableBoomerang
        "adaptive_boomerang_diagonal"
    else
        "none"
    end
    dynamics_adaptive = adaptation_class != "none"
    resolved_interval = isnothing(warmup_adaptation_interval) ?
        (t_warmup > 0 ? t_warmup / 10 : 0.0) : warmup_adaptation_interval
    return (
        requested_flow = requested_flow,
        concrete_flow_type = string(typeof(flow)),
        underlying_dynamics_type = string(typeof(underlying)),
        preconditioner_type = preconditioner,
        adaptive_scheme = adaptation_class == "adaptive_boomerang_diagonal" ?
            adaptive_scheme : "none",
        adaptation_class = adaptation_class,
        resolved_adaptation_interval = resolved_interval,
        dynamics_adaptive = dynamics_adaptive,
        initial_refresh_rate = initial_refresh_rate,
        final_refresh_rate = try Float64(PDMPSamplers.refresh_rate(flow)) catch; NaN end,
    )
end

# Streaming storage for the bridge. Online retained-phase summaries are opt-in
# through the environment, so research runners need no new arguments:
# PDMPSAMPLERS_ONLINE_BATCHES (0 = off), PDMPSAMPLERS_ONLINE_GRID_SPACING
# (physical time; NaN = no snapshots) and PDMPSAMPLERS_STREAM_EVENTS ("0" =
# write no retained-phase event chunks), PDMPSAMPLERS_ONLINE_WARMUP_BATCHES
# (0 = no warmup summaries). `T` is the absolute end time.
function _bridge_trace_storage(directory, buffer_events, T)
    batches = parse(Int, get(ENV, "PDMPSAMPLERS_ONLINE_BATCHES", "0"))
    return PDMPSamplers.StreamingTraceStorage(directory; buffer_events,
        write_events=get(ENV, "PDMPSAMPLERS_STREAM_EVENTS", "1") != "0",
        online_batches=batches,
        online_end_time=batches > 0 ? Float64(T) : NaN,
        online_grid_spacing=parse(Float64,
            get(ENV, "PDMPSAMPLERS_ONLINE_GRID_SPACING", "NaN")),
        online_warmup_batches=parse(Int,
            get(ENV, "PDMPSAMPLERS_ONLINE_WARMUP_BATCHES", "0")))
end

function _pack_result(chains::PDMPChains; requested_flow::String="",
        adaptive_scheme::String="none", initial_refresh_rate::Float64=NaN,
        warmup_adaptation_interval::Union{Nothing,Float64}=nothing,
        t_warmup::Float64=0.0)
    stats = extract_stats(chains)
    d = length(first(chains.traces[1]).position)
    provenance = _flow_provenance(chains, requested_flow, adaptive_scheme,
        initial_refresh_rate, warmup_adaptation_interval, t_warmup)
    n = length(chains.traces)
    stats["requested_flow"] = fill(provenance.requested_flow, n)
    stats["concrete_flow_type"] = fill(provenance.concrete_flow_type, n)
    stats["underlying_dynamics_type"] = fill(provenance.underlying_dynamics_type, n)
    stats["preconditioner_type"] = fill(provenance.preconditioner_type, n)
    stats["resolved_adaptive_scheme"] = fill(provenance.adaptive_scheme, n)
    stats["dynamics_adaptive"] = fill(
        provenance.dynamics_adaptive ? "true" : "false", n)
    stats["adaptation_class"] = fill(provenance.adaptation_class, n)
    stats["resolved_adaptation_interval"] = fill(
        provenance.resolved_adaptation_interval, n)
    stats["initial_refresh_rate"] = fill(provenance.initial_refresh_rate, n)
    stats["final_refresh_rate"] = fill(provenance.final_refresh_rate, n)
    # `extract_stats` uses a compact R-bridge value union (Float64 vectors,
    # String vectors, and matrices).  Store integer thread counts as numeric
    # vectors so provenance cannot make result packing fail after sampling.
    stats["julia_threads"] = fill(Float64(Threads.nthreads()), n)
    stats["blas_threads"] = fill(Float64(LinearAlgebra.BLAS.get_num_threads()), n)
    stats["openmp_threads"] = fill(get(ENV, "OMP_NUM_THREADS", "unset"), n)
    stats["mkl_threads"] = fill(get(ENV, "MKL_NUM_THREADS", "unset"), n)
    stats["openblas_threads"] = fill(get(ENV, "OPENBLAS_NUM_THREADS", "unset"), n)
    return Dict{String,Any}(
        "chains"   => chains,
        "stats"    => stats,
        "d"        => d,
        "n_chains" => PDMPSamplers.n_chains(chains),
        "flow_provenance" => Dict{String,Any}(
            "requested_flow" => provenance.requested_flow,
            "concrete_flow_type" => provenance.concrete_flow_type,
            "underlying_dynamics_type" => provenance.underlying_dynamics_type,
            "preconditioner_type" => provenance.preconditioner_type,
            "adaptive_scheme" => provenance.adaptive_scheme,
            "dynamics_adaptive" => provenance.dynamics_adaptive,
            "adaptation_class" => provenance.adaptation_class,
            "resolved_adaptation_interval" =>
                provenance.resolved_adaptation_interval,
            "initial_refresh_rate" => provenance.initial_refresh_rate,
            "final_refresh_rate" => provenance.final_refresh_rate)
    )
end

# ──────────────────────────────────────────────────────────────────────────────
# Bridge functions for R-side estimators
# ──────────────────────────────────────────────────────────────────────────────

function r_discretize(chains::PDMPChains; dt::Union{Float64, Nothing} = nothing, chain::Int = 1)
    trace = chains.traces[chain]
    if isnothing(dt)
        dt, _, _ = adaptive_dt(trace)
    end
    return Matrix(PDMPDiscretize(trace, dt))
end

r_mean(chains::PDMPChains; chain::Int = 1) = Statistics.mean(chains; chain)
r_var(chains::PDMPChains; chain::Int = 1)  = Statistics.var(chains; chain)
r_std(chains::PDMPChains; chain::Int = 1)  = Statistics.std(chains; chain)
r_cov(chains::PDMPChains; chain::Int = 1)  = Statistics.cov(chains; chain)
r_cor(chains::PDMPChains; chain::Int = 1)  = Statistics.cor(chains; chain)

function r_quantile(chains::PDMPChains, p::Float64; chain::Int = 1, coordinate::Int = -1)
    Statistics.quantile(chains, p; chain, coordinate)
end

"""
    r_quantile(chains, probabilities; chain=1, coordinate=-1)

Compute several quantiles in one position-domain sweep per coordinate.  The
scalar R bridge historically repeated the complete sweep for every requested
probability.  On long sticky traces this doubled multi-gigabyte temporary
allocations for an equal-tail interval.  Reclaim each coordinate's sweep
workspace before moving to the next coordinate so peak memory is bounded by
one marginal rather than by the full parameter dimension.
"""
function r_quantile(chains::PDMPChains, p::AbstractVector{<:Real};
        chain::Int = 1, coordinate::Int = -1)
    trace = chains.traces[chain]
    probabilities = collect(Float64, p)
    if coordinate != -1
        return Statistics.quantile(trace, probabilities; coordinate)
    end

    d = length(first(trace).position)
    result = Matrix{Float64}(undef, length(probabilities), d)
    gc_interval = max(1, tryparse(Int, get(ENV,
        "OMRF_QUANTILE_GC_INTERVAL", "16")) === nothing ? 16 :
        tryparse(Int, get(ENV, "OMRF_QUANTILE_GC_INTERVAL", "16")))
    for j in 1:d
        result[:, j] .= Statistics.quantile(trace, probabilities; coordinate=j)
        # Sweep events can be O(length(trace)) and are not reused by another
        # coordinate.  An explicit incremental collection prevents them from
        # accumulating until the enclosing all-coordinate call returns.
        (j == d || j % gc_interval == 0) && GC.gc(false)
    end
    return result
end

function r_median(chains::PDMPChains; chain::Int = 1, coordinate::Int = -1)
    Statistics.median(chains.traces[chain]; coordinate)
end

function r_cdf(chains::PDMPChains, q::Float64; chain::Int = 1, coordinate::Int)
    cdf(chains, q; chain, coordinate)
end

function r_ess(chains::PDMPChains; chain::Int = 1, n_batches::Int = 0)
    if n_batches > 0
        ess(chains; chain, n_batches)
    else
        ess(chains; chain)
    end
end

function r_summary_all(chains::PDMPChains; chain::Int = 1)
    Dict{String,Any}(
        "mean"   => Statistics.mean(chains; chain),
        "var"    => Statistics.var(chains; chain),
        "std"    => Statistics.std(chains; chain),
        "cov"    => Statistics.cov(chains; chain),
        "cor"    => Statistics.cor(chains; chain),
        "median" => Statistics.median(chains.traces[chain])
    )
end

r_inclusion_probs(chains::PDMPChains; chain::Int = 1) = inclusion_probs(chains; chain)

# ──────────────────────────────────────────────────────────────────────────────
# Transform construction and transformed estimators
# ──────────────────────────────────────────────────────────────────────────────

function _build_transforms(specs::AbstractVector)
    map(specs) do s
        t = s[:type]
        if t == "identity"
            PDMPSamplers.IdentityTransform()
        elseif t == "lower"
            PDMPSamplers.LowerBoundTransform(Float64(s[:lower]))
        elseif t == "upper"
            PDMPSamplers.UpperBoundTransform(Float64(s[:upper]))
        elseif t == "double"
            PDMPSamplers.DoubleBoundTransform(Float64(s[:lower]), Float64(s[:upper]))
        else
            error("Unknown transform type: $t")
        end
    end
end

function r_mean(chains::PDMPChains, specs::AbstractVector; chain::Int = 1)
    transforms = _build_transforms(specs)
    Statistics.mean(chains.traces[chain], transforms)
end

function r_var(chains::PDMPChains, specs::AbstractVector; chain::Int = 1)
    transforms = _build_transforms(specs)
    Statistics.var(chains.traces[chain], transforms)
end

function r_std(chains::PDMPChains, specs::AbstractVector; chain::Int = 1)
    transforms = _build_transforms(specs)
    Statistics.std(chains.traces[chain], transforms)
end

function r_quantile(chains::PDMPChains, p::Float64, specs::AbstractVector; chain::Int = 1, coordinate::Int = -1)
    transforms = _build_transforms(specs)
    Statistics.quantile(chains.traces[chain], p, transforms; coordinate)
end

function r_median(chains::PDMPChains, specs::AbstractVector; chain::Int = 1, coordinate::Int = -1)
    transforms = _build_transforms(specs)
    Statistics.median(chains.traces[chain], transforms; coordinate)
end

function r_cdf(chains::PDMPChains, q::Float64, specs::AbstractVector; chain::Int = 1, coordinate::Int)
    transforms = _build_transforms(specs)
    cdf(chains.traces[chain], q, transforms; coordinate)
end

# ──────────────────────────────────────────────────────────────────────────────
# Skeleton extraction and reconstruction for saveRDS support
# ──────────────────────────────────────────────────────────────────────────────
#
# JuliaCall only auto-converts simple Julia types to native R (scalars, plain
# Vector{Float64}, Matrix{Float64}). Complex types (Dict, Vector{Dict}) are
# returned as JuliaObject external pointers, which are invalidated by
# saveRDS/readRDS. To keep the skeleton as native R data, extraction is done
# field-by-field (one julia_call per field per chain), and reconstruction
# accepts 5 parallel flat lists — one per field.

function _dense_compact(chains::PDMPChains, chain::Int)
    trace = chains.traces[chain]
    dense = trace isa PDMPTrace ? trace : PDMPTrace(trace)
    PDMPSamplers.compact(dense)
end

function r_chain_is_factorized(chains::PDMPChains; chain::Int)
    chains.traces[chain] isa PDMPSamplers.FactorizedTrace
end

function r_chain_is_streaming(chains::PDMPChains; chain::Int)
    chains.traces[chain] isa PDMPSamplers.StreamingPDMPTrace
end

function r_chain_streaming_chunk_paths(chains::PDMPChains; chain::Int)
    trace = chains.traces[chain]
    trace isa PDMPSamplers.StreamingPDMPTrace || throw(ArgumentError(
        "requested a streaming manifest from a non-streaming trace"))
    return String[chunk.path for chunk in trace.chunks]
end

r_chain_streaming_directory(chains::PDMPChains; chain::Int) =
    String(chains.traces[chain].directory)
r_chain_streaming_chunk_count(chains::PDMPChains; chain::Int) =
    Float64(length(chains.traces[chain].chunks))
r_chain_streaming_output_bytes(chains::PDMPChains; chain::Int) =
    Float64(chains.traces[chain].chunks.total_bytes)
function r_chain_streaming_flow_cov(chains::PDMPChains; chain::Int)
    covariance = _resolved_flow_covariance(chains.traces[chain].flow)
    diagonal = diag(covariance)
    if all(i == j || iszero(covariance[i, j])
            for i in axes(covariance, 1), j in axes(covariance, 2))
        return Vector{Float64}(diagonal)
    end
    return Matrix{Float64}(covariance)
end

function r_chain_streaming_chunk_events(chains::PDMPChains; chain::Int)
    trace = chains.traces[chain]
    return Int32[chunk.n_events for chunk in trace.chunks]
end

function r_chain_streaming_chunk_bytes(chains::PDMPChains; chain::Int)
    trace = chains.traces[chain]
    return Float64[chunk.bytes for chunk in trace.chunks]
end

r_chain_streaming_physical_events(chains::PDMPChains; chain::Int) =
    Float64(chains.traces[chain].physical_events)
r_chain_streaming_computational_boundaries(chains::PDMPChains; chain::Int) =
    Float64(chains.traces[chain].computational_boundaries)
r_chain_streaming_buffer_maximum(chains::PDMPChains; chain::Int) =
    Float64(chains.traces[chain].buffer_maximum)

function r_chain_streaming_reflection_jump_summary(chains::PDMPChains;
        chain::Int)
    trace = chains.traces[chain]
    trace isa PDMPSamplers.StreamingPDMPTrace || throw(ArgumentError(
        "reflection-jump summary requires a streaming trace"))
    covariance = _resolved_flow_covariance(trace.flow)
    inverse_scale = 1.0 ./ sqrt.(diag(covariance))
    cursor = PDMPSamplers.stream_cursor(trace)
    count = 0
    total = 0.0
    maximum = 0.0
    for metadata in trace.chunks
        chunk = PDMPSamplers.load_streaming_chunk(metadata)
        @inbounds for event in eachindex(chunk.times)
            kind = chunk.kinds[event]
            jump = 0.0
            if kind == PDMPSamplers.STREAM_EVENT_REFLECT_FULL
                after = view(chunk.full_effective_velocities, :,
                    Int(chunk.full_slots[event]))
                for i in eachindex(cursor.velocity)
                    if cursor.free[i]
                        delta = (after[i] - cursor.velocity[i]) *
                            inverse_scale[i]
                        jump += delta * delta
                    end
                end
            elseif kind == PDMPSamplers.STREAM_EVENT_REFLECT_COORD
                i = Int(chunk.coordinates[event])
                delta = (chunk.velocities[event] - cursor.velocity[i]) *
                    inverse_scale[i]
                jump = delta * delta
            end
            if kind == PDMPSamplers.STREAM_EVENT_REFLECT_FULL ||
                    kind == PDMPSamplers.STREAM_EVENT_REFLECT_COORD
                count += 1
                total += jump
                maximum = max(maximum, jump)
            end
            PDMPSamplers.stream_apply_event!(cursor, chunk, event)
        end
    end
    return Float64[count, total, iszero(count) ? NaN : total / count,
        maximum]
end

r_chain_has_live_state_provenance(chains::PDMPChains; chain::Int=1) =
    length(chains.initial_states) >= chain &&
    length(chains.retained_initial_states) >= chain &&
    length(chains.terminal_states) >= chain

function r_chain_sparse_initial_time(chains::PDMPChains; chain::Int)
    Float64(chains.traces[chain].initial_state.time)
end

function r_chain_sparse_initial_position(chains::PDMPChains; chain::Int)
    Vector{Float64}(chains.traces[chain].initial_state.position)
end

function r_chain_sparse_initial_velocity(chains::PDMPChains; chain::Int)
    Vector{Float64}(chains.traces[chain].initial_state.velocity)
end

function r_chain_sparse_initial_free_mask(chains::PDMPChains; chain::Int)
    point = chains.traces[chain].initial_state
    return .!(iszero.(point.position) .& iszero.(point.velocity))
end

function r_chain_sparse_event_indices(chains::PDMPChains; chain::Int)
    Int32[e.index for e in chains.traces[chain].events]
end

function r_chain_sparse_event_times(chains::PDMPChains; chain::Int)
    Float64[e.time for e in chains.traces[chain].events]
end

function r_chain_sparse_event_positions(chains::PDMPChains; chain::Int)
    Float64[e.position for e in chains.traces[chain].events]
end

function r_chain_sparse_event_velocities(chains::PDMPChains; chain::Int)
    Float64[e.velocity for e in chains.traces[chain].events]
end


function r_chain_sparse_event_free_values(chains::PDMPChains; chain::Int)
    trace = chains.traces[chain]
    return [!(iszero(event.position) && iszero(event.velocity)) for event in trace.events]
end

function r_from_sparse_skeleton(
        initial_times_list::AbstractVector, initial_positions_list::AbstractVector,
        initial_velocities_list::AbstractVector, event_indices_list::AbstractVector,
        event_times_list::AbstractVector, event_positions_list::AbstractVector,
        event_velocities_list::AbstractVector, is_boomerang_list::AbstractVector,
        mu_list::AbstractVector)
    n = length(initial_times_list)
    traces = PDMPSamplers.FactorizedTrace[]
    sizehint!(traces, n)
    for i in 1:n
        initial_time     = Float64(initial_times_list[i])
        initial_position = Vector{Float64}(initial_positions_list[i])
        initial_velocity = Vector{Float64}(initial_velocities_list[i])
        d = length(initial_position)
        is_boom = Bool(is_boomerang_list[i])
        mu = Vector{Float64}(mu_list[i])
        flow = is_boom ? Boomerang(I(d), mu) : ZigZag(I(d), zeros(d))
        initial_state = PDMPEvent(initial_time, initial_position, initial_velocity)
        events = PDMPSamplers.FactorizedEvent.(
            Int.(event_indices_list[i]),
            Float64.(event_times_list[i]),
            Float64.(event_positions_list[i]),
            Float64.(event_velocities_list[i]),
        )
        push!(traces, _r_factorized_trace(events, flow, initial_state))
    end
    PDMPChains(traces, [PDMPSamplers.StatisticCounter() for _ in traces],
        fill(nothing, n), fill(nothing, n), fill(nothing, n))
end

function r_from_streaming_skeleton(directory_list::AbstractVector,
        chunk_count_list::AbstractVector,
        physical_events_list::AbstractVector,
        computational_boundaries_list::AbstractVector,
        buffer_maximum_list::AbstractVector,
        initial_times_list::AbstractVector,
        initial_positions_list::AbstractVector,
        initial_velocities_list::AbstractVector,
        initial_stored_velocities_list::AbstractVector,
        initial_free_masks_list::AbstractVector,
        terminal_times_list::AbstractVector,
        terminal_positions_list::AbstractVector,
        terminal_velocities_list::AbstractVector,
        terminal_stored_velocities_list::AbstractVector,
        terminal_free_masks_list::AbstractVector,
        is_boomerang_list::AbstractVector, mu_list::AbstractVector,
        flow_cov_list::AbstractVector)
    n = length(initial_times_list)
    traces = PDMPSamplers.StreamingPDMPTrace[]
    initial_states = Any[]
    terminal_states = Any[]
    sizehint!(traces, n); sizehint!(initial_states, n)
    sizehint!(terminal_states, n)
    for i in 1:n
        position0 = Vector{Float64}(initial_positions_list[i])
        d = length(position0)
        raw_covariance = flow_cov_list[i]
        covariance = raw_covariance isa AbstractVector ?
            Matrix(Diagonal(Vector{Float64}(raw_covariance))) :
            Matrix{Float64}(raw_covariance)
        precision = inv(Symmetric(covariance))
        mu = Vector{Float64}(mu_list[i])
        flow = Bool(is_boomerang_list[i]) ? Boomerang(precision, mu) :
            ZigZag(I(d), zeros(d))
        initial = PDMPSamplers.PDMPTerminalState(
            Float64(initial_times_list[i]), position0,
            Vector{Float64}(initial_velocities_list[i]),
            Vector{Float64}(initial_stored_velocities_list[i]),
            BitVector(initial_free_masks_list[i]), flow)
        terminal = PDMPSamplers.PDMPTerminalState(
            Float64(terminal_times_list[i]),
            Vector{Float64}(terminal_positions_list[i]),
            Vector{Float64}(terminal_velocities_list[i]),
            Vector{Float64}(terminal_stored_velocities_list[i]),
            BitVector(terminal_free_masks_list[i]), flow)
        directory = String(directory_list[i])
        chunk_count = Int(chunk_count_list[i])
        trace = PDMPSamplers.restore_streaming_trace(directory, chunk_count,
            flow, initial, terminal;
            physical_events=Int(physical_events_list[i]),
            computational_boundaries=Int(computational_boundaries_list[i]),
            buffer_maximum=Int(buffer_maximum_list[i]))
        push!(traces, trace); push!(initial_states, initial)
        push!(terminal_states, terminal)
    end
    return PDMPChains(traces, [PDMPSamplers.StatisticCounter() for _ in traces], initial_states,
        copy(initial_states), terminal_states)
end

@inline function _r_stream_integrate_piece!(first, second, occupancy,
        model_first, model_second, bin::Int, coordinates, edge_local,
        x, velocity, free, elapsed::Float64, flow)
    elapsed <= 0 && return nothing
    base = PDMPSamplers.underlying_flow(flow)
    is_boomerang = base isa PDMPSamplers.AnyBoomerang
    model_size = 0
    @inbounds for ii in eachindex(coordinates)
        j = coordinates[ii]
        x0 = x[j]
        if edge_local[ii] && free[j]
            occupancy[ii, bin] += elapsed
            model_size += 1
        end
        if !free[j]
            first[ii, bin] += x0 * elapsed
            second[ii, bin] += x0 * x0 * elapsed
        elseif is_boomerang
            s, c = sincos(elapsed)
            a = x0 - base.μ[j]
            v = velocity[j]
            first[ii, bin] += base.μ[j] * elapsed + a * s + v * (1 - c)
            second[ii, bin] +=
                a * a * (elapsed / 2 + sin(2elapsed) / 4) +
                v * v * (elapsed / 2 - sin(2elapsed) / 4) +
                base.μ[j]^2 * elapsed + a * v * s^2 +
                2a * base.μ[j] * s + 2v * base.μ[j] * (1 - c)
        else
            v = velocity[j]
            first[ii, bin] += x0 * elapsed + 0.5v * elapsed^2
            second[ii, bin] += x0^2 * elapsed + x0 * v * elapsed^2 +
                v^2 * elapsed^3 / 3
        end
    end
    model_first[bin] += model_size * elapsed
    model_second[bin] += model_size^2 * elapsed
    return nothing
end

function _r_stream_state_at!(out_x, out_velocity, out_free, column, x, velocity, free,
        elapsed, flow, coordinates)
    scratch_x = copy(x)
    scratch_velocity = copy(velocity)
    PDMPSamplers.stream_move!(scratch_x, scratch_velocity, free, elapsed, flow)
    @inbounds for ii in eachindex(coordinates)
        j = coordinates[ii]
        out_x[ii, column] = scratch_x[j]
        out_velocity[ii, column] = free[j] ? scratch_velocity[j] : 0.0
        out_free[ii, column] = free[j]
    end
    return nothing
end

"""One bounded chunk scan for the exact OMRF analysis contract.

`edges` and `grid_times` are relative to the retained segment origin.  No
dense event skeleton or concatenated event-time vector is constructed.
"""
function r_streaming_scientific_reduction(chains::PDMPChains,
        coordinates_raw::AbstractVector, edge_mask_raw::AbstractVector,
        edges_raw::AbstractVector, grid_times_raw;
        chain::Int=1)
    trace = chains.traces[chain]
    trace isa PDMPSamplers.StreamingPDMPTrace || throw(ArgumentError(
        "scientific streaming reduction requires a streaming trace"))
    trace.finalized || PDMPSamplers.flush_streaming_trace!(trace)
    coordinates = Int.(coordinates_raw)
    edge_mask = Bool.(edge_mask_raw)
    edge_local = edge_mask[coordinates]
    edges = Float64.(edges_raw)
    grid_times = grid_times_raw isa AbstractVector ?
        Float64.(grid_times_raw) : [Float64(grid_times_raw)]
    issorted(edges) && length(edges) >= 2 && all(diff(edges) .> 0) ||
        throw(ArgumentError("analysis edges must be strictly increasing"))
    issorted(grid_times) || throw(ArgumentError("analysis grid must be sorted"))
    initial = trace.initial_state
    terminal = trace.terminal_state
    trace_duration = terminal.t - initial.t
    tol = max(128eps(max(abs(trace_duration), 1.0)),
        1e-9 * max(abs(trace_duration), 1.0))
    abs(edges[1]) <= tol && abs(edges[end] - trace_duration) <= tol ||
        throw(ArgumentError("analysis edges must span the complete segment"))
    # The task manifest owns the intended physical duration. A sampler can
    # reach the same horizon with a few ulps of accumulated clock error.
    duration = edges[end]
    all(t -> -tol <= t <= duration + tol, grid_times) ||
        throw(ArgumentError("analysis grid lies outside the segment"))
    q = length(coordinates); nb = length(edges) - 1
    first = zeros(q, nb); second = zeros(q, nb); occupancy = zeros(q, nb)
    transitions = zeros(Int, q, nb); activations = zeros(Int, q, nb)
    deactivations = zeros(Int, q, nb)
    model_first = zeros(nb); model_second = zeros(nb)
    grid_value = Matrix{Float64}(undef, q, length(grid_times))
    grid_velocity = Matrix{Float64}(undef, q, length(grid_times))
    grid_free = Matrix{Bool}(undef, q, length(grid_times))
    cursor = PDMPSamplers.stream_cursor(trace)
    absolute_time = initial.t
    local_time = 0.0
    bin = 1
    grid_index = 1

    # Integrate one physical segment while splitting only at requested
    # analysis boundaries. Grid points strictly before its outgoing event use
    # the incoming physical state; a point on an event uses the outgoing state.
    function advance_to!(target_local)
        while local_time < target_local - tol
            stop = min(target_local, edges[bin + 1])
            elapsed = stop - local_time
            _r_stream_integrate_piece!(first, second, occupancy, model_first,
                model_second, bin, coordinates, edge_local, cursor.x,
                cursor.velocity, cursor.free, elapsed, trace.flow)
            PDMPSamplers.stream_move!(cursor.x, cursor.velocity, cursor.free,
                elapsed, trace.flow)
            local_time = stop
            if bin < nb && local_time >= edges[bin + 1] - tol
                bin += 1
            end
        end
        absolute_time = initial.t + local_time
        return nothing
    end

    for metadata in trace.chunks
        chunk = PDMPSamplers.load_streaming_chunk(metadata)
        @inbounds for event in eachindex(chunk.times)
            event_local = chunk.times[event] - initial.t
            if event_local > duration && event_local <= duration + tol
                event_local = duration
            elseif event_local > duration + tol
                throw(ArgumentError("streaming event lies beyond analysis horizon"))
            end
            while grid_index <= length(grid_times) &&
                    grid_times[grid_index] < event_local - tol
                target = grid_times[grid_index]
                _r_stream_state_at!(grid_value, grid_velocity, grid_free, grid_index,
                    cursor.x, cursor.velocity, cursor.free,
                    target - local_time, trace.flow, coordinates)
                grid_index += 1
            end
            advance_to!(event_local)
            kind = chunk.kinds[event]
            changed_coordinate = kind in (PDMPSamplers.STREAM_EVENT_FREEZE,
                PDMPSamplers.STREAM_EVENT_RELEASE) ?
                Int(chunk.coordinates[event]) : 0
            before_free = changed_coordinate == 0 ? false :
                cursor.free[changed_coordinate]
            PDMPSamplers.stream_apply_event!(cursor, chunk, event)
            transition_bin = clamp(searchsortedfirst(edges, event_local) - 1,
                1, nb)
            if changed_coordinate != 0 && before_free != cursor.free[changed_coordinate]
                for ii in eachindex(coordinates)
                    j = coordinates[ii]
                    if edge_local[ii] && j == changed_coordinate
                    transitions[ii, transition_bin] += 1
                    cursor.free[j] ? (activations[ii, transition_bin] += 1) :
                        (deactivations[ii, transition_bin] += 1)
                    end
                end
            end
            while grid_index <= length(grid_times) &&
                    abs(grid_times[grid_index] - event_local) <= tol
                for ii in eachindex(coordinates)
                    j = coordinates[ii]
                    grid_value[ii, grid_index] = cursor.x[j]
                    grid_velocity[ii, grid_index] = cursor.free[j] ?
                        cursor.velocity[j] : 0.0
                    grid_free[ii, grid_index] = cursor.free[j]
                end
                grid_index += 1
            end
        end
    end
    while grid_index <= length(grid_times)
        target = grid_times[grid_index]
        _r_stream_state_at!(grid_value, grid_velocity, grid_free, grid_index, cursor.x,
            cursor.velocity, cursor.free, target - local_time, trace.flow,
            coordinates)
        grid_index += 1
    end
    advance_to!(duration)
    replay_position = copy(cursor.x)
    replay_velocity = copy(cursor.velocity)
    replay_stored = copy(cursor.stored_velocity)
    replay_free = copy(cursor.free)
    return (edges=edges, first=first, second=second, occupancy=occupancy,
        transitions=transitions, activations=activations,
        deactivations=deactivations, model_first=model_first,
        model_second=model_second, grid_times=grid_times,
        grid_value=grid_value, grid_velocity=grid_velocity,
        grid_free=grid_free,
        replay_position=replay_position, replay_velocity=replay_velocity,
        replay_stored_velocity=replay_stored, replay_free=replay_free,
        terminal_position=copy(terminal.position),
        terminal_velocity=copy(terminal.physical_velocity),
        terminal_stored_velocity=copy(terminal.stored_frozen_velocity),
        terminal_free=copy(terminal.free), chunks=length(trace.chunks),
        physical_events=trace.physical_events,
        output_bytes=Float64(trace.chunks.total_bytes))
end

function r_streaming_replay_validation(chains::PDMPChains; chain::Int=1)
    trace = chains.traces[chain]
    trace isa PDMPSamplers.StreamingPDMPTrace || throw(ArgumentError(
        "streaming replay validation requires a streaming trace"))
    d = length(trace.initial_state.position)
    z = r_streaming_scientific_reduction(chains, collect(1:d), falses(d),
        [0.0, trace.terminal_state.t - trace.initial_state.t], Float64[];
        chain=chain)
    return (position_max_abs=maximum(abs.(z.replay_position .-
            z.terminal_position)),
        velocity_max_abs=maximum(abs.(z.replay_velocity .-
            z.terminal_velocity)),
        stored_velocity_max_abs=maximum(abs.(z.replay_stored_velocity .-
            z.terminal_stored_velocity)), mask_equal=z.replay_free == z.terminal_free)
end

# Deliberately test-only: production streaming analysis uses the bounded
# reduction above. This short-fixture materialization supplies an independent
# dense contract for regression comparisons.
function r_streaming_dense_test_reference(chains::PDMPChains; chain::Int=1)
    trace = chains.traces[chain]
    trace isa PDMPSamplers.StreamingPDMPTrace || throw(ArgumentError(
        "dense streaming test reference requires a streaming trace"))
    times = Float64[]; positions = Vector{Float64}[]
    velocities = Vector{Float64}[]; masks = BitVector[]
    iteration = iterate(trace)
    while iteration !== nothing
        push!(times, iteration[2][1])
        push!(positions, copy(iteration[2][2]))
        push!(velocities, copy(iteration[2][3]))
        push!(masks, copy(iteration[2][4].free))
        iteration = iterate(trace, iteration[2])
    end
    return (times=times, positions=reduce(hcat, positions),
        velocities=reduce(hcat, velocities), free_masks=reduce(hcat, masks),
        is_boomerang=PDMPSamplers.underlying_flow(trace.flow) isa AnyBoomerang,
        is_mutable=PDMPSamplers.underlying_flow(trace.flow) isa MutableBoomerang,
        mu=_resolved_flow_mean(trace.flow),
        flow_cov=_resolved_flow_covariance(trace.flow))
end

function _r_factorized_trace(events, flow, initial_state)
    bounds_cache = Dict{Int, Tuple{Float64, Float64}}()
    # Current constructor: the last state is marked invalid, so the trace replays
    # the events to find it when it is first needed.
    if applicable(PDMPSamplers.FactorizedTrace,
                  events, flow, initial_state, initial_state, 0.0, false,
                  bounds_cache, Ref(length(events)))
        return PDMPSamplers.FactorizedTrace(
            events, flow, initial_state, deepcopy(initial_state), 0.0, false,
            bounds_cache, Ref(length(events))
        )
    elseif applicable(PDMPSamplers.FactorizedTrace,
                      events, flow, initial_state, initial_state, false, bounds_cache)
        return PDMPSamplers.FactorizedTrace(
            events, flow, initial_state, initial_state, false, bounds_cache
        )
    else
        return PDMPSamplers.FactorizedTrace(events, flow, initial_state)
    end
end

function r_chain_times(chains::PDMPChains; chain::Int)
    Vector{Float64}(_dense_compact(chains, chain).times)
end

function r_chain_positions(chains::PDMPChains; chain::Int)
    Matrix{Float64}(_dense_compact(chains, chain).positions)
end

function r_chain_velocities(chains::PDMPChains; chain::Int)
    Matrix{Float64}(_dense_compact(chains, chain).velocities)
end

function r_chain_free_masks(chains::PDMPChains; chain::Int)
    masks = _dense_compact(chains, chain).free_masks
    masks === nothing && return falses(0, 0)
    return Matrix{Bool}(masks)
end

function r_chain_is_boomerang(chains::PDMPChains; chain::Int)
    PDMPSamplers.underlying_flow(chains.traces[chain].flow) isa AnyBoomerang
end

function r_chain_is_mutable_boomerang(chains::PDMPChains; chain::Int)
    PDMPSamplers.underlying_flow(chains.traces[chain].flow) isa MutableBoomerang
end

function r_chain_mu(chains::PDMPChains; chain::Int)
    Vector{Float64}(_resolved_flow_mean(chains.traces[chain].flow))
end

function r_chain_flow_cov(chains::PDMPChains; chain::Int)
    Matrix{Float64}(_resolved_flow_covariance(chains.traces[chain].flow))
end

_terminal(chains::PDMPChains, chain::Int) =
    PDMPSamplers.terminal_state(chains, chain)
r_chain_terminal_time(chains::PDMPChains; chain::Int) =
    Float64(_terminal(chains, chain).t)
r_chain_terminal_position(chains::PDMPChains; chain::Int) =
    Vector{Float64}(_terminal(chains, chain).position)
r_chain_terminal_physical_velocity(chains::PDMPChains; chain::Int) =
    Vector{Float64}(_terminal(chains, chain).physical_velocity)
r_chain_terminal_stored_velocity(chains::PDMPChains; chain::Int) =
    Vector{Float64}(_terminal(chains, chain).stored_frozen_velocity)
r_chain_terminal_free_mask(chains::PDMPChains; chain::Int) =
    BitVector(_terminal(chains, chain).free)
r_chain_terminal_flow_spec(chains::PDMPChains; chain::Int) =
    _flow_spec(_terminal(chains, chain).flow)

_installed_initial(chains::PDMPChains, chain::Int) =
    PDMPSamplers.initial_state(chains, chain)
r_chain_initial_time(chains::PDMPChains; chain::Int) =
    Float64(_installed_initial(chains, chain).t)
r_chain_initial_position(chains::PDMPChains; chain::Int) =
    Vector{Float64}(_installed_initial(chains, chain).position)
r_chain_initial_physical_velocity(chains::PDMPChains; chain::Int) =
    Vector{Float64}(_installed_initial(chains, chain).physical_velocity)
r_chain_initial_stored_velocity(chains::PDMPChains; chain::Int) =
    Vector{Float64}(_installed_initial(chains, chain).stored_frozen_velocity)
r_chain_initial_free_mask(chains::PDMPChains; chain::Int) =
    BitVector(_installed_initial(chains, chain).free)
r_chain_initial_flow_spec(chains::PDMPChains; chain::Int) =
    _flow_spec(_installed_initial(chains, chain).flow)

_retained_initial(chains::PDMPChains, chain::Int) =
    PDMPSamplers.retained_initial_state(chains, chain)
r_chain_retained_initial_time(chains::PDMPChains; chain::Int) =
    Float64(_retained_initial(chains, chain).t)
r_chain_retained_initial_position(chains::PDMPChains; chain::Int) =
    Vector{Float64}(_retained_initial(chains, chain).position)
r_chain_retained_initial_physical_velocity(chains::PDMPChains; chain::Int) =
    Vector{Float64}(_retained_initial(chains, chain).physical_velocity)
r_chain_retained_initial_stored_velocity(chains::PDMPChains; chain::Int) =
    Vector{Float64}(_retained_initial(chains, chain).stored_frozen_velocity)
r_chain_retained_initial_free_mask(chains::PDMPChains; chain::Int) =
    BitVector(_retained_initial(chains, chain).free)
r_chain_retained_initial_flow_spec(chains::PDMPChains; chain::Int) =
    _flow_spec(_retained_initial(chains, chain).flow)

r_chain_flow_spec(chains::PDMPChains; chain::Int) =
    _flow_spec(chains.traces[chain].flow)

function _final_lambda_ref(chains::PDMPChains, chain::Int)
    base = PDMPSamplers.underlying_flow(chains.traces[chain].flow)
    return base isa AnyBoomerang ? Float64(base.λref) : NaN
end

function r_from_skeleton(times_list::AbstractVector, positions_list::AbstractVector,
                          velocities_list::AbstractVector, is_boomerang_list::AbstractVector,
                          mu_list::AbstractVector, is_mutable_boomerang_list::AbstractVector,
                          flow_cov_list=nothing, free_masks_list=nothing)
    n = length(times_list)
    traces = PDMPTrace[]
    sizehint!(traces, n)
    for i in 1:n
        times      = Vector{Float64}(times_list[i])
        positions  = Matrix{Float64}(positions_list[i])
        velocities = Matrix{Float64}(velocities_list[i])
        d = size(positions, 1)
        is_boom    = Bool(is_boomerang_list[i])
        is_mutable = Bool(is_mutable_boomerang_list[i])
        mu         = Vector{Float64}(mu_list[i])
        covariance = if flow_cov_list === nothing || isempty(flow_cov_list[i])
            Matrix{Float64}(I, d, d)
        else
            Matrix{Float64}(flow_cov_list[i])
        end
        precision = inv(Symmetric(covariance))
        flow = if is_mutable
            MutableBoomerang(precision, mu)
        elseif is_boom
            Boomerang(precision, mu)
        else
            ZigZag(I(d), zeros(d))
        end
        masks = if free_masks_list === nothing || isempty(free_masks_list[i])
            nothing
        else
            Matrix{Bool}(free_masks_list[i])
        end
        masks === nothing || size(masks) == size(positions) ||
            throw(DimensionMismatch("free-mask dimensions must match positions"))
        push!(traces, PDMPTrace(times, positions, velocities, flow, masks))
    end
    PDMPChains(traces, [PDMPSamplers.StatisticCounter() for _ in traces],
        fill(nothing, n), fill(nothing, n), fill(nothing, n))
end

const StatsValue = Union{Vector{Float64}, Matrix{Float64}, Vector{String}}

function _ct_ess_matrix(chains::PDMPChains)
    n_chains = length(chains.traces)
    n_chains == 0 && return zeros(0, 0)

    d = length(first(chains.traces[1]).position)
    ct_ess = fill(NaN, n_chains, d)

    for i in 1:n_chains
        try
            ct_ess[i, :] .= ess(chains; chain=i)
        catch
        end
    end

    return ct_ess
end

function _counter_float(counter, name::Symbol)
    try
        return Float64(getproperty(counter, name))
    catch
        return 0.0
    end
end

function _counter_string(counter, name::Symbol)
    try
        return string(getproperty(counter, name))
    catch
        return "unavailable"
    end
end

function extract_stats(chains::PDMPChains)
    all = chains.stats
    ct_ess = try
        _ct_ess_matrix(chains)
    catch
        zeros(0, 0)
    end
    vals(name::Symbol) = Float64[_counter_float(s, name) for s in all]
    textvals(name::Symbol) = String[_counter_string(s, name) for s in all]
    ct_ess_batches = Float64[max(50, isqrt(length(trace))) for trace in chains.traces]
    return Dict{String, StatsValue}(
        "final_lambda_ref"      => Float64[_final_lambda_ref(chains, i) for i in eachindex(chains.traces)],
        "reflections_events"    => vals(:reflections_events),
        "reflections_accepted"  => vals(:reflections_accepted),
        "refreshment_events"    => vals(:refreshment_events),
        "sticky_events"         => vals(:sticky_events),
        "sticky_freezes"        => vals(:sticky_freezes),
        "sticky_unfreezes"      => vals(:sticky_unfreezes),
        "sticky_unfreeze_rejections" => vals(:sticky_unfreeze_rejections),
        "support_boundary_events" => vals(:support_boundary_events),
        "support_boundary_refresh_attempts" => vals(:support_boundary_refresh_attempts),
        "support_boundary_refresh_failures" => vals(:support_boundary_refresh_failures),
        "gradient_calls"        => vals(:∇f_calls),
        "hessian_calls"         => vals(:∇²f_calls),
        "full_gradient_calls" => vals(:full_gradient_calls),
        "prior_gradient_calls" => vals(:prior_gradient_calls),
        "residual_oracle_calls" => vals(:residual_oracle_calls),
        "fd_curvature_gradient_calls" => vals(:fd_curvature_gradient_calls),
        "elapsed_time"          => vals(:elapsed_time),
        "grid_builds"           => vals(:grid_builds),
        "grid_shrinks"          => vals(:grid_shrinks),
        "grid_grows"            => vals(:grid_grows),
        "grid_early_stops"      => vals(:grid_early_stops),
        "grid_points_evaluated" => vals(:grid_points_evaluated),
        "grid_points_skipped"   => vals(:grid_points_skipped),
        "grid_N_current"        => vals(:grid_N_current),
        "grid_schedule_samples" => vals(:grid_schedule_samples),
        "grid_N_sum"            => vals(:grid_N_sum),
        "grid_tmax_sum"         => vals(:grid_tmax_sum),
        "grid_h_sum"            => vals(:grid_h_sum),
        "grid_initial_N"        => vals(:grid_initial_N),
        "grid_final_N"          => vals(:grid_final_N),
        "grid_initial_tmax"     => vals(:grid_initial_tmax),
        "grid_final_tmax"       => vals(:grid_final_tmax),
        "grid_initial_h"        => vals(:grid_initial_h),
        "grid_final_h"          => vals(:grid_final_h),
        "grid_warmup_objective_events" => vals(:grid_warmup_objective_events),
        "grid_warmup_objective_endpoint_gradients" => vals(:grid_warmup_objective_endpoint_gradients),
        "grid_warmup_objective_acceptance_gradients" => vals(:grid_warmup_objective_acceptance_gradients),
        "grid_warmup_objective_gradients_per_event" => vals(:grid_warmup_objective_gradients_per_event),
        "grid_warmup_objective_horizon_hits" => vals(:grid_warmup_objective_horizon_hits),
        "grid_warmup_objective_horizon_rate" => vals(:grid_warmup_objective_horizon_rate),
        "grid_warmup_objective_rejections" => vals(:grid_warmup_objective_rejections),
        "grid_warmup_objective_rejection_rate" => vals(:grid_warmup_objective_rejection_rate),
        "grid_schedule_frozen" => vals(:grid_schedule_frozen),
        "curvature_backend" => textvals(:curvature_backend),
        "lazy_fallback_low_tightness" => vals(:lazy_fallback_low_tightness),
        "lazy_fallback_bound_violation" => vals(:lazy_fallback_bound_violation),
        "lazy_proposal_attempts" => vals(:lazy_proposal_attempts),
        "lazy_proposal_rejections" => vals(:lazy_proposal_rejections),
        "grid_resets_from_dynamics_adaptation" => vals(:grid_resets_from_dynamics_adaptation),
        "grid_endpoint_evaluations" => vals(:grid_endpoint_evaluations),
        "grid_endpoint_gradient_calls" => vals(:grid_endpoint_gradient_calls),
        "grid_cached_endpoint_reuses" => vals(:grid_cached_endpoint_reuses),
        "grid_acceptance_tests" => vals(:grid_acceptance_tests),
        "grid_acceptance_gradient_calls" => vals(:grid_acceptance_gradient_calls),
        "subsampling_cell_roof_proposals" => vals(:subsampling_cell_roof_proposals),
        "subsampling_aggregate_accepts" => vals(:subsampling_aggregate_accepts),
        "subsampling_subset_evaluations" => vals(:subsampling_subset_evaluations),
        "subsampling_final_reflections" => vals(:subsampling_final_reflections),
        # Exact research pipeline counters. The cumulative values are retained
        # for compatibility; the phase-local values are recorded by the Julia
        # sampler at the warmup/main boundary and are the authoritative
        # retained-phase values.
        "clock_candidates" => vals(:clock_candidates),
        "grid_bound_evaluations" => vals(:grid_bound_evaluations),
        "pointwise_screen_evaluations" => vals(:pointwise_screen_evaluations),
        "pointwise_screen_passes" => vals(:pointwise_screen_passes),
        "subset_draws" => vals(:subset_draws),
        "selected_subset_bound_evaluations" => vals(:selected_subset_bound_evaluations),
        "selected_subset_bound_passes" => vals(:selected_subset_bound_passes),
        "residual_oracle_evaluations" => vals(:residual_oracle_evaluations),
        "final_thinning_acceptances" => vals(:final_thinning_acceptances),
        "reflection_events" => vals(:reflection_events),
        "sticky_events" => vals(:sticky_events),
        "warmup_clock_candidates" => vals(:warmup_clock_candidates),
        "main_clock_candidates" => vals(:main_clock_candidates),
        "warmup_grid_bound_evaluations" => vals(:warmup_grid_bound_evaluations),
        "main_grid_bound_evaluations" => vals(:main_grid_bound_evaluations),
        "warmup_pointwise_screen_evaluations" => vals(:warmup_pointwise_screen_evaluations),
        "main_pointwise_screen_evaluations" => vals(:main_pointwise_screen_evaluations),
        "warmup_pointwise_screen_passes" => vals(:warmup_pointwise_screen_passes),
        "main_pointwise_screen_passes" => vals(:main_pointwise_screen_passes),
        "warmup_subset_draws" => vals(:warmup_subset_draws),
        "main_subset_draws" => vals(:main_subset_draws),
        "warmup_selected_subset_bound_evaluations" => vals(:warmup_selected_subset_bound_evaluations),
        "main_selected_subset_bound_evaluations" => vals(:main_selected_subset_bound_evaluations),
        "warmup_selected_subset_bound_passes" => vals(:warmup_selected_subset_bound_passes),
        "main_selected_subset_bound_passes" => vals(:main_selected_subset_bound_passes),
        "warmup_residual_oracle_evaluations" => vals(:warmup_residual_oracle_evaluations),
        "main_residual_oracle_evaluations" => vals(:main_residual_oracle_evaluations),
        "warmup_final_thinning_acceptances" => vals(:warmup_final_thinning_acceptances),
        "main_final_thinning_acceptances" => vals(:main_final_thinning_acceptances),
        "warmup_reflection_events" => vals(:warmup_reflection_events),
        "main_reflection_events" => vals(:main_reflection_events),
        "warmup_sticky_events" => vals(:warmup_sticky_events),
        "main_sticky_events" => vals(:main_sticky_events),
        "warmup_grid_schedule_builds" => vals(:warmup_grid_schedule_builds),
        "main_grid_schedule_builds" => vals(:main_grid_schedule_builds),
        "warmup_grid_schedule_searches" => vals(:warmup_grid_schedule_searches),
        "main_grid_schedule_searches" => vals(:main_grid_schedule_searches),
        "warmup_poisson_time_generations" => vals(:warmup_poisson_time_generations),
        "main_poisson_time_generations" => vals(:main_poisson_time_generations),
        "warmup_candidate_state_copies" => vals(:warmup_candidate_state_copies),
        "main_candidate_state_copies" => vals(:main_candidate_state_copies),
        "warmup_candidate_move_forwards" => vals(:warmup_candidate_move_forwards),
        "main_candidate_move_forwards" => vals(:main_candidate_move_forwards),
        "warmup_candidate_loop_iterations" => vals(:warmup_candidate_loop_iterations),
        "main_candidate_loop_iterations" => vals(:main_candidate_loop_iterations),
        "warmup_statistic_recordings" => vals(:warmup_statistic_recordings),
        "main_statistic_recordings" => vals(:main_statistic_recordings),
        "warmup_grid_bound_seconds" => vals(:warmup_grid_bound_seconds),
        "main_grid_bound_seconds" => vals(:main_grid_bound_seconds),
        "warmup_pointwise_screen_seconds" => vals(:warmup_pointwise_screen_seconds),
        "main_pointwise_screen_seconds" => vals(:main_pointwise_screen_seconds),
        "warmup_subset_draw_seconds" => vals(:warmup_subset_draw_seconds),
        "main_subset_draw_seconds" => vals(:main_subset_draw_seconds),
        "warmup_selected_subset_refinement_seconds" => vals(:warmup_selected_subset_refinement_seconds),
        "main_selected_subset_refinement_seconds" => vals(:main_selected_subset_refinement_seconds),
        "warmup_residual_oracle_seconds" => vals(:warmup_residual_oracle_seconds),
        "main_residual_oracle_seconds" => vals(:main_residual_oracle_seconds),
        "warmup_final_thinning_seconds" => vals(:warmup_final_thinning_seconds),
        "main_final_thinning_seconds" => vals(:main_final_thinning_seconds),
        "warmup_reflection_update_seconds" => vals(:warmup_reflection_update_seconds),
        "main_reflection_update_seconds" => vals(:main_reflection_update_seconds),
        "warmup_reflection_gradient_seconds" => vals(:warmup_reflection_gradient_seconds),
        "main_reflection_gradient_seconds" => vals(:main_reflection_gradient_seconds),
        "warmup_reflection_covariance_seconds" => vals(:warmup_reflection_covariance_seconds),
        "main_reflection_covariance_seconds" => vals(:main_reflection_covariance_seconds),
        "warmup_reflection_sticky_schedule_seconds" => vals(:warmup_reflection_sticky_schedule_seconds),
        "main_reflection_sticky_schedule_seconds" => vals(:main_reflection_sticky_schedule_seconds),
        "warmup_reflection_state_validation_seconds" => vals(:warmup_reflection_state_validation_seconds),
        "main_reflection_state_validation_seconds" => vals(:main_reflection_state_validation_seconds),
        "warmup_sticky_update_seconds" => vals(:warmup_sticky_update_seconds),
        "main_sticky_update_seconds" => vals(:main_sticky_update_seconds),
        "warmup_grid_schedule_build_seconds" => vals(:warmup_grid_schedule_build_seconds),
        "main_grid_schedule_build_seconds" => vals(:main_grid_schedule_build_seconds),
        "warmup_grid_schedule_search_seconds" => vals(:warmup_grid_schedule_search_seconds),
        "main_grid_schedule_search_seconds" => vals(:main_grid_schedule_search_seconds),
        "warmup_poisson_time_generation_seconds" => vals(:warmup_poisson_time_generation_seconds),
        "main_poisson_time_generation_seconds" => vals(:main_poisson_time_generation_seconds),
        "warmup_candidate_state_copy_seconds" => vals(:warmup_candidate_state_copy_seconds),
        "main_candidate_state_copy_seconds" => vals(:main_candidate_state_copy_seconds),
        "warmup_candidate_move_forward_seconds" => vals(:warmup_candidate_move_forward_seconds),
        "main_candidate_move_forward_seconds" => vals(:main_candidate_move_forward_seconds),
        "warmup_candidate_loop_overhead_seconds" => vals(:warmup_candidate_loop_overhead_seconds),
        "main_candidate_loop_overhead_seconds" => vals(:main_candidate_loop_overhead_seconds),
        "warmup_statistic_recording_seconds" => vals(:warmup_statistic_recording_seconds),
        "main_statistic_recording_seconds" => vals(:main_statistic_recording_seconds),
        "warmup_phase_elapsed_seconds" => vals(:warmup_phase_elapsed_time),
        "main_phase_elapsed_seconds" => vals(:main_phase_elapsed_time),
        "main_phase_allocated_bytes" => vals(:main_phase_allocated_bytes),
        "initialization_allocated_bytes" => vals(:initialization_allocated_bytes),
        "main_sampler_initialization_allocated_bytes" =>
            vals(:main_sampler_initialization_allocated_bytes),
        "other_main_loop_seconds" => vals(:other_main_loop_seconds),
        "grid_horizon_hits" => vals(:grid_horizon_hits),
        "grid_budget_tail_restarts" => vals(:grid_budget_tail_restarts),
        "grid_bound_violations" => vals(:grid_bound_violations),
        "sticky_inner_searches" => vals(:sticky_inner_searches),
        "sticky_inner_wins" => vals(:sticky_inner_wins),
        "sticky_inner_wasted_by_sticky" => vals(:sticky_inner_wasted_by_sticky),
        "sticky_inner_wasted_by_refresh" => vals(:sticky_inner_wasted_by_refresh),
        "sticky_all_frozen_events" => vals(:sticky_all_frozen_events),
        "warmup_events" => vals(:warmup_events),
        "main_events" => vals(:main_events),
        "final_metric_scale_min" => vals(:final_metric_scale_min),
        "final_metric_scale_max" => vals(:final_metric_scale_max),
        "adaptation_updates" => vals(:adaptation_updates),
        "preconditioner_updates" => vals(:preconditioner_updates),
        "boomerang_updates" => vals(:boomerang_updates),
        "adaptation_free_time_min" => vals(:adaptation_free_time_min),
        "adaptation_free_time_max" => vals(:adaptation_free_time_max),
        "warmup_gradient_calls" => vals(:warmup_gradient_calls),
        "main_gradient_calls" => vals(:main_gradient_calls),
        "warmup_hessian_calls" => vals(:warmup_hessian_calls),
        "main_hessian_calls" => vals(:main_hessian_calls),
        "warmup_full_gradient_calls" => vals(:warmup_full_gradient_calls),
        "main_full_gradient_calls" => vals(:main_full_gradient_calls),
        "warmup_potential_calls" => vals(:warmup_potential_calls),
        "main_potential_calls" => vals(:main_potential_calls),
        "warmup_prior_gradient_calls" => vals(:warmup_prior_gradient_calls),
        "main_prior_gradient_calls" => vals(:main_prior_gradient_calls),
        "warmup_fd_curvature_gradient_calls" => vals(:warmup_fd_curvature_gradient_calls),
        "main_fd_curvature_gradient_calls" => vals(:main_fd_curvature_gradient_calls),
        "warmup_exact_curvature_calls" => vals(:warmup_exact_curvature_calls),
        "main_exact_curvature_calls" => vals(:main_exact_curvature_calls),
        "warmup_grid_endpoint_evaluations" => vals(:warmup_grid_endpoint_evaluations),
        "main_grid_endpoint_evaluations" => vals(:main_grid_endpoint_evaluations),
        "warmup_grid_endpoint_gradient_calls" => vals(:warmup_grid_endpoint_gradient_calls),
        "main_grid_endpoint_gradient_calls" => vals(:main_grid_endpoint_gradient_calls),
        "warmup_grid_endpoint_hessian_calls" => vals(:warmup_grid_endpoint_hessian_calls),
        "main_grid_endpoint_hessian_calls" => vals(:main_grid_endpoint_hessian_calls),
        "warmup_grid_endpoint_derivative_calls" => vals(:warmup_grid_endpoint_derivative_calls),
        "main_grid_endpoint_derivative_calls" => vals(:main_grid_endpoint_derivative_calls),
        "warmup_grid_acceptance_gradient_calls" => vals(:warmup_grid_acceptance_gradient_calls),
        "main_grid_acceptance_gradient_calls" => vals(:main_grid_acceptance_gradient_calls),
        "warmup_grid_acceptance_tests" => vals(:warmup_grid_acceptance_tests),
        "main_grid_acceptance_tests" => vals(:main_grid_acceptance_tests),
        "warmup_grid_cached_endpoint_reuses" => vals(:warmup_grid_cached_endpoint_reuses),
        "main_grid_cached_endpoint_reuses" => vals(:main_grid_cached_endpoint_reuses),
        "warmup_grid_points_evaluated" => vals(:warmup_grid_points_evaluated),
        "main_grid_points_evaluated" => vals(:main_grid_points_evaluated),
        "warmup_grid_endpoint_derivative_points_loaded" => vals(:warmup_grid_endpoint_derivative_points_loaded),
        "main_grid_endpoint_derivative_points_loaded" => vals(:main_grid_endpoint_derivative_points_loaded),
        "warmup_subsampling_cell_roof_proposals" =>
            vals(:warmup_subsampling_cell_roof_proposals),
        "main_subsampling_cell_roof_proposals" =>
            vals(:main_subsampling_cell_roof_proposals),
        "warmup_subsampling_aggregate_accepts" =>
            vals(:warmup_subsampling_aggregate_accepts),
        "main_subsampling_aggregate_accepts" =>
            vals(:main_subsampling_aggregate_accepts),
        "warmup_subsampling_subset_evaluations" =>
            vals(:warmup_subsampling_subset_evaluations),
        "main_subsampling_subset_evaluations" =>
            vals(:main_subsampling_subset_evaluations),
        "warmup_subsampling_final_reflections" =>
            vals(:warmup_subsampling_final_reflections),
        "main_subsampling_final_reflections" =>
            vals(:main_subsampling_final_reflections),
        "warmup_elapsed_time" => vals(:warmup_elapsed_time),
        "main_elapsed_time" => vals(:main_elapsed_time),
        "initialization_elapsed_time" => vals(:initialization_elapsed_time),
        "warmup_phase_elapsed_time" => vals(:warmup_phase_elapsed_time),
        "main_phase_elapsed_time" => vals(:main_phase_elapsed_time),
        "transition_elapsed_time" => vals(:transition_elapsed_time),
        "warmup_adapter_finish_elapsed_time" =>
            vals(:warmup_adapter_finish_elapsed_time),
        "main_sampler_initialization_elapsed_time" =>
            vals(:main_sampler_initialization_elapsed_time),
        "algorithm_warmup_finish_elapsed_time" =>
            vals(:algorithm_warmup_finish_elapsed_time),
        "finalization_elapsed_time" => vals(:finalization_elapsed_time),
        "ct_ess"                => ct_ess,
        "ct_ess_batches"        => ct_ess_batches,
    )
end

function r_pdmp_stan(
        path_to_stan_model::String,
        path_to_stan_data::String,
        x0,
        flow_type::String,
        algorithm_type::String,
        flow_mean,
        flow_cov;
        kwargs...
    )

    model = PDMPModel(path_to_stan_model, path_to_stan_data;
        hvp=get(kwargs, :curvature_backend, "auto") != "finite_difference")
    return r_pdmp_stan(model, x0, flow_type, algorithm_type, flow_mean, flow_cov; kwargs...)
end

function r_pdmp_stan(
        model::PDMPModel,
        x0,
        flow_type::String,
        algorithm_type::String,
        flow_mean,
        flow_cov;
        theta0 = nothing,
        free0 = nothing,
        stored_velocity0 = nothing,
        c0::Float64 = 1e-2,
        grid_n::Int = 30,
        grid_t_max::Float64 = 2.0,
        curvature_backend::String = "auto",
        grid_bound::String = "value_quadratic",
        grid_curvature_bound::Union{Nothing,Float64} = nothing,
        lazy_low_tightness_threshold::Float64 = 0.1,
        lazy_max_low_tightness_rejections::Int = 3,
        lazy_max_rejections::Int = 0,
        t0::Float64 = 0.0,
        T::Float64 = 10000.0,
        t_warmup::Float64 = 0.0,
        warmup_adaptation_interval::Union{Nothing,Float64} = nothing,
        boomerang_variance_floor_fraction::Float64 = 0.0,
        sticky::Bool = false,
        can_stick = nothing,
        model_prior = nothing,
        parameter_prior = nothing,
        slab_prior = nothing,
        unc_names = String[],
        show_progress::Bool = true,
        n_chains::Int = 1,
        threaded::Bool = false,
        seed::Union{Integer, Nothing} = nothing,
        adaptive_scheme::String = "diagonal",
        refresh_rate::Float64 = 1.0,
        support_boundary_mode::String = "error",
        support_boundary_max_bisection_steps::Int = 60,
        support_boundary_time_rtol::Float64 = 1e-8,
        support_boundary_time_atol::Float64 = 1e-10,
        support_boundary_clip_fraction::Float64 = 1 - 1e-10,
        support_boundary_max_refresh_attempts::Int = 20,
        support_boundary_refresh_probe_time::Float64 = 1e-4,
        support_boundary_min_safe_time::Float64 = 1e-12,
        streaming_trace_directory::Union{Nothing,String}=nothing,
        streaming_trace_buffer_events::Int=4096
    )

    bridge_profile_start_ns = time_ns()
    d = model.d
    x0_vec = _as_float_vector(x0)
    flow_mean_vec = _as_flow_mean(flow_mean, d)
    flow_cov_mat = _as_flow_cov(flow_cov, d)
    can_stick_vec = _as_bool_vector(can_stick)
    parameter_prior_vec = isnothing(parameter_prior) ? nothing : _as_float_vector(parameter_prior)
    unc_names_vec = _as_string_vector(unc_names)
    sampling_model = _build_sticky_target(model, model_prior, slab_prior,
        can_stick_vec, unc_names_vec)

    prec = _to_precision(flow_cov_mat, d)
    flow = build_flow(flow_type, prec, flow_mean_vec; adaptive_scheme, refresh_rate)
    isfinite(boomerang_variance_floor_fraction) &&
        boomerang_variance_floor_fraction >= 0 || throw(ArgumentError(
            "boomerang_variance_floor_fraction must be finite and nonnegative"))
    boom_options = PDMPSamplers.BoomerangAdaptationOptions(
        initial_variance_floor_fraction=boomerang_variance_floor_fraction)
    boom_floor_mask = _boomerang_logscale_floor_mask(unc_names_vec, d)
    initial_refresh_rate = try Float64(PDMPSamplers.refresh_rate(flow)) catch; NaN end
    alg0 = build_algorithm(algorithm_type; c0, d, grid_n, grid_t_max, curvature_backend,
        grid_bound, grid_curvature_bound,
        lazy_low_tightness_threshold, lazy_max_low_tightness_rejections,
        lazy_max_rejections)
    alg = _build_sticky_event_algorithm(alg0, sticky, model_prior,
        parameter_prior_vec, slab_prior, can_stick_vec, flow, unc_names_vec)

    sbopts = SupportBoundaryOptions(;
        detect_boundaries = support_boundary_mode != "error",
        mode = Symbol(support_boundary_mode),
        max_bisection_steps = support_boundary_max_bisection_steps,
        time_rtol = support_boundary_time_rtol,
        time_atol = support_boundary_time_atol,
        clip_fraction = support_boundary_clip_fraction,
        max_refresh_attempts = support_boundary_max_refresh_attempts,
        refresh_probe_time = support_boundary_refresh_probe_time,
        min_safe_time = support_boundary_min_safe_time,
    )

    initial = if isnothing(theta0)
        # Preserve the package's seeded, flow-specific velocity law. Sticky
        # coordinates initialized exactly on their face start inactive; all
        # other coordinates retain the generated invariant velocity.
        velocity = PDMPSamplers.initialize_velocity(
            PDMPSamplers.make_initial_rng(seed, n_chains), flow, d)
        if sticky
            @inbounds for i in eachindex(velocity, x0_vec, can_stick_vec)
                can_stick_vec[i] && iszero(x0_vec[i]) && (velocity[i] = 0.0)
            end
        end
        SkeletonPoint(x0_vec, velocity)
    else
        SkeletonPoint(x0_vec, _as_float_vector(theta0))
    end
    adapter = PDMPSamplers.default_warmup_adapter(
        flow, sampling_model.grad, t_warmup, t0;
        warmup_adaptation_interval,
        can_stick=sticky ? can_stick_vec : nothing,
        boomerang_options=boom_options,
        boomerang_variance_floor_mask=boom_floor_mask)
    bridge_profile_before_sample_ns = time_ns()
    chains = pdmp_sample(initial, flow, sampling_model, alg, t0, T, t_warmup;
                         progress = show_progress, n_chains = n_chains, threaded = threaded,
                         seed = seed,
                         adapter = adapter,
                         warmup_stop = build_warmup_stop(t_warmup),
                         support_boundary_options = sbopts,
                         statistic_counter = PDMPSamplers.DevelStatisticCounter,
                         initial_free=isnothing(free0) ? nothing : _as_bool_vector(free0),
                         initial_stored_velocity=isnothing(stored_velocity0) ? nothing :
                            _as_float_vector(stored_velocity0),
                         trace_storage=isnothing(streaming_trace_directory) ?
                            nothing : _bridge_trace_storage(streaming_trace_directory,
                                streaming_trace_buffer_events, T))
    bridge_profile_after_sample_ns = time_ns()
    packed = _pack_result(chains; requested_flow=flow_type,
        adaptive_scheme, initial_refresh_rate,
        warmup_adaptation_interval, t_warmup)
    if get(ENV, "OMRF_PROFILE_EXCLUSIVE", "") == "1"
        packed["bridge_profile"] = Dict(
            "julia_construction_seconds" =>
                (bridge_profile_before_sample_ns - bridge_profile_start_ns) / 1e9,
            "julia_sample_seconds" =>
                (bridge_profile_after_sample_ns - bridge_profile_before_sample_ns) / 1e9,
            "julia_pack_and_counters_seconds" =>
                (time_ns() - bridge_profile_after_sample_ns) / 1e9)
    end
    return packed
end

function r_pdmp_custom(
        grad!,
        d::Integer,
        x0,
        flow_type::String,
        algorithm_type::String,
        flow_mean,
        flow_cov;
        c0::Float64 = 1e-2,
        grid_n::Int = 30,
        grid_t_max::Float64 = 2.0,
        grid_bound::String = "value_quadratic",
        grid_curvature_bound::Union{Nothing,Float64} = nothing,
        t0::Float64 = 0.0,
        T::Float64 = 10000.0,
        t_warmup::Float64 = 0.0,
        hessian = nothing,
        sticky::Bool = false,
        can_stick = nothing,
        model_prior = nothing,
        parameter_prior = nothing,
        slab_prior = nothing,
        show_progress::Bool = true,
        n_chains::Int = 1,
        threaded::Bool = false,
        seed::Union{Integer, Nothing} = nothing,
        adaptive_scheme::String = "diagonal",
        refresh_rate::Float64 = 1.0,
        support_boundary_mode::String = "error",
        support_boundary_max_bisection_steps::Int = 60,
        support_boundary_time_rtol::Float64 = 1e-8,
        support_boundary_time_atol::Float64 = 1e-10,
        support_boundary_clip_fraction::Float64 = 1 - 1e-10,
        support_boundary_max_refresh_attempts::Int = 20,
        support_boundary_refresh_probe_time::Float64 = 1e-4,
        support_boundary_min_safe_time::Float64 = 1e-12,
        streaming_trace_directory::Union{Nothing,String}=nothing,
        streaming_trace_buffer_events::Int=4096
    )

    x0_vec = _as_float_vector(x0)
    flow_mean_vec = _as_flow_mean(flow_mean, d)
    flow_cov_mat = _as_flow_cov(flow_cov, d)
    can_stick_vec = _as_bool_vector(can_stick)
    parameter_prior_vec = isnothing(parameter_prior) ? nothing : _as_float_vector(parameter_prior)

    hvp = if !isnothing(hessian)
        (out, x, v) -> begin
            hess = hessian(x)
            mul!(out, hess, v)
        end
    else
        # Centered finite-difference HVP: H(x)*v ≈ (∇f(x + εv) - ∇f(x - εv)) / (2ε)
        let ε = 1e-5, _buf1 = zeros(d), _buf2 = zeros(d)
            (out, x, v) -> begin
                @. _buf2 = x - ε * v
                grad!(_buf1, _buf2)
                @. _buf2 = x + ε * v
                grad!(out, _buf2)
                @. out = (out - _buf1) / (2ε)
            end
        end
    end
    model = if isnothing(slab_prior)
        PDMPModel(d, FullGradient(grad!), hvp)
    else
        provider = build_slab_provider(slab_prior, String[], can_stick_vec, d)
        odds = build_model_prior_odds(model_prior, PDMPSamplers.beta_indices(provider), d)
        target = DependentSlabTarget(d, grad!, provider, odds;
            posterior_hvp=PDMPSamplers.InplaceHVP(hvp, zeros(d)))
        PDMPModel(target; hvp=true)
    end

    prec = _to_precision(flow_cov_mat, d)
    flow = build_flow(flow_type, prec, flow_mean_vec; adaptive_scheme, refresh_rate)
    initial_refresh_rate = try Float64(PDMPSamplers.refresh_rate(flow)) catch; NaN end
    alg0 = build_algorithm(algorithm_type; c0, d, grid_n, grid_t_max,
        grid_bound, grid_curvature_bound)
    alg = isnothing(slab_prior) ?
        wrap_sticky(alg0, sticky, model_prior, parameter_prior_vec, can_stick_vec) :
        wrap_dependent_sticky(alg0, sticky, model_prior, slab_prior, can_stick_vec, flow)

    sbopts = SupportBoundaryOptions(;
        detect_boundaries = support_boundary_mode != "error",
        mode = Symbol(support_boundary_mode),
        max_bisection_steps = support_boundary_max_bisection_steps,
        time_rtol = support_boundary_time_rtol,
        time_atol = support_boundary_time_atol,
        clip_fraction = support_boundary_clip_fraction,
        max_refresh_attempts = support_boundary_max_refresh_attempts,
        refresh_probe_time = support_boundary_refresh_probe_time,
        min_safe_time = support_boundary_min_safe_time,
    )

    chains = pdmp_sample(x0_vec, flow, model, alg, t0, T, t_warmup;
                         progress = show_progress, n_chains = n_chains, threaded = threaded,
                         seed = seed,
                         support_boundary_options = sbopts,
                         trace_storage=isnothing(streaming_trace_directory) ?
                            nothing : _bridge_trace_storage(streaming_trace_directory,
                                streaming_trace_buffer_events, T))
    return _pack_result(chains; requested_flow=flow_type,
        adaptive_scheme, initial_refresh_rate,
        warmup_adaptation_interval=nothing, t_warmup)
end

# ──────────────────────────────────────────────────────────────────────────────
# brms backend helpers
# ──────────────────────────────────────────────────────────────────────────────

function write_cmdstan_csv(path::String, draws::Matrix{Float64}, param_names::Vector{String};
                           chain_id::Int=1, lp_values::Union{Vector{Float64}, Nothing}=nothing)
    n_samples, n_params = size(draws)
    lp_col = isnothing(lp_values) ? zeros(n_samples) : lp_values
    buf = IOBuffer(; sizehint = n_samples * (n_params + 7) * 12)
    _write_cmdstan_header(buf, n_samples, chain_id, param_names)
    for i in 1:n_samples
        print(buf, lp_col[i])
        for _ in 1:6
            write(buf, ",0.0")
        end
        @inbounds for j in 1:n_params
            write(buf, UInt8(','))
            print(buf, draws[i, j])
        end
        write(buf, UInt8('\n'))
    end
    open(path, "w") do io
        write(io, take!(buf))
    end
    return path
end

function _write_cmdstan_header(io::IO, n_samples::Int, chain_id::Int, param_names::Vector{String})
    println(io, "# model = PDMPSamplers_model")
    println(io, "# method = sample (adapt engaged=0)")
    println(io, "#   sample")
    println(io, "#     num_samples = ", n_samples)
    println(io, "#     num_warmup = 0")
    println(io, "#     save_warmup = 0")
    println(io, "#     thin = 1")
    println(io, "# id = ", chain_id)
    println(io, "# Adaptation terminated")
    println(io, "#  Elapsed Time: 0 seconds (Warm-up)")
    println(io, "#                0 seconds (Sampling)")
    println(io, "#                0 seconds (Total)")
    diag_cols = ("lp__", "accept_stat__", "stepsize__", "treedepth__",
                 "n_leapfrog__", "divergent__", "energy__")
    join(io, diag_cols, ",")
    for name in param_names
        write(io, UInt8(','))
        write(io, name)
    end
    write(io, UInt8('\n'))
end

function _draws_for_csv(chains::PDMPChains, chain_idx::Int, discretize_dt::Float64)
    trace = chains.traces[chain_idx]
    try
        if discretize_dt > 0
            return Matrix(PDMPDiscretize(trace, discretize_dt))
        else
            return adaptive_discretize(chains; chain=chain_idx)[1]
        end
    catch
        return reshape(Vector{Float64}(first(trace).position), 1, :)
    end
end

function r_constrain_and_write_csv(sm::BridgeStan.StanModel, draws_unc::Matrix{Float64},
                                   output_csv::String; chain_id::Int=1, compute_lp::Bool=false)
    param_names_c = BridgeStan.param_names(sm; include_tp=true, include_gq=true)
    rng = BridgeStan.StanRNG(sm, chain_id)
    n = size(draws_unc, 1)
    d_unc = size(draws_unc, 2)
    n_c = length(param_names_c)
    row_buf = Vector{Float64}(undef, d_unc)
    out_buf = Vector{Float64}(undef, n_c)
    buf = IOBuffer(; sizehint = n * (n_c + 7) * 12)
    _write_cmdstan_header(buf, n, chain_id, param_names_c)
    for i in 1:n
        @inbounds for j in 1:d_unc
            row_buf[j] = draws_unc[i, j]
        end
        BridgeStan.param_constrain!(sm, row_buf, out_buf; include_tp=true, include_gq=true, rng)
        lp = compute_lp ? BridgeStan.log_density(sm, row_buf) : 0.0
        print(buf, lp)
        for _ in 1:6
            write(buf, ",0.0")
        end
        @inbounds for j in 1:n_c
            write(buf, UInt8(','))
            print(buf, out_buf[j])
        end
        write(buf, UInt8('\n'))
    end
    open(output_csv, "w") do io
        write(io, take!(buf))
    end
    return output_csv
end

include(joinpath(@__DIR__, "family_subsampling_provider.jl"))
include(joinpath(@__DIR__, "stan_subsampling_provider.jl"))

function _run_subsampling_brms_subsampled(lib_path_std::String,
        data_full_file::String, data_prior_file::String,
        N::Int, m::Int, family::Symbol,
        predictor_designs::Vector{Matrix{Float64}},
        predictor_indices::Vector{Vector{Int}}, d::Int,
        offsets::Matrix{Float64}, response::Matrix{Float64}, known_se::Vector{Float64},
        observation_multipliers::Vector{Float64},
        flow_type::String, algorithm_type::String, flow_mean, flow_cov,
        output_csv::String; c0::Float64, grid_n::Int, grid_t_max::Float64,
        grid_bound::String, grid_curvature_bound::Union{Nothing,Float64},
        t0::Float64, T::Float64, t_warmup::Float64,
        adaptive_scheme::String, discretize_dt::Float64,
        show_progress::Bool, n_chains::Int, threaded::Bool,
        seed::Union{Integer,Nothing}, compute_lp::Bool,
        curvature_backend::String, sticky::Bool, can_stick, model_prior, parameter_prior,
        n_anchor_updates::Int, use_anchor_bank::Bool, bank_capacity::Int,
        use_hcv::Bool)

    algorithm_type in ("GridThinningStrategy", "ThinningStrategy") ||
        throw(ArgumentError(
            "subsampling BridgeStan subsampling requires GridThinningStrategy or ThinningStrategy"))
    all(size(design, 1) == N for design in predictor_designs) ||
        throw(DimensionMismatch("predictor designs must have N rows"))
    length(observation_multipliers) == N || throw(DimensionMismatch(
        "observation multiplier vector must have length N"))

    fmean = _as_flow_mean(flow_mean, d)
    anchor = copy(fmean)
    contexts = SubsamplingBridgeStanContext[]
    models = PDMPModel[]
    anchor_adapters = Any[]
    anchor_managers = Union{Nothing,SubsamplingBridgeAnchorManager}[]
    anchor_capacity = use_anchor_bank ? bank_capacity :
        (n_anchor_updates > 0 ? 1 : 0)
    for _ in 1:n_chains
        ctx = _new_subsampling_bss_context(lib_path_std,
            data_full_file, data_prior_file, family, predictor_designs,
            predictor_indices, d, offsets, response, known_se,
            observation_multipliers, N, m, anchor; use_hcv)
        push!(contexts, ctx)
        built = _build_subsampling_bss_model(ctx; curvature_backend, anchor_capacity)
        if anchor_capacity > 0
            model, adapter, manager = built
            adapter.update_dt = n_anchor_updates > 0 ?
                t_warmup / n_anchor_updates : Inf
            adapter.last_update = t0
            push!(models, model)
            push!(anchor_adapters, adapter)
            push!(anchor_managers, manager)
        else
            push!(models, built)
            push!(anchor_managers, nothing)
        end
    end

    prec = inv(Symmetric(_as_flow_cov(flow_cov, d)))
    flow = build_flow(flow_type, prec, fmean; adaptive_scheme)
    initial_refresh_rate = try Float64(PDMPSamplers.refresh_rate(flow)) catch; NaN end
    alg = build_algorithm(algorithm_type; c0, d, grid_n, grid_t_max,
        curvature_backend, grid_bound, grid_curvature_bound)
    alg = wrap_sticky(alg, sticky, model_prior, parameter_prior, can_stick)
    if anchor_capacity == 0
        chains = pdmp_sample(anchor, flow, models, alg, t0, T, t_warmup;
            progress=show_progress, threaded, seed,
            statistic_counter=PDMPSamplers.DevelStatisticCounter)
    else
        run_chain = function(i)
            seed_i = isnothing(seed) ? nothing : seed + i - 1
            flow_i = deepcopy(flow)
            dynamics_adapter = PDMPSamplers.default_warmup_adapter(
                flow_i, models[i].grad, t_warmup, t0;
                can_stick=sticky ? can_stick : nothing)
            adapter_i = dynamics_adapter isa PDMPSamplers.NoAdaptation ?
                anchor_adapters[i] : PDMPSamplers.SequenceAdapter(
                    (dynamics_adapter, anchor_adapters[i]))
            chain = pdmp_sample(anchor, flow_i, [models[i]], alg,
                t0, T, t_warmup; progress=show_progress && i == 1,
                adapter=adapter_i, seed=seed_i,
                statistic_counter=PDMPSamplers.DevelStatisticCounter)
            return chain.traces[1], chain.stats[1], chain.initial_states[1],
                chain.retained_initial_states[1], chain.terminal_states[1]
        end
        chain_results = if threaded && n_chains > 1
            fetch.([Threads.@spawn run_chain(i) for i in 1:n_chains])
        else
            [run_chain(i) for i in 1:n_chains]
        end
        chains = PDMPChains(getindex.(chain_results, 1), getindex.(chain_results, 2),
            getindex.(chain_results, 3), getindex.(chain_results, 4),
            getindex.(chain_results, 5))
    end

    n_ch = length(chains.traces)
    all_draws = [_draws_for_csv(chains, i, discretize_dt) for i in 1:n_ch]
    min_rows = minimum(size(draws, 1) for draws in all_draws)
    csv_paths = String[]
    for chain_idx in 1:n_ch
        draws_unc = all_draws[chain_idx][1:min_rows, :]
        csv_path = n_ch == 1 ? output_csv : replace(output_csv, r"\.csv$" => "_chain$(chain_idx).csv")
        r_constrain_and_write_csv(contexts[chain_idx].sm_full, draws_unc, csv_path;
            chain_id=chain_idx, compute_lp)
        push!(csv_paths, csv_path)
    end

    result = _pack_result(chains; requested_flow=flow_type,
        adaptive_scheme, initial_refresh_rate,
        warmup_adaptation_interval, t_warmup)
    result["csv_paths"] = csv_paths
    result["bridge_call_counts"] = [_subsampling_call_counts(ctx, manager)
        for (ctx, manager) in zip(contexts, anchor_managers)]
    result["subsampling"] = true
    return result
end

function r_pdmp_brms_subsampling(
        stan_file::String,
        data_full_file::String,
        data_prior_file::String,
        N::Integer,
        subsample_size::Integer,
        family_name::String,
        predictor_designs,
        predictor_indices,
        predictor_dimension::Integer,
        offsets,
        response,
        known_se,
        observation_multipliers,
        flow_type::String,
        algorithm_type::String,
        flow_mean,
        flow_cov,
        output_csv::String;
        c0::Float64 = 1e-2,
        grid_n::Int = 30,
        grid_t_max::Float64 = 2.0,
        grid_bound::String = "value_quadratic",
        grid_curvature_bound::Union{Nothing,Float64} = nothing,
        t0::Float64 = 0.0,
        T::Float64 = 10000.0,
        t_warmup::Float64 = 0.0,
        adaptive_scheme::String = "diagonal",
        discretize_dt::Float64 = 0.0,
        show_progress::Bool = true,
        n_chains::Int = 1,
        threaded::Bool = false,
        seed::Union{Integer, Nothing} = nothing,
        compute_lp::Bool = false,
        curvature_backend::String = "exact",
        sticky::Bool = false,
        can_stick::Union{AbstractVector{Bool}, Nothing} = nothing,
        model_prior = nothing,
        parameter_prior::Union{AbstractVector{Float64}, Nothing} = nothing,
        n_anchor_updates::Int = 0,
        use_anchor_bank::Bool = false,
        bank_capacity::Int = 20,
        use_hcv::Bool = false
    )
    lib_path_std = _compile_model(stan_file)

    N_int = Int(N)
    m = Int(subsample_size)
    designs = Matrix{Float64}[Matrix{Float64}(design)
        for design in _subsampling_collection_values(predictor_designs)]
    indices = Vector{Int}[_as_subsampling_predictor_indices(active)
        for active in _subsampling_collection_values(predictor_indices)]
    d = Int(predictor_dimension)

    return _run_subsampling_brms_subsampled(
        lib_path_std, data_full_file, data_prior_file, N_int, m,
        Symbol(family_name), designs, indices, d,
        Matrix{Float64}(offsets), Matrix{Float64}(response), Vector{Float64}(known_se),
        Vector{Float64}(observation_multipliers),
        flow_type, algorithm_type,
        flow_mean, flow_cov, output_csv;
        c0, grid_n, grid_t_max, grid_bound, grid_curvature_bound, t0, T, t_warmup,
        adaptive_scheme, discretize_dt, show_progress, n_chains, threaded,
        seed, compute_lp, curvature_backend, sticky, can_stick, model_prior, parameter_prior,
        n_anchor_updates, use_anchor_bank, bank_capacity, use_hcv)
end

function r_pdmp_stan_for_brms(
        path_to_stan_model::String,
        path_to_stan_data::String,
        flow_type::String,
        algorithm_type::String,
        flow_mean,
        flow_cov,
        output_csv::String;
        c0::Float64 = 1e-2,
        grid_n::Int = 30,
        grid_t_max::Float64 = 2.0,
        grid_bound::String = "value_quadratic",
        grid_curvature_bound::Union{Nothing,Float64} = nothing,
        t0::Float64 = 0.0,
        T::Float64 = 10000.0,
        t_warmup::Float64 = 0.0,
        adaptive_scheme::String = "diagonal",
        discretize_dt::Float64 = 0.0,
        show_progress::Bool = true,
        n_chains::Int = 1,
        threaded::Bool = false,
        seed::Union{Integer, Nothing} = nothing,
        compute_lp::Bool = false,
        curvature_backend::String = "exact",
        sticky::Bool = false,
        can_stick = nothing,
        model_prior = nothing,
        parameter_prior = nothing,
        slab_prior = nothing,
        unc_names = String[]
    )
    sm = BridgeStan.StanModel(path_to_stan_model, path_to_stan_data; warn=false)
    posterior_model = PDMPModel(sm; hvp = isnothing(slab_prior) && curvature_backend != "finite_difference")
    can_stick_vec = _as_bool_vector(can_stick)
    parameter_prior_vec = isnothing(parameter_prior) ? nothing : _as_float_vector(parameter_prior)
    unc_names_vec = _as_string_vector(unc_names)
    model = if isnothing(slab_prior)
        posterior_model
    else
        _build_dependent_slab_model(posterior_model, model_prior, slab_prior,
            can_stick_vec, unc_names_vec)
    end
    d = model.d

    fmean = _as_flow_mean(flow_mean, d)
    flow_cov_mat = _as_flow_cov(flow_cov, d)
    prec = _to_precision(flow_cov_mat, d)

    flow = build_flow(flow_type, prec, fmean; adaptive_scheme)
    initial_refresh_rate = try Float64(PDMPSamplers.refresh_rate(flow)) catch; NaN end
    alg0 = build_algorithm(algorithm_type; c0, d, grid_n, grid_t_max,
        curvature_backend, grid_bound, grid_curvature_bound)
    alg = isnothing(slab_prior) ?
        wrap_sticky(alg0, sticky, model_prior, parameter_prior_vec, can_stick_vec) :
        wrap_dependent_sticky(alg0, sticky, model_prior, slab_prior, can_stick_vec, flow, unc_names_vec)

    chains = pdmp_sample(d, flow, model, alg, t0, T, t_warmup;
                         progress = show_progress, n_chains, threaded, seed,
                         statistic_counter=PDMPSamplers.DevelStatisticCounter)

    stats = extract_stats(chains)
    n_ch = length(chains.traces)
    csv_paths = String[]
    all_draws = [_draws_for_csv(chains, i, discretize_dt) for i in 1:n_ch]
    min_rows = minimum(size(m, 1) for m in all_draws)
    for chain_idx in 1:n_ch
        draws_unc = all_draws[chain_idx][1:min_rows, :]
        csv_path = n_ch == 1 ? output_csv : replace(output_csv, r"\.csv$" => "_chain$(chain_idx).csv")
        r_constrain_and_write_csv(sm, draws_unc, csv_path; chain_id = chain_idx, compute_lp)
        push!(csv_paths, csv_path)
    end

    result = _pack_result(chains; requested_flow=flow_type,
        adaptive_scheme, initial_refresh_rate,
        warmup_adaptation_interval, t_warmup)
    result["csv_paths"] = csv_paths
    if sticky
        incl = Dict{Int,Vector{Float64}}()
        for i in 1:n_ch
            incl[i] = inclusion_probs(chains; chain=i)
        end
        result["inclusion_probs"] = incl
    end
    return result
end

function r_get_param_unc_names(path_to_stan_model::String, path_to_stan_data::String)
    sm = BridgeStan.StanModel(path_to_stan_model, path_to_stan_data; warn=false)
    return BridgeStan.param_unc_names(sm)
end

r_threading_available() = Base.Threads.nthreads() > 1

end # module PDMPSamplersRBridge

using PDMPSamplers, LinearAlgebra, BridgeStan, Random
using .PDMPSamplersRBridge

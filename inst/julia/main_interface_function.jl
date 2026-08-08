module PDMPSamplersRBridge

using PDMPSamplers, LinearAlgebra, BridgeStan, Random, SparseArrays, Statistics, Libdl

export build_flow, build_algorithm, wrap_sticky
export build_model_prior_odds, build_slab_provider, wrap_dependent_sticky
export r_discretize, r_mean, r_var, r_std, r_cov, r_cor, r_quantile, r_median, r_cdf, r_ess, r_summary_all
export r_inclusion_probs, extract_stats
export r_chain_times, r_chain_positions, r_chain_velocities, r_chain_is_boomerang, r_chain_is_mutable_boomerang, r_chain_mu
export r_from_skeleton
export r_chain_is_factorized
export r_chain_sparse_initial_time, r_chain_sparse_initial_position, r_chain_sparse_initial_velocity
export r_chain_sparse_event_indices, r_chain_sparse_event_times, r_chain_sparse_event_positions, r_chain_sparse_event_velocities
export r_from_sparse_skeleton
export r_pdmp_stan, r_pdmp_custom
export r_stan_marked_diagnostics
export write_cmdstan_csv, r_constrain_and_write_csv
export r_pdmp_brms_marked, r_pdmp_stan_for_brms
export r_get_param_unc_names
export r_threading_available

_marked_collection_values(values::Union{AbstractVector,Tuple}) = values
_marked_collection_values(values::NamedTuple) = Base.values(values)
_marked_collection_values(values::AbstractDict) = Base.values(values)
_as_marked_predictor_indices(index::Number) = [Int(index)]
_as_marked_predictor_indices(indices) = collect(Int, vec(indices))
_marked_float_vector(value::Number) = [Float64(value)]
_marked_float_vector(values) = collect(Float64, vec(values))

function build_flow(flow_type::String, prec::AbstractMatrix{Float64}, flow_mean::AbstractVector{Float64};
                    adaptive_scheme::String="diagonal")
    if flow_type == "ZigZag"
        return ZigZag(prec, flow_mean)
    elseif flow_type == "BouncyParticle"
        return BouncyParticle(prec, flow_mean)
    elseif flow_type == "Boomerang"
        return Boomerang(prec, flow_mean)
    elseif flow_type == "AdaptiveBoomerang"
        d = length(flow_mean)
        λref = parse(Float64, get(ENV, "PDMP_ADAPTIVE_BOOMERANG_LAMBDA_REF", "0.1"))
        return AdaptiveBoomerang(d; λref=λref, scheme=Symbol(adaptive_scheme))
    elseif flow_type == "PreconditionedZigZag"
        d = length(flow_mean)
        return PreconditionedZigZag(prec, flow_mean)
    elseif flow_type == "PreconditionedBPS"
        d = length(flow_mean)
        return PreconditionedBPS(prec, flow_mean)
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

function build_algorithm(algorithm_type::String; c0::Float64, d::Integer, grid_n::Int, grid_t_max::Float64,
        use_fd_hvp::Bool=false, curvature_backend::String="auto",
        post_warmup_simplify::Bool=false, grid_bound::String="constant",
        grid_curvature_bound::Union{Nothing,Float64}=nothing,
        linear_area_threshold::Float64=0.95, linear_min_area_gain::Float64=0.0,
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
            use_fd_hvp = use_fd_hvp, curvature_backend = Symbol(curvature_backend),
            post_warmup_simplify = post_warmup_simplify,
            bound = Symbol(grid_bound),
            bound_violation = bound_violation,
            curvature_bound = grid_curvature_bound,
            linear_area_threshold = linear_area_threshold,
            linear_min_area_gain = linear_min_area_gain,
            lazy_low_tightness_threshold = lazy_low_tightness_threshold,
            lazy_max_low_tightness_rejections = lazy_max_low_tightness_rejections,
            lazy_max_rejections = lazy_max_rejections)
    elseif algorithm_type == "PositiveVariationGridThinningStrategy"
        pv_max_skip_width = parse(Float64, get(ENV, "PDMP_POSITIVE_VARIATION_MAX_SKIP_WIDTH", "0.25"))
        pv_dense_cell_width = parse(Float64, get(ENV, "PDMP_POSITIVE_VARIATION_DENSE_CELL_WIDTH", "0.0"))
        pv_skip_slope_safety = parse(Float64, get(ENV, "PDMP_POSITIVE_VARIATION_SKIP_SLOPE_SAFETY", "0.0"))
        pv_use_derivative_hermite = parse(Bool, get(ENV, "PDMP_POSITIVE_VARIATION_USE_DERIVATIVE_HERMITE", "false"))
        pv_derivative_hermite_on_demand = parse(Bool, get(ENV, "PDMP_POSITIVE_VARIATION_DERIVATIVE_HERMITE_ON_DEMAND", "false"))
        pv_derivative_hermite_trigger_scale = parse(Float64, get(ENV, "PDMP_POSITIVE_VARIATION_DERIVATIVE_HERMITE_TRIGGER_SCALE", "10.0"))
        pv_validation_rtol = parse(Float64, get(ENV, "PDMP_POSITIVE_VARIATION_VALIDATION_RTOL", "0.05"))
        pv_approximate = parse(Bool, get(ENV, "PDMP_POSITIVE_VARIATION_APPROXIMATE", "false"))
        return PositiveVariationGridThinningStrategy(; N = grid_n, t_max = grid_t_max,
            validation_rtol = pv_validation_rtol,
            max_skip_width = pv_max_skip_width,
            dense_cell_width = pv_dense_cell_width,
            skip_slope_safety = pv_skip_slope_safety,
            use_derivative_hermite = pv_use_derivative_hermite,
            derivative_hermite_on_demand = pv_derivative_hermite_on_demand,
            derivative_hermite_trigger_scale = pv_derivative_hermite_trigger_scale,
            approximate = pv_approximate,
            fallback = GridThinningStrategy(; N = grid_n, t_max = grid_t_max,
                use_fd_hvp = use_fd_hvp, curvature_backend = Symbol(curvature_backend),
                post_warmup_simplify = post_warmup_simplify,
                bound = Symbol(grid_bound),
                curvature_bound = grid_curvature_bound,
                linear_area_threshold = linear_area_threshold,
                linear_min_area_gain = linear_min_area_gain,
                lazy_low_tightness_threshold = lazy_low_tightness_threshold,
                lazy_max_low_tightness_rejections = lazy_max_low_tightness_rejections,
                lazy_max_rejections = lazy_max_rejections))
    elseif algorithm_type == "VectorVariationThinningStrategy"
        vv_max_skip_width = parse(Float64, get(ENV, "PDMP_VECTOR_VARIATION_MAX_SKIP_WIDTH", "0.25"))
        vv_validation_rtol = parse(Float64, get(ENV, "PDMP_VECTOR_VARIATION_VALIDATION_RTOL", "0.05"))
        vv_max_refinement_depth = parse(Int, get(ENV, "PDMP_VECTOR_VARIATION_MAX_REFINEMENT_DEPTH", "8"))
        vv_n_min = parse(Int, get(ENV, "PDMP_VECTOR_VARIATION_N_MIN", "5"))
        vv_use_derivative_hermite = parse(Bool, get(ENV, "PDMP_VECTOR_VARIATION_USE_DERIVATIVE_HERMITE", "false"))
        vv_derivative_hermite_on_demand = parse(Bool, get(ENV, "PDMP_VECTOR_VARIATION_DERIVATIVE_HERMITE_ON_DEMAND", "false"))
        vv_derivative_hermite_trigger_scale = parse(Float64, get(ENV, "PDMP_VECTOR_VARIATION_DERIVATIVE_HERMITE_TRIGGER_SCALE", "10.0"))
        return VectorVariationThinningStrategy(; N = grid_n, t_max = grid_t_max,
            N_min = vv_n_min,
            validation_rtol = vv_validation_rtol,
            max_skip_width = vv_max_skip_width,
            max_refinement_depth = vv_max_refinement_depth,
            use_derivative_hermite = vv_use_derivative_hermite,
            derivative_hermite_on_demand = vv_derivative_hermite_on_demand,
            derivative_hermite_trigger_scale = vv_derivative_hermite_trigger_scale,
            fallback = GridThinningStrategy(; N = grid_n, t_max = grid_t_max,
                use_fd_hvp = use_fd_hvp, curvature_backend = Symbol(curvature_backend),
                post_warmup_simplify = post_warmup_simplify,
                bound = Symbol(grid_bound),
                curvature_bound = grid_curvature_bound,
                linear_area_threshold = linear_area_threshold,
                linear_min_area_gain = linear_min_area_gain,
                lazy_low_tightness_threshold = lazy_low_tightness_threshold,
                lazy_max_low_tightness_rejections = lazy_max_low_tightness_rejections,
                lazy_max_rejections = lazy_max_rejections))
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
        κ = BetaBernoulliKappa(model_prior[:a]::Float64, model_prior[:b]::Float64, parameter_prior)
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

function _as_float_sparse_matrix(spec)
    return SparseMatrixCSC(
        Int(_rget(spec, :nrow)),
        Int(_rget(spec, :ncol)),
        _as_int_vector(_rget(spec, :colptr)),
        _as_int_vector(_rget(spec, :rowval)),
        _as_float_vector(_rget(spec, :nzval)),
    )
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
    elseif _haskey(model_prior, :omega)
        omega = _as_float_vector(_rget(model_prior, :omega))
        length(omega) == m + 1 ||
            throw(DimensionMismatch("exchangeable model-size prior must have length beta dimension + 1"))
        return ExchangeableModelSizePrior(log.(omega); normalize=true)
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
    elseif t == "independent_slab_density"
        κ = _as_float_vector(_rget(slab_prior, :kappa))
        length(κ) == 1 && (κ = fill(κ[1], length(beta_idx)))
        length(κ) == length(beta_idx) ||
            throw(DimensionMismatch("independent_slab_density kappa length must be 1 or beta dimension"))
        return IndependentZeroMeanGaussianSlab(κ, beta_idx)
    elseif t == "independent_logscale_gaussian"
        log_base_scales = _as_float_vector(_rget(slab_prior, :log_base_scales))
        length(log_base_scales) == 1 && (log_base_scales = fill(log_base_scales[1], length(beta_idx)))
        length(log_base_scales) == length(beta_idx) ||
            throw(DimensionMismatch("independent_logscale_gaussian_slab log_base_scales length must be 1 or beta dimension"))
        logscale_idx = _state_indices(_rget(slab_prior, :logscale), unc_names, d, "logscale")
        length(logscale_idx) == 1 && (logscale_idx = fill(logscale_idx[1], length(beta_idx)))
        length(logscale_idx) == length(beta_idx) ||
            throw(DimensionMismatch("independent_logscale_gaussian_slab logscale length must be 1 or beta dimension"))
        isempty(intersect(beta_idx, logscale_idx)) ||
            throw(ArgumentError("independent_logscale_gaussian_slab logscale coordinates must be disjoint from slab coef coordinates"))
        any(i -> can_stick[i], unique(logscale_idx)) &&
            throw(ArgumentError("independent_logscale_gaussian_slab logscale coordinates must be non-stickable"))
        return IndependentZeroMeanLogscaleGaussianSlab(beta_idx, logscale_idx, log_base_scales)
    elseif t == "loglinear_gaussian_scale"
        log_base_scales = _as_float_vector(_rget(slab_prior, :log_base_scales))
        length(log_base_scales) == 1 &&
            (log_base_scales = fill(log_base_scales[1], length(beta_idx)))
        logscale_idx = _state_indices(_rget(slab_prior, :logscale), unc_names, d, "logscale")
        sparse_design = _rget(slab_prior, :logscale_design_sparse)
        design = isnothing(sparse_design) ?
            _as_float_matrix(_rget(slab_prior, :logscale_design)) :
            _as_float_sparse_matrix(sparse_design)
        size(design) == (length(beta_idx), length(logscale_idx)) ||
            throw(DimensionMismatch("logscale_design must have size (beta dimension, logscale dimension)"))
        isempty(intersect(beta_idx, logscale_idx)) || throw(ArgumentError(
            "loglinear_gaussian_scale_slab logscale coordinates must be disjoint from slab coef coordinates"))
        any(i -> can_stick[i], logscale_idx) && throw(ArgumentError(
            "loglinear_gaussian_scale_slab logscale coordinates must be non-stickable"))
        return LogLinearGaussianScaleSlab(beta_idx, logscale_idx,
            log_base_scales, design)
    elseif t == "global_logscale_exchangeable_gaussian"
        logscale_idx = _state_indices(_rget(slab_prior, :logscale), unc_names, d, "logscale")
        length(logscale_idx) == 1 ||
            throw(DimensionMismatch("global_logscale_exchangeable_gaussian_slab logscale length must be 1"))
        logscale_idx[1] in beta_idx &&
            throw(ArgumentError("global_logscale_exchangeable_gaussian_slab logscale coordinate must be disjoint from slab coef coordinates"))
        can_stick[logscale_idx[1]] &&
            throw(ArgumentError("global_logscale_exchangeable_gaussian_slab logscale coordinate must be non-stickable"))
        return GlobalLogscaleExchangeableGaussianSlab(beta_idx, logscale_idx[1],
            Float64(_rget(slab_prior, :u)), Float64(_rget(slab_prior, :v));
            mean=Float64(_rget(slab_prior, :mean)),
            logscale_offset=Float64(_rget(slab_prior, :logscale_offset)))
    elseif t == "callback_gaussian"
        mean_cov_r = _rget(slab_prior, :mean_cov)
        active_prior_neggrad_r = _rget(slab_prior, :active_prior_neggrad)
        mean_cov! = (mean_out, cov_out, x) -> begin
            spec = mean_cov_r(Vector{Float64}(x))
            copyto!(mean_out, _as_float_vector(_rget(spec, :mean)))
            copyto!(cov_out, _as_float_matrix(_rget(spec, :cov)))
            return nothing
        end
        active_prior_grad! = isnothing(active_prior_neggrad_r) ? nothing :
            (out, x, active) -> begin
                values = _as_float_vector(active_prior_neggrad_r(Vector{Float64}(x), Vector{Bool}(active)))
                copyto!(out, values)
                return out
            end
        return CallbackGaussianSlab(beta_idx; mean_cov!, active_prior_grad!)
    elseif t == "arbitrary_boundary"
        log_q_zero_r = _rget(slab_prior, :log_q_zero)
        active_prior_neggrad_r = _rget(slab_prior, :active_prior_neggrad)
        active_prior_neggrad! = (out, x, active) -> begin
            values = _as_float_vector(active_prior_neggrad_r(Vector{Float64}(x), Vector{Bool}(active)))
            copyto!(out, values)
            return out
        end
        log_q_zero! = (x, active, j) -> Float64(log_q_zero_r(Vector{Float64}(x), Vector{Bool}(active), Int(j)))
        return ArbitrarySlabBoundary(beta_idx; active_prior_neggrad!, log_q_zero!)
    else
        throw(ArgumentError("unknown slab_prior type: $t"))
    end
end

function wrap_dependent_sticky(alg::PDMPSamplers.PoissonTimeStrategy, sticky::Bool,
        model_prior, slab_prior, can_stick, flow_type::String,
        unc_names::AbstractVector{<:AbstractString}=String[])
    sticky || throw(ArgumentError("wrap_dependent_sticky requires sticky=true"))
    isnothing(slab_prior) && throw(ArgumentError("wrap_dependent_sticky requires slab_prior"))
    d = length(can_stick)
    provider = build_slab_provider(slab_prior, unc_names, can_stick, d)
    odds = build_model_prior_odds(model_prior, PDMPSamplers.beta_indices(provider), d)
    clock = default_aggregate_unstick_clock(provider, odds)
    return AggregateSticky(alg, clock, BitVector(can_stick))
end

function _validate_sampling_slab_prior(slab_prior)
    isnothing(slab_prior) && return nothing
    if _slab_type(slab_prior) == "callback_gaussian" && isnothing(_rget(slab_prior, :active_prior_neggrad))
        throw(ArgumentError("gaussian_scale_mixture_slab requires active_prior_neggrad when used for dependent-slab sampling"))
    end
    return nothing
end

function _build_dependent_slab_model(
        posterior_model::PDMPModel,
        model_prior,
        slab_prior,
        can_stick,
        unc_names::AbstractVector{<:AbstractString}=String[])
    _validate_sampling_slab_prior(slab_prior)
    d = posterior_model.d
    provider = build_slab_provider(slab_prior, unc_names, can_stick, d)
    odds = build_model_prior_odds(model_prior, PDMPSamplers.beta_indices(provider), d)
    target = DependentSlabTarget(d, posterior_model.grad, provider, odds)
    return PDMPModel(target; hvp=true)
end

function _to_precision(flow_cov::AbstractMatrix{Float64}, d::Int)
    if isdiag(flow_cov) && all(i -> @inbounds(flow_cov[i, i]) ≈ 1.0, 1:d)
        return Diagonal(ones(d))
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

function _cached_stan_source(path_to_stan_model::String,
        stanc_args::Vector{String}, make_args::Vector{String})
    cache_root = get(ENV, "PDMPSAMPLERSR_BRIDGESTAN_CACHE_DIR",
        joinpath(first(DEPOT_PATH), "pdmpsamplersr", "bridgestan"))
    source_bytes = read(path_to_stan_model)
    option_bytes = Vector{UInt8}(codeunits(join([stanc_args; make_args], '\0')))
    checksum = _compile_cache_checksum([source_bytes; option_bytes])
    build_dir = joinpath(cache_root, string(checksum; base=16, pad=16))
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
    source = _compile_env_enabled("PDMPSAMPLERSR_BRIDGESTAN_CACHE") ?
        _cached_stan_source(path_to_stan_model, stanc_args, make_args) : path_to_stan_model
    library = replace(source, r"\.stan$" => "_model.$(Libdl.dlext)")
    isfile(library) && return library
    return BridgeStan.compile_model(source; stanc_args, make_args)
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
    source = _compile_env_enabled("PDMPSAMPLERSR_BRIDGESTAN_CACHE") ?
        _cached_stan_source(path_to_stan_model, stanc_args, make_args) :
        path_to_stan_model
    library = replace(source, r"\.stan$" => "_model.$(Libdl.dlext)")
    isfile(library) && mtime(library) >=
        max(mtime(source), mtime(hpp_path)) && return library
    return BridgeStan.compile_model(source; stanc_args, make_args)
end

function _pack_result(chains::PDMPChains)
    stats = extract_stats(chains)
    d = length(first(chains.traces[1]).position)
    return Dict{String,Any}(
        "chains"   => chains,
        "stats"    => stats,
        "d"        => d,
        "n_chains" => PDMPSamplers.n_chains(chains)
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

function r_chain_sparse_initial_time(chains::PDMPChains; chain::Int)
    Float64(chains.traces[chain].initial_state.time)
end

function r_chain_sparse_initial_position(chains::PDMPChains; chain::Int)
    Vector{Float64}(chains.traces[chain].initial_state.position)
end

function r_chain_sparse_initial_velocity(chains::PDMPChains; chain::Int)
    Vector{Float64}(chains.traces[chain].initial_state.velocity)
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
    PDMPChains(traces, PDMPSamplers.StatisticCounter[])
end

function _r_factorized_trace(events, flow, initial_state)
    bounds_cache = Dict{Int, Tuple{Float64, Float64}}()
    if applicable(PDMPSamplers.FactorizedTrace,
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

function r_chain_is_boomerang(chains::PDMPChains; chain::Int)
    PDMPSamplers._underlying_flow(chains.traces[chain].flow) isa AnyBoomerang
end

function r_chain_is_mutable_boomerang(chains::PDMPChains; chain::Int)
    PDMPSamplers._underlying_flow(chains.traces[chain].flow) isa MutableBoomerang
end

function r_chain_mu(chains::PDMPChains; chain::Int)
    base = PDMPSamplers._underlying_flow(chains.traces[chain].flow)
    base isa AnyBoomerang ? Vector{Float64}(base.μ) : Float64[]
end

function _final_lambda_ref(chains::PDMPChains, chain::Int)
    base = PDMPSamplers._underlying_flow(chains.traces[chain].flow)
    return base isa AnyBoomerang ? Float64(base.λref) : NaN
end

function r_from_skeleton(times_list::AbstractVector, positions_list::AbstractVector,
                          velocities_list::AbstractVector, is_boomerang_list::AbstractVector,
                          mu_list::AbstractVector, is_mutable_boomerang_list::AbstractVector)
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
        flow = if is_mutable
            MutableBoomerang(I(d), mu)
        elseif is_boom
            Boomerang(I(d), mu)
        else
            ZigZag(I(d), zeros(d))
        end
        push!(traces, PDMPTrace(times, positions, velocities, flow))
    end
    PDMPChains(traces, PDMPSamplers.StatisticCounter[])
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
        "positive_variation_cells" => vals(:positive_variation_cells),
        "positive_variation_refinements" => vals(:positive_variation_refinements),
        "positive_variation_fallbacks" => vals(:positive_variation_fallbacks),
        "positive_variation_accepts" => vals(:positive_variation_accepts),
        "positive_variation_skipped_cells" => vals(:positive_variation_skipped_cells),
        "grid_resets_from_dynamics_adaptation" => vals(:grid_resets_from_dynamics_adaptation),
        "grid_endpoint_evaluations" => vals(:grid_endpoint_evaluations),
        "grid_endpoint_gradient_calls" => vals(:grid_endpoint_gradient_calls),
        "grid_cached_endpoint_reuses" => vals(:grid_cached_endpoint_reuses),
        "grid_acceptance_tests" => vals(:grid_acceptance_tests),
        "grid_acceptance_gradient_calls" => vals(:grid_acceptance_gradient_calls),
        "marked_cell_roof_proposals" => vals(:marked_cell_roof_proposals),
        "marked_aggregate_accepts" => vals(:marked_aggregate_accepts),
        "marked_subset_evaluations" => vals(:marked_subset_evaluations),
        "marked_final_reflections" => vals(:marked_final_reflections),
        "grid_horizon_hits" => vals(:grid_horizon_hits),
        "grid_budget_tail_restarts" => vals(:grid_budget_tail_restarts),
        "constant_bound_attempts" => vals(:constant_bound_attempts),
        "constant_bound_accepts" => vals(:constant_bound_accepts),
        "constant_bound_rejections" => vals(:constant_bound_rejections),
        "constant_bound_violations" => vals(:constant_bound_violations),
        "constant_bound_safety_fallbacks" => vals(:constant_bound_safety_fallbacks),
        "grid_bound_violations" => vals(:grid_bound_violations),
        "shared_node_cells" => vals(:shared_node_cells),
        "shared_node_two_point_cells" => vals(:shared_node_two_point_cells),
        "shared_node_three_point_cells" => vals(:shared_node_three_point_cells),
        "sticky_inner_searches" => vals(:sticky_inner_searches),
        "sticky_inner_wins" => vals(:sticky_inner_wins),
        "sticky_inner_wasted_by_sticky" => vals(:sticky_inner_wasted_by_sticky),
        "sticky_inner_wasted_by_refresh" => vals(:sticky_inner_wasted_by_refresh),
        "sticky_all_frozen_events" => vals(:sticky_all_frozen_events),
        "warmup_events" => vals(:warmup_events),
        "main_events" => vals(:main_events),
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
        "warmup_elapsed_time" => vals(:warmup_elapsed_time),
        "main_elapsed_time" => vals(:main_elapsed_time),
        "boomerang_interference_events" => vals(:boomerang_interference_events),
        "boomerang_target_c_share_sum" => vals(:boomerang_target_c_share_sum),
        "boomerang_target_d_share_sum" => vals(:boomerang_target_d_share_sum),
        "boomerang_nuisance_driven_target_disturbances" => vals(:boomerang_nuisance_driven_target_disturbances),
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

    use_fd_hvp = Bool(get(kwargs, :use_fd_hvp, false))
    model = PDMPModel(path_to_stan_model, path_to_stan_data; hvp=!use_fd_hvp)
    return r_pdmp_stan(model, x0, flow_type, algorithm_type, flow_mean, flow_cov; kwargs...)
end

function r_pdmp_stan(
        model::PDMPModel,
        x0,
        flow_type::String,
        algorithm_type::String,
        flow_mean,
        flow_cov;
        c0::Float64 = 1e-2,
        grid_n::Int = 30,
        grid_t_max::Float64 = 2.0,
        use_fd_hvp::Bool = false,
        curvature_backend::String = "auto",
        post_warmup_simplify::Bool = true,
        grid_bound::String = "constant",
        grid_curvature_bound::Union{Nothing,Float64} = nothing,
        linear_area_threshold::Float64 = 0.95,
        linear_min_area_gain::Float64 = 0.0,
        lazy_low_tightness_threshold::Float64 = 0.1,
        lazy_max_low_tightness_rejections::Int = 3,
        lazy_max_rejections::Int = 0,
        t0::Float64 = 0.0,
        T::Float64 = 10000.0,
        t_warmup::Float64 = 0.0,
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
        support_boundary_mode::String = "error",
        support_boundary_max_bisection_steps::Int = 60,
        support_boundary_time_rtol::Float64 = 1e-8,
        support_boundary_time_atol::Float64 = 1e-10,
        support_boundary_clip_fraction::Float64 = 1 - 1e-10,
        support_boundary_max_refresh_attempts::Int = 20,
        support_boundary_refresh_probe_time::Float64 = 1e-4,
        support_boundary_min_safe_time::Float64 = 1e-12
    )

    d = model.d
    x0_vec = _as_float_vector(x0)
    flow_mean_vec = _as_flow_mean(flow_mean, d)
    flow_cov_mat = _as_flow_cov(flow_cov, d)
    can_stick_vec = _as_bool_vector(can_stick)
    parameter_prior_vec = isnothing(parameter_prior) ? nothing : _as_float_vector(parameter_prior)
    unc_names_vec = _as_string_vector(unc_names)
    sampling_model = if isnothing(slab_prior)
        model
    else
        _build_dependent_slab_model(model, model_prior, slab_prior, can_stick_vec, unc_names_vec)
    end

    prec = _to_precision(flow_cov_mat, d)
    flow = build_flow(flow_type, prec, flow_mean_vec; adaptive_scheme)
    alg0 = build_algorithm(algorithm_type; c0, d, grid_n, grid_t_max, use_fd_hvp, curvature_backend,
        post_warmup_simplify, grid_bound, grid_curvature_bound,
        linear_area_threshold, linear_min_area_gain,
        lazy_low_tightness_threshold, lazy_max_low_tightness_rejections,
        lazy_max_rejections)
    alg = isnothing(slab_prior) ?
        wrap_sticky(alg0, sticky, model_prior, parameter_prior_vec, can_stick_vec) :
        wrap_dependent_sticky(alg0, sticky, model_prior, slab_prior, can_stick_vec, flow_type, unc_names_vec)

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

    chains = pdmp_sample(x0_vec, flow, sampling_model, alg, t0, T, t_warmup;
                         progress = show_progress, n_chains = n_chains, threaded = threaded,
                         seed = seed,
                         warmup_stop = build_warmup_stop(t_warmup),
                         support_boundary_options = sbopts,
                         statistic_counter = PDMPSamplers.DevelStatisticCounter)
    return _pack_result(chains)
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
        post_warmup_simplify::Bool = true,
        grid_bound::String = "constant",
        grid_curvature_bound::Union{Nothing,Float64} = nothing,
        linear_area_threshold::Float64 = 0.95,
        linear_min_area_gain::Float64 = 0.0,
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
        support_boundary_mode::String = "error",
        support_boundary_max_bisection_steps::Int = 60,
        support_boundary_time_rtol::Float64 = 1e-8,
        support_boundary_time_atol::Float64 = 1e-10,
        support_boundary_clip_fraction::Float64 = 1 - 1e-10,
        support_boundary_max_refresh_attempts::Int = 20,
        support_boundary_refresh_probe_time::Float64 = 1e-4,
        support_boundary_min_safe_time::Float64 = 1e-12
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
        _validate_sampling_slab_prior(slab_prior)
        provider = build_slab_provider(slab_prior, String[], can_stick_vec, d)
        odds = build_model_prior_odds(model_prior, PDMPSamplers.beta_indices(provider), d)
        target = DependentSlabTarget(d, grad!, provider, odds)
        PDMPModel(target; hvp=true)
    end

    prec = _to_precision(flow_cov_mat, d)
    flow = build_flow(flow_type, prec, flow_mean_vec; adaptive_scheme)
    alg0 = build_algorithm(algorithm_type; c0, d, grid_n, grid_t_max,
        post_warmup_simplify, grid_bound, grid_curvature_bound,
        linear_area_threshold, linear_min_area_gain)
    alg = isnothing(slab_prior) ?
        wrap_sticky(alg0, sticky, model_prior, parameter_prior_vec, can_stick_vec) :
        wrap_dependent_sticky(alg0, sticky, model_prior, slab_prior, can_stick_vec, flow_type)

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
                         support_boundary_options = sbopts)
    return _pack_result(chains)
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

include(joinpath(@__DIR__, "marked_family_provider.jl"))
include(joinpath(@__DIR__, "marked_stan_provider.jl"))

function _run_marked_brms_subsampled(lib_path_std::String,
        data_full_file::String, data_prior_file::String,
        N::Int, m::Int, family::Symbol,
        predictor_designs::Vector{Matrix{Float64}},
        predictor_indices::Vector{Vector{Int}}, d::Int,
        offsets::Matrix{Float64}, response::Matrix{Float64}, known_se::Vector{Float64},
        observation_multipliers::Vector{Float64},
        flow_type::String, algorithm_type::String, flow_mean, flow_cov,
        output_csv::String; c0::Float64, grid_n::Int, grid_t_max::Float64,
        grid_bound::String, grid_curvature_bound::Union{Nothing,Float64},
        linear_area_threshold::Float64, linear_min_area_gain::Float64,
        t0::Float64, T::Float64, t_warmup::Float64,
        adaptive_scheme::String, discretize_dt::Float64,
        show_progress::Bool, n_chains::Int, threaded::Bool,
        seed::Union{Integer,Nothing}, compute_lp::Bool,
        use_fd_hvp::Bool, post_warmup_simplify::Bool,
        sticky::Bool, can_stick, model_prior, parameter_prior,
        n_anchor_updates::Int, use_anchor_bank::Bool, bank_capacity::Int,
        use_hcv::Bool)

    algorithm_type in ("GridThinningStrategy", "ThinningStrategy") ||
        throw(ArgumentError(
            "marked BridgeStan subsampling requires GridThinningStrategy or ThinningStrategy"))
    all(size(design, 1) == N for design in predictor_designs) ||
        throw(DimensionMismatch("predictor designs must have N rows"))
    length(observation_multipliers) == N || throw(DimensionMismatch(
        "observation multiplier vector must have length N"))

    fmean = _as_flow_mean(flow_mean, d)
    anchor = copy(fmean)
    contexts = MarkedBridgeStanContext[]
    models = PDMPModel[]
    anchor_adapters = Any[]
    anchor_managers = Union{Nothing,MarkedBridgeAnchorManager}[]
    anchor_capacity = use_anchor_bank ? bank_capacity :
        (n_anchor_updates > 0 ? 1 : 0)
    for _ in 1:n_chains
        ctx = _new_marked_bss_context(lib_path_std,
            data_full_file, data_prior_file, family, predictor_designs,
            predictor_indices, d, offsets, response, known_se,
            observation_multipliers, N, m, anchor; use_hcv)
        push!(contexts, ctx)
        built = _build_marked_bss_model(ctx; use_fd_hvp, anchor_capacity)
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
    alg = build_algorithm(algorithm_type; c0, d, grid_n, grid_t_max,
        use_fd_hvp, post_warmup_simplify, grid_bound, grid_curvature_bound,
        linear_area_threshold, linear_min_area_gain)
    alg = wrap_sticky(alg, sticky, model_prior, parameter_prior, can_stick)
    if anchor_capacity == 0
        chains = pdmp_sample(anchor, flow, models, alg, t0, T, t_warmup;
            progress=show_progress, threaded, seed,
            statistic_counter=PDMPSamplers.DevelStatisticCounter)
    else
        run_chain = function(i)
            seed_i = isnothing(seed) ? nothing : seed + i - 1
            flow_i = deepcopy(flow)
            dynamics_adapter = PDMPSamplers.default_adapter(
                flow_i, models[i].grad, t_warmup / 10, t_warmup, t0)
            adapter_i = dynamics_adapter isa PDMPSamplers.NoAdaptation ?
                anchor_adapters[i] : PDMPSamplers.SequenceAdapter(
                    (dynamics_adapter, anchor_adapters[i]))
            chain = pdmp_sample(anchor, flow_i, [models[i]], alg,
                t0, T, t_warmup; progress=show_progress && i == 1,
                adapter=adapter_i, seed=seed_i,
                statistic_counter=PDMPSamplers.DevelStatisticCounter)
            return chain.traces[1], chain.stats[1]
        end
        chain_results = if threaded && n_chains > 1
            fetch.([Threads.@spawn run_chain(i) for i in 1:n_chains])
        else
            [run_chain(i) for i in 1:n_chains]
        end
        chains = PDMPChains(first.(chain_results), last.(chain_results))
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

    result = _pack_result(chains)
    result["csv_paths"] = csv_paths
    result["bridge_call_counts"] = [_marked_call_counts(ctx, manager)
        for (ctx, manager) in zip(contexts, anchor_managers)]
    result["marked_subsampling"] = true
    return result
end

function r_pdmp_brms_marked(
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
        grid_bound::String = "constant",
        grid_curvature_bound::Union{Nothing,Float64} = nothing,
        linear_area_threshold::Float64 = 0.95,
        linear_min_area_gain::Float64 = 0.0,
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
        use_fd_hvp::Bool = false,
        post_warmup_simplify::Bool = false,
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
        for design in _marked_collection_values(predictor_designs)]
    indices = Vector{Int}[_as_marked_predictor_indices(active)
        for active in _marked_collection_values(predictor_indices)]
    d = Int(predictor_dimension)

    return _run_marked_brms_subsampled(
        lib_path_std, data_full_file, data_prior_file, N_int, m,
        Symbol(family_name), designs, indices, d,
        Matrix{Float64}(offsets), Matrix{Float64}(response), Vector{Float64}(known_se),
        Vector{Float64}(observation_multipliers),
        flow_type, algorithm_type,
        flow_mean, flow_cov, output_csv;
        c0, grid_n, grid_t_max, grid_bound, grid_curvature_bound,
        linear_area_threshold, linear_min_area_gain, t0, T, t_warmup,
        adaptive_scheme, discretize_dt, show_progress, n_chains, threaded,
        seed, compute_lp, use_fd_hvp, post_warmup_simplify,
        sticky, can_stick, model_prior, parameter_prior,
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
        grid_bound::String = "constant",
        grid_curvature_bound::Union{Nothing,Float64} = nothing,
        linear_area_threshold::Float64 = 0.95,
        linear_min_area_gain::Float64 = 0.0,
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
        use_fd_hvp::Bool = false,
        post_warmup_simplify::Bool = false,
        sticky::Bool = false,
        can_stick = nothing,
        model_prior = nothing,
        parameter_prior = nothing,
        slab_prior = nothing,
        unc_names = String[]
    )
    sm = BridgeStan.StanModel(path_to_stan_model, path_to_stan_data; warn=false)
    posterior_model = PDMPModel(sm; hvp = isnothing(slab_prior) && !use_fd_hvp)
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
    alg0 = build_algorithm(algorithm_type; c0, d, grid_n, grid_t_max,
        use_fd_hvp, post_warmup_simplify, grid_bound, grid_curvature_bound,
        linear_area_threshold,
        linear_min_area_gain)
    alg = isnothing(slab_prior) ?
        wrap_sticky(alg0, sticky, model_prior, parameter_prior_vec, can_stick_vec) :
        wrap_dependent_sticky(alg0, sticky, model_prior, slab_prior, can_stick_vec, flow_type, unc_names_vec)

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

    result = _pack_result(chains)
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

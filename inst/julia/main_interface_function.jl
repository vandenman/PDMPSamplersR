module PDMPSamplersRBridge

using PDMPSamplers, LinearAlgebra, BridgeStan, Random, Statistics

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
export r_pdmp_stan, r_pdmp_custom, r_pdmp_custom_subsampled
export write_cmdstan_csv, r_constrain_and_write_csv
export r_pdmp_brms_subsampled, r_pdmp_stan_for_brms
export r_get_param_unc_names, r_stan_param_unc_num_with_header
export r_threading_available

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
        return AdaptiveBoomerang(d; scheme=Symbol(adaptive_scheme))
    elseif flow_type == "PreconditionedZigZag"
        d = length(flow_mean)
        return PreconditionedZigZag(prec, flow_mean)
    elseif flow_type == "PreconditionedBPS"
        d = length(flow_mean)
        return PreconditionedBPS(prec, flow_mean)
    else
        throw(ArgumentError("Unknown flow type: $flow_type"))
    end
end

function build_algorithm(algorithm_type::String; c0::Float64, d::Integer, grid_n::Int, grid_t_max::Float64,
        use_fd_hvp::Bool=false, post_warmup_simplify::Bool=false)
    if algorithm_type == "ThinningStrategy"
        return ThinningStrategy(GlobalBounds(c0 / d, d))
    elseif algorithm_type == "GridThinningStrategy"
        return GridThinningStrategy(; N = grid_n, t_max = grid_t_max,
            use_fd_hvp = use_fd_hvp, post_warmup_simplify = post_warmup_simplify)
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

_haskey(x, key::Symbol) = haskey(x, key) || haskey(x, String(key))
_rget(x, key::Symbol) = haskey(x, key) ? x[key] : x[String(key)]
_as_float_vector(x::Number) = [Float64(x)]
_as_float_vector(x) = Vector{Float64}(x)
_as_bool_vector(x::Bool) = Bool[x]
_as_bool_vector(x::Nothing) = nothing
_as_bool_vector(x) = Bool.(x)
_as_string_vector(x::AbstractString) = String[String(x)]
_as_string_vector(x::Nothing) = String[]
_as_string_vector(x) = String.(x)
_as_float_matrix(x::Number) = reshape([Float64(x)], 1, 1)
_as_float_matrix(x) = Matrix{Float64}(x)

function _slab_type(slab_prior)
    t = _rget(slab_prior, :type)
    return String(t)
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
        pos = findfirst(==(String(coef)), String.(unc_names))
        isnothing(pos) && throw(ArgumentError("slab coef name $(coef) was not found among unconstrained parameter names"))
        !can_stick[pos] && throw(ArgumentError("explicit slab coef must be a subset of can_stick coordinates"))
        return [Int(pos)]
    end
    if coef isa AbstractVector{<:AbstractString}
        isempty(unc_names) && throw(ArgumentError("character slab coef requires unconstrained parameter names"))
        idx = Int[]
        for name in coef
            pos = findfirst(==(String(name)), String.(unc_names))
            isnothing(pos) && throw(ArgumentError("slab coef name $(name) was not found among unconstrained parameter names"))
            push!(idx, pos)
        end
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
        pos = findfirst(==(String(value)), String.(unc_names))
        isnothing(pos) && throw(ArgumentError("$label name $(value) was not found among unconstrained parameter names"))
        return [Int(pos)]
    elseif value isa AbstractVector{<:AbstractString}
        isempty(unc_names) && throw(ArgumentError("character $label requires unconstrained parameter names"))
        idx = Int[]
        for name in value
            pos = findfirst(==(String(name)), String.(unc_names))
            isnothing(pos) && throw(ArgumentError("$label name $(name) was not found among unconstrained parameter names"))
            push!(idx, pos)
        end
        return idx
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
            return BernoulliModelPriorOdds(fill(prob[1], m))
        elseif length(prob) == m
            return BernoulliModelPriorOdds(prob)
        elseif length(prob) == d
            return BernoulliModelPriorOdds(prob[beta_indices])
        else
            throw(DimensionMismatch("Bernoulli model prior length must be 1, beta dimension, or full dimension"))
        end
    elseif _haskey(model_prior, :omega)
        omega = _as_float_vector(_rget(model_prior, :omega))
        length(omega) == m + 1 ||
            throw(DimensionMismatch("exchangeable model-size prior must have length beta dimension + 1"))
        return ExchangeableModelSizePrior(log.(omega); normalize=true)
    elseif _haskey(model_prior, :a) && _haskey(model_prior, :b)
        return BetaBernoulliModelPriorOdds(m, Float64(_rget(model_prior, :a)), Float64(_rget(model_prior, :b)))
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
    odds = build_model_prior_odds(model_prior, beta_indices(provider), d)
    clock = default_aggregate_unstick_clock(provider, odds)
    return AggregateSticky(alg, clock, BitVector(can_stick))
end

function _reject_ungated_slab_prior(slab_prior, caller::AbstractString)
    isnothing(slab_prior) && return nothing
    throw(ArgumentError("Dependent slab_prior is not yet supported for $(caller); target composition must subtract only the slab component or add back nuisance priors"))
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
        prior_model::PDMPModel,
        model_prior,
        slab_prior,
        can_stick,
        unc_names::AbstractVector{<:AbstractString}=String[])
    posterior_model.d == prior_model.d ||
        throw(DimensionMismatch("posterior and prior models must have the same unconstrained dimension"))
    _validate_sampling_slab_prior(slab_prior)
    d = posterior_model.d
    provider = build_slab_provider(slab_prior, unc_names, can_stick, d)
    odds = build_model_prior_odds(model_prior, beta_indices(provider), d)
    target = DependentSlabTarget(d, posterior_model.grad, prior_model.grad, provider, odds)
    return PDMPModel(target)
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

function _compile_model_with_header(path_to_stan_model::String, hpp_path::String)
    !endswith(path_to_stan_model, ".stan") && return path_to_stan_model
    hpp_path_for_make = replace(normpath(hpp_path), "\\" => "/")
    BridgeStan.compile_model(path_to_stan_model;
        stanc_args=["--allow-undefined"],
        make_args=["USER_HEADER=$(hpp_path_for_make)"])
end

function r_stan_param_unc_num_with_header(path_to_stan_model::String,
                                          path_to_stan_data::String,
                                          hpp_path::String)
    lib_path = _compile_model_with_header(path_to_stan_model, hpp_path)
    sm = BridgeStan.StanModel(lib_path, path_to_stan_data; warn=false)
    return Int(BridgeStan.param_unc_num(sm))
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

const StatsValue = Union{Vector{Float64}, Matrix{Float64}}

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

function extract_stats(chains::PDMPChains)
    all = chains.stats
    ct_ess = try
        _ct_ess_matrix(chains)
    catch
        zeros(0, 0)
    end
    vals(name::Symbol) = Float64[_counter_float(s, name) for s in all]
    return Dict{String, StatsValue}(
        "reflections_events"    => vals(:reflections_events),
        "reflections_accepted"  => vals(:reflections_accepted),
        "refreshment_events"    => vals(:refreshment_events),
        "sticky_events"         => vals(:sticky_events),
        "support_boundary_events" => vals(:support_boundary_events),
        "support_boundary_refresh_attempts" => vals(:support_boundary_refresh_attempts),
        "support_boundary_refresh_failures" => vals(:support_boundary_refresh_failures),
        "gradient_calls"        => vals(:∇f_calls),
        "hessian_calls"         => vals(:∇²f_calls),
        "elapsed_time"          => vals(:elapsed_time),
        "grid_builds"           => vals(:grid_builds),
        "grid_shrinks"          => vals(:grid_shrinks),
        "grid_grows"            => vals(:grid_grows),
        "grid_early_stops"      => vals(:grid_early_stops),
        "grid_points_evaluated" => vals(:grid_points_evaluated),
        "grid_points_skipped"   => vals(:grid_points_skipped),
        "grid_N_current"        => vals(:grid_N_current),
        "lazy_fallback_low_tightness" => vals(:lazy_fallback_low_tightness),
        "lazy_fallback_bound_violation" => vals(:lazy_fallback_bound_violation),
        "lazy_proposal_attempts" => vals(:lazy_proposal_attempts),
        "lazy_proposal_rejections" => vals(:lazy_proposal_rejections),
        "grid_resets_from_dynamics_adaptation" => vals(:grid_resets_from_dynamics_adaptation),
        "grid_endpoint_evaluations" => vals(:grid_endpoint_evaluations),
        "grid_cached_endpoint_reuses" => vals(:grid_cached_endpoint_reuses),
        "grid_acceptance_tests" => vals(:grid_acceptance_tests),
        "grid_acceptance_gradient_calls" => vals(:grid_acceptance_gradient_calls),
        "grid_horizon_hits" => vals(:grid_horizon_hits),
        "constant_bound_attempts" => vals(:constant_bound_attempts),
        "constant_bound_accepts" => vals(:constant_bound_accepts),
        "constant_bound_rejections" => vals(:constant_bound_rejections),
        "constant_bound_violations" => vals(:constant_bound_violations),
        "constant_bound_safety_fallbacks" => vals(:constant_bound_safety_fallbacks),
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
        "warmup_elapsed_time" => vals(:warmup_elapsed_time),
        "main_elapsed_time" => vals(:main_elapsed_time),
        "ct_ess"                => ct_ess,
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

    if haskey(kwargs, :slab_prior) && !isnothing(kwargs[:slab_prior])
        _reject_ungated_slab_prior(kwargs[:slab_prior], "r_pdmp_stan")
    end
    model = PDMPModel(path_to_stan_model, path_to_stan_data)
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
        post_warmup_simplify::Bool = true,
        t0::Float64 = 0.0,
        T::Float64 = 10000.0,
        t_warmup::Float64 = 0.0,
        sticky::Bool = false,
        can_stick = nothing,
        model_prior = nothing,
        parameter_prior = nothing,
        slab_prior = nothing,
        prior_model::Union{PDMPModel, Nothing} = nothing,
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
    _reject_ungated_slab_prior(slab_prior, "r_pdmp_stan")

    sampling_model = if isnothing(slab_prior)
        model
    else
        isnothing(prior_model) && throw(ArgumentError("r_pdmp_stan requires prior_model when slab_prior is supplied"))
        _build_dependent_slab_model(model, prior_model, model_prior, slab_prior, can_stick_vec, unc_names_vec)
    end

    prec = _to_precision(flow_cov_mat, d)
    flow = build_flow(flow_type, prec, flow_mean_vec; adaptive_scheme)
    alg0 = build_algorithm(algorithm_type; c0, d, grid_n, grid_t_max, post_warmup_simplify)
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
                         support_boundary_options = sbopts)
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
        t0::Float64 = 0.0,
        T::Float64 = 10000.0,
        t_warmup::Float64 = 0.0,
        hessian = nothing,
        sticky::Bool = false,
        can_stick = nothing,
        model_prior = nothing,
        parameter_prior = nothing,
        slab_prior = nothing,
        prior_grad! = nothing,
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
        isnothing(prior_grad!) && throw(ArgumentError("r_pdmp_custom requires slab-only prior_grad! when slab_prior is supplied"))
        _validate_sampling_slab_prior(slab_prior)
        provider = build_slab_provider(slab_prior, String[], can_stick_vec, d)
        odds = build_model_prior_odds(model_prior, beta_indices(provider), d)
        target = DependentSlabTarget(d, grad!, prior_grad!, provider, odds)
        PDMPModel(target)
    end

    prec = _to_precision(flow_cov_mat, d)
    flow = build_flow(flow_type, prec, flow_mean_vec; adaptive_scheme)
    alg0 = build_algorithm(algorithm_type; c0, d, grid_n, grid_t_max, post_warmup_simplify)
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
# Subsampled gradient bridge (R callbacks)
# ──────────────────────────────────────────────────────────────────────────────

function r_pdmp_custom_subsampled(
        grad_sub_r,
        d::Integer,
        n_obs::Integer,
        subsample_size::Integer,
        x0::AbstractVector{Float64},
        flow_type::String,
        algorithm_type::String,
        flow_mean::AbstractVector{Float64},
        flow_cov::AbstractMatrix{Float64};
        hvp_sub_r = nothing,
        grad_full_r = nothing,
        use_full_gradient_for_reflections::Bool = false,
        c0::Float64 = 1e-2,
        grid_n::Int = 30,
        grid_t_max::Float64 = 2.0,
        post_warmup_simplify::Bool = true,
        t0::Float64 = 0.0,
        T::Float64 = 10000.0,
        t_warmup::Float64 = 0.0,
        show_progress::Bool = true,
        n_chains::Int = 1,
        threaded::Bool = false,
        seed::Union{Integer, Nothing} = nothing,
        adaptive_scheme::String = "diagonal"
    )

    indices = Vector{Int}(undef, subsample_size)
    _perm = Vector{Int}(undef, n_obs)

    function resample!(nsub)
        Random.randperm!(_perm)
        for i in 1:subsample_size
            indices[i] = _perm[i]
        end
    end

    function subsampled_grad!(out, x)
        out .= grad_sub_r(x, indices)
    end

    full_grad! = if isnothing(grad_full_r)
        (out, x) -> error("No full gradient provided but it was requested")
    else
        (out, x) -> (out .= grad_full_r(x))
    end

    grad = SubsampledGradient(
        subsampled_grad!, resample!, (trace) -> nothing,
        FullGradient(full_grad!),
        subsample_size, 0, use_full_gradient_for_reflections, 0.0
    )

    hvp = if isnothing(hvp_sub_r)
        let ε = 1e-5, _buf1 = zeros(d), _buf2 = zeros(d)
            (out, x, v) -> begin
                @. _buf2 = x - ε * v
                subsampled_grad!(_buf1, _buf2)
                @. _buf2 = x + ε * v
                subsampled_grad!(out, _buf2)
                @. out = (out - _buf1) / (2ε)
            end
        end
    else
        (out, x, v) -> (out .= hvp_sub_r(x, v, indices))
    end

    model = PDMPModel(d, grad, hvp)

    prec = isempty(flow_cov) ? Diagonal(ones(d)) : _to_precision(flow_cov, d)
    fmean = isempty(flow_mean) ? zeros(d) : flow_mean
    flow = build_flow(flow_type, prec, fmean; adaptive_scheme)
    alg = build_algorithm(algorithm_type; c0, d, grid_n, grid_t_max, post_warmup_simplify)

    chains = pdmp_sample(x0, flow, model, alg, t0, T, t_warmup;
                         progress = show_progress, n_chains, threaded, seed)
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

# ──────────────────────────────────────────────────────────────────────────────
# BridgeStan subsampling for brms models (external C++ index swapping)
# ──────────────────────────────────────────────────────────────────────────────

mutable struct BridgeStanSubsamplingContext
    sm_full::BridgeStan.StanModel
    sm_sub::BridgeStan.StanModel
    sm_prior::BridgeStan.StanModel
    set_fn::Ptr{Nothing}
    idx_buf::Vector{Int32}
    correction::Vector{Float64}
    N::Int
    m::Int
    perm::Vector{Int}
    grad_buf1::Vector{Float64}
    grad_buf2::Vector{Float64}
    g_full_buf::Vector{Float64}
    hvp_buf::Vector{Float64}
    anchor::Vector{Float64}
    H_full_anchor::Matrix{Float64}
end

function _partial_shuffle!(perm::Vector{Int}, m::Int)
    n = length(perm)
    @inbounds for i in 1:m
        j = rand(i:n)
        perm[i], perm[j] = perm[j], perm[i]
    end
end

function _ccall_set_indices!(ctx::BridgeStanSubsamplingContext)
    @inbounds for i in 1:ctx.m
        ctx.idx_buf[i] = Int32(ctx.perm[i] - 1)
    end
    ccall(ctx.set_fn, Cvoid, (Ptr{Int32}, Int32), ctx.idx_buf, Int32(ctx.m))
end

function _build_bss_model(ctx::BridgeStanSubsamplingContext, n_anchor_updates::Int;
                          resample_dt::Float64=0.0, hvp_mode::String="scaled")
    s = ctx.N / ctx.m
    d = Int(BridgeStan.param_unc_num(ctx.sm_full))
    anchor_set = Ref(false)

    function neg_grad_cv!(out::Vector{Float64}, θ::Vector{Float64})
        BridgeStan.log_density_gradient!(ctx.sm_sub, θ, ctx.grad_buf1)
        BridgeStan.log_density_gradient!(ctx.sm_prior, θ, ctx.grad_buf2)
        @. out = -((1 - s) * ctx.grad_buf2 + s * ctx.grad_buf1 + ctx.correction)
        return out
    end

    function _recompute_correction!()
        BridgeStan.log_density_gradient!(ctx.sm_sub, ctx.anchor, ctx.grad_buf1)
        BridgeStan.log_density_gradient!(ctx.sm_prior, ctx.anchor, ctx.grad_buf2)
        @. ctx.correction = (s - 1) * ctx.grad_buf2 - s * ctx.grad_buf1 + ctx.g_full_buf
    end

    function resample!(nsub)
        _partial_shuffle!(ctx.perm, ctx.m)
        _ccall_set_indices!(ctx)
        if anchor_set[]
            _recompute_correction!()
        end
    end

    function update_anchor!(trace)
        Statistics.mean!(ctx.anchor, trace)
        BridgeStan.log_density_gradient!(ctx.sm_full, ctx.anchor, ctx.g_full_buf)
        _recompute_correction!()
        anchor_set[] = true
    end

    function neg_grad_full!(out::Vector{Float64}, θ::Vector{Float64})
        BridgeStan.log_density_gradient!(ctx.sm_full, θ, out)
        out .= .-out
        return out
    end

    function neg_hvp_cv!(out::Vector{Float64}, θ::Vector{Float64}, v::Vector{Float64})
        BridgeStan.log_density_hessian_vector_product!(ctx.sm_sub, θ, v, out)
        BridgeStan.log_density_hessian_vector_product!(ctx.sm_prior, θ, v, ctx.hvp_buf)
        @. out = -s * out - (1 - s) * ctx.hvp_buf
        return out
    end

    grad = SubsampledGradient(
        neg_grad_cv!, resample!, update_anchor!,
        neg_grad_full!,
        ctx.m, n_anchor_updates, true;
        resample_dt
    )
    hvp = if hvp_mode == "scaled"
        neg_hvp_cv!
    elseif hvp_mode == "none"
        nothing
    else
        error("Unknown hvp_mode: $hvp_mode. Use \"scaled\" or \"none\".")
    end
    return PDMPModel(d, grad, hvp), d
end

function _build_bss_model_hcv(ctx::BridgeStanSubsamplingContext, n_anchor_updates::Int;
                              resample_dt::Float64=0.0, hvp_mode::String="scaled",
                              use_fd_hcv::Bool=false)
    s = ctx.N / ctx.m
    d = Int(BridgeStan.param_unc_num(ctx.sm_full))
    anchor_set = Ref(false)
    hcv = HCVState(d)
    hvp_buf_hcv = zeros(d)
    hess_buf = zeros(d * d)
    hess_grad_buf = zeros(d)
    hess_buf_prior = zeros(d * d)
    grad_sub_anchor = zeros(d)
    fd_buf = zeros(d)
    h_fd = 1e-5

    function hvp_at_anchor_fd!(hvp_out::Vector{Float64}, v::Vector{Float64})
        vnorm = norm(v)
        iszero(vnorm) && (hvp_out .= 0.0; return)
        h_scaled = h_fd * max(1.0, vnorm)
        @. fd_buf = ctx.anchor + (h_scaled / vnorm) * v
        BridgeStan.log_density_gradient!(ctx.sm_sub, fd_buf, hvp_out)
        @. hvp_out = (hvp_out - grad_sub_anchor) * (vnorm / h_scaled)
    end

    function hvp_at_anchor_exact!(hvp_out::Vector{Float64}, v::Vector{Float64})
        BridgeStan.log_density_hessian_vector_product!(ctx.sm_sub, ctx.anchor, v, hvp_out)
    end

    function neg_grad_cv!(out::Vector{Float64}, θ::Vector{Float64})
        BridgeStan.log_density_gradient!(ctx.sm_sub, θ, ctx.grad_buf1)
        BridgeStan.log_density_gradient!(ctx.sm_prior, θ, ctx.grad_buf2)
        @. out = -((1 - s) * ctx.grad_buf2 + s * ctx.grad_buf1 + ctx.correction)

        if hcv.enabled
            if use_fd_hcv
                apply_hcv_correction!(out, hcv, θ, ctx.anchor, hvp_at_anchor_fd!, s, hvp_buf_hcv)
            else
                apply_hcv_correction!(out, hcv, θ, ctx.anchor, hvp_at_anchor_exact!, s, hvp_buf_hcv)
            end
        end
        return out
    end

    function _recompute_correction!()
        BridgeStan.log_density_gradient!(ctx.sm_sub, ctx.anchor, ctx.grad_buf1)
        copyto!(grad_sub_anchor, ctx.grad_buf1)
        BridgeStan.log_density_gradient!(ctx.sm_prior, ctx.anchor, ctx.grad_buf2)
        @. ctx.correction = (s - 1) * ctx.grad_buf2 - s * ctx.grad_buf1 + ctx.g_full_buf
    end

    function resample!(nsub)
        _partial_shuffle!(ctx.perm, ctx.m)
        _ccall_set_indices!(ctx)
        if anchor_set[]
            _recompute_correction!()
        end
    end

    function update_anchor!(trace)
        Statistics.mean!(ctx.anchor, trace)
        BridgeStan.log_density_hessian!(ctx.sm_full, ctx.anchor, ctx.g_full_buf, hess_buf)
        _recompute_correction!()
        anchor_set[] = true

        BridgeStan.log_density_hessian!(ctx.sm_prior, ctx.anchor, hess_grad_buf, hess_buf_prior)
        update_hcv!(hcv, hess_buf, hess_buf_prior, s, d)
    end

    function neg_grad_full!(out::Vector{Float64}, θ::Vector{Float64})
        BridgeStan.log_density_gradient!(ctx.sm_full, θ, out)
        out .= .-out
        return out
    end

    function neg_hvp_cv!(out::Vector{Float64}, θ::Vector{Float64}, v::Vector{Float64})
        BridgeStan.log_density_hessian_vector_product!(ctx.sm_sub, θ, v, out)
        BridgeStan.log_density_hessian_vector_product!(ctx.sm_prior, θ, v, ctx.hvp_buf)
        @. out = -s * out - (1 - s) * ctx.hvp_buf
        return out
    end

    grad = SubsampledGradient(
        neg_grad_cv!, resample!, update_anchor!,
        neg_grad_full!,
        ctx.m, n_anchor_updates, true;
        resample_dt
    )
    hvp = if use_fd_hcv
        nothing
    elseif hvp_mode == "scaled"
        neg_hvp_cv!
    elseif hvp_mode == "none"
        nothing
    else
        error("Unknown hvp_mode: $hvp_mode. Use \"scaled\" or \"none\".")
    end
    return PDMPModel(d, grad, hvp), d
end

function _build_bss_model_bank(ctx::BridgeStanSubsamplingContext, n_anchor_updates::Int;
                               resample_dt::Float64=0.0, hvp_mode::String="scaled",
                               bank_capacity::Int=20)
    s = ctx.N / ctx.m
    d = Int(BridgeStan.param_unc_num(ctx.sm_full))
    bank = AnchorBank(d; capacity=bank_capacity)

    function neg_grad_cv!(out::Vector{Float64}, θ::Vector{Float64})
        BridgeStan.log_density_gradient!(ctx.sm_sub, θ, ctx.grad_buf1)
        BridgeStan.log_density_gradient!(ctx.sm_prior, θ, ctx.grad_buf2)
        @. out = -((1 - s) * ctx.grad_buf2 + s * ctx.grad_buf1 + ctx.correction)
        return out
    end

    function _recompute_correction_from_active!()
        entry = active_entry(bank)
        copyto!(ctx.anchor, entry.position)
        copyto!(ctx.g_full_buf, entry.full_gradient)
        BridgeStan.log_density_gradient!(ctx.sm_sub, ctx.anchor, ctx.grad_buf1)
        BridgeStan.log_density_gradient!(ctx.sm_prior, ctx.anchor, ctx.grad_buf2)
        @. ctx.correction = (s - 1) * ctx.grad_buf2 - s * ctx.grad_buf1 + ctx.g_full_buf
    end

    function resample!(nsub)
        _partial_shuffle!(ctx.perm, ctx.m)
        _ccall_set_indices!(ctx)
        if has_active_anchor(bank)
            _recompute_correction_from_active!()
        end
    end

    function update_anchor!(trace)
        Statistics.mean!(ctx.anchor, trace)
        BridgeStan.log_density_gradient!(ctx.sm_full, ctx.anchor, ctx.g_full_buf)

        BridgeStan.log_density_gradient!(ctx.sm_sub, ctx.anchor, ctx.grad_buf1)
        BridgeStan.log_density_gradient!(ctx.sm_prior, ctx.anchor, ctx.grad_buf2)
        correction = (s - 1) .* ctx.grad_buf2 .- s .* ctx.grad_buf1 .+ ctx.g_full_buf

        add_anchor!(bank;
            position=copy(ctx.anchor),
            full_gradient=copy(ctx.g_full_buf))

        copyto!(ctx.correction, correction)
    end

    function select_fn!(x)
        has_active_anchor(bank) || return
        prev_idx = bank.active_idx
        select_nearest!(bank, x)
        if bank.active_idx != prev_idx
            _recompute_correction_from_active!()
        end
    end

    function neg_grad_full!(out::Vector{Float64}, θ::Vector{Float64})
        BridgeStan.log_density_gradient!(ctx.sm_full, θ, out)
        out .= .-out
        return out
    end

    function neg_hvp_cv!(out::Vector{Float64}, θ::Vector{Float64}, v::Vector{Float64})
        BridgeStan.log_density_hessian_vector_product!(ctx.sm_sub, θ, v, out)
        BridgeStan.log_density_hessian_vector_product!(ctx.sm_prior, θ, v, ctx.hvp_buf)
        @. out = -s * out - (1 - s) * ctx.hvp_buf
        return out
    end

    grad = SubsampledGradient(
        neg_grad_cv!, resample!, update_anchor!,
        neg_grad_full!,
        ctx.m, n_anchor_updates, true;
        resample_dt
    )
    hvp = if hvp_mode == "scaled"
        neg_hvp_cv!
    elseif hvp_mode == "none"
        nothing
    else
        error("Unknown hvp_mode: $hvp_mode. Use \"scaled\" or \"none\".")
    end

    adapter = AnchorBankAdapter(select_fn!, update_anchor!, 0.0, 0.0, true)
    return PDMPModel(d, grad, hvp), d, adapter, bank
end

function _build_bss_model_bank_hcv(ctx::BridgeStanSubsamplingContext, n_anchor_updates::Int;
                                   resample_dt::Float64=0.0, hvp_mode::String="scaled",
                                   bank_capacity::Int=20, use_fd_hcv::Bool=false)
    s = ctx.N / ctx.m
    d = Int(BridgeStan.param_unc_num(ctx.sm_full))
    bank = AnchorBank(d; capacity=bank_capacity, use_hcv=true)
    hvp_buf_hcv = zeros(d)
    hess_buf = zeros(d * d)
    hess_grad_buf = zeros(d)
    hess_buf_prior = zeros(d * d)
    grad_sub_anchor = zeros(d)
    fd_buf = zeros(d)
    h_fd = 1e-5

    function hvp_at_anchor_fd!(hvp_out::Vector{Float64}, v::Vector{Float64})
        vnorm = norm(v)
        iszero(vnorm) && (hvp_out .= 0.0; return)
        h_scaled = h_fd * max(1.0, vnorm)
        @. fd_buf = ctx.anchor + (h_scaled / vnorm) * v
        BridgeStan.log_density_gradient!(ctx.sm_sub, fd_buf, hvp_out)
        @. hvp_out = (hvp_out - grad_sub_anchor) * (vnorm / h_scaled)
    end

    function hvp_at_anchor_exact!(hvp_out::Vector{Float64}, v::Vector{Float64})
        BridgeStan.log_density_hessian_vector_product!(ctx.sm_sub, ctx.anchor, v, hvp_out)
    end

    function neg_grad_cv!(out::Vector{Float64}, θ::Vector{Float64})
        BridgeStan.log_density_gradient!(ctx.sm_sub, θ, ctx.grad_buf1)
        BridgeStan.log_density_gradient!(ctx.sm_prior, θ, ctx.grad_buf2)
        @. out = -((1 - s) * ctx.grad_buf2 + s * ctx.grad_buf1 + ctx.correction)

        if has_active_anchor(bank)
            entry = active_entry(bank)
            hcv = entry.hcv
            if hcv !== nothing && hcv.enabled
                if use_fd_hcv
                    apply_hcv_correction!(out, hcv, θ, ctx.anchor, hvp_at_anchor_fd!, s, hvp_buf_hcv)
                else
                    apply_hcv_correction!(out, hcv, θ, ctx.anchor, hvp_at_anchor_exact!, s, hvp_buf_hcv)
                end
            end
        end
        return out
    end

    function _recompute_correction_from_active!()
        entry = active_entry(bank)
        copyto!(ctx.anchor, entry.position)
        copyto!(ctx.g_full_buf, entry.full_gradient)
        BridgeStan.log_density_gradient!(ctx.sm_sub, ctx.anchor, ctx.grad_buf1)
        copyto!(grad_sub_anchor, ctx.grad_buf1)
        BridgeStan.log_density_gradient!(ctx.sm_prior, ctx.anchor, ctx.grad_buf2)
        @. ctx.correction = (s - 1) * ctx.grad_buf2 - s * ctx.grad_buf1 + ctx.g_full_buf
    end

    function resample!(nsub)
        _partial_shuffle!(ctx.perm, ctx.m)
        _ccall_set_indices!(ctx)
        if has_active_anchor(bank)
            _recompute_correction_from_active!()
        end
    end

    function update_anchor!(trace)
        Statistics.mean!(ctx.anchor, trace)
        BridgeStan.log_density_hessian!(ctx.sm_full, ctx.anchor, ctx.g_full_buf, hess_buf)
        BridgeStan.log_density_gradient!(ctx.sm_sub, ctx.anchor, ctx.grad_buf1)
        copyto!(grad_sub_anchor, ctx.grad_buf1)
        BridgeStan.log_density_gradient!(ctx.sm_prior, ctx.anchor, ctx.grad_buf2)
        correction = (s - 1) .* ctx.grad_buf2 .- s .* ctx.grad_buf1 .+ ctx.g_full_buf

        idx = add_anchor!(bank;
            position=copy(ctx.anchor),
            full_gradient=copy(ctx.g_full_buf))

        copyto!(ctx.correction, correction)

        entry = bank.entries[idx]
        if entry.hcv !== nothing
            BridgeStan.log_density_hessian!(ctx.sm_prior, ctx.anchor, hess_grad_buf, hess_buf_prior)
            update_hcv!(entry.hcv, hess_buf, hess_buf_prior, s, d)
        end
    end

    function select_fn!(x)
        has_active_anchor(bank) || return
        prev_idx = bank.active_idx
        select_nearest!(bank, x)
        if bank.active_idx != prev_idx
            _recompute_correction_from_active!()
        end
    end

    function neg_grad_full!(out::Vector{Float64}, θ::Vector{Float64})
        BridgeStan.log_density_gradient!(ctx.sm_full, θ, out)
        out .= .-out
        return out
    end

    function neg_hvp_cv!(out::Vector{Float64}, θ::Vector{Float64}, v::Vector{Float64})
        BridgeStan.log_density_hessian_vector_product!(ctx.sm_sub, θ, v, out)
        BridgeStan.log_density_hessian_vector_product!(ctx.sm_prior, θ, v, ctx.hvp_buf)
        @. out = -s * out - (1 - s) * ctx.hvp_buf
        return out
    end

    grad = SubsampledGradient(
        neg_grad_cv!, resample!, update_anchor!,
        neg_grad_full!,
        ctx.m, n_anchor_updates, true;
        resample_dt
    )
    hvp = if use_fd_hcv
        nothing
    elseif hvp_mode == "scaled"
        neg_hvp_cv!
    elseif hvp_mode == "none"
        nothing
    else
        error("Unknown hvp_mode: $hvp_mode. Use \"scaled\" or \"none\".")
    end

    adapter = AnchorBankAdapter(select_fn!, update_anchor!, 0.0, 0.0, true)
    return PDMPModel(d, grad, hvp), d, adapter, bank
end

function _resolve_set_fn(lib_path::String)
    lib = Libc.Libdl.dlopen(lib_path, Libc.Libdl.RTLD_NOLOAD | Libc.Libdl.RTLD_GLOBAL)
    Libc.Libdl.dlsym(lib, :pdmp_set_subsample_indices)
end

function _new_bss_context(lib_path_std::String, lib_path_ext::String, set_fn::Ptr{Nothing},
                          data_full_file::String, data_prior_file::String,
                          N::Int, m::Int, d::Int)
    sm_full = BridgeStan.StanModel(lib_path_std, data_full_file; warn=false)
    sm_prior = BridgeStan.StanModel(lib_path_std, data_prior_file; warn=false)
    sm_sub = BridgeStan.StanModel(lib_path_ext, data_full_file; warn=false)
    perm = collect(1:N)
    _partial_shuffle!(perm, m)
    idx_buf = Vector{Int32}(undef, m)
    ctx = BridgeStanSubsamplingContext(
        sm_full, sm_sub, sm_prior, set_fn, idx_buf,
        zeros(d), N, m, perm, zeros(d), zeros(d), zeros(d), zeros(d),
        zeros(d), zeros(d, d)
    )
    _ccall_set_indices!(ctx)
    ctx
end

function r_pdmp_brms_subsampled(
        stan_file::String,
        stan_file_ext::String,
        hpp_path::String,
        data_full_file::String,
        data_prior_file::String,
        N::Integer,
        subsample_size::Integer,
        flow_type::String,
        algorithm_type::String,
        flow_mean,
        flow_cov,
        output_csv::String;
        c0::Float64 = 1e-2,
        grid_n::Int = 30,
        grid_t_max::Float64 = 2.0,
        t0::Float64 = 0.0,
        T::Float64 = 10000.0,
        t_warmup::Float64 = 0.0,
        n_anchor_updates::Int = 10,
        adaptive_scheme::String = "diagonal",
        discretize_dt::Float64 = 0.0,
        show_progress::Bool = true,
        n_chains::Int = 1,
        threaded::Bool = false,
        seed::Union{Integer, Nothing} = nothing,
        compute_lp::Bool = false,
        resample_dt::Float64 = 0.0,
        hvp_mode::String = "scaled",
        use_hcv::Bool = false,
        use_anchor_bank::Bool = false,
        bank_capacity::Int = 20,
        use_fd_hvp::Bool = false,
        post_warmup_simplify::Bool = false,
        use_fd_hcv::Bool = false,
        sticky::Bool = false,
        can_stick::Union{AbstractVector{Bool}, Nothing} = nothing,
        model_prior = nothing,
        parameter_prior::Union{AbstractVector{Float64}, Nothing} = nothing,
        slab_prior = nothing,
        unc_names::AbstractVector{<:AbstractString} = String[]
    )
    _reject_ungated_slab_prior(slab_prior, "r_pdmp_brms_subsampled")
    lib_path_std = _compile_model_with_header(stan_file, hpp_path)
    lib_path_ext = _compile_model_with_header(stan_file_ext, hpp_path)

    N_int = Int(N)
    m = Int(subsample_size)

    sm_full = BridgeStan.StanModel(lib_path_std, data_full_file; warn=false)
    sm_prior = BridgeStan.StanModel(lib_path_std, data_prior_file; warn=false)
    sm_sub = BridgeStan.StanModel(lib_path_ext, data_full_file; warn=false)
    d = Int(BridgeStan.param_unc_num(sm_full))

    set_fn = _resolve_set_fn(lib_path_ext)

    perm = collect(1:N_int)
    _partial_shuffle!(perm, m)
    idx_buf = Vector{Int32}(undef, m)

    ctx = BridgeStanSubsamplingContext(
        sm_full, sm_sub, sm_prior, set_fn, idx_buf,
        zeros(d), N_int, m, perm, zeros(d), zeros(d), zeros(d), zeros(d),
        zeros(d), zeros(d, d)
    )
    _ccall_set_indices!(ctx)

    bank_adapter = nothing
    bank = nothing
    if use_anchor_bank && use_hcv
        model, _, bank_adapter, bank = _build_bss_model_bank_hcv(ctx, n_anchor_updates;
            resample_dt, hvp_mode, bank_capacity, use_fd_hcv)
    elseif use_anchor_bank
        model, _, bank_adapter, bank = _build_bss_model_bank(ctx, n_anchor_updates;
            resample_dt, hvp_mode, bank_capacity)
    elseif use_hcv
        model, _ = _build_bss_model_hcv(ctx, n_anchor_updates; resample_dt, hvp_mode, use_fd_hcv)
    else
        model, _ = _build_bss_model(ctx, n_anchor_updates; resample_dt, hvp_mode)
    end

    fmean = _as_flow_mean(flow_mean, d)
    prec = inv(Symmetric(_as_flow_cov(flow_cov, d)))
    flow = build_flow(flow_type, prec, fmean; adaptive_scheme)
    alg0 = build_algorithm(algorithm_type; c0, d, grid_n, grid_t_max, use_fd_hvp, post_warmup_simplify)
    alg = isnothing(slab_prior) ?
        wrap_sticky(alg0, sticky, model_prior, parameter_prior, can_stick) :
        wrap_dependent_sticky(alg0, sticky, model_prior, slab_prior, can_stick, flow_type, unc_names)

    adapter = if bank_adapter !== nothing
        PDMPSamplers.default_adapter(flow, model.grad, bank_adapter, t_warmup / 10, t_warmup, t0)
    else
        nothing
    end

    function _build_chain_model(ctx_i)
        if use_anchor_bank && use_hcv
            return _build_bss_model_bank_hcv(ctx_i, n_anchor_updates;
                resample_dt, hvp_mode, bank_capacity, use_fd_hcv)
        elseif use_anchor_bank
            return _build_bss_model_bank(ctx_i, n_anchor_updates;
                resample_dt, hvp_mode, bank_capacity)
        elseif use_hcv
            m_i, d_i = _build_bss_model_hcv(ctx_i, n_anchor_updates; resample_dt, hvp_mode, use_fd_hcv)
            return m_i, d_i, nothing, nothing
        else
            m_i, d_i = _build_bss_model(ctx_i, n_anchor_updates; resample_dt, hvp_mode)
            return m_i, d_i, nothing, nothing
        end
    end

    if !use_anchor_bank
        models = typeof(model)[model]
        for i in 2:n_chains
            ctx_i = _new_bss_context(lib_path_std, lib_path_ext, set_fn,
                                     data_full_file, data_prior_file,
                                     N_int, m, d)
            model_i, _ = _build_chain_model(ctx_i)
            push!(models, model_i)
        end
        chains = pdmp_sample(d, flow, models, alg, t0, T, t_warmup;
                             progress=show_progress, threaded, seed)
    elseif n_chains == 1
        chains = pdmp_sample(d, flow, [model], alg, t0, T, t_warmup;
                             progress=show_progress, adapter, seed)
    else
        all_traces = []
        all_stats = []
        for i in 1:n_chains
            if i == 1
                m_i, a_i = model, adapter
            else
                ctx_i = _new_bss_context(lib_path_std, lib_path_ext, set_fn,
                                         data_full_file, data_prior_file,
                                         N_int, m, d)
                m_i, _, ba_i, _ = _build_chain_model(ctx_i)
                a_i = PDMPSamplers.default_adapter(flow, m_i.grad, ba_i, t_warmup / 10, t_warmup, t0)
            end
            ch_i = pdmp_sample(d, deepcopy(flow), [m_i], alg, t0, T, t_warmup;
                               progress=(i == 1 && show_progress), adapter=a_i, seed=isnothing(seed) ? nothing : seed + i - 1)
            push!(all_traces, ch_i.traces[1])
            push!(all_stats, ch_i.stats[1])
        end
        chains = PDMPChains(all_traces, all_stats)
    end

    stats = extract_stats(chains)
    sm_constrain = sm_full
    n_ch = length(chains.traces)
    csv_paths = String[]
    all_draws = [_draws_for_csv(chains, i, discretize_dt) for i in 1:n_ch]
    min_rows = minimum(size(m, 1) for m in all_draws)
    for chain_idx in 1:n_ch
        draws_unc = all_draws[chain_idx][1:min_rows, :]
        csv_path = n_ch == 1 ? output_csv : replace(output_csv, r"\.csv$" => "_chain$(chain_idx).csv")
        r_constrain_and_write_csv(sm_constrain, draws_unc, csv_path; chain_id=chain_idx, compute_lp)
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
        prior_model::Union{PDMPModel, Nothing} = nothing,
        prior_data_file::Union{String, Nothing} = nothing,
        unc_names = String[]
    )
    sm = BridgeStan.StanModel(path_to_stan_model, path_to_stan_data; warn=false)
    _reject_ungated_slab_prior(slab_prior, "r_pdmp_stan_for_brms")
    posterior_model = PDMPModel(sm; hvp = isnothing(slab_prior) && !use_fd_hvp)
    can_stick_vec = _as_bool_vector(can_stick)
    parameter_prior_vec = isnothing(parameter_prior) ? nothing : _as_float_vector(parameter_prior)
    unc_names_vec = _as_string_vector(unc_names)
    model = if isnothing(slab_prior)
        posterior_model
    else
        base_prior_model = if isnothing(prior_model)
            isnothing(prior_data_file) &&
                throw(ArgumentError("r_pdmp_stan_for_brms requires prior_model or prior_data_file when slab_prior is supplied"))
            PDMPModel(BridgeStan.StanModel(path_to_stan_model, prior_data_file; warn=false); hvp=false)
        else
            prior_model
        end
        _build_dependent_slab_model(posterior_model, base_prior_model, model_prior, slab_prior, can_stick_vec, unc_names_vec)
    end
    d = model.d

    fmean = _as_flow_mean(flow_mean, d)
    flow_cov_mat = _as_flow_cov(flow_cov, d)
    prec = _to_precision(flow_cov_mat, d)

    flow = build_flow(flow_type, prec, fmean; adaptive_scheme)
    alg0 = build_algorithm(algorithm_type; c0, d, grid_n, grid_t_max, use_fd_hvp, post_warmup_simplify)
    alg = isnothing(slab_prior) ?
        wrap_sticky(alg0, sticky, model_prior, parameter_prior_vec, can_stick_vec) :
        wrap_dependent_sticky(alg0, sticky, model_prior, slab_prior, can_stick_vec, flow_type, unc_names_vec)

    chains = pdmp_sample(d, flow, model, alg, t0, T, t_warmup;
                         progress = show_progress, n_chains, threaded, seed)

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

using PDMPSamplers, LinearAlgebra, BridgeStan

function build_flow(flow_type::String, prec::AbstractMatrix{Float64}, flow_mean::AbstractVector{Float64})
    if flow_type == "ZigZag"
        return ZigZag(prec, flow_mean)
    elseif flow_type == "BouncyParticle"
        return BouncyParticle(prec, flow_mean)
    elseif flow_type == "Boomerang"
        return Boomerang(prec, flow_mean)
    else
        throw(ArgumentError("Unknown flow type: $flow_type"))
    end
end

function build_algorithm(algorithm_type::String; c0::Float64, d::Integer, grid_n::Int, grid_t_max::Float64)
    if algorithm_type == "ThinningStrategy"
        return ThinningStrategy(GlobalBounds(c0 / d, d))
    elseif algorithm_type == "GridThinningStrategy"
        return GridThinningStrategy(; N = grid_n, t_max = grid_t_max)
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

function discretize_chains(chains::PDMPChains; discretize_dt::Union{Float64, Nothing} = nothing)
    function _discretize_single(trace)
        if isnothing(discretize_dt)
            ts = [event.time for event in trace.events]
            dt = mean(diff(ts))
        else
            dt = discretize_dt
        end
        return Matrix(PDMPDiscretize(trace, dt))
    end
    return [_discretize_single(trace) for trace in chains.traces]
end

function _pack_result(chains::PDMPChains; discretize_dt::Union{Float64, Nothing} = nothing)
    sample_matrices = discretize_chains(chains; discretize_dt)
    stats = extract_stats(chains)
    n_chains = length(chains.traces)

    if n_chains == 1
        return Dict{String,Any}("samples" => sample_matrices[1], "stats" => stats)
    else
        chain_dicts = [Dict{String,Any}("samples" => sample_matrices[i], "stats" => stats)
                       for i in 1:n_chains]
        return chain_dicts
    end
end

function extract_stats(chains::PDMPChains)
    s = chains.stats[1]
    return Dict{String, Any}(
        "reflections_events"   => s.reflections_events,
        "reflections_accepted" => s.reflections_accepted,
        "refreshment_events"   => s.refreshment_events,
        "sticky_events"        => s.sticky_events,
        "gradient_calls"       => s.∇f_calls,
        "hessian_calls"        => s.∇²f_calls,
        "elapsed_time"         => s.elapsed_time
    )
end

function r_pdmp_stan(
        path_to_stan_model::String,
        path_to_stan_data::String,
        x0::AbstractVector{Float64},
        flow_type::String,
        algorithm_type::String,
        flow_mean::AbstractVector{Float64},
        flow_cov::AbstractMatrix{Float64};
        kwargs...
    )

    model = PDMPModel(path_to_stan_model, path_to_stan_data)
    return r_pdmp_stan(model, x0, flow_type, algorithm_type, flow_mean, flow_cov; kwargs...)
end

function r_pdmp_stan(
        model::PDMPModel,
        x0::AbstractVector{Float64},
        flow_type::String,
        algorithm_type::String,
        flow_mean::AbstractVector{Float64},
        flow_cov::AbstractMatrix{Float64};
        c0::Float64 = 1e-2,
        grid_n::Int = 30,
        grid_t_max::Float64 = 2.0,
        t0::Float64 = 0.0,
        T::Float64 = 10000.0,
        t_warmup::Float64 = 0.0,
        discretize_dt::Union{Float64, Nothing} = nothing,
        sticky::Bool = false,
        can_stick::Union{AbstractVector{Bool}, Nothing} = nothing,
        model_prior = nothing,
        parameter_prior::Union{AbstractVector{Float64}, Nothing} = nothing,
        show_progress::Bool = true,
        n_chains::Int = 1,
        threaded::Bool = false
    )

    d = model.d

    prec = inv(Symmetric(flow_cov))
    flow = build_flow(flow_type, prec, flow_mean)
    alg = build_algorithm(algorithm_type; c0, d, grid_n, grid_t_max)
    alg = wrap_sticky(alg, sticky, model_prior, parameter_prior, can_stick)

    chains = pdmp_sample(x0, flow, model, alg, t0, T, t_warmup;
                         progress = show_progress, n_chains = n_chains, threaded = threaded)
    return _pack_result(chains; discretize_dt)
end

function r_pdmp_custom(
        grad!,
        d::Integer,
        x0::AbstractVector{Float64},
        flow_type::String,
        algorithm_type::String,
        flow_mean::AbstractVector{Float64},
        flow_cov::AbstractMatrix{Float64};
        c0::Float64 = 1e-2,
        grid_n::Int = 30,
        grid_t_max::Float64 = 2.0,
        t0::Float64 = 0.0,
        T::Float64 = 10000.0,
        t_warmup::Float64 = 0.0,
        discretize_dt::Union{Float64, Nothing} = nothing,
        hessian = nothing,
        sticky::Bool = false,
        can_stick::Union{AbstractVector{Bool}, Nothing} = nothing,
        model_prior = nothing,
        parameter_prior::Union{AbstractVector{Float64}, Nothing} = nothing,
        show_progress::Bool = true,
        n_chains::Int = 1,
        threaded::Bool = false
    )

    hvp = if !isnothing(hessian)
        (out, x, v) -> begin
            hess = hessian(x)
            mul!(out, hess, v)
        end
    else
        nothing
    end
    model = PDMPModel(d, FullGradient(grad!), hvp)

    prec = inv(Symmetric(flow_cov))
    flow = build_flow(flow_type, prec, flow_mean)
    alg = build_algorithm(algorithm_type; c0, d, grid_n, grid_t_max)
    alg = wrap_sticky(alg, sticky, model_prior, parameter_prior, can_stick)

    chains = pdmp_sample(x0, flow, model, alg, t0, T, t_warmup;
                         progress = show_progress, n_chains = n_chains, threaded = threaded)
    return _pack_result(chains; discretize_dt)
end
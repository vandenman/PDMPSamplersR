# should have already been done by the calling script
using PDMPSamplers, LinearAlgebra
import BridgeStan as BS

function grad_and_hvp_from_stanmodel(path_to_stan_model, path_to_stan_data)

    # Compile Stan model
    sm = BS.StanModel(path_to_stan_model, path_to_stan_data)

    # Create log density function
    grad! = (out, x) -> begin
        BS.log_density_gradient!(sm, x, out, propto = false)
        out .*= -one(eltype(out))
        return out
    end
    # Create log density function
    hvp! = (out, x, v) -> begin
        BS.log_density_hessian_vector_product!(sm, x, v, out, propto = false)
        out .*= -one(eltype(out))
        return out
    end

    d = BS.param_num(sm)

    return grad!, hvp!, d, sm
end

function r_interface_function(
        grad!,
        d::Integer,
        x0::AbstractVector{<:Number},
        flow_type::String,
        algorithm_type::String,
        flow_mean::AbstractVector{<:Number},
        flow_cov::AbstractMatrix{<:Number},
        θ00::Union{AbstractVector{<:Number}, Nothing} = nothing,
        # for ThinningStrategy
        c0::Float64 = 1e-2,
        # for GridThinningStrategy
        grid_n::Int = 30,
        grid_t_max::Float64 = 2.0,
        # end strategies
        t0::AbstractFloat = 0.0,
        T::AbstractFloat = 10000.0,
        discretize_dt::Union{AbstractFloat, Nothing} = nothing,
        hessian! = nothing,
        hvp! = nothing,
        sticky::Bool = false,
        can_stick::Union{<:AbstractVector{<:Bool}, Nothing} = nothing,
        # TODO: this needs a type
        model_prior = nothing,
        parameter_prior::Union{<:AbstractVector{<:Real}, Nothing} = nothing,
        show_progress::Bool = true
    )

    # Create flow object
    prec = inv(Symmetric(flow_cov))
    if flow_type == "ZigZag"
        flow = ZigZag(prec, flow_mean)
    elseif flow_type == "BouncyParticle"
        flow = BouncyParticle(prec, flow_mean)
    elseif flow_type == "Boomerang"
        flow = Boomerang(prec, flow_mean)
    else
        throw(ArgumentError("Unknown flow type: $flow_type"))
    end

    # Create algorithm object
    if algorithm_type == "ThinningStrategy"
        alg = ThinningStrategy(GlobalBounds(c0 / d, d))
    elseif algorithm_type == "RootsPoissonStrategy"
        alg = RootsPoissonStrategy()
    elseif algorithm_type == "GridThinningStrategy"

        if isnothing(hvp!) && isnothing(hessian!)
            throw(ArgumentError("Either Hessian function or Hessian-vector product function must be provided for algorithm GridThinningStrategy"))
        end
        if isnothing(hvp!)
            hess! = (out, x, v) -> begin
                hess = hessian!(x)
                mul!(out, hess, v)
            end
            alg = GridThinningStrategy(; N = grid_n, t_max = grid_t_max, hvp = hess!)
        else
            alg = GridThinningStrategy(; N = grid_n, t_max = grid_t_max, hvp = hvp!)
        end

    else
        throw(ArgumentError("Unknown algorithm type: $algorithm_type"))
    end

    if sticky

        if haskey(model_prior, :prob)

            w = model_prior[:prob] ./ (1 .- model_prior[:prob])
            κ = w .* parameter_prior

        else # betabernoulli

            # A type stable approach would probably be a lot better here...
            a = model_prior[:a]::Float64
            b = model_prior[:b]::Float64
            marginal_pdfs_at_zero = parameter_prior
            κ = (i, x, γ, args...) -> begin

                # critical! do not look at sum(!iszero, x)
                # because that doesn't mean a particle is frozen.
                # instead, use γ which is state.free
                # note that this is only called whenever we're computing the time to stay frozen, so at that point
                # we have γ[i] == false, hence the assert below
                @assert !γ[i] "κ(...) was called but the coordinate to compute the freezing time for (i=$i) is not frozen!?"

                k_free = sum(γ)
                n_tot  = length(x)

                # should never happen because this should only be called when we're computing the time to stay frozen
                k_free == n_tot && throw(error("All parameters are free, cannot compute κ!"))

                # inclusion odds conditional on current sticky state
                prior_incl_odds = (a + k_free) / (b + n_tot - k_free - 1)

                if n_tot - k_free - 1 + b <= 0
                    # or error?
                    # @warn "Inclusion odds are infinite"
                    return Inf
                end

                return prior_incl_odds * marginal_pdfs_at_zero[i]
            end
        end

        alg = Sticky(alg, κ, BitVector(can_stick))

    end

    # Initial conditions
    if isnothing(θ00)
        θ0 = PDMPSamplers.initialize_velocity(flow, d)
    else
        θ0 = θ00
    end
    ξ0 = SkeletonPoint(x0, θ0)

    grad = FullGradient(grad!)

    trace, stats = pdmp_sample(ξ0, flow, grad, alg, t0, T, progress=show_progress)

    # Discretization
    if isnothing(discretize_dt)
        ts = [event.time for event in trace.events]
        dt = mean(diff(ts))
    else
        dt = discretize_dt
    end
    samples = Matrix(PDMPDiscretize(trace, dt))

    return (; samples, trace, stats)
end
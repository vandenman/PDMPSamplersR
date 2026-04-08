module PDMPSamplersRBridge

using PDMPSamplers, LinearAlgebra, BridgeStan, Random, Statistics

export build_flow, build_algorithm, wrap_sticky
export r_discretize, r_mean, r_var, r_std, r_cov, r_cor, r_quantile, r_median, r_cdf, r_ess, r_summary_all
export r_inclusion_probs, extract_stats
export r_pdmp_stan, r_pdmp_custom, r_pdmp_custom_subsampled
export write_cmdstan_csv, r_constrain_and_write_csv
export r_pdmp_brms_subsampled, r_pdmp_stan_for_brms
export r_get_param_unc_names

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

function _to_precision(flow_cov::AbstractMatrix{Float64}, d::Int)
    if isdiag(flow_cov) && all(i -> @inbounds(flow_cov[i, i]) ≈ 1.0, 1:d)
        return Diagonal(ones(d))
    end
    return inv(Symmetric(flow_cov))
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

function extract_stats(chains::PDMPChains)
    all = chains.stats
    ct_ess_per_chain = try
        [minimum(ess(chains; chain=i)) for i in 1:length(chains.traces)]
    catch
        Float64[]
    end
    ct_ess_min = isempty(ct_ess_per_chain) ? NaN : sum(ct_ess_per_chain)
    return Dict{String, Any}(
        "reflections_events"    => sum(s -> s.reflections_events, all),
        "reflections_accepted"  => sum(s -> s.reflections_accepted, all),
        "refreshment_events"    => sum(s -> s.refreshment_events, all),
        "sticky_events"         => sum(s -> s.sticky_events, all),
        "gradient_calls"        => sum(s -> s.∇f_calls, all),
        "hessian_calls"         => sum(s -> s.∇²f_calls, all),
        "elapsed_time"          => sum(s -> s.elapsed_time, all),
        "grid_builds"           => sum(s -> s.grid_builds, all),
        "grid_shrinks"          => sum(s -> s.grid_shrinks, all),
        "grid_grows"            => sum(s -> s.grid_grows, all),
        "grid_early_stops"      => sum(s -> s.grid_early_stops, all),
        "grid_points_evaluated" => sum(s -> s.grid_points_evaluated, all),
        "grid_points_skipped"   => sum(s -> s.grid_points_skipped, all),
        "grid_N_current"        => sum(s -> s.grid_N_current, all),
        "ct_ess_min"            => ct_ess_min,
        "ct_ess_per_chain"      => ct_ess_per_chain,
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
        sticky::Bool = false,
        can_stick::Union{AbstractVector{Bool}, Nothing} = nothing,
        model_prior = nothing,
        parameter_prior::Union{AbstractVector{Float64}, Nothing} = nothing,
        show_progress::Bool = true,
        n_chains::Int = 1,
        threaded::Bool = false,
        adaptive_scheme::String = "diagonal"
    )

    d = model.d

    prec = _to_precision(flow_cov, d)
    flow = build_flow(flow_type, prec, flow_mean; adaptive_scheme)
    alg0 = build_algorithm(algorithm_type; c0, d, grid_n, grid_t_max)
    alg = wrap_sticky(alg0, sticky, model_prior, parameter_prior, can_stick)

    chains = pdmp_sample(x0, flow, model, alg, t0, T, t_warmup;
                         progress = show_progress, n_chains = n_chains, threaded = threaded)
    return _pack_result(chains)
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
        hessian = nothing,
        sticky::Bool = false,
        can_stick::Union{AbstractVector{Bool}, Nothing} = nothing,
        model_prior = nothing,
        parameter_prior::Union{AbstractVector{Float64}, Nothing} = nothing,
        show_progress::Bool = true,
        n_chains::Int = 1,
        threaded::Bool = false,
        adaptive_scheme::String = "diagonal"
    )

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
    model = PDMPModel(d, FullGradient(grad!), hvp)

    prec = _to_precision(flow_cov, d)
    flow = build_flow(flow_type, prec, flow_mean; adaptive_scheme)
    alg0 = build_algorithm(algorithm_type; c0, d, grid_n, grid_t_max)
    alg = wrap_sticky(alg0, sticky, model_prior, parameter_prior, can_stick)

    chains = pdmp_sample(x0, flow, model, alg, t0, T, t_warmup;
                         progress = show_progress, n_chains = n_chains, threaded = threaded)
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
        t0::Float64 = 0.0,
        T::Float64 = 10000.0,
        t_warmup::Float64 = 0.0,
        show_progress::Bool = true,
        n_chains::Int = 1,
        threaded::Bool = false,
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
    alg = build_algorithm(algorithm_type; c0, d, grid_n, grid_t_max)

    chains = pdmp_sample(x0, flow, model, alg, t0, T, t_warmup;
                         progress = show_progress, n_chains, threaded)
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
        flow_mean::AbstractVector,
        flow_cov::AbstractMatrix,
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
        compute_lp::Bool = false,
        resample_dt::Float64 = 0.0,
        hvp_mode::String = "scaled",
        use_hcv::Bool = false,
        use_anchor_bank::Bool = false,
        bank_capacity::Int = 20,
        use_fd_hvp::Bool = false,
        post_warmup_simplify::Bool = false,
        use_fd_hcv::Bool = false
    )
    lib_path_std = if endswith(stan_file, ".stan")
        BridgeStan.compile_model(stan_file)
    else
        stan_file
    end

    lib_path_ext = if endswith(stan_file_ext, ".stan")
        hpp_path_for_make = replace(normpath(hpp_path), "\\" => "/")
        BridgeStan.compile_model(stan_file_ext;
            stanc_args=["--allow-undefined"],
            make_args=["USER_HEADER=$(hpp_path_for_make)"])
    else
        stan_file_ext
    end

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

    fmean = isempty(flow_mean) ? zeros(d) : Vector{Float64}(flow_mean)
    prec = isempty(flow_cov) ? Matrix{Float64}(I, d, d) : inv(Symmetric(Matrix{Float64}(flow_cov)))
    flow = build_flow(flow_type, prec, fmean; adaptive_scheme)
    alg = build_algorithm(algorithm_type; c0, d, grid_n, grid_t_max, use_fd_hvp, post_warmup_simplify)

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
                             progress=show_progress, threaded)
    elseif n_chains == 1
        chains = pdmp_sample(d, flow, [model], alg, t0, T, t_warmup;
                             progress=show_progress, adapter)
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
                               progress=(i == 1 && show_progress), adapter=a_i)
            push!(all_traces, ch_i.traces[1])
            push!(all_stats, ch_i.stats[1])
        end
        chains = PDMPChains(all_traces, all_stats)
    end

    stats = extract_stats(chains)
    sm_constrain = sm_full
    n_ch = length(chains.traces)
    csv_paths = String[]
    all_draws = if discretize_dt > 0
        [Matrix(PDMPDiscretize(chains.traces[i], discretize_dt)) for i in 1:n_ch]
    else
        [adaptive_discretize(chains; chain=i)[1] for i in 1:n_ch]
    end
    min_rows = minimum(size(m, 1) for m in all_draws)
    for chain_idx in 1:n_ch
        draws_unc = all_draws[chain_idx][1:min_rows, :]
        csv_path = n_ch == 1 ? output_csv : replace(output_csv, r"\.csv$" => "_chain$(chain_idx).csv")
        r_constrain_and_write_csv(sm_constrain, draws_unc, csv_path; chain_id=chain_idx, compute_lp)
        push!(csv_paths, csv_path)
    end

    return Dict{String, Any}("csv_paths" => csv_paths, "stats" => stats)
end

function r_pdmp_stan_for_brms(
        path_to_stan_model::String,
        path_to_stan_data::String,
        flow_type::String,
        algorithm_type::String,
        flow_mean::AbstractVector{Float64},
        flow_cov::AbstractMatrix{Float64},
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
        compute_lp::Bool = false,
        use_fd_hvp::Bool = false,
        post_warmup_simplify::Bool = false,
        sticky::Bool = false,
        can_stick::Union{AbstractVector{Bool}, Nothing} = nothing,
        model_prior = nothing,
        parameter_prior::Union{AbstractVector{Float64}, Nothing} = nothing
    )
    sm = BridgeStan.StanModel(path_to_stan_model, path_to_stan_data; warn=false)
    model = PDMPModel(sm; hvp = !use_fd_hvp)
    d = model.d

    fmean = isempty(flow_mean) ? zeros(d) : flow_mean
    prec = isempty(flow_cov) ? Diagonal(ones(d)) : _to_precision(flow_cov, d)

    flow = build_flow(flow_type, prec, fmean; adaptive_scheme)
    alg0 = build_algorithm(algorithm_type; c0, d, grid_n, grid_t_max, use_fd_hvp, post_warmup_simplify)
    alg = wrap_sticky(alg0, sticky, model_prior, parameter_prior, can_stick)

    chains = pdmp_sample(d, flow, model, alg, t0, T, t_warmup;
                         progress = show_progress, n_chains, threaded)

    stats = extract_stats(chains)
    n_ch = length(chains.traces)
    csv_paths = String[]
    all_draws = if discretize_dt > 0
        [Matrix(PDMPDiscretize(chains.traces[i], discretize_dt)) for i in 1:n_ch]
    else
        [adaptive_discretize(chains; chain=i)[1] for i in 1:n_ch]
    end
    min_rows = minimum(size(m, 1) for m in all_draws)
    for chain_idx in 1:n_ch
        draws_unc = all_draws[chain_idx][1:min_rows, :]
        csv_path = n_ch == 1 ? output_csv : replace(output_csv, r"\.csv$" => "_chain$(chain_idx).csv")
        r_constrain_and_write_csv(sm, draws_unc, csv_path; chain_id = chain_idx, compute_lp)
        push!(csv_paths, csv_path)
    end

    result = Dict{String, Any}("csv_paths" => csv_paths, "stats" => stats)
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

end # module PDMPSamplersRBridge

using PDMPSamplers, LinearAlgebra, BridgeStan, Random
using .PDMPSamplersRBridge
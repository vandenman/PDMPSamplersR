using Test, LinearAlgebra, Random, Statistics
using PDMPSamplers, BridgeStan

include(normpath(joinpath(@__DIR__, "..", "..", "inst", "julia",
    "main_interface_function.jl")))
using .PDMPSamplersRBridge

function omrf_person_loglik(q, x)
    P = 3
    thresholds = reshape(q[1:6], 2, 3)
    interactions = zeros(3, 3)
    edges = ((1, 2), (1, 3), (2, 3))
    for (value, (j, k)) in zip(q[7:9], edges)
        interactions[j, k] = interactions[k, j] = value
    end
    result = 0.0
    for j in 1:P
        x[j] > 0 && (result += thresholds[x[j], j])
        field = sum(x[k] * interactions[j, k] for k in 1:P)
        eta = [0.0; [thresholds[u, j] + u * field for u in 1:2]]
        maximum_eta = maximum(eta)
        result -= maximum_eta + log(sum(exp(value - maximum_eta) for value in eta))
    end
    for (j, k) in edges
        result += 2 * interactions[j, k] * x[j] * x[k]
    end
    return result
end

function finite_difference_person_neggrad(q, x)
    result = similar(q)
    step = 1e-6
    plus = copy(q)
    minus = copy(q)
    for j in eachindex(q)
        plus[j] += step
        minus[j] -= step
        result[j] = -(omrf_person_loglik(plus, x) -
            omrf_person_loglik(minus, x)) / (2step)
        plus[j] = minus[j] = q[j]
    end
    return result
end

function omrf_curvature_weights(X, seen; legacy_node_sum=false)
    N, P = size(X)
    threshold_starts = cumsum(vcat(0, seen[1:(end - 1)] .- 1))
    no_thresholds = sum(seen .- 1)
    edges = [(j, k) for j in 1:(P - 1) for k in (j + 1):P]
    weights = zeros(N)
    for n in 1:N
        person_design = zeros(no_thresholds, no_thresholds + length(edges))
        design_row = 1
        for j in 1:P
            q = seen[j] - 1
            iszero(q) && continue
            B = zeros(q, no_thresholds + length(edges))
            for u in 1:q
                B[u, threshold_starts[j] + u] = 1
                for (e, (left, right)) in enumerate(edges)
                    if left == j
                        B[u, no_thresholds + e] = u * X[n, right]
                    elseif right == j
                        B[u, no_thresholds + e] = u * X[n, left]
                    end
                end
            end
            if legacy_node_sum
                weights[n] += 0.5 * opnorm(B)^2
            else
                person_design[design_row:(design_row + q - 1), :] .= B
                design_row += q
            end
        end
        legacy_node_sum || (weights[n] = 0.5 * opnorm(person_design)^2)
    end
    return weights
end

function multivariate_normal_logdensity(x, mean, covariance)
    factor = cholesky(Symmetric(covariance))
    residual = x - mean
    return -0.5 * (length(x) * log(2pi) +
        2sum(log, diag(factor.L)) + dot(residual, factor \ residual))
end

@testset "OMRF person-level selected gradients" begin
    root = normpath(joinpath(@__DIR__, "..", ".."))
    stan = joinpath(root, "inst", "stan", "omrf", "omrf_marked.stan")
    full = joinpath(root, "tests", "stan", "omrf_marked_full.json")
    prior = joinpath(root, "tests", "stan", "omrf_marked_prior.json")
    header = joinpath(root, "inst", "stan", "pdmp_subsample.hpp")
    lib = PDMPSamplersRBridge._compile_model_with_header(stan, header)
    ctx, names = PDMPSamplersRBridge._new_stan_marked_context(lib, full, prior, 1)
    @test length(names) == 9
    q = collect(range(-0.25, 0.3; length=9))
    X = ([0, 1, 2], [2, 0, 1], [1, 2, 0])
    selected = zeros(9)
    prior_gradient = zeros(9)
    direct_sum = zeros(9)
    for n in 1:3
        PDMPSamplersRBridge._selected_gradient!(ctx, selected, q, [n])
        PDMPSamplersRBridge._clear_gradient!(ctx, :prior, prior_gradient, q)
        direct = finite_difference_person_neggrad(q, X[n])
        @test selected - prior_gradient ≈ direct atol=2e-8
        direct_sum .+= direct
    end
    full_gradient = zeros(9)
    PDMPSamplersRBridge._clear_gradient!(ctx, :full, full_gradient, q)
    PDMPSamplersRBridge._clear_gradient!(ctx, :prior, prior_gradient, q)
    @test full_gradient - prior_gradient ≈ direct_sum atol=2e-8
    @test ctx.counts.persons_evaluated == 3
end

@testset "OMRF envelope domination for initially supported flows" begin
    rng = MersenneTwister(481)
    X = [0 1 2; 2 0 1; 1 2 0]
    weights = omrf_curvature_weights(X, [3, 3, 3])
    legacy_weights = omrf_curvature_weights(
        X, [3, 3, 3]; legacy_node_sum=true)
    @test all(weights .<= legacy_weights)
    @test any(weights .< legacy_weights)
    flows = (
        ZigZag(9),
        BouncyParticle(9, 0.2),
        AdaptiveBoomerang(9; λref=0.2),
    )
    for flow in flows, repetition in 1:20
        anchor = 0.35 .* randn(rng, 9)
        position = 0.35 .* randn(rng, 9)
        velocity = 0.7 .* randn(rng, 9)
        person = rand(rng, 1:size(X, 1))
        residual = finite_difference_person_neggrad(position, view(X, person, :)) -
            finite_difference_person_neggrad(anchor, view(X, person, :))
        state = PDMPState(0.0,
            SkeletonPoint(copy(position), copy(velocity)))
        PDMPSamplers.initialize_flow_state!(state, flow)
        envelope = TrajectoryResidualEnvelope(
            reshape(weights, 1, :), anchor)
        scales = zeros(1)
        PDMPSamplers.component_scales!(scales, envelope, state, flow, 0.0)
        actual = flow isa ZigZag ?
            sum(abs.(state.ξ.θ .* residual)) :
            abs(dot(state.ξ.θ, residual))
        @test actual <= weights[person] * scales[1] * (1 + 2e-7) + 2e-8
    end
    anchor = zeros(9)
    position = collect(range(-1.1, 0.9; length=9))
    for flow in flows, person in (1, 3)
        residual = finite_difference_person_neggrad(
            position, view(X, person, :)) -
            finite_difference_person_neggrad(anchor, view(X, person, :))
        velocity = residual .+ collect(range(0.2, 0.6; length=9))
        state = PDMPState(0.0, SkeletonPoint(copy(position), velocity))
        PDMPSamplers.initialize_flow_state!(state, flow)
        envelope = TrajectoryResidualEnvelope(
            reshape(weights, 1, :), anchor)
        scales = zeros(1)
        PDMPSamplers.component_scales!(scales, envelope, state, flow, 0.0)
        actual = flow isa ZigZag ?
            sum(abs.(state.ξ.θ .* residual)) :
            abs(dot(state.ξ.θ, residual))
        @test actual <= weights[person] * scales[1] * (1 + 2e-7) + 2e-8
    end
end

@testset "OMRF marked event law and replicated posterior agreement" begin
    root = normpath(joinpath(@__DIR__, "..", ".."))
    stan = joinpath(root, "inst", "stan", "omrf", "omrf_marked.stan")
    full = joinpath(root, "tests", "stan", "omrf_marked_full.json")
    prior = joinpath(root, "tests", "stan", "omrf_marked_prior.json")
    header = joinpath(root, "inst", "stan", "pdmp_subsample.hpp")
    lib = PDMPSamplersRBridge._compile_model_with_header(stan, header)
    X = [0 1 2; 2 0 1; 1 2 0]
    weights = omrf_curvature_weights(X, [3, 3, 3])
    anchor = zeros(9)
    position = collect(range(-0.18, 0.22; length=9))
    velocity = collect(range(-0.8, 0.9; length=9))
    ctx, _ = PDMPSamplersRBridge._new_stan_marked_context(
        lib, full, prior, 1)
    oracle = PDMPSamplersRBridge.StanMarkedOracle(ctx, anchor)
    deterministic = zeros(9)
    residual = zeros(9)
    estimate_mean = zeros(9)
    oracle(deterministic, position)
    estimates = Vector{Vector{Float64}}()
    for n in axes(X, 1)
        oracle(residual, position, [n], anchor)
        estimate = deterministic + size(X, 1) * residual
        estimate_mean .+= estimate
        push!(estimates, estimate)
    end
    estimate_mean ./= size(X, 1)
    full_gradient = zeros(9)
    PDMPSamplersRBridge._clear_gradient!(ctx, :full, full_gradient, position)
    @test estimate_mean ≈ full_gradient atol=3e-10 rtol=3e-10
    for j in eachindex(velocity)
        forward = mean(max(0.0, velocity[j] * estimate[j]) for estimate in estimates)
        reverse = mean(max(0.0, -velocity[j] * estimate[j]) for estimate in estimates)
        @test forward - reverse ≈ velocity[j] * full_gradient[j] atol=3e-10
    end
    bps_forward = mean(max(0.0, dot(velocity, estimate)) for estimate in estimates)
    bps_reverse = mean(max(0.0, -dot(velocity, estimate)) for estimate in estimates)
    @test bps_forward - bps_reverse ≈ dot(velocity, full_gradient) atol=3e-10

    marked = Dict{String,Any}(
        "n_observations" => 3,
        "subsample_size" => 1,
        "anchor" => anchor,
        "residual_envelope" => Dict{String,Any}(
            "weights" => reshape(weights, 1, :),
            "growth_rates" => [0.0]))
    initial = SkeletonPoint(fill(0.05, 9), vcat(1.0, zeros(8)))
    seeds = (731, 733, 735)
    full_traces = map(seeds) do seed
        model = PDMPModel(BridgeStan.StanModel(lib, full; warn=false))
        trace, _ = pdmp_sample(
            initial, BouncyParticle(9, 0.4), model,
            GridThinningStrategy(N=12, t_max=0.5, lazy=false,
                bound_violation=:throw),
            0.0, 1_200.0; progress=false, seed)
        trace
    end
    marked_runs = map(seeds) do seed
        model, marked_ctx, _ = PDMPSamplersRBridge._build_stan_marked_model(
            lib, full, prior, marked, anchor, nothing, nothing, falses(9))
        trace, stats = pdmp_sample(
            initial, BouncyParticle(9, 0.4), model,
            GridThinningStrategy(N=12, t_max=0.5, lazy=false,
                bound_violation=:throw),
            0.0, 1_200.0; progress=false, seed,
            statistic_counter=PDMPSamplers.DevelStatisticCounter)
        (trace=trace, stats=stats, ctx=marked_ctx)
    end
    marked_traces = getproperty.(marked_runs, :trace)

    function replicated_summary(traces)
        run_means = reduce(hcat, mean.(traces))
        run_variances = reduce(hcat, diag.(cov.(traces)))
        run_ess = reduce(hcat, ess.(traces))
        R = length(traces)
        within_mcse2 = vec(sum(run_variances ./ run_ess; dims=2)) ./ R^2
        between_mcse2 = vec(var(run_means; dims=2, corrected=true)) ./ R
        return (
            mean=vec(sum(run_means; dims=2)) ./ R,
            mcse=sqrt.(max.(within_mcse2, between_mcse2)),
            minimum_ess=minimum(run_ess),
        )
    end
    full_summary = replicated_summary(full_traces)
    marked_summary = replicated_summary(marked_traces)
    combined_mcse = sqrt.(full_summary.mcse .^ 2 .+ marked_summary.mcse .^ 2)
    standardized_difference = abs.(full_summary.mean - marked_summary.mean) ./ combined_mcse
    @test full_summary.minimum_ess > 20
    @test marked_summary.minimum_ess > 20
    @test all(isfinite, standardized_difference)
    @test maximum(standardized_difference) < 4.5

    for run in marked_runs
        @test run.stats.residual_oracle_calls > 0
        @test run.ctx.counts.persons_evaluated ==
            run.ctx.m * run.ctx.counts.selected_gradient
        @test run.ctx.counts.selected_gradient ==
            2 * run.stats.residual_oracle_calls
        @test run.stats.marked_subset_evaluations ==
            run.stats.residual_oracle_calls
        @test run.stats.marked_final_reflections <=
            run.stats.marked_subset_evaluations
        @test run.stats.full_gradient_calls == 0
    end
end

@testset "dependent sticky marked target preserves nuisance prior and model odds" begin
    root = normpath(joinpath(@__DIR__, "..", ".."))
    stan = joinpath(root, "tests", "stan", "marked_gaussian_nuisance.stan")
    full = joinpath(root, "tests", "stan", "marked_gaussian_nuisance_full.json")
    prior = joinpath(root, "tests", "stan", "marked_gaussian_nuisance_prior.json")
    header = joinpath(root, "inst", "stan", "pdmp_subsample.hpp")
    lib = PDMPSamplersRBridge._compile_model_with_header(stan, header)
    marked = Dict{String,Any}(
        "n_observations" => 4,
        "subsample_size" => 2,
        "anchor" => zeros(2),
        "residual_envelope" => Dict{String,Any}(
            "weights" => fill(2.0, 1, 4),
            "growth_rates" => [0.0]))
    slab = Dict{String,Any}(
        "type" => "independent_slab_density",
        "coef" => ["beta"],
        "kappa" => [inv(2sqrt(2pi))])
    model, ctx, _ = PDMPSamplersRBridge._build_stan_marked_model(
        lib, full, prior, marked, zeros(2), slab,
        Dict("prob" => [0.5]), BitVector([true, false]))
    set_active_set!(model, BitVector([false, true]))
    x = [0.0, 0.35]
    corrected = zeros(2)
    model.grad.deterministic_gradient!(corrected, x)
    closure = zeros(2)
    residual = zeros(2)
    subsets = ([1, 2], [1, 3], [1, 4], [2, 3], [2, 4], [3, 4])
    for subset in subsets
        model.grad.residual_oracle(residual, x, subset, model.grad.anchor)
        closure .+= corrected .+ 2 .* residual
    end
    closure ./= length(subsets)
    y = [-0.4, 0.2, 1.1, 0.7]
    likelihood_gradient = sum(x[1] + x[2] - observation for observation in y)
    @test closure ≈ [likelihood_gradient,
        likelihood_gradient + (x[2] - 1) / 0.5^2] atol=2e-10
    @test abs(closure[2] - likelihood_gradient) > 0.1

    sampled = PDMPSamplersRBridge.r_pdmp_stan_marked(
        lib, full, prior, marked, [0.1, 0.8], "ZigZag",
        "GridThinningStrategy", zeros(2), Matrix{Float64}(I, 2, 2);
        T=4_000.0, grid_n=12, grid_t_max=0.5,
        sticky=true, can_stick=BitVector([true, false]),
        model_prior=Dict("prob" => [0.5]), slab_prior=slab,
        show_progress=false, seed=944)
    counters = only(sampled["marked_context_counters"])
    stats = sampled["stats"]
    @test counters["persons_evaluated"] ==
        2 * counters["selected_gradient_calls"]
    @test counters["selected_gradient_calls"] ==
        2 * stats["residual_oracle_calls"][1]
    @test counters["sampling_full_gradient_calls"] == 0
    @test counters["sampling_model_constructions"] == 0
    @test counters["sampling_data_constructions"] == 0
    @test stats["marked_subset_evaluations"][1] ==
        stats["residual_oracle_calls"][1]
    @test stats["marked_final_reflections"][1] <=
        stats["marked_subset_evaluations"][1]

    covariance_inactive = Matrix{Float64}(I, 4, 4) + 0.5^2 * ones(4, 4)
    covariance_active = Matrix{Float64}(I, 4, 4) +
        (2.0^2 + 0.5^2) * ones(4, 4)
    log_inactive = multivariate_normal_logdensity(
        y, fill(1.0, 4), covariance_inactive)
    log_active = multivariate_normal_logdensity(
        y, fill(1.0, 4), covariance_active)
    reference = 1 / (1 + exp(log_inactive - log_active))
    observed = inclusion_probs(sampled["chains"]; chain=1)[1]
    @test observed ≈ reference atol=0.08
end

@testset "binary OMRF model-size frequency matches enumerated reference" begin
    root = normpath(joinpath(@__DIR__, "..", ".."))
    stan = joinpath(root, "tests", "stan", "omrf_binary_interaction.stan")
    full = joinpath(root, "tests", "stan", "omrf_binary_interaction_full.json")
    prior = joinpath(root, "tests", "stan", "omrf_binary_interaction_prior.json")
    header = joinpath(root, "inst", "stan", "pdmp_subsample.hpp")
    lib = PDMPSamplersRBridge._compile_model_with_header(stan, header)
    X = [0 0; 0 0; 0 0; 1 1; 1 1; 1 1; 1 0; 0 1]
    weights = reshape(0.5 .* vec(sum(abs2, X; dims=2)), 1, :)
    marked = Dict{String,Any}(
        "n_observations" => 8,
        "subsample_size" => 2,
        "anchor" => [0.0],
        "residual_envelope" => Dict{String,Any}(
            "weights" => weights, "growth_rates" => [0.0]))
    prior_sd = 1.25
    slab = Dict{String,Any}(
        "type" => "independent_slab_density",
        "coef" => ["interaction"],
        "kappa" => [inv(sqrt(2pi) * prior_sd)])
    sampled = PDMPSamplersRBridge.r_pdmp_stan_marked(
        lib, full, prior, marked, [0.1], "ZigZag",
        "GridThinningStrategy", [0.0], ones(1, 1);
        T=4_000.0, grid_n=12, grid_t_max=0.5,
        sticky=true, can_stick=trues(1),
        model_prior=Dict("prob" => [0.5]), slab_prior=slab,
        show_progress=false, seed=1103)

    person_loglikelihood(beta, x1, x2) = begin
        eta1 = beta * x2
        eta2 = beta * x1
        x1 * eta1 - log1p(exp(eta1)) +
            x2 * eta2 - log1p(exp(eta2))
    end
    loglikelihood(beta) = sum(
        person_loglikelihood(beta, X[n, 1], X[n, 2])
        for n in axes(X, 1))
    grid = range(-8prior_sd, 8prior_sd; length=40_001)
    log_integrand = [loglikelihood(beta) - 0.5 * (beta / prior_sd)^2 -
        log(sqrt(2pi) * prior_sd) for beta in grid]
    shift = maximum(log_integrand)
    step = Base.step(grid)
    active_evidence = exp(shift) * step *
        (sum(exp.(log_integrand .- shift)) -
         0.5 * exp(log_integrand[1] - shift) -
         0.5 * exp(log_integrand[end] - shift))
    inactive_evidence = exp(loglikelihood(0.0))
    reference = active_evidence / (active_evidence + inactive_evidence)
    observed = inclusion_probs(sampled["chains"]; chain=1)[1]
    @test observed ≈ reference atol=0.07
    @test only(sampled["marked_context_counters"])["sampling_full_gradient_calls"] == 0
end

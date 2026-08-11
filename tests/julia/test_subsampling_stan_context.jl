using Test, PDMPSamplers, BridgeStan

include(normpath(joinpath(@__DIR__, "..", "..", "inst", "julia",
    "main_interface_function.jl")))
using .PDMPSamplersRBridge

@testset "persistent custom-Stan selected gradient closure" begin
    root = normpath(joinpath(@__DIR__, "..", ".."))
    stan = joinpath(root, "tests", "stan", "subsampling_gaussian.stan")
    full = joinpath(root, "tests", "stan", "subsampling_gaussian_full.json")
    prior = joinpath(root, "tests", "stan", "subsampling_gaussian_prior.json")
    header = joinpath(root, "inst", "stan", "pdmp_subsample.hpp")
    lib = PDMPSamplersRBridge._compile_model_with_header(stan, header)
    ctx, names = PDMPSamplersRBridge._new_stan_subsampling_context(lib, full, prior, 2)
    @test names == ["beta"]
    oracle = PDMPSamplersRBridge.StanSubsamplingOracle(ctx, [0.0], 4)
    residual = zeros(1)
    oracle(residual, [0.7], [1, 3], [0.0])
    @test residual[1] ≈ 2 * 0.7 atol=1e-10
    deterministic = zeros(1)
    oracle(deterministic, [0.7])
    @test ctx.counts.prior_gradient == 2
    estimate = deterministic[1] + (4 / 2) * residual[1]
    full_gradient = zeros(1)
    PDMPSamplersRBridge._clear_gradient!(ctx, :full, full_gradient, [0.7])
    @test estimate ≈ full_gradient[1] atol=1e-10
    @test ctx.counts.model_constructions == 3
    @test ctx.counts.data_constructions == 3
    @test ctx.counts.selected_gradient == 1
    @test ctx.counts.persons_evaluated == 2
    @test ctx.counts.anchor_cache_gradient == 4
    @test_throws ArgumentError PDMPSamplersRBridge._clear_gradient!(
        ctx, :unknown, zeros(1), [0.0])

    ctx_left, _ = PDMPSamplersRBridge._new_stan_subsampling_context(
        lib, full, prior, 2)
    ctx_right, _ = PDMPSamplersRBridge._new_stan_subsampling_context(
        lib, full, prior, 2)
    left_task = Threads.@spawn begin
        out = zeros(1)
        PDMPSamplersRBridge._selected_gradient!(ctx_left, out, [0.0], [1, 2])
        out[1]
    end
    right_task = Threads.@spawn begin
        out = zeros(1)
        PDMPSamplersRBridge._selected_gradient!(ctx_right, out, [0.0], [3, 4])
        out[1]
    end
    @test fetch(left_task) ≈ 1.0
    @test fetch(right_task) ≈ -3.0

    subsampling = Dict{String,Any}(
        "n_observations" => 4,
        "subsample_size" => 2,
        "anchor" => [0.0],
        "residual_envelope" => Dict{String,Any}(
            "weights" => ones(1, 4),
            "growth_rates" => [0.0]))
    prepared = PDMPSamplersRBridge.prepare_stan_subsampling(
        lib, full, prior, subsampling, 1)
    sampled = PDMPSamplersRBridge.r_pdmp_stan_subsampling(
        prepared, subsampling, [0.1], "ZigZag",
        "GridThinningStrategy", [0.0], ones(1, 1);
        T=0.2, grid_n=4, grid_t_max=0.1,
        show_progress=false, seed=11)
    counters = only(sampled["subsampling_context_counters"])
    @test counters["model_constructions"] == 3
    @test counters["data_constructions"] == 3
    @test counters["initialization_model_constructions"] == 3
    @test counters["initialization_data_constructions"] == 3
    @test counters["sampling_model_constructions"] == 0
    @test counters["sampling_data_constructions"] == 0
    @test counters["sampling_full_gradient_calls"] == 0
    @test counters["full_gradient_calls"] == 1
    @test counters["persons_evaluated"] ==
        2 * counters["selected_gradient_calls"]

    prepared_two = PDMPSamplersRBridge.prepare_stan_subsampling(
        lib, full, prior, subsampling, 2)
    two_chain = PDMPSamplersRBridge.r_pdmp_stan_subsampling(
        prepared_two, subsampling, [0.1], "ZigZag",
        "GridThinningStrategy", [0.0], ones(1, 1);
        T=0.05, grid_n=3, grid_t_max=0.05,
        show_progress=false, n_chains=2,
        threaded=Threads.nthreads() > 1, seed=21)
    two_counts = two_chain["subsampling_context_counters"]
    @test length(two_counts) == 2
    @test all(counter -> counter["model_constructions"] == 3, two_counts)
    @test all(counter -> counter["data_constructions"] == 3, two_counts)
    @test all(counter -> counter["sampling_model_constructions"] == 0, two_counts)
    @test all(counter -> counter["sampling_data_constructions"] == 0, two_counts)
    @test all(counter -> counter["persons_evaluated"] ==
        2 * counter["selected_gradient_calls"], two_counts)

    slab = Dict{String,Any}(
        "type" => "independent_slab_density",
        "coef" => ["beta"],
        "kappa" => [inv(2sqrt(2pi))])
    prepared_sticky = PDMPSamplersRBridge.prepare_stan_subsampling(
        lib, full, prior, subsampling, 1)
    sticky_subsampling = PDMPSamplersRBridge.r_pdmp_stan_subsampling(
        prepared_sticky, subsampling, [0.1], "ZigZag",
        "GridThinningStrategy", [0.0], ones(1, 1);
        T=0.05, grid_n=3, grid_t_max=0.05,
        sticky=true, can_stick=[true],
        model_prior=Dict("prob" => [0.5]), slab_prior=slab,
        show_progress=false, seed=31)
    @test sticky_subsampling["subsampling"] === true
end

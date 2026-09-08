using Test, LinearAlgebra, Random
using PDMPSamplers, BridgeStan

include(normpath(joinpath(@__DIR__, "..", "..", "inst", "julia",
    "main_interface_function.jl")))
using .PDMPSamplersRBridge

function write_json_data(data)
    path = tempname() * ".json"
    # This fixture only contains integer scalars, vectors, and matrices.  Keep
    # the test self-contained in the bridge project's Julia environment rather
    # than adding a JSON package solely for test data serialization.
    matrix_json(x) = "[" * join(("[" * join(string.(x[i, :]), ",") * "]"
        for i in axes(x, 1)), ",") * "]"
    vector_json(x) = "[" * join(string.(x), ",") * "]"
    N = data["N"]
    P = data["P"]
    K = data["K"]
    prior_only = data["prior_only"]
    X = data["X"]
    seen = data["seen"]
    alpha = data["prior_threshold_alpha"]
    beta = data["prior_threshold_beta"]
    slab_scale = data["slab_base_scale"]
    global_sd = data["prior_global_scale_sd"]
    node_sd = data["prior_node_scale_sd"]
    open(path, "w") do stream
        write(stream, "{\"N\":$N,\"P\":$P,")
        write(stream, "\"K\":$K,\"prior_only\":$prior_only,")
        write(stream, "\"X\":$(matrix_json(X)),")
        write(stream, "\"seen\":$(vector_json(seen)),")
        write(stream, "\"prior_threshold_alpha\":$alpha,")
        write(stream, "\"prior_threshold_beta\":$beta,")
        write(stream, "\"slab_base_scale\":$slab_scale,")
        write(stream, "\"prior_global_scale_sd\":$global_sd,")
        write(stream, "\"prior_node_scale_sd\":$node_sd}")
    end
    return path
end

function write_independent_json_data(data)
    path = tempname() * ".json"
    matrix_json(x) = "[" * join(("[" * join(string.(x[i, :]), ",") * "]"
        for i in axes(x, 1)), ",") * "]"
    vector_json(x) = "[" * join(string.(x), ",") * "]"
    N = data["N"]
    P = data["P"]
    K = data["K"]
    X = data["X"]
    seen = data["seen"]
    open(path, "w") do stream
        write(stream, "{\"N\":$N,\"P\":$P,")
        write(stream, "\"K\":$K,\"prior_only\":0,")
        write(stream, "\"X\":$(matrix_json(X)),")
        write(stream, "\"seen\":$(vector_json(seen)),")
        write(stream, "\"prior_interaction_sd\":1.0,")
        write(stream, "\"prior_threshold_alpha\":1.0,")
        write(stream, "\"prior_threshold_beta\":1.0}")
    end
    return path
end

function bridge_negative_gradient!(out, model, x)
    BridgeStan.log_density_gradient!(model, x, out)
    @. out = -out
    return out
end


@testset "production-size neutral OMRF coordinates and seeded sticky velocity" begin
    root = normpath(joinpath(@__DIR__, "..", ".."))
    P = 30
    X = zeros(Int, 2, P)
    seen = fill(4, P)
    fixtures = (
        ("omrf.stan", 525, false),
        ("omrf_node_shared.stan", 556, true),
    )
    for (stan_name, expected_d, dependent) in fixtures
        stan = joinpath(root, "inst", "stan", "omrf", stan_name)
        header = joinpath(root, "inst", "stan", "pdmp_subsample.hpp")
        lib = PDMPSamplersRBridge._compile_model_with_header(stan, header)
        data_path = if dependent
            write_json_data(Dict{String,Any}(
                "N" => 2, "P" => P, "K" => 4, "prior_only" => 0,
                "X" => X, "seen" => seen,
                "prior_threshold_alpha" => 1.0,
                "prior_threshold_beta" => 1.0,
                "slab_base_scale" => 1.0,
                "prior_global_scale_sd" => 1.0,
                "prior_node_scale_sd" => 1.0))
        else
            write_independent_json_data(Dict{String,Any}(
                "N" => 2, "P" => P, "K" => 4, "X" => X,
                "seen" => seen))
        end
        model = BridgeStan.StanModel(lib, data_path; warn=false)
        names = String.(BridgeStan.param_unc_names(model))
        d = Int(BridgeStan.param_unc_num(model))
        @test d == expected_d
        neutral = zeros(d)
        gradient = zeros(d)
        bridge_negative_gradient!(gradient, model, neutral)
        @test all(isfinite, gradient)
        interactions = findall(name -> startswith(name, "interactions_0."), names)
        @test length(interactions) == P * (P - 1) ÷ 2
        @test all(iszero, neutral[interactions])
        if dependent
            logscales = findall(name -> startswith(name, "global_scale") ||
                startswith(name, "node_scale"), names)
            @test length(logscales) == P + 1
            @test all(iszero, neutral[logscales])
        end

        precision = Diagonal(range(0.6, 1.4; length=d))
        flow = PDMPSamplersRBridge.build_flow("PreconditionedZigZag",
            precision, zeros(d))
        function initialized_velocity(seed)
            rng = PDMPSamplers._make_initial_rng(seed, 1)
            velocity = PDMPSamplers.initialize_velocity(rng, flow, d)
            velocity[interactions] .= 0.0
            return velocity
        end
        v1 = initialized_velocity(91001)
        v2 = initialized_velocity(91001)
        @test v1 == v2
        @test all(iszero, v1[interactions])
        @test all(!iszero, v1[setdiff(eachindex(v1), interactions)])
        rm(data_path; force=true)
    end
end

allocated_hvp(hvp, x, direction) = @allocated hvp(x, direction)

@testset "dependent native OMRF backend agrees with BridgeStan" begin
    root = normpath(joinpath(@__DIR__, "..", ".."))
    stan = joinpath(root, "inst", "stan", "omrf", "omrf_node_shared.stan")
    header = joinpath(root, "inst", "stan", "pdmp_subsample.hpp")
    lib = PDMPSamplersRBridge._compile_model_with_header(stan, header)
    max_composed_hvp_absolute_error = 0.0
    max_composed_hvp_relative_error = 0.0

    for (label, X) in (
            :random => rand(MersenneTwister(8101), 0:2, 7, 5),
            :hub => [0 1 0 0 0; 1 0 1 0 0; 2 1 0 1 0;
                     0 0 1 0 2; 1 0 0 2 0; 2 1 0 1 2;
                     0 2 1 0 1])
        N, P = size(X)
        seen = fill(3, P)
        data = Dict{String,Any}(
            "N" => N, "P" => P, "K" => 3, "prior_only" => 0,
            "X" => X, "seen" => seen,
            "prior_threshold_alpha" => 1.7,
            "prior_threshold_beta" => 2.3,
            "slab_base_scale" => 0.7,
            "prior_global_scale_sd" => 1.1,
            "prior_node_scale_sd" => 0.9)
        data_path = write_json_data(data)
        model = BridgeStan.StanModel(lib, data_path; warn=false)
        names = String.(BridgeStan.param_unc_names(model))
        d = Int(BridgeStan.param_unc_num(model))
        spec = Dict{String,Any}(
            "X" => X, "seen" => seen,
            "thresholds" => "thresholds_0",
            "interactions" => "interactions_0",
            "factorization" => "person_node",
            "analytic_prior" => Dict{String,Any}(
                "type" => "node_scale_gaussian",
                "threshold_alpha" => 1.7,
                "threshold_beta" => 2.3,
                "interaction_scale" => 0.7,
                "global_scale" => "global_scale",
                "node_scale" => "node_scale",
                "global_scale_sd" => 1.1,
                "node_scale_sd" => 0.9))
        native = PDMPSamplersRBridge.build_omrf_full_model(spec, names, d)
        context = PDMPSamplersRBridge._omrf_residual_context(
            spec, names, zeros(d), PDMPSamplersRBridge.StanSubsamplingCallCounts())
        norm_weights = zeros(3P, N)
        for node in 1:P, person in 1:N
            neighbour_norm = sqrt(sum(abs2, X[person, other]
                for other in 1:P if other != node))
            norm_weights[node, person] = 1
            norm_weights[P + node, person] = neighbour_norm
            norm_weights[2P + node, person] = neighbour_norm^2
        end
        compact_spec = merge(copy(spec), Dict{String,Any}(
            "type" => "omrf", "backend" => "analytic",
            "bound_type" => "node_local_norm_low_rank",
            "factorization" => "person_node",
            "factor_sampling" => "size_biased",
            "norm_weights" => norm_weights))
        compact_envelope = PDMPSamplersRBridge._build_omrf_residual_envelope(
            compact_spec, zeros(d), context)
        @test size(compact_envelope.weights) == (3P, N)
        @test sum(length, compact_envelope.component_blocks) == 3P
        harmonic_flow = PDMPSamplersRBridge.build_flow("AdaptiveBoomerang",
            Matrix(Diagonal(range(0.7, 1.3; length=d))), zeros(d);
            adaptive_scheme="diagonal", refresh_rate=0.1)
        harmonic_rng = MersenneTwister(8099 + Int(label === :hub))
        harmonic_state = PDMPState(0.0, SkeletonPoint(0.05randn(harmonic_rng, d),
            PDMPSamplers.initialize_velocity(harmonic_rng, harmonic_flow, d)))
        PDMPSamplers.initialize_flow_state!(harmonic_state, harmonic_flow)
        PDMPSamplers.component_cell_scales!(compact_envelope.cell_scales,
            compact_envelope, harmonic_state, harmonic_flow, -0.2, 0.4)
        @test all(isfinite, compact_envelope.cell_scales)
        @test @allocated(PDMPSamplers.component_cell_scales!(
            compact_envelope.cell_scales, compact_envelope, harmonic_state,
            harmonic_flow, -0.2, 0.4)) == 0
        @test native.grad.f isa PDMPSamplersRBridge.OMRFFullGradientProvider
        @test native.hvp isa PDMPSamplers.InplaceHVP
        @test occursin("OMRFFullHVPProvider", string(typeof(native.hvp.f)))

        rng = MersenneTwister(8102 + Int(label === :hub))
        for state_number in 1:4
            x = randn(rng, d) .* 0.25
            direction = randn(rng, d) .* 0.3
            native_gradient = zeros(d)
            bridge_gradient = zeros(d)
            compute_gradient!(native.grad, x, native_gradient)
            bridge_negative_gradient!(bridge_gradient, model, x)
            @test native_gradient ≈ bridge_gradient atol=4e-8 rtol=4e-8

            native_hvp = copy(native.hvp(x, direction))
            step = 2e-6
            plus = zeros(d)
            minus = zeros(d)
            bridge_negative_gradient!(plus, model, x .+ step .* direction)
            bridge_negative_gradient!(minus, model, x .- step .* direction)
            bridge_hvp = (plus - minus) / (2step)
            @test native_hvp ≈ bridge_hvp atol=3e-6 rtol=3e-6

            # The direct deterministic interval certificate is checked against
            # the actual full native gradient along linear and harmonic cells.
            # This is an all-active comparator; inactive-slab correction is
            # exercised separately through the selected subsampling model.
            precision = Diagonal(collect(range(0.55, 1.75; length=d)))
            flows = (
                PDMPSamplersRBridge.build_flow(
                    "PreconditionedZigZag", precision, fill(0.13, d)),
                PDMPSamplersRBridge.build_flow(
                    "PreconditionedBPS", precision, fill(0.13, d);
                    refresh_rate=0.1),
                PDMPSamplersRBridge.build_flow(
                    "AdaptiveBoomerang", precision, fill(0.13, d);
                    adaptive_scheme="diagonal", refresh_rate=0.1),
            )
            for flow in flows
                velocity = PDMPSamplers.initialize_velocity(rng, flow, d)
                state = StickyPDMPState(0.0,
                    SkeletonPoint(copy(x), velocity), trues(d))
                bound = PDMPSamplers.deterministic_rate_cell_bound(
                    native.grad, state, flow, 0.0, 0.35)
                @test isfinite(bound) && bound >= 0
                for t in range(0.0, 0.35; length=65)
                    point = copy(state)
                    PDMPSamplers.move_forward_time!(point, t, flow)
                    native.grad.f(native_gradient, point.ξ.x)
                    actual_flow = flow isa PDMPSamplers.PreconditionedDynamics ?
                        flow.dynamics : flow
                    if actual_flow isa PDMPSamplers.AnyBoomerang
                        native_gradient .-= actual_flow.Γ *
                            (point.ξ.x - actual_flow.μ)
                    end
                    rate = actual_flow isa PDMPSamplers.ZigZag ?
                        sum(max.(0.0, point.ξ.θ .* native_gradient)) :
                        max(0.0, dot(point.ξ.θ, native_gradient))
                    @test rate <= bound
                end
                PDMPSamplers.deterministic_rate_cell_bound(
                    native.grad, state, flow, 0.0, 0.35)
                @test @allocated(PDMPSamplers.deterministic_rate_cell_bound(
                    native.grad, state, flow, 0.0, 0.35)) == 0
            end

            # Zeroing one or several interactions is a partially frozen state,
            # not a changed all-active target. The native/BridgeStan equality
            # must still hold before the active-slab wrapper is applied.
            interaction_indices = context.interaction_indices
            for frozen_count in (1, 2)
                frozen_x = copy(x)
                frozen_x[interaction_indices[1:frozen_count]] .= 0.0
                compute_gradient!(native.grad, frozen_x, native_gradient)
                bridge_negative_gradient!(bridge_gradient, model, frozen_x)
                @test native_gradient ≈ bridge_gradient atol=4e-8 rtol=4e-8
            end
        end

        edge_pairs = [(a, b) for a in 1:(P - 1) for b in (a + 1):P]
        slab_design = zeros(length(edge_pairs), 1 + P)
        slab_design[:, 1] .= 1.0
        for (edge, (left, right)) in enumerate(edge_pairs)
            slab_design[edge, 1 + left] = 0.5
            slab_design[edge, 1 + right] = 0.5
        end
        slab = Dict{String,Any}(
            "type" => "loglinear_gaussian_scale",
            "log_base_scales" => fill(log(0.7), length(context.interaction_indices)),
            "logscale" => ["global_scale", "node_scale"],
            "logscale_design" => Dict{String,Any}(
                "storage" => "dense", "dims" => [length(edge_pairs), 1 + P],
                "x" => vec(slab_design)),
            "coef" => "interactions_0")
        model_prior = Dict{String,Any}("a" => 1.0, "b" => 4.0)
        can_stick = falses(d)
        can_stick[context.interaction_indices] .= true

        # The dependent subsampling deterministic target must expose the same
        # exact fixed-stratum HVP composition.  A finite-difference fallback
        # here previously under-estimated one Zig--Zag deterministic grid
        # cell in the production-size seed-42003 diagnostic.
        prior_data = copy(data)
        prior_data["prior_only"] = 1
        prior_data_path = write_json_data(prior_data)
        sub_ctx, sub_names = PDMPSamplersRBridge._new_stan_subsampling_context(
            lib, data_path, prior_data_path, 1)
        sub_spec = Dict{String,Any}(
            "n_observations" => N * P,
            "subsample_size" => 1,
            "residual_envelope" => compact_spec)
        sub_model, _, _ = PDMPSamplersRBridge._build_stan_subsampling_model(
            sub_ctx, sub_names, sub_spec, zeros(d), slab, can_stick)
        @test sub_model.hvp isa PDMPSamplers.InplaceHVP
        @test occursin("StanSubsamplingSlabDeterministic",
            string(typeof(sub_model.hvp.f)))
        sub_free = trues(d)
        sub_free[context.interaction_indices[1:2]] .= false
        PDMPSamplers.set_active_set!(sub_model, sub_free)
        @test sub_model.hvp.f.oracle.active_free == sub_free
        sub_x = randn(rng, d) .* 0.2
        sub_direction = randn(rng, d) .* 0.3
        sub_plus = zeros(d); sub_minus = zeros(d)
        compute_gradient!(sub_model.grad,
            sub_x .+ 2e-6 .* sub_direction, sub_plus)
        compute_gradient!(sub_model.grad,
            sub_x .- 2e-6 .* sub_direction, sub_minus)
        sub_hvp = copy(sub_model.hvp(sub_x, sub_direction))
        @test sub_hvp ≈ (sub_plus - sub_minus) / (4e-6) atol=4e-6 rtol=4e-6
        sub_model.hvp(sub_x, sub_direction)
        @test allocated_hvp(sub_model.hvp, sub_x, sub_direction) == 0

        # Regression for the fixed-stratum cache defect: before the oracle
        # mask was synchronized, this certificate still used the all-active
        # node-scale slab and could under-bound the wrapped active target.
        sub_flow = PDMPSamplersRBridge.build_flow(
            "AdaptiveBoomerang", Diagonal(collect(range(0.7, 1.6; length=d))),
            fill(0.11, d); adaptive_scheme="diagonal", refresh_rate=0.1)
        sub_state = StickyPDMPState(0.0,
            SkeletonPoint(copy(sub_x),
                PDMPSamplers.initialize_velocity(rng, sub_flow, d)),
            BitVector(sub_free))
        sub_state.ξ.x[.!sub_free] .= 0.0
        sub_state.ξ.θ[.!sub_free] .= 0.0
        sub_bound = PDMPSamplers.deterministic_rate_cell_bound(
            sub_model.hvp.f, sub_state, sub_flow, 0.0, 0.4)
        raw_target = sub_model.grad.deterministic_gradient!
        stats_counter = PDMPSamplers.DevelStatisticCounter()
        counted_model = PDMPSamplers.with_stats(sub_model, stats_counter)
        counted_target = counted_model.grad.deterministic_gradient!
        @test PDMPSamplers.has_direct_deterministic_rate_cell_bound(raw_target)
        @test PDMPSamplers.has_direct_deterministic_rate_cell_bound(counted_target)
        @test PDMPSamplers.deterministic_rate_cell_bound(
            counted_target, sub_state, sub_flow, 0.0, 0.4) == sub_bound

        # Normal grid initialization wraps the gradient and may independently
        # choose HVP, VHV, or finite-difference rate derivatives.  Certificate
        # selection must continue to resolve from the gradient in every case.
        grid_cache = (; ∇ϕx=zeros(d))
        gradient_provider = PDMPSamplers.GradientProvider(
            copy(sub_state.ξ.θ), sub_flow, counted_model.grad, grid_cache)
        derivative_providers = (
            PDMPSamplers.GradHVPProvider(gradient_provider, counted_model.hvp),
            PDMPSamplers.VHVProvider(gradient_provider,
                (x, v, w) -> 0.0, zeros(d)),
            PDMPSamplers.FiniteDiffVHV(gradient_provider,
                zeros(d), zeros(d), zeros(d), stats_counter),
        )
        for provider in derivative_providers
            @test PDMPSamplers.has_direct_deterministic_rate_cell_bound(provider)
            @test PDMPSamplers.deterministic_rate_cell_bound(
                provider, sub_state, sub_flow, 0.0, 0.4) == sub_bound
        end
        sub_gradient = zeros(d)
        for t in range(0.0, 0.4; length=101)
            point = copy(sub_state)
            PDMPSamplers.move_forward_time!(point, t, sub_flow)
            compute_gradient!(sub_model.grad, point.ξ.x, sub_gradient)
            @test max(0.0, dot(point.ξ.θ, sub_gradient)) <=
                sub_bound + 2e-10 * max(1.0, sub_bound)
        end
        rm(prior_data_path; force=true)

        wrapped = PDMPSamplersRBridge._build_dependent_slab_model(
            native, model_prior, slab, can_stick, names)
        @test occursin("_DependentSlabAnalyticHVP", string(typeof(wrapped.hvp.f)))
        @test !occursin("FiniteDiff", string(typeof(wrapped.hvp.f)))
        all_free = trues(d)
        PDMPSamplers.set_active_set!(wrapped, all_free)
        wrapped_gradient = zeros(d)
        native_gradient = zeros(d)
        x = randn(rng, d) .* 0.2
        compute_gradient!(native.grad, x, native_gradient)
        compute_gradient!(wrapped.grad, x, wrapped_gradient)
        @test wrapped_gradient ≈ native_gradient atol=5e-8 rtol=5e-8

        frozen = copy(all_free)
        frozen[context.interaction_indices[1]] = false
        PDMPSamplers.set_active_set!(wrapped, frozen)
        compute_gradient!(wrapped.grad, x, wrapped_gradient)
        @test abs(wrapped_gradient[context.interaction_indices[1]] -
            native_gradient[context.interaction_indices[1]]) > 1e-8
        @test any(abs.(wrapped_gradient[context.node_scale_indices] -
            native_gradient[context.node_scale_indices]) .> 1e-8)

        # The convex-endpoint BPS certificate must also dispatch through the
        # active-slab wrapper.  Its prior interval uses the fixed stratum, not
        # the all-active analytic prior embedded in the native posterior.
        dependent_precision = Diagonal(collect(range(0.55, 1.75; length=d)))
        dependent_bps = PDMPSamplersRBridge.build_flow(
            "PreconditionedBPS", dependent_precision, fill(0.13, d);
            refresh_rate=0.1)
        dependent_state = StickyPDMPState(0.0,
            SkeletonPoint(copy(x), PDMPSamplers.initialize_velocity(
                rng, dependent_bps, d)), BitVector(frozen))
        dependent_state.ξ.x[.!dependent_state.free] .= 0.0
        dependent_state.ξ.θ[.!dependent_state.free] .= 0.0
        PDMPSamplers.set_active_set!(wrapped, dependent_state.free)
        dependent_bound = PDMPSamplers.deterministic_rate_cell_bound(
            wrapped.grad, dependent_state, dependent_bps, 0.0, 0.35)
        for t in range(0.0, 0.35; length=65)
            point = copy(dependent_state)
            PDMPSamplers.move_forward_time!(point, t, dependent_bps)
            compute_gradient!(wrapped.grad, point.ξ.x, wrapped_gradient)
            @test max(0.0, dot(point.ξ.θ, wrapped_gradient)) <=
                dependent_bound
        end
        PDMPSamplers.deterministic_rate_cell_bound(
            wrapped.grad, dependent_state, dependent_bps, 0.0, 0.35)
        @test @allocated(PDMPSamplers.deterministic_rate_cell_bound(
            wrapped.grad, dependent_state, dependent_bps, 0.0, 0.35)) == 0

        PDMPSamplers.set_active_set!(wrapped, all_free)
        direction = randn(rng, d)
        wrapped_hvp = copy(wrapped.hvp(x, direction))
        @test all(isfinite, wrapped_hvp)
        @test wrapped_hvp ≈ native.hvp(x, direction) atol=5e-9 rtol=5e-9

        # The active set is fixed between sticky events, so the compositional
        # HVP must differentiate the corresponding fixed-stratum gradient.
        # Exercise empty, sparse, and dense interaction masks through the
        # actual R-to-Julia construction path.
        for active_count in (0, min(3, length(context.interaction_indices)),
                length(context.interaction_indices))
            free = trues(d)
            free[context.interaction_indices] .= false
            active_count > 0 &&
                (free[context.interaction_indices[1:active_count]] .= true)
            PDMPSamplers.set_active_set!(wrapped, free)
            step = 2e-6
            plus = zeros(d)
            minus = zeros(d)
            compute_gradient!(wrapped.grad, x .+ step .* direction, plus)
            compute_gradient!(wrapped.grad, x .- step .* direction, minus)
            finite_difference = (plus - minus) / (2step)
            analytic = copy(wrapped.hvp(x, direction))
            absolute_error = maximum(abs.(analytic .- finite_difference))
            relative_error = absolute_error / max(maximum(abs.(finite_difference)), eps())
            max_composed_hvp_absolute_error = max(
                max_composed_hvp_absolute_error, absolute_error)
            max_composed_hvp_relative_error = max(
                max_composed_hvp_relative_error, relative_error)
            @test analytic ≈ finite_difference atol=4e-6 rtol=4e-6
        end
        # Prewarming must leave the composed dependent HVP allocation-free.
        PDMPSamplers.set_active_set!(wrapped, all_free)
        wrapped.hvp(x, direction)
        @test allocated_hvp(wrapped.hvp, x, direction) == 0
        rm(data_path; force=true)
    end
    @info "dependent OMRF composed HVP error" max_composed_hvp_absolute_error max_composed_hvp_relative_error
end

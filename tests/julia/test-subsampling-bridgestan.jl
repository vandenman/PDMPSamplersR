using Test, LinearAlgebra, BridgeStan, PDMPSamplers

include(joinpath(@__DIR__, "..", "..", "inst", "julia", "main_interface_function.jl"))
const RB = PDMPSamplersRBridge

function _new_test_context(lib_standard, data_path, prior_path,
        likelihood_design, offsets, N, m, d, anchor;
        observation_multipliers=ones(N), kwargs...)
    active = findall(j -> any(!iszero, @view likelihood_design[:, j]), 1:d)
    designs = [Matrix(@view likelihood_design[:, active])]
    indices = [active]
    return RB._new_subsampling_bss_context(lib_standard, data_path, prior_path,
        :bernoulli, designs, indices, d, reshape(copy(offsets), N, 1),
        zeros(N, 1), Float64[], observation_multipliers, N, m, anchor; kwargs...)
end

@testset "analytic subsampling Bernoulli-logit bridge" begin
    mktempdir() do dir
        standard_path = joinpath(dir, "logistic_standard.stan")
        data_path = joinpath(dir, "data.json")
        prior_path = joinpath(dir, "prior.json")

        write(standard_path, """
        data {
          int<lower=1> N;
          int<lower=1> K;
          matrix[N, K] X;
          vector[N] offsets;
          array[N] int<lower=0, upper=1> Y;
          int<lower=0, upper=1> prior_only;
        }
        parameters {
          vector[K] beta;
          real<lower=0> nuisance_scale;
        }
        model {
          beta ~ normal(0, 2);
          nuisance_scale ~ lognormal(0, 0.7);
          target += -0.01 * nuisance_scale^4;
          if (!prior_only) Y ~ bernoulli_logit(X * beta + offsets);
        }
        """)

        write(data_path, """{
          "N": 4, "K": 2,
          "X": [[1,-1.5],[1,-0.5],[1,0.75],[1,2]],
          "offsets": [0.2,-0.1,0.3,0], "Y": [0,1,0,1],
          "prior_only": 0
        }""")
        write(prior_path, replace(read(data_path, String),
            "\"prior_only\": 0" => "\"prior_only\": 1"))

        X = [1.0 -1.5; 1.0 -0.5; 1.0 0.75; 1.0 2.0]
        offsets = [0.2, -0.1, 0.3, 0.0]
        likelihood_design = hcat(X, zeros(size(X, 1)))
        anchor = zeros(3)
        lib_standard = RB._compile_model(standard_path)
        ctx = _new_test_context(lib_standard,
            data_path, prior_path, likelihood_design, offsets, 4, 2, 3, anchor)
        @test ctx.predictor_operator_norms ≈
              vec(sqrt.(sum(abs2, likelihood_design; dims=2)))
        model = RB._build_subsampling_bss_model(ctx)
        @test ctx.likelihood_hessian === nothing

        theta = [0.35, -0.2, 0.4]
        velocity = [0.7, -1.1, 0.3]
        deterministic = zeros(3)
        residual = zeros(3)
        estimate_mean = zeros(3)
        subsets = ([1, 2], [1, 3], [1, 4], [2, 3], [2, 4], [3, 4])
        weights = vec(sum(abs2, likelihood_design; dims=2)) ./ 4
        for subset in subsets
            model.grad.deterministic_gradient!(deterministic, theta)
            model.grad.residual_oracle(residual, theta, subset, anchor)
            estimate_mean .+= deterministic .+ 2 .* residual

            eta_x = X[subset, :] * theta[1:2] .+ offsets[collect(subset)]
            eta_a = X[subset, :] * anchor[1:2] .+ offsets[collect(subset)]
            direct = zeros(3)
            direct[1:2] .= X[subset, :]' * (RB._subsampling_logistic.(eta_x) .-
                                             RB._subsampling_logistic.(eta_a))
            @test residual ≈ direct atol=1e-12 rtol=1e-12

            residual_bound = 2 * sum(weights[collect(subset)]) *
                norm(velocity) * norm(theta - anchor)
            @test abs(dot(velocity, 2 .* residual)) <= residual_bound * (1 + 1e-12)
        end
        estimate_mean ./= length(subsets)

        full_log_gradient = zeros(3)
        BridgeStan.log_density_gradient!(ctx.sm_full, theta, full_log_gradient)
        @test estimate_mean ≈ -full_log_gradient atol=1e-10 rtol=1e-10
        @test ctx.calls.analytic_residual == length(subsets)
        @test ctx.calls.full_gradient == 1
        model.grad.residual_oracle(residual, theta, subsets[1], anchor)
        @test @allocated(model.grad.residual_oracle(
            residual, theta, subsets[1], anchor)) == 0

        hcv_ctx = _new_test_context(lib_standard,
            data_path, prior_path, likelihood_design, offsets,
            4, 2, 3, anchor; use_hcv=true, hcv_damping=2.5)
        hcv_model = RB._build_subsampling_bss_model(hcv_ctx)
        @test hcv_ctx.likelihood_hessian isa Matrix{Float64}
        @test hcv_model.grad.envelope.weights[1, :] ≈ weights
        @test hcv_model.grad.envelope.weights[2, :] ≈
            vec(sum(abs2, likelihood_design; dims=2)).^(3 / 2) ./ (12sqrt(3.0))
        fill!(estimate_mean, 0.0)
        for subset in subsets
            hcv_model.grad.deterministic_gradient!(deterministic, theta)
            hcv_model.grad.residual_oracle(residual, theta, subset, anchor)
            estimate_mean .+= deterministic .+ 2 .* residual

            direct = zeros(3)
            δ = theta - anchor
            α = hcv_ctx.hcv_damping / (hcv_ctx.hcv_damping + dot(δ, δ))
            for i in subset
                eta_x = dot(X[i, :], theta[1:2]) + offsets[i]
                eta_a = dot(X[i, :], anchor[1:2]) + offsets[i]
                p_a = RB._subsampling_logistic(eta_a)
                coefficient = RB._subsampling_logistic(eta_x) - p_a -
                    α * p_a * (1 - p_a) * (eta_x - eta_a)
                direct[1:2] .+= X[i, :] .* coefficient
            end
            @test residual ≈ direct atol=1e-12 rtol=1e-12

            hcv_velocity = theta - anchor
            state = PDMPState(0.0,
                SkeletonPoint(copy(anchor), copy(hcv_velocity)))
            flow = BouncyParticle(3, 0.0)
            PDMPSamplers.initialize_flow_state!(state, flow)
            point_scales = zeros(2)
            PDMPSamplers.component_scales!(point_scales,
                hcv_model.grad.envelope, state, flow, 1.0)
            subset_bound = 2 * dot(point_scales,
                vec(sum(hcv_model.grad.envelope.weights[:, collect(subset)]; dims=2)))
            @test abs(dot(hcv_velocity, 2 .* residual)) <=
                subset_bound * (1 + 1e-12) + 1e-12
        end
        estimate_mean ./= length(subsets)
        BridgeStan.log_density_gradient!(hcv_ctx.sm_full, theta, full_log_gradient)
        @test estimate_mean ≈ -full_log_gradient atol=1e-10 rtol=1e-10
        @test hcv_ctx.calls.full_gradient == 1
        analytic_hvp = zeros(3)
        hcv_model.grad.deterministic_hvp!(analytic_hvp, theta, velocity)
        epsilon = 1e-5
        gradient_plus = zeros(3); gradient_minus = zeros(3)
        hcv_model.grad.deterministic_gradient!(
            gradient_plus, theta .+ epsilon .* velocity)
        hcv_model.grad.deterministic_gradient!(
            gradient_minus, theta .- epsilon .* velocity)
        finite_difference_hvp = (gradient_plus - gradient_minus) / (2epsilon)
        @test analytic_hvp ≈ finite_difference_hvp atol=1e-7 rtol=1e-7

        bank_ctx = _new_test_context(lib_standard,
            data_path, prior_path, likelihood_design, offsets, 4, 2, 3, anchor)
        bank_model, _, manager = RB._build_subsampling_bss_model(
            bank_ctx; anchor_capacity=2)
        second_anchor = [0.2, -0.15, 0.1]
        prepared = RB._prepare_subsampling_anchor(bank_ctx, second_anchor)
        @test prepared.likelihood_hessian === nothing
        idx = RB._insert_subsampling_anchor!(manager, prepared)
        manager.active_idx = idx
        full_calls_before_activation = bank_ctx.calls.full_gradient
        PDMPSamplers.refresh_anchor!(bank_model.grad, second_anchor)
        @test bank_ctx.calls.full_gradient == full_calls_before_activation
        @test bank_ctx.calls.anchor_activations == 1
        @test bank_model.grad.anchor == second_anchor
        @test bank_model.grad.envelope.component_scales!.anchor == second_anchor

        fill!(estimate_mean, 0.0)
        for subset in subsets
            bank_model.grad.deterministic_gradient!(deterministic, theta)
            bank_model.grad.residual_oracle(
                residual, theta, subset, bank_model.grad.anchor)
            estimate_mean .+= deterministic .+ 2 .* residual
        end
        estimate_mean ./= length(subsets)
        BridgeStan.log_density_gradient!(bank_ctx.sm_full, theta, full_log_gradient)
        @test estimate_mean ≈ -full_log_gradient atol=1e-10 rtol=1e-10

        full_calls_before_selection = bank_ctx.calls.full_gradient
        RB._select_subsampling_anchor!(manager, bank_model.grad, anchor)
        @test bank_model.grad.anchor == anchor
        @test bank_ctx.calls.full_gradient == full_calls_before_selection
        @test bank_ctx.calls.anchor_activations == 2

        # Replace an inactive LRU slot: the active provider must not change.
        third_anchor = [-0.3, 0.25, -0.05]
        third = RB._prepare_subsampling_anchor(bank_ctx, third_anchor)
        manager.entries[1].age = 0
        manager.entries[2].age = 10
        activations_before_inactive_replacement = bank_ctx.calls.anchor_activations
        @test RB._store_prepared_subsampling_anchor!(
            manager, bank_model.grad, third) == 2
        @test bank_model.grad.anchor == anchor
        @test bank_ctx.calls.anchor_activations ==
            activations_before_inactive_replacement
        @test manager.entries[2].state.anchor == third_anchor

        # Select the replacement, then replace that active LRU slot. Active
        # replacement must atomically install the newly prepared state.
        RB._select_subsampling_anchor!(manager, bank_model.grad, third_anchor)
        @test bank_model.grad.anchor == third_anchor
        fourth_anchor = [0.45, 0.1, -0.2]
        fourth = RB._prepare_subsampling_anchor(bank_ctx, fourth_anchor)
        manager.entries[1].age = 0
        manager.entries[2].age = 10
        @test RB._store_prepared_subsampling_anchor!(
            manager, bank_model.grad, fourth) == 2
        @test manager.active_idx == 2
        @test bank_model.grad.anchor == fourth_anchor
        @test manager.entries[2].state.anchor == fourth_anchor
        @test bank_ctx.calls.full_gradient ==
            1 + bank_ctx.calls.anchor_preparations

        for flow in (
                BouncyParticle(3), ZigZag(3),
                PreconditionedBPS(3; scale=[0.7, 1.0, 1.4]),
                PreconditionedZigZag(3; scale=[0.7, 1.0, 1.4]),
                DensePreconditionedBPS(3), DensePreconditionedZigZag(3),
                Boomerang(Matrix{Float64}(I, 3, 3), zeros(3), 0.1),
                AdaptiveBoomerang(3; λref=0.1))
            for analytic_hcv in (false, true)
                flow_ctx = _new_test_context(lib_standard,
                    data_path, prior_path, likelihood_design, offsets,
                    4, 2, 3, anchor; use_hcv=analytic_hcv)
                flow_model = RB._build_subsampling_bss_model(flow_ctx)
                chains = pdmp_sample(anchor, flow, [flow_model],
                    GridThinningStrategy(N=12, t_max=0.5, lazy=false),
                    0.0, 0.5, 0.0; progress=false, seed=72,
                    statistic_counter=PDMPSamplers.DevelStatisticCounter)
                @test length(chains.traces) == 1
                @test flow_ctx.calls.full_gradient == 1
                @test flow_ctx.calls.analytic_residual ==
                      chains.stats[1].residual_oracle_calls
                @test chains.stats[1].subsampling_cell_roof_proposals >=
                      chains.stats[1].subsampling_aggregate_accepts
                @test chains.stats[1].subsampling_aggregate_accepts ==
                      chains.stats[1].subsampling_subset_evaluations ==
                      chains.stats[1].residual_oracle_calls
                @test chains.stats[1].subsampling_final_reflections <=
                      chains.stats[1].subsampling_subset_evaluations
            end
        end

        multiplier = [0.5, 2.0, 0.0, 1.5]
        weighted_ctx = _new_test_context(lib_standard,
            data_path, prior_path, likelihood_design, offsets,
            4, 2, 3, anchor; observation_multipliers=multiplier)
        weighted_model = RB._build_subsampling_bss_model(weighted_ctx)
        weighted_model.grad.residual_oracle(residual, theta, (1, 2), anchor)
        eta_x = X[1:2, :] * theta[1:2] .+ offsets[1:2]
        eta_a = X[1:2, :] * anchor[1:2] .+ offsets[1:2]
        direct = zeros(3)
        direct[1:2] .= X[1:2, :]' * (multiplier[1:2] .*
            (RB._subsampling_logistic.(eta_x) .- RB._subsampling_logistic.(eta_a)))
        @test residual ≈ direct atol=1e-12 rtol=1e-12

        threaded_contexts = [
            _new_test_context(lib_standard,
                data_path, prior_path, likelihood_design, offsets, 4, 2, 3, anchor)
            for _ in 1:2
        ]
        threaded_models = [RB._build_subsampling_bss_model(ctx_i) for ctx_i in threaded_contexts]
        threaded = pdmp_sample(anchor, BouncyParticle(3), threaded_models,
            GridThinningStrategy(N=8, t_max=0.25, lazy=false),
            0.0, 0.25, 0.0; progress=false, threaded=true, seed=91)
        @test length(threaded.traces) == 2
        @test all(ctx_i.calls.full_gradient == 1 for ctx_i in threaded_contexts)
    end
end

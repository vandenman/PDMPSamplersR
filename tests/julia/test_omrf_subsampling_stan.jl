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

function omrf_person_loglik_variable(q, x, seen)
    P = length(seen)
    no_thresholds = sum(seen .- 1)
    threshold_starts = cumsum(vcat(1, seen[1:(end - 1)] .- 1))
    edges = [(left, right) for left in 1:(P - 1) for right in (left + 1):P]
    result = 0.0
    for node in 1:P
        field = 0.0
        for (edge, (left, right)) in enumerate(edges)
            if left == node
                field += x[right] * q[no_thresholds + edge]
            elseif right == node
                field += x[left] * q[no_thresholds + edge]
            end
        end
        start = threshold_starts[node]
        x[node] > 0 && (result += q[start + x[node] - 1])
        logits = [0.0; [q[start + u - 1] + u * field
            for u in 1:(seen[node] - 1)]]
        maximum_logit = maximum(logits)
        result -= maximum_logit +
            log(sum(exp(value - maximum_logit) for value in logits))
    end
    for (edge, (left, right)) in enumerate(edges)
        result += 2q[no_thresholds + edge] * x[left] * x[right]
    end
    return result
end

function finite_difference_person_neggrad_variable(q, x, seen)
    result = similar(q)
    step = 1e-6
    plus = copy(q)
    minus = copy(q)
    for j in eachindex(q)
        plus[j] += step
        minus[j] -= step
        result[j] = -(omrf_person_loglik_variable(plus, x, seen) -
            omrf_person_loglik_variable(minus, x, seen)) / (2step)
        plus[j] = minus[j] = q[j]
    end
    return result
end

@testset "OMRF analytic residual supports variable category counts" begin
    seen = [2, 3, 4]
    X = [0 1 3; 1 2 0]
    no_thresholds = sum(seen .- 1)
    no_edges = 3
    d = no_thresholds + no_edges + 2
    names = vcat(
        ["thresholds_0.$j" for j in 1:no_thresholds],
        ["interactions_0.$j" for j in 1:no_edges],
        ["global_scale", "node_scale.1"])
    envelope_spec = Dict{String,Any}(
        "X" => X,
        "seen" => seen,
        "thresholds" => "thresholds_0",
        "interactions" => "interactions_0")
    counts = PDMPSamplersRBridge.StanSubsamplingCallCounts()
    anchor = collect(range(-0.22, 0.19; length=d))
    position = collect(range(0.31, -0.27; length=d))
    context = PDMPSamplersRBridge._omrf_residual_context(
        envelope_spec, names, anchor, counts)
    residual = zeros(d)
    for person in axes(X, 1)
        PDMPSamplersRBridge.omrf_residual!(
            context, residual, position, [person])
        expected = finite_difference_person_neggrad_variable(
            position[1:(no_thresholds + no_edges)], view(X, person, :), seen) -
            finite_difference_person_neggrad_variable(
                anchor[1:(no_thresholds + no_edges)], view(X, person, :), seen)
        @test residual[1:(no_thresholds + no_edges)] ≈ expected atol=2e-9
        @test residual[(no_thresholds + no_edges + 1):end] == zeros(2)
    end

    person_sum = zeros(d)
    PDMPSamplersRBridge.omrf_residual!(
        context, person_sum, position, collect(axes(X, 1)))
    factor_spec = copy(envelope_spec)
    factor_spec["factorization"] = "person_node"
    factor_context = PDMPSamplersRBridge._omrf_residual_context(
        factor_spec, names, anchor,
        PDMPSamplersRBridge.StanSubsamplingCallCounts())
    factor_sum = zeros(d)
    PDMPSamplersRBridge.omrf_residual!(factor_context, factor_sum, position,
        collect(1:(size(X, 1) * size(X, 2))))
    @test factor_sum ≈ person_sum atol=2e-14 rtol=2e-14
    supports = PDMPSamplersRBridge._omrf_node_supports(factor_context)
    factor_residual = zeros(d)
    for factor in 1:(size(X, 1) * size(X, 2))
        PDMPSamplersRBridge.omrf_residual!(
            factor_context, factor_residual, position, [factor])
        node = cld(factor, size(X, 1))
        outside = setdiff(1:d, supports[node])
        @test factor_residual[outside] == zeros(length(outside))
    end
    @test factor_context.counts.analytic_factors ==
        2 * size(X, 1) * size(X, 2)
    @test factor_context.counts.analytic_node_conditionals ==
        2 * size(X, 1) * size(X, 2)
end

function omrf_curvature_weights(X, seen; legacy_node_sum=false,
        node_local=false)
    N, P = size(X)
    threshold_starts = cumsum(vcat(0, seen[1:(end - 1)] .- 1))
    no_thresholds = sum(seen .- 1)
    edges = [(j, k) for j in 1:(P - 1) for k in (j + 1):P]
    weights = zeros(N)
    node_weights = zeros(P, N)
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
            node_weights[j, n] = 0.5 * opnorm(B)^2
            if legacy_node_sum
                weights[n] += 0.5 * opnorm(B)^2
            else
                person_design[design_row:(design_row + q - 1), :] .= B
                design_row += q
            end
        end
        legacy_node_sum || (weights[n] = 0.5 * opnorm(person_design)^2)
    end
    return node_local ? node_weights : weights
end

function multivariate_normal_logdensity(x, mean, covariance)
    factor = cholesky(Symmetric(covariance))
    residual = x - mean
    return -0.5 * (length(x) * log(2pi) +
        2sum(log, diag(factor.L)) + dot(residual, factor \ residual))
end

@testset "OMRF person-level selected gradients" begin
    root = normpath(joinpath(@__DIR__, "..", ".."))
    stan = joinpath(root, "inst", "stan", "omrf", "omrf.stan")
    full = joinpath(root, "tests", "stan", "omrf_subsampling_full.json")
    prior = joinpath(root, "tests", "stan", "omrf_subsampling_prior.json")
    header = joinpath(root, "inst", "stan", "pdmp_subsample.hpp")
    lib = PDMPSamplersRBridge._compile_model_with_header(stan, header)
    ctx, names = PDMPSamplersRBridge._new_stan_subsampling_context(lib, full, prior, 1)
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

    anchor = collect(range(-0.12, 0.16; length=9))
    position = collect(range(0.21, -0.17; length=9))
    envelope_spec = Dict{String,Any}(
        "type" => "omrf",
        "backend" => "analytic",
        "X" => [0 1 2; 2 0 1; 1 2 0],
        "seen" => [3, 3, 3],
        "thresholds" => "thresholds_0",
        "interactions" => "interactions_0",
        "weights" => reshape(omrf_curvature_weights(
            [0 1 2; 2 0 1; 1 2 0], [3, 3, 3]), 1, :),
        "node_weights" => omrf_curvature_weights(
            [0 1 2; 2 0 1; 1 2 0], [3, 3, 3]; node_local=true),
        "growth_rates" => [0.0])
    analytic = PDMPSamplersRBridge.OMRFAnalyticSubsamplingOracle(
        ctx, envelope_spec, names, anchor)
    analytic_residual = zeros(9)
    selected_position = zeros(9)
    selected_anchor = zeros(9)
    prior_position = zeros(9)
    prior_anchor = zeros(9)
    for subset in ([1], [2], [3])
        analytic(analytic_residual, position, subset, anchor)
        PDMPSamplersRBridge._stan_subset_gradient!(ctx,
            selected_position, position, subset;
            require_configured_size=false)
        PDMPSamplersRBridge._stan_subset_gradient!(ctx,
            selected_anchor, anchor, subset;
            require_configured_size=false)
        PDMPSamplersRBridge._clear_gradient!(ctx, :prior,
            prior_position, position)
        PDMPSamplersRBridge._clear_gradient!(ctx, :prior,
            prior_anchor, anchor)
        @test analytic_residual ≈
            (selected_position - prior_position) -
            (selected_anchor - prior_anchor) atol=3e-10 rtol=3e-10
    end
    analytic(analytic_residual, anchor, [1], anchor)
    @test iszero(norm(analytic_residual))
    @test ctx.counts.analytic_residual == 4
    @test ctx.counts.analytic_persons == 4
    allocation_subset = [2]
    analytic(analytic_residual, position, allocation_subset, anchor)
    @test @allocated(analytic(
        analytic_residual, position, allocation_subset, anchor)) == 0

    analytic_prior_spec = copy(envelope_spec)
    analytic_prior_spec["analytic_prior"] = Dict{String,Any}(
        "type" => "gaussian", "threshold_alpha" => 1.0,
        "threshold_beta" => 1.0, "interaction_scale" => 1.0)
    analytic_prior_ctx, analytic_prior_names =
        PDMPSamplersRBridge._new_stan_subsampling_context(
            lib, full, prior, 1)
    analytic_prior_oracle =
        PDMPSamplersRBridge.OMRFAnalyticSubsamplingOracle(
            analytic_prior_ctx, analytic_prior_spec,
            analytic_prior_names, anchor)
    analytic_prior_gradient = copy(
        PDMPSamplersRBridge._prior_gradient_at_x!(
            analytic_prior_oracle, position))
    stan_prior_gradient = zeros(9)
    PDMPSamplersRBridge._clear_gradient!(
        analytic_prior_ctx, :prior, stan_prior_gradient, position)
    @test analytic_prior_gradient ≈ stan_prior_gradient atol=2e-12 rtol=2e-12
    @test analytic_prior_oracle.analytic_prior !== nothing

    full_native_model = PDMPSamplersRBridge.build_omrf_full_model(
        analytic_prior_spec, analytic_prior_names, 9)
    native_gradient = zeros(9)
    compute_gradient!(full_native_model.grad, position, native_gradient)
    stan_full_gradient = zeros(9)
    PDMPSamplersRBridge._clear_gradient!(
        analytic_prior_ctx, :full, stan_full_gradient, position)
    @test native_gradient ≈ stan_full_gradient atol=3e-10 rtol=3e-10

    direction = collect(range(-0.3, 0.4; length=9))
    native_hvp = copy(full_native_model.hvp(position, direction))
    step = 2e-6
    gradient_plus = zeros(9)
    gradient_minus = zeros(9)
    compute_gradient!(full_native_model.grad,
        position .+ step .* direction, gradient_plus)
    compute_gradient!(full_native_model.grad,
        position .- step .* direction, gradient_minus)
    @test native_hvp ≈
        (gradient_plus - gradient_minus) / (2step) atol=2e-9 rtol=2e-9
    compute_gradient!(full_native_model.grad, position, native_gradient)
    @test @allocated(compute_gradient!(
        full_native_model.grad, position, native_gradient)) == 0
    full_native_model.hvp(position, direction)
    @test @allocated(full_native_model.hvp(position, direction)) == 0

    cauchy_spec = copy(analytic_prior_spec)
    cauchy_spec["analytic_prior"] = Dict{String,Any}(
        "type" => "cauchy", "threshold_alpha" => 1.7,
        "threshold_beta" => 2.3, "interaction_scale" => 0.45)
    cauchy_model = PDMPSamplersRBridge.build_omrf_full_model(
        cauchy_spec, analytic_prior_names, 9)
    cauchy_hvp = copy(cauchy_model.hvp(position, direction))
    compute_gradient!(cauchy_model.grad,
        position .+ step .* direction, gradient_plus)
    compute_gradient!(cauchy_model.grad,
        position .- step .* direction, gradient_minus)
    @test cauchy_hvp ≈
        (gradient_plus - gradient_minus) / (2step) atol=2e-9 rtol=2e-9

    unsupported_spec = copy(analytic_prior_spec)
    unsupported_names = [analytic_prior_names; "nuisance"]
    @test_throws ArgumentError PDMPSamplersRBridge.build_omrf_full_model(
        unsupported_spec, unsupported_names, 10)

    model_ctx, model_names = PDMPSamplersRBridge._new_stan_subsampling_context(
        lib, full, prior, 1)
    subsampling = Dict{String,Any}(
        "n_observations" => 3,
        "subsample_size" => 1,
        "residual_envelope" => envelope_spec)
    model, _, _ = PDMPSamplersRBridge._build_stan_subsampling_model(
        model_ctx, model_names, subsampling, anchor, nothing, falses(9))
    @test model.grad.residual_oracle isa
        PDMPSamplersRBridge.OMRFAnalyticSubsamplingOracle
    @test model.grad.envelope.component_scales! isa
        PDMPSamplersRBridge.OMRFNodeLocalScales
    @test size(model.grad.envelope.weights) == (3, 3)
    @test model_ctx.counts.selected_gradient == 0
    @test model_ctx.counts.anchor_cache_gradient == 0

    factor_spec = copy(envelope_spec)
    factor_spec["factorization"] = "person_node"
    factor_weights = zeros(3, 9)
    for person in 1:3, node in 1:3
        factor_weights[node, (node - 1) * 3 + person] =
            envelope_spec["node_weights"][node, person]
    end
    factor_spec["weights"] = factor_weights
    factor_spec["node_weights"] = factor_weights
    factor_subsampling = Dict{String,Any}(
        "n_observations" => 9,
        "subsample_size" => 3,
        "residual_envelope" => factor_spec)
    factor_ctx, factor_names =
        PDMPSamplersRBridge._new_stan_subsampling_context(
            lib, full, prior, 3)
    factor_model, _, _ = PDMPSamplersRBridge._build_stan_subsampling_model(
        factor_ctx, factor_names, factor_subsampling, anchor, nothing, falses(9))
    @test factor_model.grad.residual_oracle.residual_context.factorization ===
        :person_node
    @test size(factor_model.grad.envelope.weights) == (3, 9)
    factor_deterministic = zeros(9)
    factor_residual_sum = zeros(9)
    factor_residual = zeros(9)
    factor_model.grad.residual_oracle(factor_deterministic, position)
    for factor in 1:9
        factor_model.grad.residual_oracle(
            factor_residual, position, [factor], anchor)
        factor_residual_sum .+= factor_residual
    end
    exact_factor_full = zeros(9)
    PDMPSamplersRBridge._clear_gradient!(
        factor_ctx, :full, exact_factor_full, position)
    @test factor_deterministic + factor_residual_sum ≈
        exact_factor_full atol=3e-10 rtol=3e-10
    factor_model.grad.residual_oracle(
        factor_residual, position, [1], anchor)
    sparse_gradient = copy(factor_deterministic)
    expected_sparse_gradient = factor_deterministic + 9 .* factor_residual
    sparse_state = PDMPState(0.0,
        SkeletonPoint(copy(position), fill(1.0, 9)))
    deterministic_rate = PDMPSamplers.λ(
        sparse_state, factor_deterministic, ZigZag(9))
    sparse_rate = PDMPSamplers.subsampling_candidate_rate!(
        factor_model.grad.residual_oracle, sparse_state, sparse_gradient,
        factor_residual, 9.0, ZigZag(9), deterministic_rate, [1])
    @test sparse_gradient ≈ expected_sparse_gradient atol=2e-14 rtol=2e-14
    @test sparse_rate ≈ PDMPSamplers.λ(
        sparse_state, expected_sparse_gradient, ZigZag(9)) atol=2e-14 rtol=2e-14

    hcv_spec = copy(envelope_spec)
    hcv_spec["use_hcv"] = true
    hcv_spec["hcv_damping"] = 9.0
    hcv_spec["hcv_remainder_weights"] =
        0.75 .* (2 .* hcv_spec["node_weights"]) .^ 1.5
    hcv_ctx, hcv_names = PDMPSamplersRBridge._new_stan_subsampling_context(
        lib, full, prior, 1)
    hcv_oracle = PDMPSamplersRBridge.OMRFAnalyticSubsamplingOracle(
        hcv_ctx, hcv_spec, hcv_names, anchor)
    hcv_envelope = PDMPSamplersRBridge._build_stan_residual_envelope(
        hcv_spec, anchor, hcv_oracle)
    @test hcv_envelope.component_scales! isa
        PDMPSamplersRBridge.OMRFDampedHCVNodeLocalScales
    @test size(hcv_envelope.weights) == (6, 3)
    @test hcv_oracle.residual_context.likelihood_hessian !== nothing
    hcv_deterministic = zeros(9)
    hcv_residual_sum = zeros(9)
    hcv_residual = zeros(9)
    hcv_oracle(hcv_deterministic, position)
    for person in 1:3
        hcv_oracle(hcv_residual, position, [person], anchor)
        hcv_residual_sum .+= hcv_residual
    end
    exact_full = zeros(9)
    PDMPSamplersRBridge._clear_gradient!(hcv_ctx, :full, exact_full, position)
    @test hcv_deterministic + hcv_residual_sum ≈ exact_full atol=3e-10 rtol=3e-10
    hcv_oracle(hcv_residual, anchor, [1], anchor)
    @test iszero(norm(hcv_residual))

    staged_spec = copy(hcv_spec)
    staged_spec["hcv_after_warmup"] = true
    staged_subsampling = Dict{String,Any}(
        "n_observations" => 3,
        "subsample_size" => 1,
        "residual_envelope" => staged_spec)
    staged_ctx, staged_names =
        PDMPSamplersRBridge._new_stan_subsampling_context(
            lib, full, prior, 1)
    staged_model, _, _, staged_adapter, staged_manager =
        PDMPSamplersRBridge._build_stan_subsampling_model(
            staged_ctx, staged_names, staged_subsampling, anchor, nothing,
            falses(9); anchor_capacity=2)
    @test !staged_model.grad.residual_oracle.residual_context.hcv_active
    staged_anchor = position .+ 0.001
    staged_prepared = PDMPSamplersRBridge._prepare_omrf_anchor(
        staged_manager, staged_anchor)
    PDMPSamplersRBridge._insert_omrf_anchor!(
        staged_manager, staged_prepared)
    staged_scales = zeros(6)
    staged_state = PDMPState(0.0,
        SkeletonPoint(copy(position), fill(1.0, 9)))
    staged_model.grad.envelope.component_scales!(
        staged_scales, staged_state, BouncyParticle(9, 0.0), 0.0)
    @test all(iszero, staged_scales[4:6])
    staged_deterministic = zeros(9)
    staged_residual_sum = zeros(9)
    staged_residual = zeros(9)
    staged_model.grad.residual_oracle(staged_deterministic, position)
    for person in 1:3
        staged_model.grad.residual_oracle(
            staged_residual, position, [person], anchor)
        staged_residual_sum .+= staged_residual
    end
    @test staged_deterministic + staged_residual_sum ≈
        exact_full atol=3e-10 rtol=3e-10
    @test PDMPSamplers.finish_warmup!(staged_adapter, staged_state,
        BouncyParticle(9, 0.0), staged_model.grad, nothing, nothing)
    @test staged_model.grad.residual_oracle.residual_context.hcv_active
    @test staged_manager.hcv_active
    @test staged_model.grad.residual_oracle.anchor == staged_anchor
    staged_model.grad.envelope.component_scales!(
        staged_scales, staged_state, BouncyParticle(9, 0.0), 0.0)
    @test any(!iszero, staged_scales[4:6])
    fill!(staged_residual_sum, 0.0)
    staged_model.grad.residual_oracle(staged_deterministic, position)
    for person in 1:3
        staged_model.grad.residual_oracle(
            staged_residual, position, [person], anchor)
        staged_residual_sum .+= staged_residual
    end
    @test staged_deterministic + staged_residual_sum ≈
        exact_full atol=3e-10 rtol=3e-10

    bank_ctx, bank_names = PDMPSamplersRBridge._new_stan_subsampling_context(
        lib, full, prior, 1)
    managed_model, _, _, _, manager =
        PDMPSamplersRBridge._build_stan_subsampling_model(
            bank_ctx, bank_names, subsampling, anchor, nothing, falses(9);
            anchor_capacity=2)
    new_anchor = anchor .+ collect(range(-0.04, 0.05; length=9))
    prepared_anchor = PDMPSamplersRBridge._prepare_omrf_anchor(
        manager, new_anchor)
    new_idx = PDMPSamplersRBridge._insert_omrf_anchor!(
        manager, prepared_anchor)
    manager.active_idx = new_idx
    full_calls_before_activation = bank_ctx.counts.full_gradient
    prior_calls_before_activation = bank_ctx.counts.prior_gradient
    PDMPSamplers.refresh_anchor!(managed_model.grad, new_anchor)
    @test manager.oracle.anchor == new_anchor
    @test manager.activations == 1
    @test manager.preparations == 1
    @test manager.preparation_seconds > 0
    @test length(manager.entries) == 2
    @test bank_ctx.counts.full_gradient == full_calls_before_activation
    @test bank_ctx.counts.prior_gradient == prior_calls_before_activation
    @test managed_model.grad.envelope.component_scales!.anchor ===
        manager.oracle.anchor
    manager.oracle(analytic_residual, new_anchor, [1], new_anchor)
    @test iszero(norm(analytic_residual))
    refresh_position = new_anchor .+ 0.1
    manager.main_refresh_distance = 0.01
    preparations_before_main = manager.preparations
    PDMPSamplersRBridge._select_omrf_anchor!(
        manager, managed_model.grad, refresh_position; phase=:main)
    @test manager.preparations == preparations_before_main + 1
    @test manager.main_refreshes == 1
    @test manager.main_activations == 1
    @test manager.oracle.anchor == refresh_position
    @test managed_model.grad.anchor == refresh_position

    node_weights = envelope_spec["node_weights"]
    pattern_spec = copy(envelope_spec)
    pattern_spec["bound_type"] = "pattern_local_range"
    pattern_spec["structural_groups"] =
        permutedims(reshape(collect(1:9), 3, 3))
    pattern_spec["n_structural_components"] = 9
    pattern_spec["component_nodes"] = repeat(collect(1:3); inner=3)
    pattern_spec["component_covariates"] = repeat(
        [0 1 2; 2 0 1; 1 2 0], 3, 1)
    pattern_envelope = PDMPSamplersRBridge._build_omrf_residual_envelope(
        pattern_spec, anchor, analytic.residual_context)
    @test pattern_envelope isa PDMPSamplers.GroupedResidualEnvelope
    covariate_spec = copy(analytic_prior_spec)
    covariate_spec["bound_type"] = "covariate_local_expansion"
    covariate_weights = zeros(27, 3)
    destination = 1
    for node in 1:3
        neighbours = analytic.residual_context.incident_neighbours[node]
        for person in 1:3
            covariate_weights[destination, person] = 1.0
            for block in 0:1, k in eachindex(neighbours)
                covariate_weights[destination + block * length(neighbours) + k,
                    person] = analytic.residual_context.X[
                        person, neighbours[k]]
            end
            cross_start = destination + 2length(neighbours) + 1
            for k in eachindex(neighbours), l in eachindex(neighbours)
                covariate_weights[cross_start +
                    (k - 1) * length(neighbours) + l - 1, person] =
                    analytic.residual_context.X[person, neighbours[k]] *
                    analytic.residual_context.X[person, neighbours[l]]
            end
        end
        destination += (length(neighbours) + 1)^2
    end
    covariate_spec["covariate_weights"] = covariate_weights
    covariate_envelope = PDMPSamplersRBridge._build_omrf_residual_envelope(
        covariate_spec, anchor, analytic.residual_context)
    @test covariate_envelope.component_scales! isa
        PDMPSamplersRBridge.OMRFCovariateLocalScales
    covariate_subsampling = Dict{String,Any}(
        "n_observations" => 3,
        "subsample_size" => 1,
        "residual_envelope" => covariate_spec)
    covariate_bank_ctx, covariate_bank_names =
        PDMPSamplersRBridge._new_stan_subsampling_context(
            lib, full, prior, 1)
    covariate_model, _, _, _, covariate_manager =
        PDMPSamplersRBridge._build_stan_subsampling_model(
            covariate_bank_ctx, covariate_bank_names,
            covariate_subsampling, anchor, nothing, falses(9);
            anchor_capacity=2)
    caller_anchor = copy(anchor)
    covariate_prepared = PDMPSamplersRBridge._prepare_omrf_anchor(
        covariate_manager, anchor .+ 0.01)
    @test covariate_prepared.envelope.weights ===
        covariate_model.grad.envelope.weights
    @test length(covariate_prepared.envelope.component_scales!.
        coordinate_offsets) == length(anchor)
    first_inserted_idx = PDMPSamplersRBridge._insert_omrf_anchor!(
        covariate_manager, covariate_prepared)
    covariate_manager.active_idx = first_inserted_idx
    PDMPSamplers.refresh_anchor!(covariate_model.grad,
        covariate_prepared.anchor)

    # The first replacement creates an inactive scratch state. The following
    # preparation must recycle its arrays without mutating the active anchor,
    # and the fused native pass must agree with BridgeStan's exact gradient.
    replacement_anchor = anchor .- 0.02
    replacement = PDMPSamplersRBridge._prepare_omrf_anchor(
        covariate_manager, replacement_anchor)
    PDMPSamplersRBridge._insert_omrf_anchor!(
        covariate_manager, replacement)
    recycled_context = covariate_manager.recycled_state.residual_context
    recycled_probabilities = recycled_context.anchor_probabilities
    recycled_envelope = covariate_manager.recycled_state.envelope
    recycled_scales = recycled_envelope.component_scales!
    active_before_recycling = copy(covariate_manager.oracle.anchor)
    recycled_anchor = anchor .+ 0.03
    recycled_before = covariate_manager.recycled_preparations
    recycled_prepared = PDMPSamplersRBridge._prepare_omrf_anchor(
        covariate_manager, recycled_anchor)
    @test recycled_prepared.residual_context === recycled_context
    @test recycled_prepared.residual_context.anchor_probabilities ===
        recycled_probabilities
    @test recycled_prepared.envelope.component_scales! === recycled_scales
    @test recycled_prepared.envelope.component_scales!.anchor ==
        recycled_prepared.anchor
    @test covariate_manager.oracle.anchor == active_before_recycling
    @test anchor == caller_anchor
    @test covariate_manager.recycled_preparations == recycled_before + 1
    recycled_stan_full = zeros(9)
    PDMPSamplersRBridge._clear_gradient!(covariate_bank_ctx, :full,
        recycled_stan_full, recycled_anchor)
    @test recycled_prepared.full_anchor ≈
        recycled_stan_full atol=3e-10 rtol=3e-10
    recycled_prior = zeros(9)
    PDMPSamplersRBridge._clear_gradient!(covariate_bank_ctx, :prior,
        recycled_prior, recycled_anchor)
    @test recycled_prepared.prior_anchor ≈
        recycled_prior atol=2e-12 rtol=2e-12
    covariate_hcv_spec = copy(covariate_spec)
    covariate_hcv_spec["use_hcv"] = true
    covariate_hcv_spec["hcv_damping"] = 9.0
    covariate_hcv_spec["hcv_remainder_weights"] =
        0.75 .* (2 .* covariate_hcv_spec["node_weights"]) .^ 1.5
    covariate_hcv_subsampling = copy(covariate_subsampling)
    covariate_hcv_subsampling["residual_envelope"] = covariate_hcv_spec
    covariate_hcv_ctx, covariate_hcv_names =
        PDMPSamplersRBridge._new_stan_subsampling_context(
            lib, full, prior, 1)
    _, _, _, _, covariate_hcv_manager =
        PDMPSamplersRBridge._build_stan_subsampling_model(
            covariate_hcv_ctx, covariate_hcv_names,
            covariate_hcv_subsampling, anchor, nothing, falses(9);
            anchor_capacity=2)
    covariate_hcv_prepared =
        PDMPSamplersRBridge._prepare_omrf_anchor(
            covariate_hcv_manager, anchor .+ 0.01)
    @test covariate_hcv_prepared.envelope.component_scales! isa
        PDMPSamplersRBridge.OMRFDampedHCVNodeLocalScales
    velocity = collect(range(-0.75, 0.85; length=9))
    for flow in (ZigZag(9), BouncyParticle(9, 0.2), Boomerang(9))
        state = PDMPState(0.0,
            SkeletonPoint(copy(position), copy(velocity)))
        PDMPSamplers.initialize_flow_state!(state, flow)
        cell_scales = similar(model.grad.envelope.cell_scales)
        PDMPSamplers.component_cell_scales!(cell_scales,
            model.grad.envelope, state, flow, 0.05, 0.35)
        pattern_cell_scales = similar(pattern_envelope.cell_scales)
        PDMPSamplers.component_cell_scales!(pattern_cell_scales,
            pattern_envelope, state, flow, 0.05, 0.35)
        covariate_cell_scales = similar(covariate_envelope.cell_scales)
        PDMPSamplers.component_cell_scales!(covariate_cell_scales,
            covariate_envelope, state, flow, 0.05, 0.35)
        for person in 1:3, t in range(0.05, 0.35; length=9)
            candidate = move_forward_time(state, t, flow)
            residual_at_t = finite_difference_person_neggrad(
                candidate.ξ.x, X[person]) -
                finite_difference_person_neggrad(anchor, X[person])
            actual = PDMPSamplers.λ(candidate, residual_at_t, flow) +
                PDMPSamplers.λ(candidate, -residual_at_t, flow)
            bound = dot(cell_scales, @view node_weights[:, person])
            @test actual <= bound * (1 + 2e-7) + 2e-8
            pattern_bound = sum(partition -> pattern_cell_scales[
                    pattern_envelope.groups[partition, person]],
                axes(pattern_envelope.groups, 1))
            @test actual <= pattern_bound * (1 + 2e-7) + 2e-8
            covariate_bound = dot(covariate_cell_scales,
                @view covariate_weights[:, person])
            @test actual <= covariate_bound * (1 + 2e-7) + 2e-8

            factor_node = mod1(person + 1, 3)
            factor = (factor_node - 1) * 3 + person
            factor_model.grad.residual_oracle(
                factor_residual, candidate.ξ.x, [factor], anchor)
            factor_bound = dot(cell_scales,
                @view factor_model.grad.envelope.weights[:, factor])
            factor_actual = PDMPSamplers.λ(candidate, factor_residual, flow) +
                PDMPSamplers.λ(candidate, -factor_residual, flow)
            @test factor_actual <= factor_bound * (1 + 2e-7) + 2e-8

            hcv_oracle(hcv_residual, candidate.ξ.x, [person], anchor)
            hcv_cell_scales = similar(hcv_envelope.cell_scales)
            PDMPSamplers.component_cell_scales!(hcv_cell_scales,
                hcv_envelope, state, flow, 0.05, 0.35)
            hcv_actual = PDMPSamplers.λ(candidate, hcv_residual, flow) +
                PDMPSamplers.λ(candidate, -hcv_residual, flow)
            hcv_bound = dot(hcv_cell_scales,
                @view hcv_envelope.weights[:, person])
            @test hcv_actual <= hcv_bound * (1 + 2e-7) + 2e-8
        end
        PDMPSamplers.component_cell_scales!(pattern_cell_scales,
            pattern_envelope, state, flow, 0.05, 0.35)
        @test @allocated(PDMPSamplers.component_cell_scales!(
            pattern_cell_scales, pattern_envelope, state, flow,
            0.05, 0.35)) == 0
        covariate_workspace = covariate_envelope.component_scales!
        workspace_ids = objectid.((covariate_workspace.coordinate_offsets,
            covariate_workspace.coordinate_cosines,
            covariate_workspace.coordinate_sines))
        PDMPSamplers.component_scales!(covariate_envelope.scales,
            covariate_envelope, state, flow, 0.2)
        PDMPSamplers.component_cell_scales!(covariate_cell_scales,
            covariate_envelope, state, flow, 0.05, 0.35)
        @test objectid.((covariate_workspace.coordinate_offsets,
            covariate_workspace.coordinate_cosines,
            covariate_workspace.coordinate_sines)) == workspace_ids
        @test @allocated(PDMPSamplers.component_scales!(
            covariate_envelope.scales, covariate_envelope,
            state, flow, 0.2)) == 0
        @test @allocated(PDMPSamplers.component_cell_scales!(
            covariate_cell_scales, covariate_envelope, state, flow,
            0.05, 0.35)) == 0
        if flow isa Union{BouncyParticle,PDMPSamplers.AnyBoomerang}
            candidate = move_forward_time(state, 0.2, flow)
            PDMPSamplers.component_scales!(model.grad.envelope.scales,
                model.grad.envelope, candidate, flow, 0.0)
            PDMPSamplers.component_scales!(pattern_envelope.scales,
                pattern_envelope, candidate, flow, 0.0)
            for person in 1:3
                norm_bound = PDMPSamplers.observation_residual_bound(
                    model.grad.envelope, person)
                pattern_bound =
                    PDMPSamplersRBridge._omrf_person_pattern_bound_at(
                        analytic.residual_context, anchor, candidate, flow,
                        person, 0.0)
                @test pattern_bound ≈
                    PDMPSamplers.observation_residual_bound(
                        pattern_envelope, person) atol=2e-14 rtol=2e-14
                @test pattern_bound <= norm_bound * (1 + 2e-14) + 2e-14
                wrapped = PDMPSamplers.WithResidualStats(analytic, nothing)
                tightened = PDMPSamplers.subsampling_residual_subset_bound(
                    wrapped, candidate, flow, 0.0, norm_bound,
                    [person], 1.0)
                @test tightened ≈ pattern_bound atol=2e-14 rtol=2e-14
            end
            one_person = [1]
            wrapped = PDMPSamplers.WithResidualStats(analytic, nothing)
            norm_bound = PDMPSamplers.observation_residual_bound(
                model.grad.envelope, 1)
            PDMPSamplers.subsampling_residual_subset_bound(wrapped, candidate,
                flow, 0.0, norm_bound, one_person, 1.0)
            @test PDMPSamplers.subsampling_residual_subset_bound(
                wrapped, candidate, flow, 0.0, norm_bound,
                one_person, 1.0) <= norm_bound
        end
    end
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

@testset "OMRF subsampling event law and replicated posterior agreement" begin
    root = normpath(joinpath(@__DIR__, "..", ".."))
    stan = joinpath(root, "inst", "stan", "omrf", "omrf.stan")
    full = joinpath(root, "tests", "stan", "omrf_subsampling_full.json")
    prior = joinpath(root, "tests", "stan", "omrf_subsampling_prior.json")
    header = joinpath(root, "inst", "stan", "pdmp_subsample.hpp")
    lib = PDMPSamplersRBridge._compile_model_with_header(stan, header)
    X = [0 1 2; 2 0 1; 1 2 0]
    weights = omrf_curvature_weights(X, [3, 3, 3])
    anchor = zeros(9)
    position = collect(range(-0.18, 0.22; length=9))
    velocity = collect(range(-0.8, 0.9; length=9))
    ctx, _ = PDMPSamplersRBridge._new_stan_subsampling_context(
        lib, full, prior, 1)
    oracle = PDMPSamplersRBridge.StanSubsamplingOracle(ctx, anchor, 3)
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

    subsampling = Dict{String,Any}(
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
    subsampling_runs = map(seeds) do seed
        prepared = PDMPSamplersRBridge.prepare_stan_subsampling(
            lib, full, prior, subsampling, 1)
        model, subsampling_ctx, _ = PDMPSamplersRBridge._build_stan_subsampling_model(
            only(prepared.contexts), prepared.unc_names, subsampling,
            anchor, nothing, falses(9))
        trace, stats = pdmp_sample(
            initial, BouncyParticle(9, 0.4), model,
            GridThinningStrategy(N=12, t_max=0.5, lazy=false,
                bound_violation=:throw),
            0.0, 1_200.0; progress=false, seed,
            statistic_counter=PDMPSamplers.DevelStatisticCounter)
        (trace=trace, stats=stats, ctx=subsampling_ctx)
    end
    subsampling_traces = getproperty.(subsampling_runs, :trace)

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
    subsampling_summary = replicated_summary(subsampling_traces)
    combined_mcse = sqrt.(full_summary.mcse .^ 2 .+ subsampling_summary.mcse .^ 2)
    standardized_difference = abs.(full_summary.mean - subsampling_summary.mean) ./ combined_mcse
    @test full_summary.minimum_ess > 20
    @test subsampling_summary.minimum_ess > 20
    @test all(isfinite, standardized_difference)
    @test maximum(standardized_difference) < 4.5

    for run in subsampling_runs
        @test run.stats.residual_oracle_calls > 0
        @test run.ctx.counts.persons_evaluated ==
            run.ctx.m * run.ctx.counts.selected_gradient
        @test run.ctx.counts.selected_gradient ==
            run.stats.residual_oracle_calls
        @test run.stats.subsampling_subset_evaluations ==
            run.stats.residual_oracle_calls
        @test run.stats.subsampling_final_reflections <=
            run.stats.subsampling_subset_evaluations
        @test run.stats.full_gradient_calls == 0
    end
end

@testset "dependent sticky subsampling target preserves nuisance prior and model odds" begin
    root = normpath(joinpath(@__DIR__, "..", ".."))
    stan = joinpath(root, "tests", "stan", "subsampling_gaussian_nuisance.stan")
    full = joinpath(root, "tests", "stan", "subsampling_gaussian_nuisance_full.json")
    prior = joinpath(root, "tests", "stan", "subsampling_gaussian_nuisance_prior.json")
    header = joinpath(root, "inst", "stan", "pdmp_subsample.hpp")
    lib = PDMPSamplersRBridge._compile_model_with_header(stan, header)
    subsampling = Dict{String,Any}(
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
    prepared = PDMPSamplersRBridge.prepare_stan_subsampling(
        lib, full, prior, subsampling, 1)
    model, ctx, _ = PDMPSamplersRBridge._build_stan_subsampling_model(
        only(prepared.contexts), prepared.unc_names, subsampling,
        zeros(2), slab, BitVector([true, false]))
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

    prepared = PDMPSamplersRBridge.prepare_stan_subsampling(
        lib, full, prior, subsampling, 1)
    sampled = PDMPSamplersRBridge.r_pdmp_stan_subsampling(
        prepared, subsampling, [0.1, 0.8], "ZigZag",
        "GridThinningStrategy", zeros(2), Matrix{Float64}(I, 2, 2);
        T=4_000.0, grid_n=12, grid_t_max=0.5,
        sticky=true, can_stick=BitVector([true, false]),
        model_prior=Dict("prob" => [0.5]), slab_prior=slab,
        show_progress=false, seed=944)
    counters = only(sampled["subsampling_context_counters"])
    stats = sampled["stats"]
    @test counters["persons_evaluated"] ==
        2 * counters["selected_gradient_calls"]
    @test counters["selected_gradient_calls"] ==
        stats["residual_oracle_calls"][1]
    @test counters["sampling_full_gradient_calls"] == 0
    @test counters["sampling_model_constructions"] == 0
    @test counters["sampling_data_constructions"] == 0
    @test stats["subsampling_subset_evaluations"][1] ==
        stats["residual_oracle_calls"][1]
    @test stats["subsampling_final_reflections"][1] <=
        stats["subsampling_subset_evaluations"][1]

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
    subsampling = Dict{String,Any}(
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
    prepared = PDMPSamplersRBridge.prepare_stan_subsampling(
        lib, full, prior, subsampling, 1)
    sampled = PDMPSamplersRBridge.r_pdmp_stan_subsampling(
        prepared, subsampling, [0.1], "ZigZag",
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
    @test only(sampled["subsampling_context_counters"])["sampling_full_gradient_calls"] == 0
end

using Test, LinearAlgebra, Random, Statistics
using PDMPSamplers, BridgeStan

include(normpath(joinpath(@__DIR__, "..", "..", "inst", "julia",
    "main_interface_function.jl")))
using .PDMPSamplersRBridge

@testset "OMRF observed-threshold setup counts" begin
    X = [0 1; 1 0; 2 1]
    observed = PDMPSamplersRBridge._omrf_observed_threshold_counts(
        X, [1, 3], 3)
    @test observed == [1.0, 1.0, 2.0]
end

@testset "independent-prior interval uses parameter mask indices" begin
    prior = PDMPSamplersRBridge.OMRFIndependentAnalyticPrior(
        [1, 2], [3, 4, 5], 1.0, 1.0, :gaussian, 1.0)
    active = BitVector([true, true, false, true, true])
    xlo = fill(0.5, 5)
    xhi = copy(xlo)
    glo = zeros(5)
    ghi = zeros(5)
    PDMPSamplersRBridge._omrf_prior_intervals!(
        glo, ghi, prior, xlo, xhi, active)
    @test glo[3] == ghi[3] == 0.0
    @test glo[4] <= 0.5 <= ghi[4]
    @test glo[5] <= 0.5 <= ghi[5]
end

function omrf_allocation_next_event!(rng, model, flow, alg, state, cache, stats)
    PDMPSamplers.next_event_time(rng, model, flow, alg, state, cache, stats,
        0.25, false, :horizon_hit)
    return nothing
end

function omrf_allocation_signed_cell!(out, envelope, state, flow, left, right)
    PDMPSamplers.component_cell_scales!(out, envelope, state, flow,
        left, right)
    return nothing
end

function omrf_allocation_residual_bound!(fixture)
    PDMPSamplers.subsampling_residual_subset_bound(
        fixture.oracle, fixture.state, fixture.gradient, fixture.flow,
        fixture.D, fixture.M, fixture.subset, fixture.scale)
    return nothing
end

function omrf_allocation_subset_bound!(fixture)
    PDMPSamplers.subsampling_subset_bound(
        fixture.oracle, fixture.state, fixture.gradient, fixture.flow,
        fixture.D, fixture.M, fixture.subset, fixture.scale)
    return nothing
end

function omrf_allocation_anchor_select!(fixture)
    PDMPSamplersRBridge._select_omrf_anchor!(
        fixture.manager, fixture.cv, fixture.position; phase=:main)
    return nothing
end

@testset "Beta-Bernoulli sticky wrapper excludes always-active parameters" begin
    inner = PDMPSamplers.GridThinningStrategy()
    mask = BitVector([false, true, true, false, true])
    wrapped = PDMPSamplersRBridge.wrap_sticky(inner, true,
        Dict(:a => 1.0, :b => 4.0), fill(0.5, 5), mask)
    free = BitVector([true, true, false, true, true])
    # Only the three stickable coordinates define the model size; two are
    # active.  The always-active coordinates at positions one and four do not
    # enter the beta-binomial odds.
    @test wrapped.κ(3, zeros(5), free) ==
        (1 + 2) / (4 + 3 - 2 - 1) * 0.5
end

# Residual of one person-node factor of an undamped-HCV analytic OMRF oracle.
function omrf_factor_residual!(out, oracle, x, person, node)
    context = oracle.residual_context
    fill!(out, 0.0)
    @. context.displacement = x - oracle.anchor
    PDMPSamplersRBridge._omrf_node_residual!(
        context, out, x, person, node, 1.0)
    return out
end

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

function omrf_person_hessian(q, x, seen)
    P = length(seen)
    starts = cumsum(vcat(1, seen[1:(end - 1)] .- 1))
    n_thresholds = sum(seen .- 1)
    edges = [(j, k) for j in 1:(P - 1) for k in (j + 1):P]
    H = zeros(length(q), length(q))
    for node in 1:P
        B = zeros(seen[node] - 1, length(q))
        for u in axes(B, 1)
            B[u, starts[node] + u - 1] = 1
            for (e, (j, k)) in enumerate(edges)
                if node == j
                    B[u, n_thresholds + e] = u * x[k]
                elseif node == k
                    B[u, n_thresholds + e] = u * x[j]
                end
            end
        end
        logits = [0.0; B * q]
        probabilities = exp.(logits .- maximum(logits))
        probabilities ./= sum(probabilities)
        p = probabilities[2:end]
        H .+= B' * (Diagonal(p) - p * p') * B
    end
    return H
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
    context.displacement .= position - anchor
    expected_support = [1, 7, 8, 2, 3, 9, 4, 5, 6]
    @test context.likelihood_support == expected_support
    @test all(<=(no_thresholds + no_edges), context.likelihood_support)
    PDMPSamplersRBridge.omrf_residual!(context, residual, position, Int[])
    @test iszero(residual)
    @test context.likelihood_support == expected_support
    for person in axes(X, 1)
        PDMPSamplersRBridge.omrf_residual!(
            context, residual, position, [person])
        expected = finite_difference_person_neggrad_variable(
            position[1:(no_thresholds + no_edges)], view(X, person, :), seen) -
            finite_difference_person_neggrad_variable(
                anchor[1:(no_thresholds + no_edges)], view(X, person, :), seen)
        H = omrf_person_hessian(anchor[1:length(expected)],
            view(X, person, :), seen)
        expected -= H * (position[1:length(expected)] - anchor[1:length(expected)])
        @test residual[1:(no_thresholds + no_edges)] ≈ expected atol=2e-9
        @test residual[(no_thresholds + no_edges + 1):end] == zeros(2)
    end

end

@testset "node-scale analytic prior includes transformed scale terms" begin
    seen = [2, 2, 2]
    no_thresholds = sum(seen .- 1)
    no_edges = 3
    d = no_thresholds + no_edges + 1 + length(seen)
    names = vcat(
        ["thresholds_0.$j" for j in 1:no_thresholds],
        ["interactions_0.$j" for j in 1:no_edges],
        ["global_scale"], ["node_scale.$j" for j in 1:length(seen)])
    spec = Dict{String,Any}(
        "X" => [0 1 0; 1 0 1],
        "seen" => seen,
        "thresholds" => "thresholds_0",
        "interactions" => "interactions_0",
        "analytic_prior" => Dict{String,Any}(
            "type" => "node_scale_gaussian",
            "threshold_alpha" => 1.0,
            "threshold_beta" => 1.0,
            "interaction_scale" => 0.7,
            "global_scale" => "global_scale",
            "node_scale" => "node_scale",
            "global_scale_sd" => 1.2,
            "node_scale_sd" => 0.8))
    context = PDMPSamplersRBridge._omrf_residual_context(
        spec, names, zeros(d), PDMPSamplersRBridge.StanSubsamplingCallCounts())
    @test context.global_scale_index == no_thresholds + no_edges + 1
    @test context.node_scale_indices == (no_thresholds + no_edges + 2):d

    position = [0.2, -0.1, 0.3, 0.4, -0.2, 0.1, 0.5, -0.3, 0.25, 0.1]
    gradient = zeros(d)
    prior = PDMPSamplersRBridge._omrf_analytic_prior(spec, context)
    PDMPSamplersRBridge._omrf_prior_gradient!(gradient, prior, position)
    global_scale = exp(position[7])
    node_scales = exp.(position[8:10])
    expected_global = global_scale^2 / 1.2^2 - 1
    expected_nodes = node_scales .^ 2 ./ 0.8^2 .- 1
    edge = 1
    for left in 1:2, right in (left + 1):3
        scale2 = 0.7^2 * global_scale^2 * node_scales[left] * node_scales[right]
        ratio = position[3 + edge]^2 / scale2
        expected_global += 1 - ratio
        expected_nodes[left] += 0.5 * (1 - ratio)
        expected_nodes[right] += 0.5 * (1 - ratio)
        @test gradient[3 + edge] ≈ position[3 + edge] / scale2
        edge += 1
    end
    @test gradient[7] ≈ expected_global
    @test gradient[8:10] ≈ expected_nodes

    # The node-scale prior HVP must agree with a finite difference of the
    # analytic prior gradient.  In particular, the global log scale enters
    # edge log variances with coefficient two, while each node log scale enters
    # with coefficient one.
    direction = collect(range(-0.27, 0.36; length=d))
    pure_hvp = zeros(d)
    PDMPSamplersRBridge._omrf_prior_hvp!(pure_hvp, prior, position, direction)
    step = 2e-6
    gradient_plus = zeros(d)
    gradient_minus = zeros(d)
    PDMPSamplersRBridge._omrf_prior_gradient!(
        gradient_plus, prior, position .+ step .* direction)
    PDMPSamplersRBridge._omrf_prior_gradient!(
        gradient_minus, prior, position .- step .* direction)
    @test pure_hvp ≈ (gradient_plus - gradient_minus) / (2step) atol=2e-8 rtol=2e-8

    # Check symmetry independently of the finite-difference direction.
    direction_two = collect(range(0.41, -0.22; length=d))
    pure_hvp_two = zeros(d)
    PDMPSamplersRBridge._omrf_prior_hvp!(
        pure_hvp_two, prior, position, direction_two)
    @test dot(direction_two, pure_hvp) ≈ dot(direction, pure_hvp_two) atol=2e-10 rtol=2e-10

end

@testset "Boomerang log-scale floor name resolution is exact" begin
    names = ["thresholds_0.1", "global_scale", "node_scale.1",
        "node_scale_extra", "interactions_0.1"]
    @test PDMPSamplersRBridge._boomerang_logscale_floor_mask(names, 5) ==
        BitVector([false, true, true, false, false])
    @test_throws DimensionMismatch PDMPSamplersRBridge._boomerang_logscale_floor_mask(
        names, 4)
end

@testset "retained full-data invariant reports every work category" begin
    clean = [Dict{String,Any}(
        "anchor_main_preparations" => 0,
        "main_native_full_data_scans" => 0)]
    @test isnothing(PDMPSamplersRBridge._enforce_no_main_full_data!(
        [0], clean))
    dirty = [Dict{String,Any}(
        "anchor_main_preparations" => 2,
        "main_native_full_data_scans" => 3)]
    error = try
        PDMPSamplersRBridge._enforce_no_main_full_data!([1], dirty)
        nothing
    catch exception
        exception
    end
    @test error isa ErrorException
    message = sprint(showerror, error)
    @test contains(message, "main_full_gradient_calls=1")
    @test contains(message, "anchor_main_preparations=2")
    @test contains(message, "main_native_full_data_scans=3")
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
        H = omrf_person_hessian(anchor, X[only(subset)], [3, 3, 3])
        @test analytic_residual ≈
            (selected_position - prior_position) -
            (selected_anchor - prior_anchor) - H * (position - anchor) atol=3e-10 rtol=3e-10
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

    stan_full_gradient = zeros(9)
    PDMPSamplersRBridge._clear_gradient!(
        analytic_prior_ctx, :full, stan_full_gradient, position)
    direction = collect(range(-0.3, 0.4; length=9))
    step = 2e-6

    # The undamped global HCV splits the full gradient exactly into its
    # deterministic part and the residual summed over every person.
    hcv_spec = copy(analytic_prior_spec)
    hcv_oracle = PDMPSamplersRBridge.OMRFAnalyticSubsamplingOracle(
        analytic_prior_ctx, hcv_spec, analytic_prior_names, anchor)
    hcv_deterministic = zeros(9)
    hcv_residual = zeros(9)
    hcv_oracle(hcv_deterministic, position)
    hcv_oracle(hcv_residual, position, [1, 2, 3], anchor)
    @test hcv_deterministic + hcv_residual ≈ stan_full_gradient atol=3e-10 rtol=3e-10

    @testset "categorical HCV one-eighth remainder" begin
        rng = Xoshiro(8192)
        for q in (1, 2, 5, 12), scale in (1e-8, 0.1, 1.0, 8.0), _ in 1:100
            logits = scale .* randn(rng, q)
            delta = scale .* randn(rng, q)
            velocity = scale .* randn(rng, q)
            remainder, bound =
                PDMPSamplersRBridge._omrf_categorical_hcv_remainder(
                    logits, delta, velocity)
            @test abs(remainder) <= bound + 2e-12 * (1 + bound)
        end
        adversarial = (
            ([0.0], [20.0], [-30.0]),
            ([-40.0, 40.0], [40.0, -40.0], [60.0, -60.0]),
            ([100.0, -100.0, 0.0], [-200.0, 200.0, 1e-10],
                [1e3, -1e3, 0.0]),
        )
        for (logits, delta, velocity) in adversarial
            remainder, bound =
                PDMPSamplersRBridge._omrf_categorical_hcv_remainder(
                    logits, delta, velocity)
            @test abs(remainder) <= bound + 2e-11 * (1 + bound)
        end
    end

    @testset "direct factor-roof person HCV cell certificate" begin
        direct_spec = copy(analytic_prior_spec)
        direct_spec["exact_selected_refinement"] = true
        direct_spec["weights"] = zeros(1, 3)
        direct_oracle = PDMPSamplersRBridge.OMRFAnalyticSubsamplingOracle(
            analytic_prior_ctx, direct_spec, analytic_prior_names, anchor)
        direct_envelope = PDMPSamplersRBridge._build_stan_residual_envelope(
            direct_spec, direct_oracle.anchor, direct_oracle)
        @test direct_envelope isa
            PDMPSamplersRBridge.OMRFDirectPersonNodeHCVEnvelope

        velocity = collect(range(-0.8, 0.9; length=9))
        direct_state = PDMPState(0.0,
            SkeletonPoint(copy(position), copy(velocity)))
        direct_flow = BouncyParticle(9)
        PDMPSamplers.initialize_flow_state!(direct_state, direct_flow)
        left, right = 0.0, 0.3
        PDMPSamplers.component_cell_scales!(direct_envelope.cell_scales,
            direct_envelope, direct_state, direct_flow, left, right)
        # The certificate table is indexed by person-node factor
        # (node - 1) * N + person even though persons are sampled.
        direct_residual = zeros(9)
        maximum_violation = 0.0
        for factor in 1:9, t in range(left, right; length=65)
            point = position .+ t .* velocity
            node0, person0 = divrem(factor - 1, 3)
            omrf_factor_residual!(direct_residual, direct_oracle, point,
                person0 + 1, node0 + 1)
            exact = abs(dot(velocity, direct_residual))
            bound = PDMPSamplers.observation_residual_bound(
                direct_envelope, factor)
            maximum_violation = max(maximum_violation, exact - bound)
        end
        @test maximum_violation <= 2e-12

        # Factor roof, person gradient: the certified person-node table is
        # summed per person, one person is sampled, and the sum of all of
        # that person's node residuals is evaluated.  The size-biased person
        # law is C_i/sum(C_i).
        hybrid_spec = direct_spec
        hybrid_oracle = PDMPSamplersRBridge.OMRFAnalyticSubsamplingOracle(
            analytic_prior_ctx, hybrid_spec, analytic_prior_names, anchor)
        hybrid_envelope = PDMPSamplersRBridge._build_stan_residual_envelope(
            hybrid_spec, hybrid_oracle.anchor, hybrid_oracle)
        @test PDMPSamplers.n_observations(hybrid_envelope) == 3
        PDMPSamplers.component_cell_scales!(hybrid_envelope.cell_scales,
            hybrid_envelope, direct_state, direct_flow, left, right)
        PDMPSamplers.prepare_residual_sampling!(
            hybrid_envelope, direct_state, direct_flow, 0.0)
        @test hybrid_envelope.cell_mass[hybrid_envelope.active_cell] ==
            direct_envelope.cell_mass[direct_envelope.active_cell]

        # Research S1a keeps the whole-person estimator and mark law while
        # summing signed node endpoints before taking an absolute value.
        signed_spec = copy(hybrid_spec)
        signed_spec["person_signed_roof"] = true
        signed_envelope = PDMPSamplersRBridge._build_stan_residual_envelope(
            signed_spec, hybrid_oracle.anchor, hybrid_oracle)
        PDMPSamplers.component_cell_scales!(signed_envelope.cell_scales,
            signed_envelope, direct_state, direct_flow, left, right)
        signed_B = signed_envelope.cell_mass[signed_envelope.active_cell]
        baseline_B = hybrid_envelope.cell_mass[hybrid_envelope.active_cell]
        @test signed_B <= baseline_B + 1e-12
        for person in 1:3, t in range(left, right; length=129)
            point = position .+ t .* velocity
            hybrid_oracle(direct_residual, point, [person], anchor)
            exact = abs(dot(velocity, direct_residual))
            certified = PDMPSamplersRBridge._omrf_direct_person_cell_bound(
                signed_envelope, person)
            @test exact <= certified + 2e-11 * (1 + certified)
        end
        @test sum(PDMPSamplersRBridge._omrf_direct_person_cell_bound(
            signed_envelope, person) for person in 1:3) ≈ signed_B
        signed_det_now = zeros(9)
        hybrid_oracle(signed_det_now, position)
        signed_D = PDMPSamplers.λ(direct_state,
            signed_det_now, direct_flow)
        signed_only_spec = copy(signed_spec)
        signed_only_spec["signed_only_person_roof"] = true
        signed_only = PDMPSamplersRBridge._build_stan_residual_envelope(
            signed_only_spec, hybrid_oracle.anchor, hybrid_oracle)
        PDMPSamplers.component_cell_scales!(signed_only.cell_scales,
            signed_only, direct_state, direct_flow, left, right)
        signed_only_B = signed_only.cell_mass[signed_only.active_cell]
        @test signed_only_B >= signed_B - 1e-12
        @test_throws ArgumentError PDMPSamplers.observation_residual_bound(
            signed_only, 1)
        for person in 1:3, t in range(left, right; length=129)
            hybrid_oracle(direct_residual, position .+ t .* velocity,
                [person], anchor)
            certified = PDMPSamplersRBridge._omrf_direct_person_cell_bound(
                signed_only, person)
            @test abs(dot(velocity, direct_residual)) <=
                certified + 2e-11 * (1 + certified)
        end
        signed_only_cv = PDMPSamplers.SubsampledControlVariate(
            hybrid_oracle, hybrid_oracle, signed_only, anchor, 1)
        PDMPSamplers.prepare_residual_sampling!(
            signed_only, direct_state, direct_flow, left)
        mark_rng = Xoshiro(0x5292)
        for _ in 1:100
            bound = PDMPSamplers.draw_subset!(
                mark_rng, signed_only_cv, signed_D, signed_only_B)
            person = signed_only_cv.subset[1]
            @test 1 <= person <= 3
            @test bound >= signed_D +
                3PDMPSamplersRBridge._omrf_direct_person_cell_bound(
                    signed_only, person) - 1e-12
        end
        cancelling = PDMPSamplersRBridge._build_stan_residual_envelope(
            signed_spec, hybrid_oracle.anchor, hybrid_oracle)
        cancelling.cell_left[1] = 0.0
        cancelling.cell_right[1] = 0.1
        cancelling.factor_workspace[1] = 2.0
        cancelling.factor_workspace[4] = 2.0
        cancelling.person_left_signed[1] = 0.0  # +2 and -2 nodes
        cancelling.person_right_signed[1] = 0.0
        cancelling.person_left_abs[1] = 4.0
        cancelling.person_right_abs[1] = 4.0
        PDMPSamplersRBridge._omrf_direct_finish_person_cell!(
            cancelling, 1, true)
        @test cancelling.cell_mass[1] < 1e-12
        @test cancelling.person_workspace[1] == 4.0
        signed_deterministic = zeros(9)
        hybrid_oracle(signed_deterministic, position)
        signed_D = PDMPSamplers.λ(direct_state,
            signed_deterministic, direct_flow)
        for person in 1:3
            hybrid_oracle(direct_residual, position, [person], anchor)
            selected_gradient = copy(signed_deterministic)
            actual = PDMPSamplers.subsampling_candidate_rate!(
                hybrid_oracle, direct_state, selected_gradient,
                direct_residual, 3.0, direct_flow, signed_D, [person])
            C = signed_D + 3 *
                PDMPSamplersRBridge._omrf_direct_person_cell_bound(
                    signed_envelope, person)
            @test actual <= C + 2e-11 * (1 + C)
            proposal_probability = C / (3 * (signed_D + signed_B))
            @test isapprox((signed_D + signed_B) * proposal_probability *
                actual / C, actual / 3; rtol=8eps(), atol=8eps())
        end

        hybrid_cv = PDMPSamplers.SubsampledControlVariate(
            hybrid_oracle, hybrid_oracle, hybrid_envelope, anchor, 1)
        hybrid_rng = Xoshiro(0xfa7700)
        person_cell_bounds = zeros(3)
        for person in 1:3, node in 1:3
            factor = (node - 1) * 3 + person
            person_cell_bounds[person] +=
                PDMPSamplers.observation_residual_bound(
                    hybrid_envelope, factor)
        end
        draws = zeros(Int, 3)
        B = hybrid_envelope.cell_mass[hybrid_envelope.active_cell]
        for _ in 1:100_000
            selected_roof = PDMPSamplers.draw_subset!(
                hybrid_rng, hybrid_cv, 0.0, B)
            person = hybrid_cv.subset[1]
            draws[person] += 1
            @test isapprox(selected_roof, 3person_cell_bounds[person];
                rtol=2e-14, atol=2e-14)
        end
        @test maximum(abs.(draws ./ sum(draws) .-
            person_cell_bounds ./ sum(person_cell_bounds))) < 0.01

        hybrid_residual = zeros(9)
        factor_residual = zeros(9)
        for person in 1:3
            hybrid_oracle(hybrid_residual, position, [person], anchor)
            fill!(factor_residual, 0.0)
            for node in 1:3
                omrf_factor_residual!(direct_residual, direct_oracle,
                    position, person, node)
                factor_residual .+= direct_residual
            end
            @test hybrid_residual ≈ factor_residual atol=2e-13 rtol=2e-13
        end
        for flow in (BouncyParticle(9), ZigZag(9))
            state = PDMPState(0.0,
                SkeletonPoint(copy(position), copy(velocity)))
            PDMPSamplers.initialize_flow_state!(state, flow)
            PDMPSamplers.component_cell_scales!(hybrid_envelope.cell_scales,
                hybrid_envelope, state, flow, left, right)
            PDMPSamplers.prepare_residual_sampling!(
                hybrid_envelope, state, flow, 0.0)
            deterministic = zeros(9)
            hybrid_oracle(deterministic, position)
            D = PDMPSamplers.λ(state, deterministic, flow)
            B = hybrid_envelope.cell_mass[hybrid_envelope.active_cell]
            for person in 1:3
                person_bound = sum(PDMPSamplers.observation_residual_bound(
                    hybrid_envelope, (node - 1) * 3 + person)
                    for node in 1:3)
                hybrid_oracle(hybrid_residual, position, [person], anchor)
                selected_gradient = copy(deterministic)
                actual = PDMPSamplers.subsampling_candidate_rate!(
                    hybrid_oracle, state, selected_gradient, hybrid_residual,
                    3.0, flow, D, [person])
                selected_roof = D + 3person_bound
                @test actual <= selected_roof + 2e-11
                proposal_probability = selected_roof / (3 * (D + B))
                @test isapprox((D + B) * proposal_probability *
                    actual / selected_roof, actual / 3;
                    rtol=8eps(), atol=8eps())
            end
        end

        # A dominant non-reference anchor category can make
        # `1 - sum(p_nonreference)` round to zero even though its reference
        # probability is needed after a large displacement.  The direct HCV
        # cache must retain the stably normalized reference mass itself.
        extreme_anchor = copy(anchor)
        extreme_anchor[1] = 50.0
        extreme_anchor[2] = 100.0
        extreme_oracle = PDMPSamplersRBridge.OMRFAnalyticSubsamplingOracle(
            analytic_prior_ctx, direct_spec, analytic_prior_names,
            extreme_anchor)
        extreme_envelope = PDMPSamplersRBridge._build_stan_residual_envelope(
            direct_spec, extreme_oracle.anchor, extreme_oracle)
        extreme_context = extreme_oracle.residual_context
        @test 1.0 - sum(view(extreme_context.anchor_probabilities,
            1:2, 1)) == 0.0
        @test extreme_context.anchor_reference_probabilities[1, 1] > 0.0
        extreme_position = copy(extreme_anchor)
        extreme_position[1] -= 60.0
        extreme_position[2] -= 120.0
        extreme_velocity = collect(range(-1.1, 0.7; length=9))
        extreme_state = PDMPState(0.0,
            SkeletonPoint(extreme_position, extreme_velocity))
        PDMPSamplers.initialize_flow_state!(extreme_state, direct_flow)
        PDMPSamplers.component_cell_scales!(extreme_envelope.cell_scales,
            extreme_envelope, extreme_state, direct_flow, 0.0, 0.05)
        extreme_residual = zeros(9)
        omrf_factor_residual!(extreme_residual, extreme_oracle,
            extreme_position, 1, 1)
        @test abs(dot(extreme_velocity, extreme_residual)) <=
            PDMPSamplers.observation_residual_bound(extreme_envelope, 1) +
            2e-12

        # Affine path caching only changes floating-point association. The
        # outward-protected cached result agrees with the independent scalar
        # construction to a few ulps and remains a valid certificate.
        for factor in 1:9
            node0, person0 = divrem(factor - 1, 3)
            scalar = PDMPSamplersRBridge._omrf_direct_hcv_factor_cell_bound(
                direct_envelope, direct_state, direct_flow,
                person0 + 1, node0 + 1, left, right)
            @test isapprox(scalar, direct_envelope.factor_workspace[factor];
                rtol=64eps(), atol=4eps())
            @test scalar <= PDMPSamplers.observation_residual_bound(
                direct_envelope, factor) + 4eps()
        end

        # Random factor/cell checks exceed the requested 20,000 comparisons.
        rng = Xoshiro(0xd1ce7)
        for _ in 1:20_001
            factor = rand(rng, 1:9)
            node0, person0 = divrem(factor - 1, 3)
            random_position = anchor .+ 0.5 .* randn(rng, 9)
            random_velocity = randn(rng, 9)
            random_state = PDMPState(0.0,
                SkeletonPoint(random_position, random_velocity))
            l = -0.2 + 0.4rand(rng)
            width = 0.4rand(rng)
            r = l + width
            bound = PDMPSamplersRBridge._omrf_direct_hcv_factor_cell_bound(
                direct_envelope, random_state, direct_flow,
                person0 + 1, node0 + 1, l, r)
            t = l + rand(rng) * width
            point = random_position .+ t .* random_velocity
            omrf_factor_residual!(direct_residual, direct_oracle, point,
                person0 + 1, node0 + 1)
            @test abs(dot(random_velocity, direct_residual)) <=
                bound + 3e-12 * (1 + bound)
        end

        # Mathematical adversaries: equal/dominant logits, zero velocity or
        # displacement, anchor-crossing cells, and opposite endpoint signs.
        adversarial_cells = (
            (zeros(3), zeros(3), zeros(3), 0.5),
            ([50.0, -50.0, -50.0], [-1.0, 0.5, 2.0],
                [3.0, -2.0, 1.0], 0.7),
            (zeros(3), [-0.3, 0.2, 0.1], [1.0, -1.0, 0.5], 0.6),
            ([1.0, -1.0, 0.0], zeros(3), [2.0, -3.0, 1.0], 0.4),
        )
        for (logits, delta, z, width) in adversarial_cells
            bound = PDMPSamplersRBridge._omrf_direct_categorical_cell_bound(
                logits, delta, z, width)
            for t in range(0.0, width; length=257)
                point_delta = delta .+ t .* z
                h, _ = PDMPSamplersRBridge._omrf_categorical_hcv_remainder(
                    logits, point_delta, z)
                @test abs(h) <= bound + 3e-12 * (1 + bound)
            end
        end

        # Inactive interactions are represented by exactly zero physical
        # velocity and remain covered, including a cell crossing the anchor.
        inactive_velocity = copy(velocity)
        inactive_velocity[direct_oracle.residual_context.interaction_indices[1]] = 0.0
        inactive_state = StickyPDMPState(Ref(0.0),
            SkeletonPoint(copy(anchor), inactive_velocity),
            trues(9))
        inactive_state.free[
            direct_oracle.residual_context.interaction_indices[1]] = false
        PDMPSamplers.component_cell_scales!(direct_envelope.cell_scales,
            direct_envelope, inactive_state, direct_flow, -0.2, 0.2)
        @test all(isfinite, view(direct_envelope.cell_cumulative, :, 1))

        # Cell reconstruction and cached marked selection are allocation-free
        # once buffers and method specializations have been prepared.
        direct_state2 = PDMPState(1.0,
            SkeletonPoint(copy(position), copy(velocity)))
        PDMPSamplers.component_cell_scales!(direct_envelope.cell_scales,
            direct_envelope, direct_state2, direct_flow, 0.0, 0.1)
        first_path_preparations =
            direct_envelope.context.counts.direct_path_preparations
        first_scans = direct_envelope.context.counts.direct_factor_bound_scans
        PDMPSamplers.component_cell_scales!(direct_envelope.cell_scales,
            direct_envelope, direct_state2, direct_flow, 0.0, 0.1)
        @test direct_envelope.context.counts.direct_factor_bound_scans ==
            first_scans

        # Adjacent cells reuse every left endpoint. Disabling reuse and
        # evaluating the same cell independently gives bit-identical bounds
        # and cumulative masses because both paths use the cached affine
        # representation.
        PDMPSamplers.component_cell_scales!(direct_envelope.cell_scales,
            direct_envelope, direct_state2, direct_flow, 0.1, 0.2)
        @test direct_envelope.context.counts.direct_adjacent_endpoint_reuses >= 9
        adjacent_bounds = copy(view(direct_envelope.cell_cumulative, :, 2))
        direct_envelope.last_boundary_valid = false
        PDMPSamplersRBridge._omrf_direct_hcv_fill_cell_cached!(
            direct_envelope, 0.1, 0.2, 3, false)
        @test adjacent_bounds == view(direct_envelope.cell_cumulative, :, 3)

        # Noncontiguous and out-of-order requests fall back to two endpoint
        # evaluations; a subsequent adjacent extension resumes reuse.
        before_left = direct_envelope.context.counts.direct_left_endpoint_evaluations
        PDMPSamplers.component_cell_scales!(direct_envelope.cell_scales,
            direct_envelope, direct_state2, direct_flow, 0.4, 0.5)
        @test direct_envelope.context.counts.direct_left_endpoint_evaluations ==
            before_left + 9
        before_reuse = direct_envelope.context.counts.direct_adjacent_endpoint_reuses
        PDMPSamplers.component_cell_scales!(direct_envelope.cell_scales,
            direct_envelope, direct_state2, direct_flow, 0.5, 0.6)
        @test direct_envelope.context.counts.direct_adjacent_endpoint_reuses ==
            before_reuse + 9
        before_left = direct_envelope.context.counts.direct_left_endpoint_evaluations
        PDMPSamplers.component_cell_scales!(direct_envelope.cell_scales,
            direct_envelope, direct_state2, direct_flow, 0.25, 0.3)
        @test direct_envelope.context.counts.direct_left_endpoint_evaluations ==
            before_left + 9

        # Every property that changes the physical affine flight invalidates
        # the cache explicitly, including active-set and anchor generations.
        @test PDMPSamplersRBridge._omrf_direct_path_matches(
            direct_envelope, direct_state2)
        direct_state2.ξ.θ[1] += 0.125
        @test !PDMPSamplersRBridge._omrf_direct_path_matches(
            direct_envelope, direct_state2)
        direct_state2.ξ.θ[1] -= 0.125
        direct_state2.ξ.x[1] += 0.125
        @test !PDMPSamplersRBridge._omrf_direct_path_matches(
            direct_envelope, direct_state2)
        direct_state2.ξ.x[1] -= 0.125
        direct_envelope.context.active_generation += 1
        @test !PDMPSamplersRBridge._omrf_direct_path_matches(
            direct_envelope, direct_state2)
        direct_envelope.context.active_generation -= 1
        direct_envelope.anchor[1] += 0.125
        @test !PDMPSamplersRBridge._omrf_direct_path_matches(
            direct_envelope, direct_state2)
        direct_envelope.anchor[1] -= 0.125
        direct_state2.t[] += 0.125
        @test !PDMPSamplersRBridge._omrf_direct_path_matches(
            direct_envelope, direct_state2)
        direct_state2.t[] -= 0.125

        # Path preparation, adjacent extension, and noncontiguous fallback
        # all reuse their setup buffers and allocate no steady-state bytes.
        PDMPSamplersRBridge._omrf_direct_prepare_path!(
            direct_envelope, direct_state2, direct_flow)
        GC.gc()
        @test @allocated(PDMPSamplersRBridge._omrf_direct_prepare_path!(
            direct_envelope, direct_state2, direct_flow)) == 0
        direct_envelope.cell_count = 0
        direct_envelope.last_boundary_valid = false
        PDMPSamplersRBridge._omrf_direct_hcv_fill_cell_cached!(
            direct_envelope, 0.0, 0.1, 1, true)
        GC.gc()
        @test @allocated(PDMPSamplersRBridge._omrf_direct_hcv_fill_cell_cached!(
            direct_envelope, 0.1, 0.2, 2, true)) == 0
        direct_envelope.last_boundary_valid = false
        GC.gc()
        @test @allocated(PDMPSamplersRBridge._omrf_direct_hcv_fill_cell_cached!(
            direct_envelope, 0.3, 0.4, 3, true)) == 0
        direct_envelope.cell_count = 3
        @test direct_envelope.context.counts.direct_path_preparations >=
            first_path_preparations
        candidate = copy(direct_state2)
        move_forward_time!(candidate, 0.15, direct_flow)
        PDMPSamplers.screening_residual_bound(
            direct_envelope, candidate, direct_flow, 0.0)
        direct_cv = SubsampledControlVariate(x -> x, direct_oracle,
            direct_envelope, anchor, 1)
        draw_rng = Xoshiro(123)
        PDMPSamplers.draw_subset!(draw_rng, direct_cv, 0.2,
            direct_envelope.scales[1])
        GC.gc()
        @test @allocated(PDMPSamplers.draw_subset!(draw_rng, direct_cv, 0.2,
            direct_envelope.scales[1])) == 0

        direct_model_context, direct_model_names =
            PDMPSamplersRBridge._new_stan_subsampling_context(
                lib, full, prior, 1)
        direct_model, _, _ =
            PDMPSamplersRBridge._build_stan_subsampling_model(
                direct_model_context, direct_model_names,
                Dict{String,Any}("n_observations" => 3,
                    "subsample_size" => 1,
                    "residual_envelope" => direct_spec),
                anchor, nothing, falses(9); anchor_capacity=1)
        # The deterministic HCV likelihood product is prepared once per
        # physical segment.  Its affine/harmonic evaluation must reproduce
        # the ordinary dense H*(x-anchor) path for every supported dynamic.
        for cached_flow in (BouncyParticle(9, 0.0), ZigZag(9), Boomerang(9))
            cached_state = PDMPState(0.0,
                SkeletonPoint(copy(position), copy(velocity)))
            PDMPSamplers.initialize_flow_state!(cached_state, cached_flow)
            cached_envelope = direct_model.grad.envelope
            PDMPSamplersRBridge._omrf_direct_prepare_path!(
                cached_envelope, cached_state, cached_flow)
            cached_candidate = copy(cached_state)
            move_forward_time!(cached_candidate, 0.071, cached_flow)
            dense_gradient = zeros(9)
            cached_gradient = zeros(9)
            direct_model.grad.deterministic_gradient!(
                dense_gradient, cached_candidate.ξ.x)
            PDMPSamplers.deterministic_gradient!(cached_gradient,
                direct_model.grad, cached_candidate, cached_flow)
            @test cached_gradient ≈ dense_gradient atol=5e-12 rtol=5e-12
            if cached_flow isa Union{BouncyParticle,PDMPSamplers.AnyBoomerang}
                corrected_dense = copy(dense_gradient)
                correction_cache = (; z=zeros(9))
                PDMPSamplers.correct_gradient!(corrected_dense,
                    cached_candidate.ξ.x, cached_candidate.ξ.θ,
                    cached_flow, correction_cache)
                deferred = PDMPSamplers.subsampling_deterministic_signed_rate(
                    direct_model.grad.residual_oracle, direct_model.grad,
                    cached_candidate, cached_flow)
                @test deferred.active
                signed = deferred.signed
                @test isapprox(signed,
                    dot(cached_candidate.ξ.θ, corrected_dense);
                    atol=5e-12, rtol=5e-12)
                marked_residual = zeros(9)
                selected_mark = [1]
                direct_model.grad.residual_oracle(marked_residual,
                    cached_candidate.ξ.x, selected_mark, anchor)
                scale = 3.0
                deferred_rate =
                    PDMPSamplers.subsampling_deferred_candidate_rate(
                        direct_model.grad.residual_oracle, cached_candidate,
                        marked_residual, scale, cached_flow, signed,
                        selected_mark)
                axpy!(scale, marked_residual, corrected_dense)
                @test isapprox(deferred_rate, PDMPSamplers.λ(
                    cached_candidate, corrected_dense, cached_flow);
                    atol=5e-12, rtol=5e-12)
            end
        end
        direct_event_rng = Xoshiro(0xd1ec7)
        direct_event_flow = BouncyParticle(9, 0.0)
        direct_event_state, direct_event_model, direct_event_alg,
                direct_event_cache, direct_event_stats =
            PDMPSamplers.initialize_state(direct_event_rng,
                direct_event_flow, direct_model,
                GridThinningStrategy(N=2, N_min=2, t_max=0.2,
                    bound_violation=:throw), 0.0,
                SkeletonPoint(copy(position), copy(velocity)))
        direct_event_alg.schedule_frozen[] = true
        for _ in 1:4
            omrf_allocation_next_event!(direct_event_rng,
                direct_event_model, direct_event_flow, direct_event_alg,
                direct_event_state, direct_event_cache,
                direct_event_stats)
        end
        GC.gc()
        @test @allocated(omrf_allocation_next_event!(direct_event_rng,
            direct_event_model, direct_event_flow, direct_event_alg,
            direct_event_state, direct_event_cache,
            direct_event_stats)) == 0
    end

    # Exact selected refinement: the subset bound evaluates the selected
    # person's residual once, caches it for the candidate rate, and the
    # cached residual equals a direct oracle evaluation.
    exact_spec = copy(analytic_prior_spec)
    exact_spec["exact_selected_refinement"] = true
    exact_spec["weights"] = zeros(1, 3)
    exact_ctx, exact_names =
        PDMPSamplersRBridge._new_stan_subsampling_context(
            lib, full, prior, 1)
    exact_model, _, _ =
        PDMPSamplersRBridge._build_stan_subsampling_model(
            exact_ctx, exact_names,
            Dict{String,Any}(
                "n_observations" => 3,
                "subsample_size" => 1,
                "residual_envelope" => exact_spec),
            anchor, nothing, falses(9))
    exact_oracle = exact_model.grad.residual_oracle
    @test exact_oracle.exact_selected_refinement
    exact_state = PDMPState(0.0,
        SkeletonPoint(copy(position), fill(1.0, 9)))
    exact_gradient = zeros(9)
    compute_gradient!(exact_model.grad, position, exact_gradient)
    exact_flow = ZigZag(9)
    exact_D = PDMPSamplers.λ(exact_state, exact_gradient, exact_flow)
    exact_M = exact_D + 1e6
    exact_subset = [1]
    exact_scale = 3.0
    exact_bound = PDMPSamplers.subsampling_subset_bound(
        exact_oracle, exact_state, exact_gradient, exact_flow,
        exact_D, exact_M, exact_subset, exact_scale)
    cached_residual = zeros(9)
    @test @allocated(PDMPSamplers.subsampling_subset_bound(
        exact_oracle, exact_state, exact_gradient, exact_flow,
        exact_D, exact_M, exact_subset, exact_scale)) == 0
    @test @allocated(PDMPSamplers.subsampling_cached_residual!(
        exact_oracle, cached_residual, exact_state, exact_subset, anchor)) == 0
    PDMPSamplers.subsampling_subset_bound(
        exact_oracle, exact_state, exact_gradient, exact_flow,
        exact_D, exact_M, exact_subset, exact_scale)
    @test PDMPSamplers.subsampling_cached_residual!(
        exact_oracle, cached_residual, exact_state, exact_subset, anchor)
    @test !PDMPSamplers.subsampling_cached_residual!(
        exact_oracle, cached_residual, exact_state, exact_subset, anchor)
    PDMPSamplers.begin_subsampling_mark!(exact_oracle)
    @test !exact_oracle.cached_residual_valid
    exact_oracle(cached_residual, exact_state.ξ.x, exact_subset, anchor)
    expected_exact_gradient = exact_gradient + exact_scale .* cached_residual
    candidate_gradient = copy(exact_gradient)
    exact_rate = PDMPSamplers.subsampling_candidate_rate!(
        exact_oracle, exact_state, candidate_gradient, cached_residual,
        exact_scale, exact_flow, exact_D, exact_subset)
    @test candidate_gradient ≈ expected_exact_gradient atol=2e-14 rtol=2e-14
    @test exact_rate <= exact_bound + 2e-12
    exact_residual = zeros(9)
    exact_oracle(exact_residual, position, exact_subset, anchor)
    @test cached_residual ≈ exact_residual atol=2e-14 rtol=2e-14
    @test !exact_oracle.cached_residual_valid
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
        trace, _ = only(pdmp_sample(
            initial, BouncyParticle(9, 0.4), model,
            GridThinningStrategy(N=12, t_max=0.5, bound_violation=:throw),
            0.0, 1_200.0; progress=false, seed))
        trace
    end
    subsampling_runs = map(seeds) do seed
        prepared = PDMPSamplersRBridge.prepare_stan_subsampling(
            lib, full, prior, subsampling, 1)
        model, subsampling_ctx, _ = PDMPSamplersRBridge._build_stan_subsampling_model(
            only(prepared.contexts), prepared.unc_names, subsampling,
            anchor, nothing, falses(9))
        trace, stats = only(pdmp_sample(
            initial, BouncyParticle(9, 0.4), model,
            GridThinningStrategy(N=12, t_max=0.5, bound_violation=:throw),
            0.0, 1_200.0; progress=false, seed,
            statistic_counter=PDMPSamplers.DevelStatisticCounter))
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
        "type" => "dense_gaussian",
        "coef" => ["beta"],
        "mean" => [0.0],
        "cov" => fill(4.0, 1, 1))
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
        "type" => "dense_gaussian",
        "coef" => ["interaction"],
        "mean" => [0.0],
        "cov" => fill(prior_sd^2, 1, 1))
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

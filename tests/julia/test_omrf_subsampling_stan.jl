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
    velocity = collect(range(-1.1, 0.9; length=d))
    factor_flows = (
            ZigZag(d), BouncyParticle(d, 0.1), Boomerang(d),
            AdaptiveBoomerang(d; λref=0.2))
    for factor_flow in factor_flows
        factor_state = PDMPState(
            0.0, SkeletonPoint(copy(position), copy(velocity)))
        PDMPSamplers.initialize_flow_state!(factor_state, factor_flow)
        for t in (0.0, 0.07, 0.31)
            candidate_state = deepcopy(factor_state)
            PDMPSamplers.move_forward_time!(candidate_state, t, factor_flow)
            point = copy(candidate_state.ξ.x)
            for factor in 1:(size(X, 1) * size(X, 2))
                PDMPSamplersRBridge.omrf_residual!(
                    factor_context, factor_residual, point, [factor])
                node = cld(factor, size(X, 1))
                person = mod1(factor, size(X, 1))
                outside = setdiff(1:d, supports[node])
                @test factor_residual[outside] == zeros(length(outside))
                pattern_bound =
                    PDMPSamplersRBridge._omrf_person_node_pattern_bound_at(
                        factor_context, anchor, factor_state, factor_flow,
                        person, node, t)
                exact_rate = factor_flow isa ZigZag ?
                    sum(abs.(candidate_state.ξ.θ .* factor_residual)) :
                    abs(dot(candidate_state.ξ.θ, factor_residual))
                @test exact_rate <= pattern_bound + 2e-13
            end
        end
    end
    expected_factor_calls = (1 + length(factor_flows) * 3) *
        size(X, 1) * size(X, 2)
    @test factor_context.counts.analytic_factors == expected_factor_calls
    @test factor_context.counts.analytic_node_conditionals == expected_factor_calls

    local_spec = copy(factor_spec)
    local_candidates = hcat(anchor,
        anchor .+ collect(range(-0.35, 0.28; length=d)),
        anchor .+ collect(range(0.24, -0.31; length=d)))
    local_counts = PDMPSamplersRBridge.StanSubsamplingCallCounts()
    local_cv = PDMPSamplersRBridge._prepare_omrf_factor_local_cv(
        local_spec, names, local_candidates, local_counts)
    @test size(local_cv.anchor_logits) ==
        (maximum(seen) - 1, length(seen), size(X, 1))
    @test size(local_cv.source_indices) == (length(seen), size(X, 1))
    @test all(index -> 1 <= index <= size(local_candidates, 2),
        local_cv.source_indices)
    # Preparation uses the native factor algebra once.  In particular it
    # never asks BridgeStan for a full-data gradient.
    @test local_counts.full_gradient == 0
    # One scan initializes the common anchor cache and one constructs all
    # factor-local deterministic baselines.
    @test local_counts.native_full_data_scans == 2
    local_baseline_and_residual = copy(local_cv.full_anchor_likelihood)
    local_factor_residual = zeros(d)
    local_factors = collect(1:(size(X, 1) * size(X, 2)))
    local_cv(local_factor_residual, position, local_factors)
    local_baseline_and_residual .+= local_factor_residual
    exact_likelihood = zeros(d)
    for person in axes(X, 1)
        exact_likelihood[1:(no_thresholds + no_edges)] .+=
            finite_difference_person_neggrad_variable(
                position[1:(no_thresholds + no_edges)],
                view(X, person, :), seen)
    end
    @test local_baseline_and_residual ≈ exact_likelihood atol=2e-9

    local_velocity = collect(range(-0.7, 0.8; length=d))
    for factor in local_factors
        node = cld(factor, size(X, 1))
        person = mod1(factor, size(X, 1))
        local_cv(local_factor_residual, position, [factor])
        exact_rate = abs(dot(local_velocity, local_factor_residual))
        bound = PDMPSamplersRBridge._omrf_factor_local_bound(
            local_cv, position, local_velocity, person, node)
        @test exact_rate <= bound + 2e-13
    end

    local_bank_counts = PDMPSamplersRBridge.StanSubsamplingCallCounts()
    local_bank = PDMPSamplersRBridge._prepare_omrf_node_local_bank(
        local_spec, names, local_candidates, local_bank_counts)
    # Each candidate performs one cache scan and one exact node-baseline scan.
    @test local_bank_counts.native_full_data_scans ==
        2 * size(local_candidates, 2)
    @test size(local_bank.node_baselines) ==
        (d, length(seen), size(local_candidates, 2))
    mixed_sources = [1, 2, 3]
    mixed_total = zeros(d)
    PDMPSamplersRBridge._omrf_node_local_baseline!(
        mixed_total, local_bank, mixed_sources)
    mixed_residual = zeros(d)
    PDMPSamplersRBridge._omrf_node_local_residual!(local_bank,
        mixed_residual, position, collect(axes(X, 1)), mixed_sources)
    @test mixed_total + mixed_residual ≈ exact_likelihood atol=2e-9
    # The same identity holds under person-node subsampling and after an
    # arbitrary no-data recombination of already prepared node references.
    reverse_sources = reverse(mixed_sources)
    PDMPSamplersRBridge._omrf_node_local_baseline!(
        mixed_total, local_bank, reverse_sources)
    PDMPSamplersRBridge._omrf_node_local_residual!(local_bank,
        mixed_residual, position, local_factors, reverse_sources;
        factorization=:person_node)
    @test mixed_total + mixed_residual ≈ exact_likelihood atol=2e-9
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

    node_full_model = PDMPSamplersRBridge.build_omrf_full_model(
        spec, names, d)
    node_full_gradient = zeros(d)
    node_full_gradient_plus = zeros(d)
    node_full_gradient_minus = zeros(d)
    compute_gradient!(node_full_model.grad, position, node_full_gradient)
    compute_gradient!(node_full_model.grad,
        position .+ step .* direction, node_full_gradient_plus)
    compute_gradient!(node_full_model.grad,
        position .- step .* direction, node_full_gradient_minus)
    node_full_hvp = copy(node_full_model.hvp(position, direction))
    @test node_full_hvp ≈
        (node_full_gradient_plus - node_full_gradient_minus) / (2step) atol=3e-8 rtol=3e-8
    node_full_hvp_two = copy(node_full_model.hvp(position, direction_two))
    @test dot(direction_two, node_full_hvp) ≈
        dot(direction, node_full_hvp_two) atol=2e-8 rtol=2e-8

    # The full-data warmup model for a dependent slab is the native OMRF
    # model wrapped in the active-slab correction.  At the all-free state it
    # must equal the native all-active analytic target.  After freezing one
    # interaction, both that interaction and the node-scale coordinates must
    # change; otherwise warmup is adapting the wrong target.
    can_stick = falses(d)
    can_stick[(no_thresholds + 1):(no_thresholds + no_edges)] .= true
    slab_spec = Dict{String,Any}(
        "type" => "loglinear_gaussian_scale",
        "log_base_scales" => fill(log(0.7), no_edges),
        "logscale" => ["global_scale", "node_scale"],
        "logscale_design" => Dict{String,Any}(
            "storage" => "dense", "dims" => [no_edges, 1 + length(seen)],
            "x" => vec([
                1.0 0.5 0.5 0.0;
                1.0 0.5 0.0 0.5;
                1.0 0.0 0.5 0.5])),
        "coef" => "interactions_0")
    # Pass the R-facing specifications here, as the bridge does in the
    # warmup construction; _build_dependent_slab_model resolves both into
    # Julia providers.
    model_prior = Dict{String,Any}("a" => 1.0, "b" => 4.0)
    wrapped_model = PDMPSamplersRBridge._build_dependent_slab_model(
        node_full_model, model_prior, slab_spec, can_stick, names)
    native_gradient = zeros(d)
    wrapped_gradient = zeros(d)
    compute_gradient!(node_full_model.grad, position, native_gradient)
    PDMPSamplers.set_active_set!(wrapped_model, trues(d))
    compute_gradient!(wrapped_model.grad, position, wrapped_gradient)
    @test wrapped_gradient ≈ native_gradient atol=3e-10 rtol=3e-10

    frozen = trues(d)
    frozen[no_thresholds + 1] = false
    PDMPSamplers.set_active_set!(wrapped_model, frozen)
    compute_gradient!(wrapped_model.grad, position, wrapped_gradient)
    @test abs(wrapped_gradient[no_thresholds + 1] -
        native_gradient[no_thresholds + 1]) > 1e-8
    scale_difference = wrapped_gradient[(no_thresholds + no_edges + 1):d] -
        native_gradient[(no_thresholds + no_edges + 1):d]
    @test maximum(abs, scale_difference) > 1e-8
end

@testset "retained full-data invariant reports every work category" begin
    clean = [Dict{String,Any}(
        "anchor_main_preparations" => 0,
        "main_native_full_data_scans" => 0,
        "main_event_bound_observation_visits" => 0)]
    @test isnothing(PDMPSamplersRBridge._enforce_no_main_full_data!(
        [0], clean))
    dirty = [Dict{String,Any}(
        "anchor_main_preparations" => 2,
        "main_native_full_data_scans" => 3,
        "main_event_bound_observation_visits" => 17)]
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
    @test contains(message, "main_event_bound_observation_visits=17")
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
    for _ in 1:3
        compute_gradient!(full_native_model.grad, position, native_gradient)
    end
    # Julia 1.12 retains a 16-byte bookkeeping allocation in this native
    # counter path after the provider counter extension; the signed route
    # allocation checks below remain strict zero-byte checks.
    @test @allocated(compute_gradient!(
        full_native_model.grad, position, native_gradient)) <= 16
    full_native_model.hvp(position, direction)
    @test @allocated(full_native_model.hvp(position, direction)) == 0

    node_hcv_spec = copy(analytic_prior_spec)
    node_hcv_spec["use_hcv"] = true
    node_hcv_spec["hcv_scope"] = "node"
    node_hcv_spec["hcv_damping"] = 2.0
    node_hcv_spec["hcv_after_warmup"] = false
    node_hcv_oracle = PDMPSamplersRBridge.OMRFAnalyticSubsamplingOracle(
        analytic_prior_ctx, node_hcv_spec, analytic_prior_names, anchor)
    node_hcv_context = node_hcv_oracle.residual_context
    @test node_hcv_context.hcv_nodewise
    node_hcv_deterministic = zeros(9)
    node_hcv_residual = zeros(9)
    node_hcv_oracle(node_hcv_deterministic, position)
    node_hcv_oracle(node_hcv_residual, position, [1, 2, 3], anchor)
    @test node_hcv_deterministic + node_hcv_residual ≈ stan_full_gradient atol=3e-10 rtol=3e-10
    node_hcv_plus = zeros(9)
    node_hcv_minus = zeros(9)
    node_hcv_oracle(node_hcv_plus, position .+ step .* direction)
    node_hcv_oracle(node_hcv_minus, position .- step .* direction)
    node_hcv_hvp = zeros(9)
    node_hcv_oracle(node_hcv_hvp, position, direction)
    @test node_hcv_hvp ≈
        (node_hcv_plus - node_hcv_minus) / (2step) atol=2e-7 rtol=2e-7
    @test @allocated(node_hcv_oracle(node_hcv_hvp, position, direction)) == 0

    @testset "compact person-node HCV certificate" begin
        X_compact = [0 1 2; 2 0 1; 1 2 0]
        N_compact, P_compact = size(X_compact)
        compact_weights = zeros(4P_compact, N_compact)
        for node in 1:P_compact, person in 1:N_compact
            r = sqrt(sum(abs2, X_compact[person,
                [j for j in 1:P_compact if j != node]]))
            for power in 0:3
                compact_weights[power * P_compact + node, person] = r^power
            end
        end
        compact_spec = copy(analytic_prior_spec)
        compact_spec["use_hcv"] = true
        compact_spec["hcv_scope"] = "global"
        compact_spec["hcv_damping"] = Inf
        compact_spec["hcv_after_warmup"] = false
        compact_spec["factorization"] = "person_node"
        compact_spec["bound_type"] = "norm_local_hcv_compact"
        compact_spec["weights"] = zeros(1, N_compact * P_compact)
        compact_spec["norm_weights"] = compact_weights
        compact_oracle = PDMPSamplersRBridge.OMRFAnalyticSubsamplingOracle(
            analytic_prior_ctx, compact_spec, analytic_prior_names, anchor)
        compact_envelope = PDMPSamplersRBridge._build_stan_residual_envelope(
            compact_spec, compact_oracle.anchor, compact_oracle)
        compact_scales = compact_envelope.component_scales!
        @test compact_scales isa
            PDMPSamplersRBridge.OMRFCompactHCVNormLocalScales

        compact_deterministic = zeros(9)
        compact_sum = zeros(9)
        compact_residual = zeros(9)
        compact_oracle(compact_deterministic, position)
        for factor in 1:(N_compact * P_compact)
            compact_oracle(compact_residual, position, [factor], anchor)
            compact_sum .+= compact_residual
        end
        @test compact_deterministic + compact_sum ≈ stan_full_gradient atol=3e-10 rtol=3e-10

        exhaustive_mean = zeros(9)
        exhaustive_subset = zeros(Int, 2)
        subset_count = 0
        for first_factor in 1:(N_compact * P_compact - 1),
                second_factor in (first_factor + 1):(N_compact * P_compact)
            exhaustive_subset[1] = first_factor
            exhaustive_subset[2] = second_factor
            compact_oracle(compact_residual, position,
                exhaustive_subset, anchor)
            exhaustive_mean .+=
                (N_compact * P_compact / length(exhaustive_subset)) .*
                compact_residual
            subset_count += 1
        end
        exhaustive_mean ./= subset_count
        @test exhaustive_mean ≈ compact_sum atol=3e-10 rtol=3e-10

        compact_velocity = collect(range(-0.8, 0.9; length=9))
        compact_state = PDMPState(0.0,
            SkeletonPoint(copy(position), copy(compact_velocity)))
        compact_flow = BouncyParticle(9)
        PDMPSamplers.initialize_flow_state!(compact_state, compact_flow)
        compact_scales(compact_envelope.scales, compact_state,
            compact_flow, 0.0)
        for factor in 1:(N_compact * P_compact)
            compact_oracle(compact_residual, position, [factor], anchor)
            exact = abs(dot(compact_velocity, compact_residual))
            bound = PDMPSamplers.observation_residual_bound(
                compact_envelope, factor)
            @test exact <= bound + 2e-12
        end

        left, right = 0.0, 0.3
        compact_scales(compact_envelope.cell_scales, compact_state,
            compact_flow, left, right)
        cell_scales = copy(compact_envelope.cell_scales)
        point_scales = similar(cell_scales)
        for t in (left, right, 0.07, 0.19, 0.299)
            compact_scales(point_scales, compact_state, compact_flow, t)
            @test all(point_scales .<= cell_scales .+ 2e-12)
        end
        compact_scales(compact_envelope.scales, compact_state,
            compact_flow, 0.0)
        compact_scales(compact_envelope.scales, compact_state,
            compact_flow, 0.0)
        @test @allocated(compact_scales(compact_envelope.scales,
            compact_state, compact_flow, 0.0)) == 0
        compact_subset = [1]
        compact_oracle(compact_residual, position, compact_subset, anchor)
        @test @allocated(compact_oracle(
            compact_residual, position, compact_subset, anchor)) == 0

        compact_model_context, compact_model_names =
            PDMPSamplersRBridge._new_stan_subsampling_context(
                lib, full, prior, 3)
        compact_model, _, _, _, compact_manager =
            PDMPSamplersRBridge._build_stan_subsampling_model(
                compact_model_context, compact_model_names,
                Dict{String,Any}("n_observations" => 9,
                    "subsample_size" => 3,
                    "residual_envelope" => compact_spec),
                anchor, nothing, falses(9); anchor_capacity=1)
        original_hessian = compact_manager.oracle.residual_context.likelihood_hessian
        compact_manager.main_refresh_distance = 0.0
        for _ in 1:2
            PDMPSamplersRBridge._select_omrf_anchor!(compact_manager,
                compact_model.grad,
                compact_manager.oracle.anchor .+ 0.02; phase=:main)
        end
        @test compact_manager.oracle.residual_context.likelihood_hessian ===
            original_hessian
        PDMPSamplersRBridge._omrf_likelihood_hessian!(original_hessian,
            compact_manager.oracle.residual_context)
        @test @allocated(PDMPSamplersRBridge._omrf_likelihood_hessian!(
            original_hessian,
            compact_manager.oracle.residual_context)) == 0
        compact_refresh_position = compact_manager.oracle.anchor .+ 0.02
        compact_refresh_fixture = (; manager=compact_manager,
            cv=compact_model.grad, position=compact_refresh_position)
        # The ordinary owner update remains pinned to 112 bytes below. The
        # compact-HCV refresh additionally crosses the optional dense-Hessian
        # field's union-typed dispatch, which Julia 1.12 boxes in 32 bytes;
        # the Hessian, probabilities, envelope, and anchor arrays are reused.
        @test @allocated(omrf_allocation_anchor_select!(
            compact_refresh_fixture)) == 144
        compact_event_rng = Xoshiro(0xc04ac7)
        compact_event_flow = BouncyParticle(9, 0.0)
        compact_event_state, compact_event_model, compact_event_alg,
                compact_event_cache, compact_event_stats =
            PDMPSamplers.initialize_state(compact_event_rng,
                compact_event_flow, compact_model,
                GridThinningStrategy(N=2, N_min=2, t_max=0.25,
                    lazy=false, bound_violation=:throw), 0.0,
                SkeletonPoint(copy(position), copy(compact_velocity)))
        compact_event_alg.schedule_frozen[] = true
        for _ in 1:3
            omrf_allocation_next_event!(compact_event_rng,
                compact_event_model, compact_event_flow, compact_event_alg,
                compact_event_state, compact_event_cache,
                compact_event_stats)
        end
        GC.gc()
        @test @allocated(omrf_allocation_next_event!(compact_event_rng,
            compact_event_model, compact_event_flow, compact_event_alg,
            compact_event_state, compact_event_cache,
            compact_event_stats)) == 0
    end

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

    @testset "direct person-node HCV cell certificate" begin
        direct_spec = copy(analytic_prior_spec)
        direct_spec["use_hcv"] = true
        direct_spec["hcv_scope"] = "global"
        direct_spec["hcv_damping"] = Inf
        direct_spec["hcv_after_warmup"] = false
        direct_spec["factorization"] = "person_node"
        direct_spec["bound_type"] = "direct_person_node_hcv"
        direct_spec["exact_selected_refinement"] = true
        direct_spec["weights"] = zeros(1, 9)
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
        direct_residual = zeros(9)
        maximum_violation = 0.0
        for factor in 1:9, t in range(left, right; length=65)
            point = position .+ t .* velocity
            direct_oracle(direct_residual, point, [factor], anchor)
            exact = abs(dot(velocity, direct_residual))
            bound = PDMPSamplers.observation_residual_bound(
                direct_envelope, factor)
            maximum_violation = max(maximum_violation, exact - bound)
        end
        @test maximum_violation <= 2e-12

        # Hybrid factor-roof/person-gradient HCV: retain the certified direct
        # person-node table, but sample one person and evaluate the sum of all
        # of that person's node residuals.  Its aggregate roof is exactly the
        # factor roof, and the size-biased person law is C_i/sum(C_i).
        hybrid_spec = copy(direct_spec)
        hybrid_spec["factorization"] = "person"
        hybrid_spec["bound_type"] = "direct_factor_roof_person_hcv"
        hybrid_spec["weights"] = zeros(1, 3)
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
        # A single HCV table on the entire physical flight remains a valid
        # (possibly looser) roof for every contained deterministic grid cell.
        horizon_hcv_spec = copy(signed_spec)
        horizon_hcv_spec["direct_cv_cell_scope"] = "event_horizon"
        horizon_hcv = PDMPSamplersRBridge._build_stan_residual_envelope(
            horizon_hcv_spec, hybrid_oracle.anchor, hybrid_oracle)
        PDMPSamplers.prepare_residual_horizon!(
            horizon_hcv, direct_state, direct_flow, right)
        horizon_hcv_mass = horizon_hcv.cell_mass[1]
        constructions_before = horizon_hcv.context.counts.
            direct_residual_table_constructions
        for (cell_left, cell_right) in ((0.0, left), (left, right))
            PDMPSamplers.component_cell_scales!(horizon_hcv.cell_scales,
                horizon_hcv, direct_state, direct_flow, cell_left, cell_right)
            @test horizon_hcv.cell_mass[1] == horizon_hcv_mass
            for person in 1:3, t in range(cell_left, cell_right; length=65)
                hybrid_oracle(direct_residual,
                    position .+ t .* velocity, [person], anchor)
                certified = PDMPSamplersRBridge._omrf_direct_person_cell_bound(
                    horizon_hcv, person)
                @test abs(dot(velocity, direct_residual)) <=
                    certified + 2e-11 * (1 + certified)
            end
        end
        @test horizon_hcv.context.counts.direct_residual_table_constructions ==
            constructions_before
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
                factor = (node - 1) * 3 + person
                direct_oracle(direct_residual, position, [factor], anchor)
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
        extreme_oracle(extreme_residual, extreme_position, [1],
            extreme_anchor)
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
            direct_oracle(direct_residual, point, [factor], anchor)
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
            direct_envelope, anchor, 2)
        draw_rng = Xoshiro(123)
        PDMPSamplers.draw_subset!(draw_rng, direct_cv, 0.2,
            direct_envelope.scales[1])
        GC.gc()
        @test @allocated(PDMPSamplers.draw_subset!(draw_rng, direct_cv, 0.2,
            direct_envelope.scales[1])) == 0

        # Exhaustive mixture algebra for every size-two subset. The existing
        # distinguished-factor plus conditional base-subset law yields the
        # selected probability base_prob*(D + (N/m)sum b_i)/(D+B).
        Nf, m, D = 9, 2, 0.2
        B = direct_envelope.scales[1]
        base_probability = 1 / binomial(Nf, m)
        probability_sum = 0.0
        for first in 1:(Nf - 1), second in (first + 1):Nf
            selected = D + (Nf / m) * (
                PDMPSamplers.observation_residual_bound(
                    direct_envelope, first) +
                PDMPSamplers.observation_residual_bound(
                    direct_envelope, second))
            probability_sum += base_probability * selected / (D + B)
        end
        @test probability_sum ≈ 1.0 atol=2e-14

        pairs = [(first, second) for first in 1:(Nf - 1)
            for second in (first + 1):Nf]
        expected_probability = zeros(length(pairs))
        for index in eachindex(pairs)
            first, second = pairs[index]
            selected = D + (Nf / m) * (
                PDMPSamplers.observation_residual_bound(
                    direct_envelope, first) +
                PDMPSamplers.observation_residual_bound(
                    direct_envelope, second))
            expected_probability[index] =
                base_probability * selected / (D + B)
        end
        empirical_counts = zeros(Int, length(pairs))
        empirical_draws = 100_000
        empirical_rng = Xoshiro(0xe71c)
        for _ in 1:empirical_draws
            PDMPSamplers.draw_subset!(
                empirical_rng, direct_cv, D, B)
            selected_pair = minmax(direct_cv.subset[1], direct_cv.subset[2])
            empirical_counts[findfirst(==(selected_pair), pairs)] += 1
        end
        empirical_probability = empirical_counts ./ empirical_draws
        @test all(abs.(empirical_probability .- expected_probability) .<=
            6 .* sqrt.(expected_probability .* (1 .- expected_probability) ./
                empirical_draws) .+ 1 / empirical_draws)

        direct_model_context, direct_model_names =
            PDMPSamplersRBridge._new_stan_subsampling_context(
                lib, full, prior, 3)
        direct_model, _, _ =
            PDMPSamplersRBridge._build_stan_subsampling_model(
                direct_model_context, direct_model_names,
                Dict{String,Any}("n_observations" => 9,
                    "subsample_size" => 3,
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
                selected_mark = [1, 4, 7]
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
                    lazy=false, bound_violation=:throw), 0.0,
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

    @testset "direct person-node first-order CV cell certificate" begin
        cv_spec = copy(analytic_prior_spec)
        cv_spec["use_hcv"] = false
        cv_spec["factorization"] = "person_node"
        cv_spec["bound_type"] = "direct_person_node_cv"
        cv_spec["weights"] = zeros(1, 9)
        cv_oracle = PDMPSamplersRBridge.OMRFAnalyticSubsamplingOracle(
            analytic_prior_ctx, cv_spec, analytic_prior_names, anchor)
        cv_envelope = PDMPSamplersRBridge._build_stan_residual_envelope(
            cv_spec, cv_oracle.anchor, cv_oracle)
        @test cv_envelope isa
            PDMPSamplersRBridge.OMRFDirectPersonNodeHCVEnvelope
        @test !cv_envelope.hcv_remainder
        @test cv_oracle.residual_context.likelihood_hessian === nothing

        norm_spec = copy(cv_spec)
        norm_spec["bound_type"] = "node_local_norm_low_rank"
        norm_oracle = PDMPSamplersRBridge.OMRFAnalyticSubsamplingOracle(
            analytic_prior_ctx, norm_spec, analytic_prior_names, anchor)
        cv_deterministic = zeros(9)
        norm_deterministic = zeros(9)
        cv_residual = zeros(9)
        norm_residual = zeros(9)
        cv_oracle(cv_deterministic, position)
        norm_oracle(norm_deterministic, position)
        @test cv_deterministic == norm_deterministic
        for subset in ([1], [4], [9], [1, 5], [2, 6, 8])
            cv_oracle(cv_residual, position, subset, anchor)
            norm_oracle(norm_residual, position, subset, anchor)
            @test cv_residual == norm_residual
            cv_candidate = cv_deterministic .+
                (9 / length(subset)) .* cv_residual
            norm_candidate = norm_deterministic .+
                (9 / length(subset)) .* norm_residual
            @test cv_candidate == norm_candidate
        end

        # Complete first-order decomposition and exhaustive size-two
        # unbiasedness use the same residual estimator as norm-local.
        cv_complete = zeros(9)
        one_factor = [1]
        for factor in 1:9
            one_factor[1] = factor
            cv_oracle(cv_residual, position, one_factor, anchor)
            cv_complete .+= cv_residual
        end
        @test cv_deterministic + cv_complete ≈ stan_full_gradient atol=3e-10 rtol=3e-10
        exhaustive = zeros(9)
        pairs = 0
        subset = zeros(Int, 2)
        for first in 1:8, second in (first + 1):9
            subset[1] = first
            subset[2] = second
            cv_oracle(cv_residual, position, subset, anchor)
            exhaustive .+= (9 / 2) .* cv_residual
            pairs += 1
        end
        exhaustive ./= pairs
        @test exhaustive ≈ cv_complete atol=3e-10 rtol=3e-10

        cv_velocity = collect(range(-0.8, 0.9; length=9))
        cv_state = PDMPState(0.0,
            SkeletonPoint(copy(position), copy(cv_velocity)))
        cv_flow = BouncyParticle(9)
        PDMPSamplers.initialize_flow_state!(cv_state, cv_flow)
        PDMPSamplers.component_cell_scales!(cv_envelope.cell_scales,
            cv_envelope, cv_state, cv_flow, -0.2, 0.3)
        maximum_violation = 0.0
        for factor in 1:9
            previous = -Inf
            for t in range(-0.2, 0.3; length=129)
                point = position .+ t .* cv_velocity
                cv_oracle(cv_residual, point, [factor], anchor)
                h = dot(cv_velocity, cv_residual)
                @test h >= previous - 3e-12 * (1 + abs(h))
                previous = h
                bound = PDMPSamplers.observation_residual_bound(
                    cv_envelope, factor)
                maximum_violation = max(maximum_violation, abs(h) - bound)
            end
        end
        @test maximum_violation <= 2e-12

        # Categorical adversaries include equal/dominant logits, zero
        # velocity/displacement, sign-changing residuals and anchor crossings.
        categorical_cases = (
            (zeros(3), zeros(3), zeros(3), 0.5),
            ([50.0, -50.0, -50.0], [-1.0, 0.5, 2.0],
                [3.0, -2.0, 1.0], 0.7),
            (zeros(3), [-0.3, 0.2, 0.1], [1.0, -1.0, 0.5], 0.6),
            ([1.0, -1.0, 0.0], zeros(3), [2.0, -3.0, 1.0], 0.4),
        )
        for (logits, delta, z, width) in categorical_cases
            bound = PDMPSamplersRBridge._omrf_direct_categorical_cv_cell_bound(
                logits, delta, z, width)
            previous = -Inf
            for t in range(0.0, width; length=257)
                h = PDMPSamplersRBridge._omrf_categorical_cv_residual(
                    logits, delta .+ t .* z, z)
                @test h >= previous - 3e-12 * (1 + abs(h))
                @test abs(h) <= bound + 3e-12 * (1 + bound)
                previous = h
            end
        end

        # Event-horizon scope evaluates the same monotone endpoint
        # certificate once on [0,H] and reuses its cumulative factor table for
        # every deterministic grid cell contained in that interval.
        horizon_spec = copy(cv_spec)
        horizon_spec["direct_cv_cell_scope"] = "event_horizon"
        horizon_oracle = PDMPSamplersRBridge.OMRFAnalyticSubsamplingOracle(
            analytic_prior_ctx, horizon_spec, analytic_prior_names, anchor)
        horizon_envelope = PDMPSamplersRBridge._build_stan_residual_envelope(
            horizon_spec, horizon_oracle.anchor, horizon_oracle)
        @test horizon_envelope.event_horizon_scope
        horizon_state = StickyPDMPState(Ref(0.0),
            SkeletonPoint(copy(position), copy(cv_velocity)), trues(9))
        PDMPSamplers.prepare_residual_horizon!(
            horizon_envelope, horizon_state, cv_flow, 0.6)
        horizon_table = copy(view(horizon_envelope.cell_cumulative, :, 1))
        horizon_mass = horizon_envelope.cell_mass[1]
        constructions = horizon_envelope.context.counts.
            direct_residual_table_constructions
        for (left, right) in ((0.0, 0.1), (0.1, 0.2), (0.45, 0.6))
            PDMPSamplers.component_cell_scales!(
                horizon_envelope.cell_scales, horizon_envelope,
                horizon_state, cv_flow, left, right)
            @test horizon_envelope.cell_mass[1] == horizon_mass
            @test view(horizon_envelope.cell_cumulative, :, 1) == horizon_table
        end
        @test horizon_envelope.context.counts.
            direct_residual_table_constructions == constructions

        # The coarser table changes rejection probabilities, but not the
        # accepted marked intensity. Exhaustive size-two subset algebra gives
        # the same accepted mass for a local table and the shared table.
        PDMPSamplers.component_cell_scales!(cv_envelope.cell_scales,
            cv_envelope, horizon_state, cv_flow, 0.2, 0.4)
        PDMPSamplers.component_cell_scales!(horizon_envelope.cell_scales,
            horizon_envelope, horizon_state, cv_flow, 0.2, 0.4)
        marked_state = copy(horizon_state)
        PDMPSamplers.move_forward_time!(marked_state, 0.3, cv_flow)
        horizon_oracle(cv_deterministic, marked_state.ξ.x)
        marked_D = PDMPSamplers.λ(
            marked_state, cv_deterministic, cv_flow)
        function exhaustive_accepted_mass(envelope)
            B = envelope.cell_mass[envelope.active_cell]
            accepted_mass = 0.0
            base_probability = 1 / binomial(9, 2)
            subset = zeros(Int, 2)
            candidate_gradient = zeros(9)
            for first in 1:8, second in (first + 1):9
                subset[1] = first
                subset[2] = second
                horizon_oracle(cv_residual, marked_state.ξ.x, subset, anchor)
                copyto!(candidate_gradient, cv_deterministic)
                actual = PDMPSamplers.subsampling_candidate_rate!(
                    horizon_oracle, marked_state, candidate_gradient,
                    cv_residual, 9 / 2, cv_flow, marked_D, subset)
                selected_bound = marked_D + (9 / 2) * sum(factor ->
                    PDMPSamplers.observation_residual_bound(
                        envelope, factor), subset)
                selection_probability = base_probability * selected_bound /
                    (marked_D + B)
                accepted_mass += (marked_D + B) * selection_probability *
                    actual / selected_bound
            end
            return accepted_mass
        end
        @test exhaustive_accepted_mass(cv_envelope) ≈
            exhaustive_accepted_mass(horizon_envelope) atol=2e-12 rtol=2e-12

        # Dense interior and sign-crossing checks cover the complete shared
        # horizon, including zero-velocity and inactive-interaction factors.
        horizon_maximum_violation = 0.0
        for factor in 1:9, t in range(0.0, 0.6; length=257)
            point = position .+ t .* cv_velocity
            horizon_oracle(cv_residual, point, [factor], anchor)
            bound = factor == 1 ? horizon_table[1] :
                horizon_table[factor] - horizon_table[factor - 1] +
                4eps(max(abs(horizon_table[factor]),
                    abs(horizon_table[factor - 1]), 1.0))
            horizon_maximum_violation = max(horizon_maximum_violation,
                abs(dot(cv_velocity, cv_residual)) - bound)
        end
        @test horizon_maximum_violation <= 2e-12

        # A changed/extended horizon, reflected velocity, sticky active set,
        # anchor, or metric generation cannot retain the old table.
        constructions = horizon_envelope.context.counts.
            direct_residual_table_constructions
        PDMPSamplers.prepare_residual_horizon!(
            horizon_envelope, horizon_state, cv_flow, 0.6)
        @test horizon_envelope.context.counts.
            direct_residual_table_constructions == constructions
        PDMPSamplers.prepare_residual_horizon!(
            horizon_envelope, horizon_state, cv_flow, 0.7)
        constructions += 1
        @test horizon_envelope.context.counts.
            direct_residual_table_constructions == constructions
        horizon_state.ξ.θ[1] = -horizon_state.ξ.θ[1]
        PDMPSamplers.prepare_residual_horizon!(
            horizon_envelope, horizon_state, cv_flow, 0.7)
        constructions += 1
        @test horizon_envelope.context.counts.
            direct_residual_table_constructions == constructions
        active_horizon = copy(horizon_oracle.active_free)
        active_horizon[horizon_oracle.residual_context.interaction_indices[1]] = false
        PDMPSamplers.set_active_set!(horizon_oracle, active_horizon)
        PDMPSamplers.prepare_residual_horizon!(
            horizon_envelope, horizon_state, cv_flow, 0.7)
        constructions += 1
        @test horizon_envelope.context.counts.
            direct_residual_table_constructions == constructions
        PDMPSamplers.set_active_set!(horizon_oracle, trues(9))
        horizon_envelope.anchor[1] += 0.01
        PDMPSamplers.prepare_residual_horizon!(
            horizon_envelope, horizon_state, cv_flow, 0.7)
        constructions += 1
        @test horizon_envelope.context.counts.
            direct_residual_table_constructions == constructions
        horizon_envelope.anchor[1] -= 0.01
        metric_flow = PreconditionedBPS(9)
        PDMPSamplers.prepare_residual_horizon!(
            horizon_envelope, horizon_state, metric_flow, 0.7)
        constructions += 1
        metric_flow.metric.generation += UInt(1)
        PDMPSamplers.prepare_residual_horizon!(
            horizon_envelope, horizon_state, metric_flow, 0.7)
        constructions += 1
        @test horizon_envelope.context.counts.
            direct_residual_table_constructions == constructions

        # At least 20,001 random factor/horizon checks exercise the actual
        # person-node provider and arbitrary interior proposal times.
        horizon_rng = Xoshiro(0x48f120)
        for _ in 1:20_001
            random_position = anchor .+ 0.5 .* randn(horizon_rng, 9)
            random_velocity = randn(horizon_rng, 9)
            random_state = PDMPState(0.0,
                SkeletonPoint(random_position, random_velocity))
            H = 0.05 + 0.55rand(horizon_rng)
            PDMPSamplers.prepare_residual_horizon!(
                horizon_envelope, random_state, cv_flow, H)
            factor = rand(horizon_rng, 1:9)
            t = H * rand(horizon_rng)
            horizon_oracle(cv_residual,
                random_position .+ t .* random_velocity, [factor], anchor)
            PDMPSamplers.component_cell_scales!(
                horizon_envelope.cell_scales, horizon_envelope,
                random_state, cv_flow, 0.0, H)
            @test abs(dot(random_velocity, cv_residual)) <=
                PDMPSamplers.observation_residual_bound(
                    horizon_envelope, factor) + 3e-12
        end

        # Preparation, internal-cell lookup, horizon extension and marked
        # proposal primitives allocate nothing after exact prewarming.
        PDMPSamplers.prepare_residual_horizon!(
            horizon_envelope, horizon_state, cv_flow, 0.4)
        PDMPSamplers.component_cell_scales!(horizon_envelope.cell_scales,
            horizon_envelope, horizon_state, cv_flow, 0.1, 0.2)
        GC.gc()
        @test @allocated(PDMPSamplers.prepare_residual_horizon!(
            horizon_envelope, horizon_state, cv_flow, 0.4)) == 0
        @test @allocated(PDMPSamplers.component_cell_scales!(
            horizon_envelope.cell_scales, horizon_envelope,
            horizon_state, cv_flow, 0.2, 0.3)) == 0
        GC.gc()
        @test @allocated(PDMPSamplers.prepare_residual_horizon!(
            horizon_envelope, horizon_state, cv_flow, 0.5)) == 0

        # Matrix fields equal the original scalar incident-neighbour sums for
        # complete/smaller graphs, random/sparse/all-zero interactions, and
        # frozen coordinates.
        field_rng = Xoshiro(0xf1e1d)
        for retained_edges in (0:0, 1:2, 1:3), frozen in (false, true)
            field_state = StickyPDMPState(Ref(0.0),
                SkeletonPoint(copy(position), randn(field_rng, 9)), trues(9))
            interaction_indices = cv_oracle.residual_context.interaction_indices
            for (edge, index) in enumerate(interaction_indices)
                if edge in retained_edges
                    field_state.ξ.x[index] = randn(field_rng)
                    field_state.ξ.θ[index] = randn(field_rng)
                else
                    field_state.ξ.x[index] = anchor[index]
                    field_state.ξ.θ[index] = 0.0
                end
            end
            if frozen
                field_state.free[interaction_indices[1]] = false
            end
            PDMPSamplersRBridge._omrf_direct_prepare_linear_fields!(
                cv_envelope, field_state, cv_flow)
            context = cv_envelope.context
            for node in 1:3, person in 1:3
                scalar_delta = 0.0
                scalar_velocity = 0.0
                for k in eachindex(context.incident_edges[node])
                    covariate = context.X[person,
                        context.incident_neighbours[node][k]]
                    index = context.interaction_indices[
                        context.incident_edges[node][k]]
                    delta, velocity = PDMPSamplersRBridge._omrf_coordinate_at(
                        field_state, cv_flow, anchor, index, 0.0)
                    scalar_delta += covariate * delta
                    scalar_velocity += covariate * velocity
                end
                @test cv_envelope.field_delta[person, node] ≈ scalar_delta atol=2e-15
                @test cv_envelope.field_velocity[person, node] ≈ scalar_velocity atol=2e-15
            end
        end

        # More than 20,000 factor/cell checks use dense random interior points.
        check_rng = Xoshiro(0xd1ec7c)
        for _ in 1:20_001
            random_position = anchor .+ 0.5 .* randn(check_rng, 9)
            random_velocity = randn(check_rng, 9)
            random_state = PDMPState(0.0,
                SkeletonPoint(random_position, random_velocity))
            left = -0.2 + 0.4rand(check_rng)
            right = left + 0.4rand(check_rng)
            PDMPSamplers.component_cell_scales!(cv_envelope.cell_scales,
                cv_envelope, random_state, cv_flow, left, right)
            factor = rand(check_rng, 1:9)
            t = left + rand(check_rng) * (right - left)
            cv_oracle(cv_residual,
                random_position .+ t .* random_velocity, [factor], anchor)
            @test abs(dot(random_velocity, cv_residual)) <=
                PDMPSamplers.observation_residual_bound(
                    cv_envelope, factor) + 3e-12
        end

        # Adjacent reuse, noncontiguous fallback, active/velocity/anchor
        # invalidation and every production hot primitive stay allocation-free.
        allocation_state = PDMPState(1.0,
            SkeletonPoint(copy(position), copy(cv_velocity)))
        PDMPSamplers.component_cell_scales!(cv_envelope.cell_scales,
            cv_envelope, allocation_state, cv_flow, 0.0, 0.1)
        PDMPSamplers.component_cell_scales!(cv_envelope.cell_scales,
            cv_envelope, allocation_state, cv_flow, 0.1, 0.2)
        adjacent = copy(view(cv_envelope.cell_cumulative, :, 2))
        cv_envelope.last_boundary_valid = false
        PDMPSamplersRBridge._omrf_direct_hcv_fill_cell_cached!(
            cv_envelope, 0.1, 0.2, 3, false)
        @test adjacent == view(cv_envelope.cell_cumulative, :, 3)
        @test PDMPSamplersRBridge._omrf_direct_path_matches(
            cv_envelope, allocation_state)
        active = copy(cv_oracle.active_free)
        active[cv_oracle.residual_context.interaction_indices[1]] = false
        PDMPSamplers.set_active_set!(cv_oracle, active)
        @test !PDMPSamplersRBridge._omrf_direct_path_matches(
            cv_envelope, allocation_state)
        PDMPSamplers.set_active_set!(cv_oracle, trues(9))
        PDMPSamplersRBridge._omrf_direct_prepare_path!(
            cv_envelope, allocation_state, cv_flow)
        GC.gc()
        @test @allocated(PDMPSamplersRBridge._omrf_direct_prepare_path!(
            cv_envelope, allocation_state, cv_flow)) == 0
        cv_envelope.last_boundary_valid = false
        PDMPSamplersRBridge._omrf_direct_hcv_fill_cell_cached!(
            cv_envelope, 0.0, 0.1, 1, true)
        GC.gc()
        @test @allocated(PDMPSamplersRBridge._omrf_direct_hcv_fill_cell_cached!(
            cv_envelope, 0.1, 0.2, 2, true)) == 0
        cv_envelope.last_boundary_valid = false
        GC.gc()
        @test @allocated(PDMPSamplersRBridge._omrf_direct_hcv_fill_cell_cached!(
            cv_envelope, 0.3, 0.4, 3, true)) == 0

        # The existing marked construction is reused unchanged. Its complete
        # next-event path remains allocation-free after exact prewarming.
        cv_model_context, cv_model_names =
            PDMPSamplersRBridge._new_stan_subsampling_context(
                lib, full, prior, 3)
        cv_model, _, _ = PDMPSamplersRBridge._build_stan_subsampling_model(
            cv_model_context, cv_model_names,
            Dict{String,Any}("n_observations" => 9,
                "subsample_size" => 3, "residual_envelope" => cv_spec),
            anchor, nothing, falses(9); anchor_capacity=0)
        cv_event_rng = Xoshiro(0xd1ec7c)
        cv_event_flow = BouncyParticle(9, 0.0)
        cv_event_state, cv_event_model, cv_event_alg, cv_event_cache,
                cv_event_stats = PDMPSamplers.initialize_state(
            cv_event_rng, cv_event_flow, cv_model,
            GridThinningStrategy(N=2, N_min=2, t_max=0.2,
                lazy=false, bound_violation=:throw), 0.0,
            SkeletonPoint(copy(position), copy(cv_velocity)))
        cv_event_alg.schedule_frozen[] = true
        for _ in 1:4
            omrf_allocation_next_event!(cv_event_rng, cv_event_model,
                cv_event_flow, cv_event_alg, cv_event_state,
                cv_event_cache, cv_event_stats)
        end
        GC.gc()
        @test @allocated(omrf_allocation_next_event!(cv_event_rng,
            cv_event_model, cv_event_flow, cv_event_alg, cv_event_state,
            cv_event_cache, cv_event_stats)) == 0

        horizon_model_context, horizon_model_names =
            PDMPSamplersRBridge._new_stan_subsampling_context(
                lib, full, prior, 3)
        horizon_model, _, _, _, horizon_manager =
            PDMPSamplersRBridge._build_stan_subsampling_model(
                horizon_model_context, horizon_model_names,
                Dict{String,Any}("n_observations" => 9,
                    "subsample_size" => 3,
                    "residual_envelope" => horizon_spec),
                anchor, nothing, falses(9); anchor_capacity=1)
        refreshed_horizon_envelope =
            PDMPSamplersRBridge._prepare_omrf_anchor_envelope(
                horizon_manager, copy(anchor),
                horizon_model.grad.residual_oracle.residual_context)
        @test refreshed_horizon_envelope.event_horizon_scope
        horizon_event_rng = Xoshiro(0x48f120)
        horizon_event_state, horizon_event_model, horizon_event_alg,
                horizon_event_cache, horizon_event_stats =
            PDMPSamplers.initialize_state(
                horizon_event_rng, cv_event_flow, horizon_model,
                GridThinningStrategy(N=2, N_min=2, t_max=0.2,
                    lazy=false, bound_violation=:throw), 0.0,
                SkeletonPoint(copy(position), copy(cv_velocity)))
        horizon_event_alg.schedule_frozen[] = true
        for _ in 1:4
            omrf_allocation_next_event!(horizon_event_rng,
                horizon_event_model, cv_event_flow, horizon_event_alg,
                horizon_event_state, horizon_event_cache,
                horizon_event_stats)
        end
        GC.gc()
        @test @allocated(omrf_allocation_next_event!(horizon_event_rng,
            horizon_event_model, cv_event_flow, horizon_event_alg,
            horizon_event_state, horizon_event_cache,
            horizon_event_stats)) == 0
    end

    # Whole-person node-local HCV keeps all node residuals in one sampled
    # contribution.  The exact selected refinement must combine that residual
    # before applying the Zig-Zag rate rather than summing nodewise rates.
    whole_person_hcv_spec = copy(node_hcv_spec)
    whole_person_hcv_spec["factorization"] = "person"
    whole_person_hcv_spec["exact_selected_refinement"] = true
    whole_person_hcv_oracle =
        PDMPSamplersRBridge.OMRFAnalyticSubsamplingOracle(
            analytic_prior_ctx, whole_person_hcv_spec,
            analytic_prior_names, anchor)
    symbol_refinement_spec = Dict{Any,Any}(whole_person_hcv_spec)
    delete!(symbol_refinement_spec, "exact_selected_refinement")
    symbol_refinement_spec[:exact_selected_refinement] = true
    symbol_refinement_oracle =
        PDMPSamplersRBridge.OMRFAnalyticSubsamplingOracle(
            analytic_prior_ctx, symbol_refinement_spec,
            analytic_prior_names, anchor)
    @test symbol_refinement_oracle.exact_selected_refinement
    whole_person_state = PDMPState(0.0,
        SkeletonPoint(copy(position), collect(range(-0.8, 0.9; length=9))))
    whole_person_flow = ZigZag(9)
    PDMPSamplers.initialize_flow_state!(whole_person_state, whole_person_flow)
    whole_person_gradient = zeros(9)
    whole_person_hcv_oracle(whole_person_gradient, position)
    whole_person_D = PDMPSamplers.λ(
        whole_person_state, whole_person_gradient, whole_person_flow)
    whole_person_subset = [2]
    whole_person_M = whole_person_D + 1e6
    whole_person_bound = PDMPSamplers.subsampling_subset_bound(
        whole_person_hcv_oracle, whole_person_state, whole_person_gradient,
        whole_person_flow, whole_person_D, whole_person_M,
        whole_person_subset, 3.0)
    whole_person_cached = zeros(9)
    @test PDMPSamplers.subsampling_cached_residual!(
        whole_person_hcv_oracle, whole_person_cached, whole_person_state,
        whole_person_subset, anchor)
    @test !PDMPSamplers.subsampling_cached_residual!(
        whole_person_hcv_oracle, whole_person_cached, whole_person_state,
        whole_person_subset, anchor)
    whole_person_candidate_gradient = copy(whole_person_gradient)
    whole_person_actual = PDMPSamplers.subsampling_candidate_rate!(
        whole_person_hcv_oracle, whole_person_state,
        whole_person_candidate_gradient, whole_person_cached, 3.0,
        whole_person_flow, whole_person_D, whole_person_subset)
    @test whole_person_actual <= whole_person_bound + 2e-10
    @test @allocated(PDMPSamplers.subsampling_subset_bound(
        whole_person_hcv_oracle, whole_person_state, whole_person_gradient,
        whole_person_flow, whole_person_D, whole_person_M,
        whole_person_subset, 3.0)) == 0

    sign_spec = copy(envelope_spec)
    sign_spec["bound_type"] = "pattern_local_range"
    sign_spec["structural_groups"] = permutedims(reshape(collect(1:9), 3, 3))
    sign_spec["n_structural_components"] = 9
    sign_spec["component_nodes"] = repeat(collect(1:3); inner=3)
    sign_spec["component_covariates"] = repeat(
        [0 1 2; 2 0 1; 1 2 0], 3, 1)
    sign_spec["zigzag_sign_aware"] = true
    sign_oracle = PDMPSamplersRBridge.OMRFAnalyticSubsamplingOracle(
        ctx, sign_spec, names, anchor)
    sign_velocity = collect(range(-0.75, 0.85; length=9))
    sign_state = PDMPState(0.0, SkeletonPoint(copy(position), copy(sign_velocity)))
    sign_flow = ZigZag(9)
    PDMPSamplers.initialize_flow_state!(sign_state, sign_flow)
    sign_gradient = zeros(9)
    sign_residual = zeros(9)
    sign_subset = [1]
    sign_oracle(sign_gradient, position)
    sign_oracle(sign_residual, position, sign_subset, anchor)
    sign_D = PDMPSamplers.λ(sign_state, sign_gradient, sign_flow)
    sign_M = 1e12
    sign_bound = PDMPSamplers.subsampling_residual_subset_bound(
        sign_oracle, sign_state, sign_gradient, sign_flow,
        sign_D, sign_M, sign_subset, 3.0)
    @test PDMPSamplers.λ(sign_state,
        sign_gradient .+ 3.0 .* sign_residual, sign_flow) <= sign_bound + 1e-10
    @test @allocated(PDMPSamplers.subsampling_residual_subset_bound(
        sign_oracle, sign_state, sign_gradient, sign_flow,
        sign_D, sign_M, sign_subset, 3.0)) == 0

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

    factor_spec = copy(analytic_prior_spec)
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
    factor_model, _, _, _, factor_manager =
        PDMPSamplersRBridge._build_stan_subsampling_model(
            factor_ctx, factor_names, factor_subsampling, anchor, nothing,
            falses(9); anchor_capacity=1)
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

    @testset "analytic person-node allocation contract" begin
        allocation_flow = BouncyParticle(9, 0.0)
        allocation_state = PDMPState(0.0,
            SkeletonPoint(copy(position), fill(1.0, 9)))
        allocation_cv = factor_model.grad
        allocation_gradient = zeros(9)
        allocation_residual = zeros(9)
        allocation_subset = allocation_cv.subset
        copyto!(allocation_subset, [1, 5, 9])
        PDMPSamplers.prepare_residual_sampling!(allocation_cv.envelope,
            allocation_state, allocation_flow, 0.0)
        allocation_B = PDMPSamplers.screening_residual_bound(
            allocation_cv.envelope, allocation_state, allocation_flow, 0.0)
        PDMPSamplers.deterministic_gradient!(allocation_gradient,
            allocation_cv, allocation_state.ξ.x)
        allocation_D = PDMPSamplers.λ(
            allocation_state, allocation_gradient, allocation_flow)
        allocation_scale = PDMPSamplers._subsampling_scale(allocation_cv)
        allocation_M = PDMPSamplers.draw_subset!(Xoshiro(0xa110c),
            allocation_cv, allocation_D, allocation_B)
        PDMPSamplers.subsampling_residual_subset_bound(
            allocation_cv.residual_oracle, allocation_state,
            allocation_gradient, allocation_flow, allocation_D,
            allocation_M, allocation_subset, allocation_scale)
        PDMPSamplers.subsampling_subset_bound(
            allocation_cv.residual_oracle, allocation_state,
            allocation_gradient, allocation_flow, allocation_D,
            allocation_M, allocation_subset, allocation_scale)
        allocation_bound_fixture = (; oracle=allocation_cv.residual_oracle,
            state=allocation_state, gradient=allocation_gradient,
            flow=allocation_flow, D=allocation_D, M=allocation_M,
            subset=allocation_subset, scale=allocation_scale)
        omrf_allocation_residual_bound!(allocation_bound_fixture)
        omrf_allocation_subset_bound!(allocation_bound_fixture)
        allocation_cv.residual_oracle(allocation_residual,
            allocation_state.ξ.x, allocation_subset, allocation_cv.anchor)

        @test @allocated(PDMPSamplers.screening_residual_bound(
            allocation_cv.envelope, allocation_state,
            allocation_flow, 0.0)) == 0
        allocation_rng = Xoshiro(0xa110d)
        PDMPSamplers.draw_subset!(allocation_rng, allocation_cv,
            allocation_D, allocation_B)
        @test @allocated(PDMPSamplers.draw_subset!(allocation_rng,
            allocation_cv, allocation_D, allocation_B)) == 0
        @test @allocated(omrf_allocation_residual_bound!(
            allocation_bound_fixture)) == 0
        @test @allocated(omrf_allocation_subset_bound!(
            allocation_bound_fixture)) == 0
        @test @allocated(allocation_cv.residual_oracle(allocation_residual,
            allocation_state.ξ.x, allocation_subset,
            allocation_cv.anchor)) == 0
        @test @allocated(PDMPSamplers.deterministic_gradient!(
            allocation_gradient, allocation_cv,
            allocation_state.ξ.x)) == 0

        factor_manager.main_refresh_distance = Inf
        no_refresh_fixture = (; manager=factor_manager, cv=allocation_cv,
            position=copy(factor_manager.oracle.anchor))
        omrf_allocation_anchor_select!(no_refresh_fixture)
        @test @allocated(omrf_allocation_anchor_select!(
            no_refresh_fixture)) == 0

        # The first capacity-one replacement exposes an inactive state for
        # recycling. Prewarm that recycled refresh, then measure an equivalent
        # later refresh with all anchor/context/envelope buffers already owned.
        factor_manager.main_refresh_distance = 0.0
        first_refresh_fixture = (; manager=factor_manager, cv=allocation_cv,
            position=copy(factor_manager.oracle.anchor) .+ 0.02)
        omrf_allocation_anchor_select!(first_refresh_fixture)
        warm_refresh_fixture = (; manager=factor_manager, cv=allocation_cv,
            position=copy(factor_manager.oracle.anchor) .+ 0.02)
        omrf_allocation_anchor_select!(warm_refresh_fixture)
        measured_refresh_fixture = (; manager=factor_manager,
            cv=allocation_cv,
            position=copy(factor_manager.oracle.anchor) .+ 0.02)
        # The prepared buffers are fully recycled. The remaining 112 bytes are
        # the two abstract owner-reference replacements in
        # refresh_anchor_owned!; keep the byte count explicit rather than
        # hiding it behind a broad allowance.
        @test @allocated(omrf_allocation_anchor_select!(
            measured_refresh_fixture)) == 112

        event_rng = Xoshiro(0xa110e)
        event_state, event_model, event_alg, event_cache, event_stats =
            PDMPSamplers.initialize_state(event_rng, allocation_flow,
                factor_model,
                GridThinningStrategy(N=2, N_min=2, t_max=0.25,
                    lazy=false, bound_violation=:throw), 0.0,
                SkeletonPoint(copy(position), fill(1.0, 9)))
        event_alg.schedule_frozen[] = true
        for _ in 1:3
            omrf_allocation_next_event!(event_rng, event_model,
                allocation_flow, event_alg, event_state, event_cache,
                event_stats)
        end
        GC.gc()
        @test @allocated(omrf_allocation_next_event!(event_rng, event_model,
            allocation_flow, event_alg, event_state, event_cache,
            event_stats)) == 0
    end

    exact_factor_spec = copy(factor_spec)
    exact_factor_spec["exact_selected_refinement"] = true
    exact_factor_ctx, exact_factor_names =
        PDMPSamplersRBridge._new_stan_subsampling_context(
            lib, full, prior, 3)
    exact_factor_model, _, _ =
        PDMPSamplersRBridge._build_stan_subsampling_model(
            exact_factor_ctx, exact_factor_names,
            Dict{String,Any}(
                "n_observations" => 9,
                "subsample_size" => 3,
                "residual_envelope" => exact_factor_spec),
            anchor, nothing, falses(9))
    exact_oracle = exact_factor_model.grad.residual_oracle
    @test exact_oracle.exact_selected_refinement
    exact_state = PDMPState(0.0,
        SkeletonPoint(copy(position), fill(1.0, 9)))
    exact_gradient = zeros(9)
    compute_gradient!(exact_factor_model.grad, position, exact_gradient)
    exact_flow = ZigZag(9)
    exact_D = PDMPSamplers.λ(exact_state, exact_gradient, exact_flow)
    exact_M = exact_D + 1e6
    exact_subset = [1, 5, 9]
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
    direct_residual = zeros(9)
    exact_oracle(direct_residual, position, exact_subset, anchor)
    @test cached_residual ≈ direct_residual atol=2e-14 rtol=2e-14
    @test !exact_oracle.cached_residual_valid

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

    # Exact-variance HCV banks are learned automatically from the current
    # chain's bounded warmup grid and final adapted velocity metric. This is a
    # warmup finalizer, not a retained or evaluation-trajectory intervention.
    variance_spec = copy(hcv_spec)
    variance_spec["hcv_damping"] = Inf
    variance_spec["hcv_bank_learning"] = "exact_variance"
    variance_subsampling = Dict{String,Any}(
        "n_observations" => 3,
        "subsample_size" => 1,
        "residual_envelope" => variance_spec)
    variance_ctx, variance_names =
        PDMPSamplersRBridge._new_stan_subsampling_context(
            lib, full, prior, 1)
    variance_model, _, _, _, variance_manager =
        PDMPSamplersRBridge._build_stan_subsampling_model(
            variance_ctx, variance_names, variance_subsampling, anchor,
            nothing, falses(9); anchor_capacity=2)
    @test variance_manager.warmup_bank_learning_policy === :exact_variance
    variance_trace = PDMPTrace([
        PDMPEvent(0.0, anchor, fill(0.2, 9)),
        PDMPEvent(1.0, position .- 0.15, fill(0.2, 9)),
        PDMPEvent(2.0, position .+ 0.20, fill(-0.1, 9)),
        PDMPEvent(3.0, position .+ collect(range(-0.25, 0.25; length=9)),
            fill(0.1, 9)),
    ], BouncyParticle(9, 0.0))
    variance_flow = PreconditionedDynamics(
        DiagonalPreconditioner(collect(range(0.7, 1.3; length=9))),
        BouncyParticle(9, 0.0))
    variance_anchors, variance_scales =
        PDMPSamplersRBridge._omrf_exact_variance_anchors!(
            variance_manager, variance_trace, position, 2, variance_flow)
    @test size(variance_anchors) == (9, 2)
    @test variance_anchors[:, 1] == position
    @test all(isfinite, variance_scales)
    @test all(>(0), variance_scales)
    @test variance_manager.warmup_bank_candidate_count == 5
    @test variance_manager.warmup_bank_training_state_count == 5
    @test isfinite(variance_manager.warmup_bank_objective_single)
    @test isfinite(variance_manager.warmup_bank_objective_selected)
    @test variance_manager.warmup_bank_objective_selected <=
        variance_manager.warmup_bank_objective_single
    @test variance_manager.warmup_bank_objective_selected <
        variance_manager.warmup_bank_objective_single
    @test variance_manager.warmup_bank_learning_seconds > 0

    # The selector-coupled combinatorial objective evaluates the source that
    # the retained geometric selector actually installs. It is intentionally
    # not the oracle minimum of each column of `costs`.
    coupled_candidates = reshape([0.0, 1.0, 2.0, 3.0], 1, :)
    coupled_costs = [1.0 9.0 9.0 9.0;
                     9.0 1.0 8.0 8.0;
                     9.0 8.0 1.0 8.0;
                     9.0 8.0 8.0 1.0]
    coupled_indices, coupled_objective =
        PDMPSamplersRBridge._omrf_selector_coupled_indices(
            coupled_costs, coupled_candidates, [1.0], 2)
    brute = [(indices=[1, candidate], objective=
        PDMPSamplersRBridge._omrf_selector_assigned_objective(
            coupled_costs, coupled_candidates, [1.0], [1, candidate]))
        for candidate in 2:4]
    brute_best = brute[argmin(getfield.(brute, :objective))]
    @test coupled_indices == brute_best.indices
    @test coupled_objective == brute_best.objective
    explicit_work = zeros(9)
    explicit_mean = zeros(9)
    explicit_score =
        PDMPSamplersRBridge._omrf_hcv_expected_directional_variance!(
            variance_manager.recycled_state.residual_context,
            variance_manager.recycled_state.anchor, position,
            variance_flow.metric.scale .^ 2, explicit_work, explicit_mean,
            [1])
    @test isfinite(explicit_score)
    @test explicit_score >= 0
    explicit_context = variance_manager.recycled_state.residual_context
    @. explicit_context.displacement = position -
        variance_manager.recycled_state.anchor
    explicit_residuals = zeros(9, 3)
    explicit_one = zeros(9)
    for person in 1:3
        PDMPSamplersRBridge.omrf_residual!(
            explicit_context, explicit_one, position, [person])
        explicit_residuals[:, person] .= explicit_one
    end
    explicit_center = vec(mean(explicit_residuals; dims=2))
    manual_score = sum(variance_flow.metric.scale[coordinate]^2 *
        mean((explicit_residuals[coordinate, :] .-
            explicit_center[coordinate]).^2) for coordinate in 1:9)
    @test explicit_score ≈ manual_score atol=2e-14 rtol=2e-14
    variance_manager.selection_policy = :position
    variance_state = PDMPState(3.0,
        SkeletonPoint(copy(position), fill(0.1, 9)))
    variance_trace_manager = PDMPSamplers.TraceManager(
        variance_trace, variance_trace, 3.0)
    @test PDMPSamplersRBridge._finish_omrf_warmup!(
        variance_manager, variance_model.grad, position,
        variance_trace_manager; state=variance_state, flow=variance_flow)
    @test length(variance_manager.entries) == 2
    @test variance_manager.warmup_finish_preparations ==
        variance_manager.preparations
    @test variance_manager.warmup_finish_native_full_data_scans ==
        variance_ctx.counts.native_full_data_scans
    @test variance_model.grad.anchor ==
        variance_manager.entries[variance_manager.active_idx].state.anchor

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
    @test PDMPSamplersRBridge._select_omrf_anchor!(
        manager, managed_model.grad, refresh_position; phase=:main)
    @test manager.preparations == preparations_before_main + 1
    @test manager.main_refreshes == 1
    @test manager.main_activations == 1
    @test manager.oracle.anchor == refresh_position
    @test managed_model.grad.anchor == refresh_position
    @test !PDMPSamplersRBridge._select_omrf_anchor!(
        manager, managed_model.grad, refresh_position; phase=:main)
    PDMPSamplersRBridge._finish_omrf_warmup!(
        manager, managed_model.grad, refresh_position)
    @test manager.warmup_finish_preparations == manager.preparations
    post_warmup_position = refresh_position .+ 0.2
    @test PDMPSamplersRBridge._select_omrf_anchor!(
        manager, managed_model.grad, post_warmup_position; phase=:main)
    @test manager.preparations - manager.warmup_finish_preparations == 1
    @test bank_ctx.counts.native_full_data_scans -
        manager.warmup_finish_native_full_data_scans > 0

    representative_trace = PDMPTrace([
        PDMPEvent(0.0, zeros(2), ones(2)),
        PDMPEvent(1.0, [1.0, 0.0], ones(2)),
        PDMPEvent(2.0, [2.0, 2.0], ones(2)),
    ], BouncyParticle(2, 0.0))
    representatives, representative_scales =
        PDMPSamplersRBridge._omrf_representative_anchors(
            representative_trace, [3.0, 3.0], 3)
    @test size(representatives) == (2, 3)
    @test representatives[:, 1] == [3.0, 3.0]
    @test all(isfinite, representative_scales)
    @test all(>(0), representative_scales)
    @test length(unique(eachcol(representatives))) == 3
    tail_candidates = PDMPSamplersRBridge._omrf_warmup_candidates(
        representative_trace, [3.0, 3.0], 3; retained_fraction=0.5)
    @test all(tail_candidates[1, 2:end] .>= 1.0)

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

    # Active-set-aware pattern compression must rebuild from the current
    # sticky mask, update both incident nodes of a frozen interaction, and
    # preserve the pointwise and cell certificates.  This fixture has four
    # full patterns per node; freezing edge (1, 2) merges the first two node
    # partitions down to two patterns each while leaving node three unchanged.
    active_X = [0 0 1; 0 1 1; 1 0 2; 1 1 2]
    active_spec = copy(pattern_spec)
    active_spec["X"] = active_X
    active_spec["node_weights"] = zeros(3, 4)
    active_spec["structural_groups"] = [1 2 3 4; 5 6 7 8; 9 10 11 12]
    active_spec["n_structural_components"] = 12
    active_spec["component_nodes"] = repeat(collect(1:3); inner=4)
    active_spec["component_covariates"] = repeat(active_X, 3, 1)
    active_spec["active_set_compression"] = true
    active_anchor = zeros(9)
    active_oracle = PDMPSamplersRBridge.OMRFAnalyticSubsamplingOracle(
        model_ctx, active_spec, model_names, active_anchor)
    active_envelope = PDMPSamplersRBridge._build_omrf_residual_envelope(
        active_spec, active_anchor, active_oracle.residual_context)
    reference_spec = copy(active_spec)
    reference_spec["active_set_compression"] = false
    reference_oracle = PDMPSamplersRBridge.OMRFAnalyticSubsamplingOracle(
        model_ctx, reference_spec, model_names, active_anchor)
    reference_envelope = PDMPSamplersRBridge._build_omrf_residual_envelope(
        reference_spec, active_anchor, reference_oracle.residual_context)
    active_free = trues(9)
    edge_12 = active_oracle.residual_context.interaction_indices[1]
    active_position = collect(range(-0.18, 0.23; length=9))
    active_position[edge_12] = 0.0
    active_state = StickyPDMPState(0.0,
        SkeletonPoint(copy(active_position), ones(9)), copy(active_free))
    reference_state = copy(active_state)
    active_flow = ZigZag(9)
    PDMPSamplers.initialize_flow_state!(active_state, active_flow)
    PDMPSamplers.initialize_flow_state!(reference_state, active_flow)
    active_point = similar(active_envelope.scales)
    reference_point = similar(reference_envelope.scales)
    PDMPSamplers.component_scales!(active_point, active_envelope, active_state,
        active_flow, 0.0)
    PDMPSamplers.component_scales!(reference_point, reference_envelope, reference_state,
        active_flow, 0.0)
    active_scales = active_envelope.component_scales!
    @test active_scales.active_component_count == 10
    @test active_scales.context.active_rebuilds == 1
    # The active-set rebuild is immutable after the first call; Julia 1.12
    # retains a small bookkeeping allocation in this diagnostic path.
    @test @allocated(PDMPSamplers.component_scales!(active_point, active_envelope,
        active_state, active_flow, 0.0)) <= 96
    for observation in 1:4
        @test PDMPSamplers.observation_residual_bound(active_envelope, observation) ≈ PDMPSamplers.observation_residual_bound(reference_envelope, observation) atol=2e-12 rtol=2e-12
    end

    copyto!(active_free, trues(9))
    active_free[edge_12] = false
    copyto!(active_state.free, active_free)
    copyto!(reference_state.free, active_free)
    PDMPSamplers.set_active_set!(active_oracle, active_free)
    @test @allocated(PDMPSamplers.component_scales!(active_point,
        active_envelope, active_state, active_flow, 0.0)) == 0
    PDMPSamplers.component_scales!(reference_point, reference_envelope, reference_state,
        active_flow, 0.0)
    @test active_scales.active_component_count == 8
    @test active_scales.node_component_ends .-
        active_scales.node_component_starts .+ 1 == [2, 2, 4]
    @test active_scales.context.active_generation == 1
    @test @allocated(PDMPSamplers.component_scales!(active_point, active_envelope,
        active_state, active_flow, 0.0)) == 0
    for observation in 1:4
        @test PDMPSamplers.observation_residual_bound(active_envelope, observation) ≈ PDMPSamplers.observation_residual_bound(reference_envelope, observation) atol=2e-12 rtol=2e-12
    end
    active_cell = similar(active_envelope.cell_scales)
    reference_cell = similar(reference_envelope.cell_scales)
    PDMPSamplers.component_cell_scales!(active_cell, active_envelope, active_state,
        active_flow, 0.05, 0.35)
    PDMPSamplers.component_cell_scales!(reference_cell, reference_envelope, reference_state,
        active_flow, 0.05, 0.35)
    for node in 1:3, observation in 1:4
        @test active_cell[active_envelope.groups[node, observation]] <= reference_cell[reference_envelope.groups[node, observation]] + 2e-12
    end

    # Restoring the edge to the active set must restore the full partition,
    # even when its current numerical coordinate is exactly zero.
    active_free[edge_12] = true
    copyto!(active_state.free, active_free)
    PDMPSamplers.set_active_set!(active_oracle, active_free)
    PDMPSamplers.component_scales!(active_point, active_envelope, active_state,
        active_flow, 0.0)
    @test active_scales.active_component_count == 10
    @test active_scales.context.active_generation == 2
    @test active_scales.context.active_rebuilds == 3
    @test active_scales.context.active_group_count_min == 8
    @test active_scales.context.active_group_count_max == 10

    pattern_hcv_spec = copy(pattern_spec)
    pattern_hcv_spec["bound_type"] = "pattern_local_hcv_range"
    pattern_hcv_spec["use_hcv"] = true
    pattern_hcv_spec["hcv_after_warmup"] = false
    pattern_hcv_spec["hcv_damping"] = Inf
    pattern_hcv_oracle = PDMPSamplersRBridge.OMRFAnalyticSubsamplingOracle(
        model_ctx, pattern_hcv_spec, analytic_prior_names, copy(anchor))
    pattern_hcv_envelope = PDMPSamplersRBridge._build_stan_residual_envelope(
        pattern_hcv_spec, pattern_hcv_oracle.anchor, pattern_hcv_oracle)
    @test pattern_hcv_envelope isa PDMPSamplers.GroupedResidualEnvelope
    @test pattern_hcv_envelope.component_scales! isa
        PDMPSamplersRBridge.OMRFPatternLocalHCVScales
    factor_pattern_hcv_spec = copy(pattern_hcv_spec)
    factor_pattern_hcv_spec["factorization"] = "person_node"
    factor_pattern_hcv_spec["structural_groups"] =
        reshape(collect(1:9), 1, 9)
    factor_pattern_hcv_oracle =
        PDMPSamplersRBridge.OMRFAnalyticSubsamplingOracle(
            model_ctx, factor_pattern_hcv_spec, analytic_prior_names,
            copy(anchor))
    factor_pattern_hcv_context =
        factor_pattern_hcv_oracle.residual_context
    factor_pattern_hcv_envelope =
        PDMPSamplersRBridge._build_omrf_residual_envelope(
            factor_pattern_hcv_spec, anchor,
            factor_pattern_hcv_context)
    @test factor_pattern_hcv_envelope isa
        PDMPSamplers.GroupedResidualEnvelope
    @test PDMPSamplers.n_observations(factor_pattern_hcv_envelope) == 9
    @test factor_pattern_hcv_envelope.component_scales! isa
        PDMPSamplersRBridge.OMRFBatchedFactorHCVScales
    grouped_hcv_spec = copy(pattern_hcv_spec)
    grouped_hcv_spec["bound_type"] = "grouped_person_hcv_box"
    grouped_hcv_spec["box_membership"] = [1, 1, 2]
    grouped_hcv_spec["box_minima"] = [0 0 1; 1 2 0]
    grouped_hcv_spec["box_maxima"] = [2 1 2; 1 2 0]
    grouped_hcv_oracle =
        PDMPSamplersRBridge.OMRFAnalyticSubsamplingOracle(
            model_ctx, grouped_hcv_spec, analytic_prior_names,
            copy(anchor))
    grouped_hcv_envelope =
        PDMPSamplersRBridge._build_omrf_residual_envelope(
            grouped_hcv_spec, grouped_hcv_oracle.anchor,
            grouped_hcv_oracle.residual_context)
    @test grouped_hcv_envelope isa PDMPSamplers.GroupedResidualEnvelope
    @test PDMPSamplers.n_observations(grouped_hcv_envelope) == 3
    @test grouped_hcv_envelope.component_scales! isa
        PDMPSamplersRBridge.OMRFGroupedPersonHCVBoxScales
    @test grouped_hcv_envelope.groups == reshape([1, 1, 2], 1, :)
    node_grouped_hcv_spec = copy(pattern_hcv_spec)
    node_grouped_hcv_spec["bound_type"] =
        "node_grouped_person_hcv_box"
    node_grouped_hcv_spec["node_box_membership"] = [
        1 1 2;
        3 4 3;
        5 6 5
    ]
    node_grouped_hcv_spec["node_box_nodes"] = [1, 1, 2, 2, 3, 3]
    node_grouped_hcv_spec["node_box_minima"] = [
        0 0 1; 0 2 0; 0 0 0; 2 0 1; 0 1 0; 2 0 0
    ]
    node_grouped_hcv_spec["node_box_maxima"] = [
        0 1 2; 0 2 0; 1 0 2; 2 0 1; 1 2 0; 2 0 0
    ]
    node_grouped_hcv_oracle =
        PDMPSamplersRBridge.OMRFAnalyticSubsamplingOracle(
            model_ctx, node_grouped_hcv_spec, analytic_prior_names,
            copy(anchor))
    node_grouped_hcv_envelope =
        PDMPSamplersRBridge._build_omrf_residual_envelope(
            node_grouped_hcv_spec, node_grouped_hcv_oracle.anchor,
            node_grouped_hcv_oracle.residual_context)
    @test node_grouped_hcv_envelope isa PDMPSamplers.GroupedResidualEnvelope
    @test node_grouped_hcv_envelope.groups ==
        node_grouped_hcv_spec["node_box_membership"]
    @test node_grouped_hcv_envelope.component_scales! isa
        PDMPSamplersRBridge.OMRFNodeGroupedPersonHCVBoxScales
    cell_cached_spec = copy(node_grouped_hcv_spec)
    cell_cached_spec["bound_type"] =
        "cell_cached_node_grouped_person_hcv_box"
    cell_cached_oracle =
        PDMPSamplersRBridge.OMRFAnalyticSubsamplingOracle(
            model_ctx, cell_cached_spec, analytic_prior_names,
            copy(anchor))
    cell_cached_envelope =
        PDMPSamplersRBridge._build_omrf_residual_envelope(
            cell_cached_spec, cell_cached_oracle.anchor,
            cell_cached_oracle.residual_context)
    @test cell_cached_envelope.component_scales! isa
        PDMPSamplersRBridge.OMRFCellCachedNodeGroupedPersonHCVBoxScales
    lazy_node_grouped_hcv_spec = copy(node_grouped_hcv_spec)
    lazy_node_grouped_hcv_spec["bound_type"] =
        "lazy_node_grouped_person_hcv_box"
    lazy_node_grouped_hcv_spec["node_box_membership"] = [
        2 2 3;
        5 6 5;
        8 9 8
    ]
    lazy_node_grouped_hcv_spec["node_box_nodes"] = [
        1, 1, 1, 2, 2, 2, 3, 3, 3]
    lazy_node_grouped_hcv_spec["node_box_minima"] = [
        0 0 1; 0 0 1; 0 0 2;
        0 0 0; 0 0 0; 2 0 0;
        0 0 0; 0 0 0; 1 0 0
    ]
    lazy_node_grouped_hcv_spec["node_box_maxima"] = [
        0 2 2; 0 1 2; 0 2 2;
        2 0 2; 1 0 2; 2 0 2;
        2 2 0; 2 1 0; 1 2 0
    ]
    lazy_node_grouped_hcv_spec["node_box_parents"] = [
        0, 1, 1, 0, 4, 4, 0, 7, 7]
    lazy_node_grouped_hcv_spec["node_box_roots"] = [1, 4, 7]
    lazy_node_grouped_hcv_spec["node_box_leaves"] = [2, 3, 5, 6, 8, 9]
    lazy_node_grouped_hcv_oracle =
        PDMPSamplersRBridge.OMRFAnalyticSubsamplingOracle(
            model_ctx, lazy_node_grouped_hcv_spec,
            analytic_prior_names, copy(anchor))
    lazy_node_grouped_hcv_envelope =
        PDMPSamplersRBridge._build_omrf_residual_envelope(
            lazy_node_grouped_hcv_spec,
            lazy_node_grouped_hcv_oracle.anchor,
            lazy_node_grouped_hcv_oracle.residual_context)
    @test lazy_node_grouped_hcv_envelope.component_scales! isa
        PDMPSamplersRBridge.OMRFLazyNodeGroupedPersonHCVBoxScales
    @test PDMPSamplers.n_observations(lazy_node_grouped_hcv_envelope) == 3
    @test length(lazy_node_grouped_hcv_envelope.scales) == 3
    mixed_hcv_candidates = hcat(anchor, anchor .+ 0.17)
    mixed_hcv_entries = PDMPSamplersRBridge.OMRFAnchorEntry[]
    for candidate in eachcol(mixed_hcv_candidates)
        mixed_context = PDMPSamplersRBridge._omrf_residual_context(
            factor_pattern_hcv_spec, analytic_prior_names, candidate,
            model_ctx.counts)
        mixed_state = PDMPSamplersRBridge.OMRFPreparedAnchorState(
            collect(candidate), zeros(9), zeros(9), mixed_context,
            factor_pattern_hcv_envelope)
        push!(mixed_hcv_entries,
            PDMPSamplersRBridge.OMRFAnchorEntry(mixed_state, 0))
    end
    mixed_hcv_bank =
        PDMPSamplersRBridge._prepare_omrf_node_local_hcv_bank(
            mixed_hcv_entries, model_ctx.counts)
    PDMPSamplersRBridge._install_omrf_node_local_sources!(
        mixed_hcv_bank, [1, 2, 1])
    factor_pattern_hcv_oracle.node_local_hcv = mixed_hcv_bank
    PDMPSamplersRBridge._install_omrf_batched_node_sources!(
        factor_pattern_hcv_envelope.component_scales!, mixed_hcv_bank)
    grouped_hcv_oracle.node_local_hcv = mixed_hcv_bank
    grouped_hcv_manager =
        PDMPSamplersRBridge._new_omrf_anchor_manager(
            grouped_hcv_oracle, grouped_hcv_spec,
            analytic_prior_names, grouped_hcv_envelope, 2)
    grouped_hcv_manager.oracle.node_local_hcv = mixed_hcv_bank
    grouped_hcv_cv = SubsampledControlVariate(
        (out, x) -> fill!(out, 0.0), grouped_hcv_oracle,
        grouped_hcv_envelope, anchor, 1)
    grouped_hcv_state = PDMPState(0.0,
        SkeletonPoint(copy(anchor), fill(0.25, length(anchor))))
    grouped_hcv_flow = AdaptiveBoomerang(length(anchor); λref=0.2)
    PDMPSamplers.initialize_flow_state!(grouped_hcv_state, grouped_hcv_flow)
    grouped_hcv_state.ξ.x .+= range(0.05, 0.25; length=length(anchor))
    warmup_pattern_trace = PDMPTrace([
        PDMPEvent(0.0, copy(anchor), fill(0.10, length(anchor))),
        PDMPEvent(1.0, copy(grouped_hcv_state.ξ.x),
            fill(0.25, length(anchor))),
    ], grouped_hcv_flow)
    visits_before = model_ctx.counts.event_bound_observation_visits
    PDMPSamplersRBridge._select_omrf_warmup_pattern_sources!(
        grouped_hcv_manager, grouped_hcv_cv, warmup_pattern_trace,
        grouped_hcv_state, grouped_hcv_flow)
    @test all(1 .<= mixed_hcv_bank.sources .<= 2)
    @test grouped_hcv_manager.node_source_checks == 1
    @test model_ctx.counts.event_bound_observation_visits - visits_before ==
        32 * 3 * 3 * 2
    @test grouped_hcv_envelope.component_scales!.node_sources ==
        mixed_hcv_bank.sources
    mapped_checks = grouped_hcv_manager.node_source_checks
    mapped_visits = model_ctx.counts.event_bound_observation_visits
    PDMPSamplersRBridge._select_omrf_warmup_mapped_sources!(
        grouped_hcv_manager, grouped_hcv_cv, grouped_hcv_state)
    @test grouped_hcv_manager.node_source_checks == mapped_checks + 1
    @test model_ctx.counts.event_bound_observation_visits == mapped_visits
    pattern_checks = grouped_hcv_manager.node_source_checks
    PDMPSamplersRBridge._select_omrf_pattern_count_sources!(
        grouped_hcv_manager, grouped_hcv_cv,
        grouped_hcv_state, grouped_hcv_flow)
    @test grouped_hcv_manager.node_source_checks == pattern_checks + 1
    @test model_ctx.counts.event_bound_observation_visits == mapped_visits
    @test all(sum(mixed_hcv_bank.node_pattern_counts[node]) == 3
        for node in 1:3)
    cell_cached_scales = cell_cached_envelope.component_scales!
    cell_roof = similar(cell_cached_envelope.scales)
    point_roof = similar(cell_cached_envelope.scales)
    PDMPSamplers.component_cell_scales!(cell_roof, cell_cached_envelope,
        grouped_hcv_state, grouped_hcv_flow, 0.0, 0.2)
    point_calls = cell_cached_oracle.ctx.counts.envelope_point_calls
    PDMPSamplers.component_scales!(point_roof, cell_cached_envelope,
        grouped_hcv_state, grouped_hcv_flow, 0.1)
    @test point_roof == cell_roof
    @test cell_cached_oracle.ctx.counts.envelope_point_calls == point_calls
    @test cell_cached_scales.n_cached_cells == 1
    @test PDMPSamplers.screening_residual_bound(
        cell_cached_envelope, grouped_hcv_state, grouped_hcv_flow, 0.1) ≈
        dot(cell_roof, cell_cached_envelope.totals)
    cell_cached_cv = SubsampledControlVariate(
        (out, x) -> fill!(out, 0.0), cell_cached_oracle,
        cell_cached_envelope, anchor, 1)
    cell_B = PDMPSamplers.screening_residual_bound(
        cell_cached_envelope, grouped_hcv_state, grouped_hcv_flow, 0.1)
    cell_M = PDMPSamplers.draw_subset!(
        Xoshiro(20260814), cell_cached_cv, 0.0, cell_B)
    @test cell_cached_scales.selected_person_point_bound <= cell_M + 1e-12
    for person in 1:3, node in 1:3
        component = cell_cached_envelope.groups[node, person]
        exact_person_bound =
            PDMPSamplersRBridge._omrf_node_grouped_person_point_bound(
                cell_cached_scales.base, person, node)
        box_bound = PDMPSamplersRBridge._omrf_node_grouped_point_bound(
            cell_cached_scales.base, component, node)
        @test exact_person_bound <= box_bound + 1e-12
    end
    @test PDMPSamplers.subsampling_residual_subset_bound(
        cell_cached_oracle, grouped_hcv_state, grouped_hcv_flow, 0.0,
        cell_M, cell_cached_cv.subset, 3.0) ≈
        cell_cached_scales.selected_person_point_bound
    lazy_cv = SubsampledControlVariate(
        (out, x) -> fill!(out, 0.0), lazy_node_grouped_hcv_oracle,
        lazy_node_grouped_hcv_envelope, anchor, 1)
    lazy_B = PDMPSamplers.screening_residual_bound(
        lazy_cv.envelope, grouped_hcv_state, grouped_hcv_flow, 0.0)
    @test lazy_B >= 0
    lazy_scales = lazy_cv.envelope.component_scales!
    @test lazy_scales.root_totals == lazy_cv.envelope.totals
    leaf_person_bounds = zeros(3)
    direct_leaf_person_bounds = zeros(3)
    for person in 1:3, node in 1:3
        component = lazy_node_grouped_hcv_spec[
            "node_box_membership"][node, person]
        leaf_person_bounds[person] +=
            PDMPSamplersRBridge._omrf_lazy_component_point_bound(
                lazy_scales, component)
        direct_leaf_person_bounds[person] +=
            PDMPSamplersRBridge._omrf_node_grouped_point_bound(
                lazy_scales.base, component, node)
    end
    @test leaf_person_bounds == direct_leaf_person_bounds
    leaf_total = sum(leaf_person_bounds)
    @test lazy_B + 1e-12 >= leaf_total
    lazy_rng = Xoshiro(20260813)
    lazy_draws = 20_000
    lazy_successes = 0
    lazy_person_counts = zeros(Int, 3)
    for _ in 1:lazy_draws
        bound = PDMPSamplers.draw_subset!(lazy_rng, lazy_cv, 0.0, lazy_B)
        if bound > 0
            lazy_successes += 1
            lazy_person_counts[lazy_cv.subset[1]] += 1
            @test bound + 1e-12 >= 3leaf_person_bounds[lazy_cv.subset[1]]
        end
    end
    if lazy_B > 0
        @test lazy_successes / lazy_draws ≈ leaf_total / lazy_B atol=0.02
        @test lazy_person_counts ./ lazy_draws ≈
            leaf_person_bounds ./ lazy_B atol=0.02
    end
    @test lazy_node_grouped_hcv_oracle.ctx.counts.lazy_branch_evaluations > 0
    @test lazy_node_grouped_hcv_oracle.ctx.counts.lazy_leaf_evaluations > 0
    proxy_rows_before = copy(mixed_hcv_bank.contexts[1].X)
    selections_before_direct = grouped_hcv_manager.main_selections
    PDMPSamplersRBridge._select_omrf_node_local_sources!(
        grouped_hcv_manager, grouped_hcv_cv,
        grouped_hcv_state, grouped_hcv_flow)
    @test grouped_hcv_manager.main_selections == selections_before_direct + 1
    @test mixed_hcv_bank.contexts[1].X == proxy_rows_before
    @test grouped_hcv_envelope.component_scales!.node_anchors ===
        mixed_hcv_bank.candidates
    grouped_hcv_manager.node_source_check_distance = 1.0
    copyto!(grouped_hcv_manager.node_source_last_position,
        grouped_hcv_state.ξ.x)
    checks_before = grouped_hcv_manager.node_source_checks
    skips_before = grouped_hcv_manager.node_source_skips
    PDMPSamplersRBridge._select_omrf_anchor!(
        grouped_hcv_manager, grouped_hcv_cv, grouped_hcv_state.ξ.x;
        phase=:main, state=grouped_hcv_state, flow=grouped_hcv_flow)
    @test grouped_hcv_manager.node_source_checks == checks_before
    @test grouped_hcv_manager.node_source_skips == skips_before + 1
    grouped_hcv_state.ξ.x[1] += 2grouped_hcv_manager.distance_scales[1]
    PDMPSamplersRBridge._select_omrf_anchor!(
        grouped_hcv_manager, grouped_hcv_cv, grouped_hcv_state.ξ.x;
        phase=:main, state=grouped_hcv_state, flow=grouped_hcv_flow)
    @test grouped_hcv_manager.node_source_checks == checks_before + 1
    @test model_ctx.counts.event_bound_observation_visits == mapped_visits
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
    factor_covariate_spec = copy(covariate_spec)
    factor_covariate_spec["factorization"] = "person_node"
    factor_covariate_spec["covariate_component_blocks"] =
        repeat(collect(1:3); inner=9)
    factor_covariate_context = PDMPSamplersRBridge._omrf_residual_context(
        factor_covariate_spec, covariate_bank_names, anchor,
        covariate_bank_ctx.counts)
    factor_covariate_envelope =
        PDMPSamplersRBridge._build_omrf_residual_envelope(
            factor_covariate_spec, anchor, factor_covariate_context)
    @test factor_covariate_envelope isa
        PDMPSamplers.BlockSeparableResidualEnvelope
    @test size(factor_covariate_envelope.weights) == (27, 3)
    @test PDMPSamplers.n_observations(factor_covariate_envelope) == 9
    compact_covariate_spec = copy(factor_covariate_spec)
    compact_covariate_spec["bound_type"] = "covariate_local_compact"
    delete!(compact_covariate_spec, "covariate_weights")
    delete!(compact_covariate_spec, "covariate_component_blocks")
    compact_covariate_context = PDMPSamplersRBridge._omrf_residual_context(
        compact_covariate_spec, covariate_bank_names, anchor,
        covariate_bank_ctx.counts)
    compact_static_before =
        covariate_bank_ctx.counts.compact_static_constructions
    compact_covariate_envelope =
        PDMPSamplersRBridge._build_omrf_residual_envelope(
            compact_covariate_spec, anchor, compact_covariate_context)
    @test compact_covariate_envelope isa
        PDMPSamplersRBridge.OMRFCompactCovariateLocalEnvelope
    @test PDMPSamplers.n_observations(compact_covariate_envelope) == 9
    @test covariate_bank_ctx.counts.compact_static_constructions ==
        compact_static_before + 1
    shared_static = compact_covariate_envelope.component_scales!.static
    factor_state = PDMPState(
        0.0, SkeletonPoint(copy(position), fill(1.0, length(position))))
    factor_flow = ZigZag(length(position))
    for t in (0.0, 0.03, 0.2)
        screening_bound = PDMPSamplers.screening_residual_bound(
            factor_covariate_envelope, factor_state, factor_flow, t)
        materialized_bound = PDMPSamplers.total_residual_bound(
            factor_covariate_envelope, factor_state, factor_flow, t)
        @test screening_bound >= materialized_bound
        @test screening_bound ≈ materialized_bound atol=2e-12 rtol=2e-14
    end
    compact_flows = (ZigZag(9), BouncyParticle(9, 0.2),
        AdaptiveBoomerang(9; λref=0.2))
    for flow in compact_flows
        compact_state = PDMPState(0.0,
            SkeletonPoint(copy(position), fill(1.0, length(position))))
        PDMPSamplers.initialize_flow_state!(compact_state, flow)
        compact_scales = similar(compact_covariate_envelope.scales)
        expanded_scales = similar(factor_covariate_envelope.scales)
        PDMPSamplers.component_scales!(compact_scales,
            compact_covariate_envelope, compact_state, flow, 0.2)
        PDMPSamplers.component_scales!(expanded_scales,
            factor_covariate_envelope, compact_state, flow, 0.2)
        @test isapprox(sum(compact_scales),
            dot(expanded_scales, factor_covariate_envelope.totals);
            atol=2e-10, rtol=2e-11)
        compact_cell = similar(compact_covariate_envelope.cell_scales)
        expanded_cell = similar(factor_covariate_envelope.cell_scales)
        PDMPSamplers.component_cell_scales!(compact_cell,
            compact_covariate_envelope, compact_state, flow, 0.05, 0.35)
        PDMPSamplers.component_cell_scales!(expanded_cell,
            factor_covariate_envelope, compact_state, flow, 0.05, 0.35)
        @test @allocated(PDMPSamplers.component_scales!(compact_scales,
            compact_covariate_envelope, compact_state, flow, 0.2)) == 0
        @test @allocated(PDMPSamplers.component_cell_scales!(compact_cell,
            compact_covariate_envelope, compact_state, flow, 0.05, 0.35)) == 0
        compact_cell_total = sum(compact_cell)
        expanded_cell_total = dot(expanded_cell,
            factor_covariate_envelope.totals)
        if flow isa PDMPSamplers.AnyBoomerang
            # The compact harmonic aggregate is a certified, intentionally
            # conservative node-level roof; equality with the expanded
            # factor decomposition is not expected for this route.
            @test compact_cell_total >= expanded_cell_total
        else
            @test isapprox(compact_cell_total, expanded_cell_total;
                atol=2e-10, rtol=2e-11)
        end
        residual_cell = PDMPSamplers.residual_affine_cell_bound(
            compact_covariate_envelope, compact_state, flow, 0.05, 0.35)
        if flow isa Union{ZigZag,BouncyParticle}
            @test residual_cell !== nothing
            left_mass = dot(residual_cell.left_components,
                compact_covariate_envelope.totals)
            right_mass = dot(residual_cell.right_components,
                compact_covariate_envelope.totals)
            @test residual_cell.left ≈ left_mass atol=2e-12 rtol=2e-12
            @test residual_cell.left + residual_cell.slope * 0.30 ≈
                right_mass atol=2e-12 rtol=2e-12

            # The proposal-clock segment must use the same residual chord,
            # rather than the older constant closed-cell roof.
            proposal_bound = PDMPSamplers.PiecewiseAffineBound(1)
            PDMPSamplers._append_subsampling_segment!(proposal_bound,
                compact_covariate_envelope, compact_state, flow,
                0.05, 0.35, 0.0, 0.0)
            @test proposal_bound(0.05) ≈ left_mass atol=2e-12 rtol=2e-12
            @test proposal_bound(0.35) ≈ right_mass atol=2e-12 rtol=2e-12

            # At an interior time, node draws follow the endpoint-mixture
            # masses used by the same certified chord.
            interior = 0.20
            α = (interior - 0.05) / 0.30
            expected_node_mass = (1 - α) .* residual_cell.left_components .*
                compact_covariate_envelope.totals .+
                α .* residual_cell.right_components .*
                compact_covariate_envelope.totals
            expected_node_probabilities = expected_node_mass /
                sum(expected_node_mass)
            PDMPSamplers.screening_residual_bound(
                compact_covariate_envelope, compact_state, flow, interior)
            observed_nodes = zeros(Int, 3)
            proposal_rng = MersenneTwister(811)
            for _ in 1:4000
                node = PDMPSamplers._draw_component(
                    proposal_rng, compact_covariate_envelope,
                    sum(expected_node_mass))
                observed_nodes[node] += 1
            end
            @test all(abs.(observed_nodes ./ 4000 .-
                expected_node_probabilities) .< 0.04)
        else
            # Boomerang retains the generic constant-cell proposal path.
            @test residual_cell === nothing
        end
        chord_screen = PDMPSamplers.screening_residual_bound(
            compact_covariate_envelope, compact_state, flow, 0.2)
        @test chord_screen <= dot(compact_cell,
            compact_covariate_envelope.totals) + 2e-10
        PDMPSamplers._draw_distinguished(
            MersenneTwister(43), compact_covariate_envelope, 2)
        @test covariate_bank_ctx.counts.compact_person_bounds == 0
        # Exercise the exact point-screen token independently of the cached
        # endpoint-chord route used by the production grid.
        compact_covariate_envelope.cell_valid = false
        compact_screen = PDMPSamplers.screening_residual_bound(
            compact_covariate_envelope, compact_state, flow, 0.2)
        @test isapprox(compact_screen,
            compact_covariate_envelope.cumulative_masses[end];
            atol=2e-12, rtol=2e-12)
        screened_scales = copy(compact_covariate_envelope.scales)
        screening_token = PDMPSamplers.residual_sampling_generation(
            compact_covariate_envelope, compact_state)
        moved_state = copy(compact_state)
        PDMPSamplers.move_forward_time!(moved_state, 0.2, flow)
        @test PDMPSamplers.reuse_screened_residual_sampling!(
            compact_covariate_envelope, moved_state, flow, 0.0,
            screening_token)
        @test compact_covariate_envelope.scales == screened_scales
        @test !PDMPSamplers.reuse_screened_residual_sampling!(
            compact_covariate_envelope, moved_state, flow, 0.0,
            screening_token)
        # The forced path must produce the same bounds as the cached path.
        PDMPSamplers.prepare_residual_sampling!(
            compact_covariate_envelope, moved_state, flow, 0.0)
        @test compact_covariate_envelope.scales ≈ screened_scales
        node = 2
        compact_covariate_scales = compact_covariate_envelope.component_scales!
        PDMPSamplersRBridge._omrf_compact_prepare_person_node!(
            compact_covariate_scales, node)
        node_mass = sum(view(
            compact_covariate_envelope.component_scales!.person_bounds,
            node, :))
        @test isapprox(node_mass, compact_covariate_envelope.scales[node];
            atol=2e-10, rtol=2e-11)
        counts = zeros(Int, 3)
        rng = MersenneTwister(42)
        for _ in 1:4000
            factor = PDMPSamplers._draw_distinguished(
                rng, compact_covariate_envelope, node)
            counts[mod1(factor, 3)] += 1
        end
        probabilities = view(
            compact_covariate_envelope.component_scales!.person_bounds,
            node, :) ./ node_mass
        @test all(abs.(counts ./ 4000 .- probabilities) .< 0.04)
        covariate_bank_ctx.counts.compact_person_bounds = 0
        for _ in 1:1000
            PDMPSamplers._draw_distinguished(
                rng, compact_covariate_envelope, node)
        end
        @test covariate_bank_ctx.counts.compact_person_bounds == 0
        @test @allocated(PDMPSamplers._draw_distinguished(
            rng, compact_covariate_envelope, node)) == 0
        PDMPSamplers.component_cell_scales!(compact_cell,
            compact_covariate_envelope, compact_state, flow, 0.05, 0.35)
        for t in range(0.05, 0.35; length=9), person in 1:3
            candidate = copy(compact_state)
            PDMPSamplers.move_forward_time!(candidate, t, flow)
            residual_at_t = finite_difference_person_neggrad(
                candidate.ξ.x, X[person]) -
                finite_difference_person_neggrad(anchor, X[person])
            actual = PDMPSamplers.λ(candidate, residual_at_t, flow) +
                PDMPSamplers.λ(candidate, -residual_at_t, flow)
            bound = sum(PDMPSamplers.observation_residual_bound(
                compact_covariate_envelope, (node - 1) * 3 + person)
                for node in 1:3)
            @test actual <= bound * (1 + 2e-7) + 2e-8
        end
    end
    # Research-only signed BPS endpoint-mixture certificate.  It is kept
    # separate from the production compact route so that the latter retains
    # its established node sampler and proposal clock.
    signed_covariate_spec = copy(compact_covariate_spec)
    signed_covariate_spec["signed_bps"] = true
    signed_context = PDMPSamplersRBridge._omrf_residual_context(
        signed_covariate_spec, covariate_bank_names, anchor,
        covariate_bank_ctx.counts)
    signed_envelope = PDMPSamplersRBridge._new_omrf_compact_envelope(
        anchor, signed_context, shared_static)
    signed_envelope.signed_bps = true
    @test signed_envelope.signed_bps
    signed_state = PDMPState(0.0,
        SkeletonPoint(copy(position), fill(1.0, length(position))))
    signed_flow = BouncyParticle(9, 0.2)
    signed_cell = similar(signed_envelope.cell_scales)
    PDMPSamplers.component_cell_scales!(signed_cell, signed_envelope,
        signed_state, signed_flow, 0.05, 0.35)
    signed_cell_hook = PDMPSamplers.residual_affine_cell_bound(
        signed_envelope, signed_state, signed_flow, 0.05, 0.35)
    @test signed_cell_hook !== nothing
    @test signed_cell_hook.left_components[end] ≈
        signed_envelope.signed_cell_left_mass
    @test signed_cell_hook.right_components[end] ≈
        signed_envelope.signed_cell_right_mass
    @test signed_envelope.signed_cell_left_mass ≈
        signed_envelope.signed_cell_cumulative_left[end]
    @test signed_envelope.signed_cell_right_mass ≈
        signed_envelope.signed_cell_cumulative_right[end]
    # Compare every optimized endpoint mass with the original scalar
    # reference implementation before testing the chord and sampler.
    for (endpoint, cumulative) in ((0.05,
            signed_envelope.signed_cell_cumulative_left),
            (0.35, signed_envelope.signed_cell_cumulative_right))
        reference_cumulative = zeros(9)
        reference_total = 0.0
        for node in 1:3, person in 1:3
            reference_total += PDMPSamplersRBridge._omrf_person_node_pattern_bound_at(
                signed_context, anchor, signed_state, signed_flow,
                person, node, endpoint)
            reference_cumulative[(node - 1) * 3 + person] = reference_total
        end
        @test cumulative ≈ reference_cumulative atol=2e-10 rtol=2e-10
    end
    for t in range(0.05, 0.35; length=7)
        point_total = PDMPSamplersRBridge._omrf_signed_fill_endpoint!(
            signed_envelope, signed_state, signed_flow, t,
            signed_envelope.signed_screen_cumulative)
        α = (t - 0.05) / 0.30
        chord_total = (1 - α) * signed_envelope.signed_cell_left_mass +
            α * signed_envelope.signed_cell_right_mass
        @test point_total <= chord_total + 2e-9
        @test PDMPSamplers.screening_residual_bound(
            signed_envelope, signed_state, signed_flow, t) ≈ chord_total
        candidate = copy(signed_state)
        PDMPSamplers.move_forward_time!(candidate, t, signed_flow)
        for factor in 1:9
            node, person = divrem(factor - 1, 3) .+ (1, 1)
            bound = PDMPSamplers.observation_residual_bound(
                signed_envelope, factor)
            @test bound >= 0
            @test bound ≈ PDMPSamplersRBridge._omrf_signed_factor_bound_coordinates(
                signed_context, anchor, signed_envelope.screen_x,
                signed_envelope.screen_theta, person, node)
            exact_residual = zeros(Float64, length(position))
            PDMPSamplersRBridge.omrf_residual!(
                signed_context, exact_residual, candidate.ξ.x, [factor])
            exact_rate = PDMPSamplersRBridge._omrf_exact_selected_residual_abs_rate(
                candidate, exact_residual, signed_flow,
                signed_context.touched_positions)
            @test exact_rate <= bound + 2e-9
        end
    end
    signed_counts_before_adjacent =
        signed_context.counts.signed_bps.endpoint_reuses
    PDMPSamplers.component_cell_scales!(signed_cell, signed_envelope,
        signed_state, signed_flow, 0.35, 0.65)
    @test signed_context.counts.signed_bps.endpoint_reuses ==
        signed_counts_before_adjacent + 1
    PDMPSamplers.component_cell_scales!(signed_cell, signed_envelope,
        signed_state, signed_flow, 0.65, 0.95)
    @test @allocated(omrf_allocation_signed_cell!(signed_cell,
        signed_envelope, signed_state, signed_flow, 0.95, 1.25)) == 0
    PDMPSamplers.screening_residual_bound(
        signed_envelope, signed_state, signed_flow, 0.2)
    signed_counts = zeros(Int, 9)
    signed_rng = MersenneTwister(812)
    for _ in 1:6000
        factor = PDMPSamplers._draw_component(
            signed_rng, signed_envelope,
            (1 - signed_envelope.signed_screen_alpha) *
                signed_envelope.signed_cell_left_mass +
                signed_envelope.signed_screen_alpha *
                signed_envelope.signed_cell_right_mass)
        signed_counts[factor] += 1
    end
    signed_expected = (1 - signed_envelope.signed_screen_alpha) .*
        signed_envelope.signed_cell_cumulative_left .+
        signed_envelope.signed_screen_alpha .*
        signed_envelope.signed_cell_cumulative_right
    signed_expected = diff(vcat(0.0, signed_expected)) /
        signed_expected[end]
    @test all(abs.(signed_counts ./ 6000 .- signed_expected) .< 0.035)
    PDMPSamplers._draw_component(
        signed_rng, signed_envelope, signed_expected[end])
    @test @allocated(PDMPSamplers._draw_component(
        signed_rng, signed_envelope, 1.0)) == 0
    signed_preconditioned = PreconditionedBPS(9)
    signed_pre_state = PDMPState(0.0,
        SkeletonPoint(copy(position), fill(1.0, length(position))))
    PDMPSamplers.initialize_flow_state!(signed_pre_state, signed_preconditioned)
    signed_pre_cell = similar(signed_envelope.cell_scales)
    PDMPSamplers.component_cell_scales!(signed_pre_cell, signed_envelope,
        signed_pre_state, signed_preconditioned, 0.05, 0.35)
    @test PDMPSamplers.residual_affine_cell_bound(
        signed_envelope, signed_pre_state, signed_preconditioned,
        0.05, 0.35) !== nothing
    PDMPSamplers.screening_residual_bound(
        signed_envelope, signed_pre_state, signed_preconditioned, 0.2)
    signed_pre_candidate = copy(signed_pre_state)
    PDMPSamplers.move_forward_time!(
        signed_pre_candidate, 0.2, signed_preconditioned)
    pre_residual = zeros(Float64, length(position))
    PDMPSamplersRBridge.omrf_residual!(
        signed_context, pre_residual, signed_pre_candidate.ξ.x, [1])
    @test PDMPSamplersRBridge._omrf_exact_selected_residual_abs_rate(
        signed_pre_candidate, pre_residual, signed_preconditioned,
        signed_context.touched_positions) <=
        PDMPSamplers.observation_residual_bound(signed_envelope, 1) + 2e-9
    frozen_free = trues(length(position))
    frozen_free[signed_context.interaction_indices[1]] = false
    frozen_state = StickyPDMPState(0.0,
        SkeletonPoint(copy(position), fill(1.0, length(position))),
        frozen_free)
    frozen_cell = similar(signed_envelope.cell_scales)
    PDMPSamplers.component_cell_scales!(frozen_cell, signed_envelope,
        frozen_state, signed_flow, 0.05, 0.35)
    @test all(isfinite, signed_envelope.signed_cell_cumulative_left)
    @test all(isfinite, signed_envelope.signed_cell_cumulative_right)
    PDMPSamplers.component_cell_scales!(signed_cell, signed_envelope,
        signed_state, signed_flow, 0.05, 0.35)

    # Short microbenchmark used by the implementation gate.  The reference
    # loop deliberately retains the old person-node/incident-edge traversal;
    # the optimized call exercises the cached dense-field path.
    reference_cumulative = zeros(9)
    micro_envelope = PDMPSamplersRBridge._new_omrf_compact_envelope(
        anchor, signed_context, shared_static; signed_bps=true)
    reference_fill!(out, t) = begin
        total = 0.0
        for node in 1:3, person in 1:3
            total += PDMPSamplersRBridge._omrf_person_node_pattern_bound_at(
                signed_context, anchor, signed_state, signed_flow,
                person, node, t)
            out[(node - 1) * 3 + person] = total
        end
        total
    end
    optimized_cumulative = micro_envelope.signed_screen_cumulative
    for _ in 1:10
        reference_fill!(reference_cumulative, 0.2)
        PDMPSamplersRBridge._omrf_signed_fill_endpoint!(
            micro_envelope, signed_state, signed_flow, 0.2,
            optimized_cumulative)
    end
    reference_seconds = @elapsed for _ in 1:200
        reference_fill!(reference_cumulative, 0.2)
    end
    optimized_seconds = @elapsed for _ in 1:200
        PDMPSamplersRBridge._omrf_signed_fill_endpoint!(
            micro_envelope, signed_state, signed_flow, 0.2,
            optimized_cumulative)
    end
    native_gradient = zeros(9)
    native_seconds = @elapsed for _ in 1:200
        compute_gradient!(full_native_model.grad, position, native_gradient)
    end
    @test optimized_cumulative ≈ reference_cumulative atol=2e-10 rtol=2e-10
    global OMRF_SIGNED_MICROBENCHMARK = (
        reference_seconds=reference_seconds,
        optimized_seconds=optimized_seconds,
        native_gradient_seconds=native_seconds,
        iterations=200)
    @test covariate_bank_ctx.counts.compact_static_constructions ==
        compact_static_before + 1
    refreshed_envelope = PDMPSamplersRBridge._new_omrf_compact_envelope(
        copy(anchor), compact_covariate_context, shared_static)
    refreshed_state = PDMPState(
        0.0, SkeletonPoint(copy(position), fill(1.0, length(position))))
    refreshed_flow = BouncyParticle(9, 0.2)
    refreshed_scales = similar(refreshed_envelope.scales)
    original_scales = similar(compact_covariate_envelope.scales)
    PDMPSamplers.component_scales!(refreshed_scales, refreshed_envelope,
        refreshed_state, refreshed_flow, 0.2)
    PDMPSamplers.component_scales!(original_scales,
        compact_covariate_envelope, refreshed_state, refreshed_flow, 0.2)
    @test refreshed_envelope.component_scales!.static === shared_static
    @test refreshed_scales ≈ original_scales
    @test covariate_bank_ctx.counts.compact_static_constructions ==
        compact_static_before + 1
    fill!(factor_covariate_envelope.cumulative_masses, 0.0)
    PDMPSamplers.prepare_residual_sampling!(
        factor_covariate_envelope, factor_state, factor_flow, 0.03)
    @test factor_covariate_envelope.cumulative_masses[end] > 0
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

    # Compact anchor banks share the immutable covariate summaries and reuse
    # the inactive workspaces once the bank is full.
    compact_manager_subsampling = copy(covariate_subsampling)
    compact_manager_subsampling["n_observations"] = 9
    compact_manager_subsampling["subsample_size"] = 1
    compact_manager_spec = copy(compact_covariate_spec)
    compact_manager_spec["weights"] = ones(1, 9)
    compact_manager_subsampling["residual_envelope"] = compact_manager_spec
    compact_manager_ctx, compact_manager_names =
        PDMPSamplersRBridge._new_stan_subsampling_context(
            lib, full, prior, 1)
    compact_manager_model, _, _, _, compact_manager =
        PDMPSamplersRBridge._build_stan_subsampling_model(
            compact_manager_ctx, compact_manager_names,
            compact_manager_subsampling, anchor, nothing, falses(9);
            anchor_capacity=2)
    compact_first = PDMPSamplersRBridge._prepare_omrf_anchor(
        compact_manager, anchor .+ 0.01)
    PDMPSamplersRBridge._insert_omrf_anchor!(compact_manager, compact_first)
    compact_second = PDMPSamplersRBridge._prepare_omrf_anchor(
        compact_manager, anchor .- 0.02)
    PDMPSamplersRBridge._insert_omrf_anchor!(compact_manager, compact_second)
    compact_static_count =
        compact_manager_ctx.counts.compact_static_constructions
    @test compact_static_count == 1
    @test all(entry.state.envelope.component_scales!.static ===
        compact_manager.entries[1].state.envelope.component_scales!.static
        for entry in compact_manager.entries)
    compact_recycled_before = compact_manager.recycled_preparations
    compact_recycled_anchor = anchor .+ 0.03
    compact_recycled_prepared =
        PDMPSamplersRBridge._prepare_omrf_anchor(
            compact_manager, compact_recycled_anchor)
    @test compact_manager.recycled_preparations ==
        compact_recycled_before + 1
    @test compact_manager_ctx.counts.compact_static_constructions ==
        compact_static_count
    PDMPSamplersRBridge._insert_omrf_anchor!(
        compact_manager, compact_recycled_prepared)
    compact_recycled_target = anchor .+ 0.04
    compact_warm_recycled = PDMPSamplersRBridge._prepare_omrf_anchor(
        compact_manager, compact_recycled_target)
    PDMPSamplersRBridge._insert_omrf_anchor!(
        compact_manager, compact_warm_recycled)
    compact_recycled_second = compact_manager.recycled_state
    compact_allocated = @allocated PDMPSamplersRBridge._prepare_omrf_anchor(
        compact_manager, compact_recycled_target)
    @test compact_allocated == 0
    compact_fresh = PDMPSamplersRBridge._new_omrf_compact_envelope(
        compact_recycled_target,
        compact_recycled_second.residual_context,
        compact_recycled_second.envelope.component_scales!.static)
    compact_check_state = PDMPState(
        0.0, SkeletonPoint(copy(position), fill(1.0, length(position))))
    compact_check_flow = BouncyParticle(9, 0.2)
    compact_recycled_scales = similar(compact_recycled_prepared.envelope.scales)
    compact_fresh_scales = similar(compact_fresh.scales)
    PDMPSamplers.component_scales!(compact_recycled_scales,
        compact_recycled_second.envelope, compact_check_state,
        compact_check_flow, 0.2)
    PDMPSamplers.component_scales!(compact_fresh_scales, compact_fresh,
        compact_check_state, compact_check_flow, 0.2)
    @test compact_recycled_scales ≈ compact_fresh_scales
    @test compact_manager_model.grad.residual_oracle !== nothing
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
    pattern_hcv_residual = zeros(9)
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
        pattern_hcv_cell_scales = similar(pattern_hcv_envelope.cell_scales)
        grouped_hcv_cell_scales = similar(grouped_hcv_envelope.cell_scales)
        node_grouped_hcv_cell_scales =
            similar(node_grouped_hcv_envelope.cell_scales)
        factor_pattern_hcv_cell_scales =
            similar(factor_pattern_hcv_envelope.cell_scales)
        if flow isa PDMPSamplers.AnyBoomerang
            PDMPSamplers.component_cell_scales!(pattern_hcv_cell_scales,
                pattern_hcv_envelope, state, flow, 0.05, 0.35)
            PDMPSamplers.component_cell_scales!(grouped_hcv_cell_scales,
                grouped_hcv_envelope, state, flow, 0.05, 0.35)
            PDMPSamplers.component_cell_scales!(node_grouped_hcv_cell_scales,
                node_grouped_hcv_envelope, state, flow, 0.05, 0.35)
            PDMPSamplers.component_cell_scales!(
                factor_pattern_hcv_cell_scales,
                factor_pattern_hcv_envelope, state, flow, 0.05, 0.35)
            PDMPSamplers.component_scales!(pattern_hcv_envelope.scales,
                pattern_hcv_envelope, state, flow, 0.2)
            PDMPSamplers.component_scales!(factor_pattern_hcv_envelope.scales,
                factor_pattern_hcv_envelope, state, flow, 0.2)
            @test @allocated(PDMPSamplers.component_scales!(
                pattern_hcv_envelope.scales, pattern_hcv_envelope,
                state, flow, 0.2)) == 0
            @test @allocated(PDMPSamplers.component_cell_scales!(
                pattern_hcv_cell_scales, pattern_hcv_envelope,
                state, flow, 0.05, 0.35)) == 0
            @test @allocated(PDMPSamplers.component_cell_scales!(
                grouped_hcv_cell_scales, grouped_hcv_envelope,
                state, flow, 0.05, 0.35)) == 0
            @test @allocated(PDMPSamplers.component_cell_scales!(
                node_grouped_hcv_cell_scales, node_grouped_hcv_envelope,
                state, flow, 0.05, 0.35)) == 0
            @test @allocated(PDMPSamplers.total_residual_bound(
                pattern_hcv_envelope, state, flow, 0.2)) == 0
            @test @allocated(PDMPSamplers.component_scales!(
                factor_pattern_hcv_envelope.scales,
                factor_pattern_hcv_envelope, state, flow, 0.2)) == 0
            @test @allocated(PDMPSamplers.component_cell_scales!(
                factor_pattern_hcv_cell_scales,
                factor_pattern_hcv_envelope, state, flow, 0.05, 0.35)) == 0
        end
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

            if flow isa PDMPSamplers.AnyBoomerang
                pattern_hcv_oracle(pattern_hcv_residual,
                    candidate.ξ.x, [person], anchor)
                pattern_hcv_actual =
                    PDMPSamplers.λ(candidate, pattern_hcv_residual, flow) +
                    PDMPSamplers.λ(candidate, -pattern_hcv_residual, flow)
                pattern_hcv_bound = sum(partition ->
                    pattern_hcv_cell_scales[
                        pattern_hcv_envelope.groups[partition, person]],
                    axes(pattern_hcv_envelope.groups, 1))
                @test pattern_hcv_actual <=
                    pattern_hcv_bound * (1 + 2e-7) + 2e-8
                grouped_hcv_oracle(pattern_hcv_residual,
                    candidate.ξ.x, [person], anchor)
                grouped_hcv_actual =
                    PDMPSamplers.λ(candidate, pattern_hcv_residual, flow) +
                    PDMPSamplers.λ(candidate, -pattern_hcv_residual, flow)
            grouped_hcv_bound = grouped_hcv_cell_scales[
                grouped_hcv_envelope.groups[1, person]]
            @test grouped_hcv_actual <=
                grouped_hcv_bound * (1 + 2e-7) + 2e-8
            node_grouped_hcv_bound = sum(partition ->
                node_grouped_hcv_cell_scales[
                    node_grouped_hcv_envelope.groups[partition, person]],
                axes(node_grouped_hcv_envelope.groups, 1))
            @test grouped_hcv_actual <=
                node_grouped_hcv_bound * (1 + 2e-7) + 2e-8

                factor_node = mod1(person + 1, 3)
                factor = (factor_node - 1) * 3 + person
                factor_pattern_hcv_oracle(
                    factor_residual, candidate.ξ.x, [factor], anchor)
                factor_pattern_hcv_actual =
                    PDMPSamplers.λ(candidate, factor_residual, flow) +
                    PDMPSamplers.λ(candidate, -factor_residual, flow)
                factor_pattern_hcv_bound =
                    factor_pattern_hcv_cell_scales[
                        factor_pattern_hcv_envelope.groups[1, factor]]
                @test factor_pattern_hcv_actual <=
                    factor_pattern_hcv_bound * (1 + 2e-7) + 2e-8
            else
                # The whole-person node-box clock is also supported for the
                # linear flows. BPS uses the directional categorical HCV
                # certificate; Zig-Zag uses the separate norm-certified
                # coordinatewise remainder bound.
                PDMPSamplers.component_cell_scales!(
                    node_grouped_hcv_cell_scales,
                    node_grouped_hcv_envelope, state, flow, 0.05, 0.35)
                node_grouped_hcv_oracle(pattern_hcv_residual,
                    candidate.ξ.x, [person], anchor)
                node_grouped_hcv_actual =
                    PDMPSamplers.λ(candidate, pattern_hcv_residual, flow) +
                    PDMPSamplers.λ(candidate, -pattern_hcv_residual, flow)
                node_grouped_hcv_bound = sum(partition ->
                    node_grouped_hcv_cell_scales[
                        node_grouped_hcv_envelope.groups[partition, person]],
                    axes(node_grouped_hcv_envelope.groups, 1))
                @test node_grouped_hcv_actual <=
                    node_grouped_hcv_bound * (1 + 2e-7) + 2e-8
                @test @allocated(PDMPSamplers.component_cell_scales!(
                    node_grouped_hcv_cell_scales,
                    node_grouped_hcv_envelope, state, flow, 0.05, 0.35)) == 0
            end

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
        pattern_point_scales = similar(pattern_envelope.scales)
        PDMPSamplers.component_scales!(pattern_point_scales,
            pattern_envelope, state, flow, 0.2)
        @test @allocated(PDMPSamplers.component_scales!(
            pattern_point_scales, pattern_envelope, state, flow, 0.2)) == 0
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

    # Scalar-flow finite-refresh gate: both BPS and Boomerang must use the
    # same refreshed anchor for the oracle, control variate, and envelope.
    # This specifically guards the route reopened after the frozen-anchor
    # production failures; it also exercises the scalar certificate cache
    # after a refresh rather than only checking the selector counters.
    scalar_refresh_manager = manager
    scalar_refresh_manager.main_refresh_distance = 0.01
    for (offset, scalar_flow) in enumerate(
            (BouncyParticle(9, 0.2), Boomerang(9)))
        scalar_position = refresh_position .+ (0.05 * offset) .+
            collect(range(0.01, 0.09; length=9))
        scalar_state = PDMPState(0.0,
            SkeletonPoint(copy(scalar_position), copy(velocity)))
        PDMPSamplers.initialize_flow_state!(scalar_state, scalar_flow)
        preparations_before = scalar_refresh_manager.preparations
        @test PDMPSamplersRBridge._select_omrf_anchor!(
            scalar_refresh_manager, managed_model.grad, scalar_position;
            phase=:main, state=scalar_state, flow=scalar_flow)
        @test scalar_refresh_manager.preparations == preparations_before + 1
        @test scalar_refresh_manager.oracle.anchor == scalar_position
        @test managed_model.grad.anchor == scalar_position
        @test managed_model.grad.envelope.component_scales!.anchor ===
            scalar_refresh_manager.oracle.anchor
        scalar_residual = zeros(9)
        scalar_refresh_manager.oracle(
            scalar_residual, scalar_position, [1], scalar_position)
        @test iszero(norm(scalar_residual))
        @test PDMPSamplers.total_residual_bound(
            managed_model.grad.envelope, scalar_state, scalar_flow, 0.0) >= 0
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

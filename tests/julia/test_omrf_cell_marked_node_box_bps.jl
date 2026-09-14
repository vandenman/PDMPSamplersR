using Test, Random, LinearAlgebra, PDMPSamplers, BridgeStan
include(normpath(joinpath(@__DIR__, "..", "..", "inst", "julia",
    "main_interface_function.jl")))
using .PDMPSamplersRBridge

@inline _bps_cell_marked_draw!(fixture) = PDMPSamplers.draw_subset!(
    fixture.rng, fixture.cv, fixture.D, fixture.B)
@inline _bps_cell_marked_selected_bound!(fixture) =
    PDMPSamplers.subsampling_residual_subset_bound(fixture.oracle,
        fixture.state, fixture.gradient, fixture.flow, fixture.D, fixture.M,
        fixture.subset, fixture.scale)
@inline function _bps_integrated_rejection_path!(fixture)
    PDMPSamplers.next_event_time(fixture.rng, fixture.model, fixture.flow,
        fixture.alg, fixture.state, fixture.cache, fixture.stats,
        0.2, false, :horizon_hit)
    return nothing
end

function _bps_node_box_spec(X, seen; categorical_residual_cap=false)
    N, P = size(X)
    membership = zeros(Int, P, N)
    nodes = Int[]
    minima = Vector{Vector{Float64}}()
    maxima = Vector{Vector{Float64}}()
    offset = 0
    for node in 1:P
        local_membership = [1 + mod(person - 1, min(2, N)) for person in 1:N]
        for group in 1:maximum(local_membership)
            people = findall(==(group), local_membership)
            lower = zeros(P); upper = zeros(P)
            for neighbour in 1:P
                neighbour == node && continue
                lower[neighbour] = minimum(X[people, neighbour])
                upper[neighbour] = maximum(X[people, neighbour])
            end
            push!(nodes, node); push!(minima, lower); push!(maxima, upper)
        end
        membership[node, :] .= offset .+ local_membership
        offset += maximum(local_membership)
    end
    spec = Dict{String,Any}(
        "type"=>"omrf", "backend"=>"analytic", "X"=>X, "seen"=>seen,
        "thresholds"=>"thresholds_0", "interactions"=>"interactions_0",
        "bound_type"=>"cell_marked_node_box_first_order",
        "factorization"=>"person_node", "factor_sampling"=>"size_biased",
        "prior_backend"=>"analytic", "use_hcv"=>false,
        "analytic_prior"=>Dict{String,Any}("type"=>"gaussian",
            "threshold_alpha"=>1.0, "threshold_beta"=>1.0,
            "interaction_scale"=>1.0),
        "exact_selected_refinement"=>true,
        "node_box_membership"=>membership, "node_box_nodes"=>nodes,
        "node_box_minima"=>reduce(vcat, permutedims.(minima)),
        "node_box_maxima"=>reduce(vcat, permutedims.(maxima)))
    categorical_residual_cap &&
        (spec["categorical_residual_cap"] = true)
    return spec
end

@testset "OMRF cell-marked node-box linear BPS" begin
    categorical_bound =
        PDMPSamplersRBridge._omrf_categorical_residual_bound
    @test categorical_bound(2.0, 3.0, true) == nextfloat(1.5)
    @test categorical_bound(8.0, 3.0, true) == nextfloat(3.0)
    @test categorical_bound(8.0, 3.0, true) <=
        categorical_bound(8.0, 3.0, false)
    @test categorical_bound(0.0, 3.0, true) == 0.0
    @test categorical_bound(3.0, 0.0, true) == 0.0
    # Adding a common reference-category shift cannot change either range.
    a = [-1.0, 0.0, 2.0]; b = [-3.0, 1.0, 4.0]
    shifted_a = a .+ 17.0; shifted_b = b .- 9.0
    @test categorical_bound(maximum(a)-minimum(a),
        maximum(b)-minimum(b), true) ==
        categorical_bound(maximum(shifted_a)-minimum(shifted_a),
            maximum(shifted_b)-minimum(shifted_b), true)
    root = normpath(joinpath(@__DIR__, "..", ".."))
    lib = PDMPSamplersRBridge._compile_model_with_header(
        joinpath(root, "inst", "stan", "omrf", "omrf.stan"),
        joinpath(root, "inst", "stan", "pdmp_subsample.hpp"))
    full = joinpath(root, "tests", "stan", "omrf_subsampling_full.json")
    prior = joinpath(root, "tests", "stan", "omrf_subsampling_prior.json")
    ctx, names = PDMPSamplersRBridge._new_stan_subsampling_context(
        lib, full, prior, 1)
    X = [0 1 2; 2 0 1; 1 2 0]
    seen = [3, 3, 3]; N, P = size(X)
    spec = _bps_node_box_spec(X, seen)
    anchor = zeros(9)
    oracle = PDMPSamplersRBridge.OMRFAnalyticSubsamplingOracle(
        ctx, spec, names, anchor)
    envelope = PDMPSamplersRBridge._build_omrf_residual_envelope(
        spec, oracle.anchor, oracle.residual_context)
    @test Base.summarysize(envelope) < 1_000_000
    @test PDMPSamplers.n_observations(envelope) == N * P

    scale = exp.(range(-0.7, 0.6; length=9))
    flow = PDMPSamplers.PreconditionedBPS(9; refresh_rate=0.7, scale=scale)
    rng = Xoshiro(48031)
    state = StickyPDMPState(0.0,
        SkeletonPoint(0.25randn(rng, 9),
            PDMPSamplers.initialize_velocity(rng, flow, 9)), trues(9))
    selected = zeros(9)
    maximum_exact_point_violation = -Inf
    maximum_point_cell_violation = -Inf
    maximum_interior_violation = -Inf

    # The deterministic first-order CV is the cached likelihood gradient plus
    # the explicit prior.  Its direct interval roof is independent of N after
    # anchor preparation and must survive frozen interactions and anchor
    # crossings without touching the residual marked law.
    PDMPSamplers.set_active_set!(oracle, state.free)
    deterministic = zeros(9)
    prior_gradient = zeros(9)
    oracle(deterministic, state.ξ.x)
    PDMPSamplersRBridge._omrf_prior_gradient!(
        prior_gradient, oracle.analytic_prior, state.ξ.x)
    @test deterministic .- prior_gradient ≈
        oracle.full_anchor .- oracle.prior_anchor atol=2e-12 rtol=2e-12
    direct_bound = PDMPSamplers.deterministic_rate_cell_bound(
        oracle, state, flow, -0.25, 0.55)
    for t in range(-0.25, 0.55; length=129)
        candidate = move_forward_time(state, t, flow)
        oracle(deterministic, candidate.ξ.x)
        @test PDMPSamplers.λ(candidate, deterministic, flow) <= direct_bound
    end
    PDMPSamplers.deterministic_rate_cell_bound(
        oracle, state, flow, -0.25, 0.55)
    @test @allocated(PDMPSamplers.deterministic_rate_cell_bound(
        oracle, state, flow, -0.25, 0.55)) == 0

    # The full-gradient BPS comparator exploits convexity of the likelihood
    # along a linear path: its signed rate is largest at the right endpoint.
    # Only the prior needs coordinate intervals.  The older general
    # observation-wise interval construction remains a valid, looser oracle
    # for comparison and for nonlinear/coordinatewise dynamics.
    full_native = PDMPSamplersRBridge.build_omrf_full_model(spec, names, 9)
    full_gradient = zeros(9)
    full_bound = PDMPSamplers.deterministic_rate_cell_bound(
        full_native.grad, state, flow, -0.25, 0.55)
    general_full_bound = PDMPSamplersRBridge._omrf_full_deterministic_cell_bound_general(
        full_native.grad.f, nothing, state, flow, -0.25, 0.55)
    previous_likelihood_rate = -Inf
    for t in range(-0.25, 0.55; length=129)
        candidate = move_forward_time(state, t, flow)
        PDMPSamplers.compute_gradient!(
            full_native.grad, candidate.ξ.x, full_gradient)
        @test PDMPSamplers.λ(candidate, full_gradient, flow) <= full_bound
        fill!(prior_gradient, 0.0)
        PDMPSamplersRBridge._omrf_prior_gradient!(
            prior_gradient, full_native.grad.f.prior, candidate.ξ.x)
        likelihood_rate = dot(candidate.ξ.θ,
            full_gradient .- prior_gradient)
        if isfinite(previous_likelihood_rate)
            @test likelihood_rate >= previous_likelihood_rate -
                2e-11 * max(1.0, abs(likelihood_rate),
                    abs(previous_likelihood_rate))
        end
        previous_likelihood_rate = likelihood_rate
    end
    @test full_bound <= general_full_bound
    PDMPSamplers.deterministic_rate_cell_bound(
        full_native.grad, state, flow, -0.25, 0.55)
    @test @allocated(PDMPSamplers.deterministic_rate_cell_bound(
        full_native.grad, state, flow, -0.25, 0.55)) == 0

    # Active-set invalidation is explicit.  A frozen interaction has constant
    # position and zero physical velocity throughout the certified cell.
    frozen = oracle.residual_context.interaction_indices[1]
    state.free[frozen] = false
    state.ξ.x[frozen] = 0.0
    state.ξ.θ[frozen] = 0.0
    PDMPSamplers.set_active_set!(oracle, state.free)
    frozen_bound = PDMPSamplers.deterministic_rate_cell_bound(
        oracle, state, flow, -0.25, 0.55)
    for t in range(-0.25, 0.55; length=65)
        candidate = move_forward_time(state, t, flow)
        oracle(deterministic, candidate.ξ.x)
        @test PDMPSamplers.λ(candidate, deterministic, flow) <= frozen_bound
    end
    state.free .= true
    PDMPSamplers.set_active_set!(oracle, state.free)

    # A separately prepared anchor owns distinct constant and interval
    # workspaces; changing anchors cannot reuse the old deterministic cache.
    shifted_anchor = fill(0.07, 9)
    shifted_oracle = PDMPSamplersRBridge.OMRFAnalyticSubsamplingOracle(
        ctx, spec, names, shifted_anchor)
    PDMPSamplers.set_active_set!(shifted_oracle, state.free)
    shifted_bound = PDMPSamplers.deterministic_rate_cell_bound(
        shifted_oracle, state, flow, -0.25, 0.55)
    for t in range(-0.25, 0.55; length=65)
        candidate = move_forward_time(state, t, flow)
        shifted_oracle(deterministic, candidate.ξ.x)
        @test PDMPSamplers.λ(candidate, deterministic, flow) <= shifted_bound
    end

    # A Cauchy gradient is nonmonotone.  Intervals crossing both stationary
    # points +/-s must include the exact extrema +/-1/s.
    cauchy = PDMPSamplersRBridge.OMRFIndependentAnalyticPrior(
        collect(1:3), collect(4:6), 1.0, 1.0, :cauchy, 0.7)
    c_lo = fill(-1.4, 9); c_hi = fill(1.4, 9)
    c_glo = zeros(9); c_ghi = zeros(9)
    PDMPSamplersRBridge._omrf_prior_intervals!(
        c_glo, c_ghi, cauchy, c_lo, c_hi)
    @test all(c_glo[4:6] .<= -inv(0.7))
    @test all(c_ghi[4:6] .>= inv(0.7))

    # Linear endpoint construction is exact coordinate-wise; randomized
    # checks include anchor crossings, inactive interactions, and interior
    # times. The factor certificate is the categorical covariance range bound.
    for iteration in 1:20_001
        state.ξ.x .= 0.5randn(rng, 9)
        state.ξ.θ .= PDMPSamplers.initialize_velocity(rng, flow, 9)
        state.free .= true
        if iteration % 5 == 0
            frozen = oracle.residual_context.interaction_indices[
                1 + mod(iteration, length(oracle.residual_context.interaction_indices))]
            state.free[frozen] = false
            state.ξ.x[frozen] = 0.0
            state.ξ.θ[frozen] = 0.0
        end
        left = iteration % 7 == 0 ? -0.5 : rand(rng)
        right = left + 0.001 + 0.8rand(rng)
        envelope.component_scales!.n_cached_cells = 0
        PDMPSamplers.component_cell_scales!(envelope.cell_scales,
            envelope, state, flow, left, right)
        factor = rand(rng, 1:(N * P))
        component = envelope.groups[1, factor]
        for t in (left, right, left + rand(rng) * (right-left))
            B = PDMPSamplers.screening_residual_bound(
                envelope, state, flow, t)
            @test B == envelope.component_scales!.cell_total_masses[
                envelope.component_scales!.active_cell]
            point = PDMPSamplersRBridge._omrf_cell_marked_factor_point_bound(
                envelope.component_scales!, factor)
            cell = envelope.scales[component]
            candidate = move_forward_time(state, t, flow)
            fill!(selected, 0.0)
            oracle(selected, candidate.ξ.x, [factor], oracle.anchor)
            exact = PDMPSamplersRBridge._omrf_exact_selected_residual_abs_rate(
                candidate, selected, flow,
                oracle.residual_context.touched_positions)
            maximum_exact_point_violation = max(
                maximum_exact_point_violation, exact-point)
            maximum_point_cell_violation = max(
                maximum_point_cell_violation, point-cell)
            maximum_interior_violation = max(maximum_interior_violation,
                exact-cell)
            @test exact <= point * (1 + 2e-7) + 2e-8
            @test point <= cell * (1 + 64eps(Float64))
        end
    end
    @test maximum_exact_point_violation <= 2e-8
    @test maximum_point_cell_violation <= 0.0
    @test maximum_interior_violation <= 2e-8

    # The capped construction encloses points for mixed node references and
    # frozen coordinates, and never exceeds the corresponding legacy roof.
    scales = envelope.component_scales!
    frozen_state = deepcopy(state)
    frozen_index = oracle.residual_context.interaction_indices[1]
    frozen_state.free[frozen_index] = false
    frozen_state.ξ.θ[frozen_index] = 0.0
    scales.base.categorical_residual_cap = false
    scales.n_cached_cells = 0
    PDMPSamplers.component_cell_scales!(envelope.cell_scales,
        envelope, frozen_state, flow, 0.0, 0.4)
    legacy_cell = copy(envelope.cell_scales)
    scales.base.categorical_residual_cap = true
    scales.n_cached_cells = 0
    PDMPSamplers.component_cell_scales!(envelope.cell_scales,
        envelope, frozen_state, flow, 0.0, 0.4)
    @test all(envelope.cell_scales .<= legacy_cell .+ 1e-12)
    for factor in 1:(N * P), t in (0.0, 0.2, 0.4)
        PDMPSamplers.screening_residual_bound(
            envelope, frozen_state, flow, t)
        point = PDMPSamplersRBridge._omrf_cell_marked_factor_point_bound(
            scales, factor)
        component = envelope.groups[1, factor]
        candidate = move_forward_time(frozen_state, t, flow)
        fill!(selected, 0.0)
        oracle(selected, candidate.ξ.x, [factor], oracle.anchor)
        exact = PDMPSamplersRBridge._omrf_exact_selected_residual_abs_rate(
            candidate, selected, flow,
            oracle.residual_context.touched_positions)
        @test exact <= point * (1 + 2e-7) + 2e-8
        @test point <= envelope.scales[component] * (1 + 64eps(Float64))
    end
    scales.base.categorical_residual_cap = false
    scales.n_cached_cells = 0

    # Enumerated and Monte Carlo component/factor balance, including unequal
    # component sizes and zero component mass.
    envelope.component_scales!.base.categorical_residual_cap = true
    envelope.component_scales!.n_cached_cells = 0
    PDMPSamplers.component_cell_scales!(envelope.cell_scales,
        envelope, state, flow, 0.0, 0.4)
    B = PDMPSamplers.screening_residual_bound(envelope, state, flow, 0.2)
    cv = SubsampledControlVariate((out, x)->fill!(out, 0.0), oracle,
        envelope, anchor, 1)
    expected = [envelope.scales[envelope.groups[1, factor]] / B
        for factor in 1:(N * P)]
    counts = zeros(Int, N * P)
    draw_rng = Xoshiro(48032)
    for _ in 1:100_000
        PDMPSamplers.draw_subset!(draw_rng, cv, 0.0, B)
        counts[cv.subset[1]] += 1
    end
    @test maximum(abs.(counts ./ sum(counts) .- expected)) < 0.007
    for factor in 1:(N * P)
        component = envelope.groups[1, factor]
        nc = envelope.totals[component]
        cell = envelope.scales[component]
        point = PDMPSamplersRBridge._omrf_cell_marked_factor_point_bound(
            envelope.component_scales!, factor)
        @test B * (nc*cell/B) * inv(nc) * (point/cell) ≈ point
    end
    @test 4.0 * (0.0 / (4.0 + B)) == 0.0

    # Deterministic and residual mixture branches retain the same cell roof.
    marks = envelope.component_scales!.component_marks
    PDMPSamplers.draw_subset!(draw_rng, cv, 1e300, B)
    @test envelope.component_scales!.component_marks == marks
    PDMPSamplers.draw_subset!(draw_rng, cv, 0.0, B)
    @test envelope.component_scales!.component_marks == marks + 1

    fill!(deterministic, 0.0)
    residual = zeros(9)
    candidate = move_forward_time(state, 0.2, flow)
    oracle(deterministic, candidate.ξ.x)
    deterministic_rate = PDMPSamplers.λ(candidate, deterministic, flow)
    marked_rate = 0.0
    recovered_rate = 0.0
    population = N * P
    for factor in 1:population
        fill!(residual, 0.0)
        oracle(residual, candidate.ξ.x, [factor], oracle.anchor)
        exact_mark = PDMPSamplers.λ(candidate,
            deterministic .+ population .* residual, flow)
        component = envelope.groups[1, factor]
        selected_roof = deterministic_rate + population *
            envelope.scales[component]
        probability = PDMPSamplers._selected_subset_probability(
            cv, deterministic_rate, B, selected_roof)
        marked_rate += exact_mark / population
        recovered_rate += (deterministic_rate + B) * probability *
            exact_mark / selected_roof
    end
    @test recovered_rate ≈ marked_rate atol=2e-12 rtol=2e-12
    envelope.component_scales!.base.categorical_residual_cap = false

    # Sticky diagonal-preconditioned BPS reflection reverses the active
    # directional rate, preserves metric kinetic energy, and leaves frozen
    # coordinates untouched.
    reflection_state = deepcopy(state)
    reflection_state.free[oracle.residual_context.interaction_indices[1]] = false
    gradient = randn(rng, 9)
    before = copy(reflection_state.ξ.θ)
    free = reflection_state.free
    before_rate = dot(before[free], gradient[free])
    before_energy = sum(abs2(before[i]/scale[i]) for i in eachindex(before) if free[i])
    cache = (; z=zeros(9))
    PDMPSamplers.reflect!(rng, reflection_state, gradient, flow, cache)
    @test dot(reflection_state.ξ.θ[free], gradient[free]) ≈ -before_rate
    @test sum(abs2(reflection_state.ξ.θ[i]/scale[i])
        for i in eachindex(before) if free[i]) ≈ before_energy
    @test reflection_state.ξ.θ[.!free] == before[.!free]

    # Cell construction, cached lookup, component marking, selected point
    # evaluation, exact refinement, and reflection are allocation-free.
    envelope.component_scales!.n_cached_cells = 0
    PDMPSamplers.component_cell_scales!(envelope.cell_scales,
        envelope, state, flow, 0.0, 0.3)
    B = PDMPSamplers.screening_residual_bound(envelope, state, flow, 0.1)
    PDMPSamplers.draw_subset!(draw_rng, cv, 0.0, B)
    oracle(deterministic, state.ξ.x)
    D = PDMPSamplers.λ(state, deterministic, flow)
    selected_roof = D + (N * P) * PDMPSamplers.observation_residual_bound(
        envelope, cv.subset[1])
    PDMPSamplers.subsampling_residual_subset_bound(oracle, state,
        deterministic, flow, D, selected_roof, cv.subset, N * P)
    @test @allocated(PDMPSamplers.component_cell_scales!(
        envelope.cell_scales, envelope, state, flow, 0.3, 0.6)) == 0
    @test @allocated(PDMPSamplers.screening_residual_bound(
        envelope, state, flow, 0.4)) == 0
    draw_fixture = (; rng=draw_rng, cv, D=0.0, B)
    selected_fixture = (; oracle, state, gradient=deterministic, flow, D,
        M=selected_roof, subset=cv.subset, scale=N * P)
    _bps_cell_marked_draw!(draw_fixture)
    _bps_cell_marked_selected_bound!(selected_fixture)
    @test @allocated(_bps_cell_marked_draw!(draw_fixture)) == 0
    @test @allocated(_bps_cell_marked_selected_bound!(selected_fixture)) == 0
    @test @allocated(PDMPSamplers.reflect!(rng, reflection_state,
        gradient, flow, cache)) == 0

    # Integrated wrapper regression: use the real statistics-wrapped CV and a
    # StickyPDMPState.  After all branches are prepared, the rejection-heavy
    # event search must not allocate per proposal.
    event_strategy = PDMPSamplers.Sticky(
        PDMPSamplers.GridThinningStrategy(N=8, N_min=8, t_max=0.2,
            lazy=false, bound_violation=:throw), fill(Inf, 9))
    event_rng = Xoshiro(48039)
    event_state, event_model, event_alg, event_cache, event_stats =
        PDMPSamplers.initialize_state(event_rng, flow,
            PDMPSamplers.PDMPModel(9, cv), event_strategy, 0.0, state.ξ;
            statistic_counter=PDMPSamplers.DevelStatisticCounter,
            initial_free=state.free, initial_stored_velocity=zeros(9))
    integrated = (; rng=event_rng, model=event_model, flow,
        alg=event_alg.inner_alg_state, state=event_state,
        cache=event_cache, stats=event_stats)
    for _ in 1:5
        _bps_integrated_rejection_path!(integrated)
    end
    before_candidates = PDMPSamplers._get_counter_clock_candidates(event_stats)
    GC.gc()
    integrated_bytes = @allocated _bps_integrated_rejection_path!(integrated)
    after_candidates = PDMPSamplers._get_counter_clock_candidates(event_stats)
    @test after_candidates > before_candidates
    @test integrated_bytes == 0
    println("bps_cell_marked_exact_point_max_violation=",
        maximum_exact_point_violation)
    println("bps_cell_marked_point_cell_max_violation=",
        maximum_point_cell_violation)
end

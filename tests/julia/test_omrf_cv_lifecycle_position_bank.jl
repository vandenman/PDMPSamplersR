using Test, Random, LinearAlgebra, PDMPSamplers, BridgeStan
include(normpath(joinpath(@__DIR__, "..", "..", "inst", "julia",
    "main_interface_function.jl")))
using .PDMPSamplersRBridge

function audit_spec(X, seen)
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
    Dict{String,Any}(
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
end

function live_cv(grad)
    grad isa SubsampledControlVariate && return grad
    hasproperty(grad, :f) && return live_cv(getproperty(grad, :f))
    hasproperty(grad, :posterior_grad!) &&
        return live_cv(getproperty(grad, :posterior_grad!))
    return nothing
end

@testset "OMRF live CV lifecycle and position-selected bank" begin
    root = normpath(joinpath(@__DIR__, "..", ".."))
    lib = PDMPSamplersRBridge._compile_model_with_header(
        joinpath(root, "inst", "stan", "omrf", "omrf.stan"),
        joinpath(root, "inst", "stan", "pdmp_subsample.hpp"))
    full = joinpath(root, "tests", "stan", "omrf_subsampling_full.json")
    prior = joinpath(root, "tests", "stan", "omrf_subsampling_prior.json")
    X = [0 1 2; 2 0 1; 1 2 0]
    spec = audit_spec(X, [3, 3, 3])
    sub = Dict{String,Any}("n_observations"=>9, "subsample_size"=>1,
        "anchor"=>zeros(9), "residual_envelope"=>spec)

    function built(capacity=1, m=1)
        ctx, names = PDMPSamplersRBridge._new_stan_subsampling_context(
            lib, full, prior, m)
        local_sub = copy(sub)
        local_sub["subsample_size"] = m
        PDMPSamplersRBridge._build_stan_subsampling_model(
            ctx, names, local_sub, zeros(9), nothing, falses(9);
            anchor_capacity=capacity)
    end

    function entry_case(; separate=false, warmup=0.04, updates=true)
        model, _, names, adapter, manager = built(1)
        adapter.update_dt = updates ? 0.01 : Inf
        flow = BouncyParticle(9, 0.5)
        alg = GridThinningStrategy(; N=4, t_max=0.02, lazy=false,
            bound_violation=:throw)
        observed = Ref(false)
        PDMPSamplers.set_phase_first_proposal_observer!(
            (phase, state, live_model, flow, alg, cache) -> begin
                phase === :main || return
                cv = live_cv(live_model.grad)
                @test cv !== nothing
                @test cv.anchor == manager.oracle.anchor
                @test PDMPSamplers._trajectory_scale_anchor(cv.envelope) ==
                    manager.oracle.anchor
                @test cv.residual_oracle.f.anchor == manager.oracle.anchor
                observed[] = true
            end)
        warm_model = separate ?
            PDMPSamplersRBridge.build_omrf_full_model(spec, names, 9) : nothing
        try
            pdmp_sample(SkeletonPoint(fill(0.15, 9), fill(0.4, 9)), flow,
                model, alg, 0.0, warmup + 0.02, warmup;
                adapter, warmup_model=warm_model,
                warmup_algorithm=separate ? deepcopy(alg) : nothing,
                progress=false, seed=7301,
                statistic_counter=PDMPSamplers.DevelStatisticCounter)
        finally
            PDMPSamplers.set_phase_first_proposal_observer!(nothing)
        end
        @test observed[]
        @test manager.oracle.anchor != zeros(9)
    end

    # Actual sampler-entry regressions: separate full warmup, ordinary
    # subsampled warmup, and a zero-warmup restart.  Capacity one forces bank
    # replacement without an index change.
    entry_case(separate=true, warmup=0.04, updates=true)
    entry_case(separate=false, warmup=0.04, updates=true)
    entry_case(separate=false, warmup=0.0, updates=false)

    # Position-bank finalization is explicit and independent of the adaptive
    # route's historical `main_refresh_distance == Inf` sentinel, including
    # the separate full-gradient warmup lifecycle.
    fixed_model, _, fixed_names, fixed_adapter, fixed_manager = built(4)
    fixed_manager.selection_policy = :position
    fixed_manager.main_refresh_distance = 0.2
    fixed_adapter.update_dt = Inf
    fixed_warm_model = PDMPSamplersRBridge.build_omrf_full_model(
        spec, fixed_names, 9)
    pdmp_sample(SkeletonPoint(fill(0.15, 9), fill(0.4, 9)),
        BouncyParticle(9, 0.5), fixed_model,
        GridThinningStrategy(; N=4, t_max=0.02, lazy=false,
            bound_violation=:throw), 0.0, 0.06, 0.04;
        adapter=fixed_adapter, warmup_model=fixed_warm_model,
        warmup_algorithm=GridThinningStrategy(; N=4, t_max=0.02,
            lazy=false, bound_violation=:throw),
        progress=false, seed=7302,
        statistic_counter=PDMPSamplers.DevelStatisticCounter)
    @test fixed_manager.capacity == 4
    @test 1 <= length(fixed_manager.entries) <= 4
    @test all(isfinite, fixed_manager.distance_scales)
    @test all(>(0), fixed_manager.distance_scales)
    @test fixed_manager.preparations ==
        fixed_manager.warmup_finish_preparations

    # Opt-in protection: explicitly requesting the existing adaptive policy
    # produces the same seeded event stream as leaving the new feature off.
    adaptive_model_a, _, _, adaptive_adapter_a, _ = built(1)
    adaptive_model_b, _, _, adaptive_adapter_b, adaptive_manager_b = built(1)
    adaptive_manager_b.selection_policy = :adaptive
    function adaptive_path(model, adapter)
        pdmp_sample(SkeletonPoint(fill(0.12, 9), fill(0.3, 9)),
            BouncyParticle(9, 0.5), model,
            GridThinningStrategy(; N=4, t_max=0.02, lazy=false,
                bound_violation=:throw), 0.0, 0.03, 0.0;
            adapter, progress=false, seed=7303,
            statistic_counter=PDMPSamplers.DevelStatisticCounter).traces[1]
    end
    adaptive_trace_a = adaptive_path(adaptive_model_a, adaptive_adapter_a)
    adaptive_trace_b = adaptive_path(adaptive_model_b, adaptive_adapter_b)
    @test adaptive_trace_a.times == adaptive_trace_b.times
    @test adaptive_trace_a.positions == adaptive_trace_b.positions
    @test adaptive_trace_a.velocities == adaptive_trace_b.velocities

    # Direct finalization regression for the same-index replacement branch.
    model, _, _, _, manager = built(1)
    stats = PDMPSamplers.DevelStatisticCounter()
    wrapped = PDMPSamplers.with_stats(model, stats)
    forced = fill(0.31, 9)
    prepared = PDMPSamplersRBridge._prepare_omrf_anchor(manager, forced)
    @test PDMPSamplersRBridge._insert_omrf_anchor!(manager, prepared) == 1
    manager.active_idx = 1
    PDMPSamplersRBridge._install_omrf_anchor!(manager, prepared)
    @test PDMPSamplers._trajectory_scale_anchor(wrapped.grad.envelope) != forced
    @test PDMPSamplersRBridge._finish_omrf_warmup!(
        manager, wrapped.grad, forced)
    @test wrapped.grad.anchor == forced
    @test PDMPSamplers._trajectory_scale_anchor(wrapped.grad.envelope) == forced

    # Exact weighted Voronoi crossings.  The bank is prepared once, then the
    # position policy performs no preparation or replacement.
    model, _, _, adapter, manager = built(3)
    cv = model.grad
    anchors = ([0.0; zeros(8)], [2.0; zeros(8)], [-2.0; zeros(8)])
    for a in anchors[2:3]
        PDMPSamplersRBridge._insert_omrf_anchor!(manager,
            PDMPSamplersRBridge._prepare_omrf_anchor(manager, a))
    end
    manager.selection_policy = :position
    manager.distance_scales .= 1.0
    manager.active_idx = 1
    PDMPSamplers.refresh_anchor_owned!(cv, manager.entries[1].state.anchor)
    preparations = manager.preparations
    state = StickyPDMPState(0.0,
        SkeletonPoint(zeros(9), [1.0; zeros(8)]), trues(9))
    state.stored_velocity[4] = -0.7
    # Statistics wrapping copies the envelope.  The real adapter horizon is
    # the pre-proposal hook and must attach its counters before event search.
    cv.envelope.component_scales!.base.categorical_residual_cap = true
    preproposal = PDMPSamplers.with_stats(
        model, PDMPSamplers.DevelStatisticCounter())
    preproposal_scales = preproposal.grad.envelope.component_scales!
    @test preproposal_scales.base.context !== manager.oracle.residual_context
    @test preproposal_scales.base.categorical_residual_cap
    PDMPSamplers.adaptation_horizon(adapter, state,
        BouncyParticle(9, 0.0), preproposal.grad, :main, stats)
    @test preproposal_scales.base.context === manager.oracle.residual_context
    @test preproposal_scales.base.categorical_residual_cap
    cv.envelope.component_scales!.base.categorical_residual_cap = false
    before = (copy(state.ξ.x), copy(state.ξ.θ), copy(state.free),
        copy(state.stored_velocity))
    horizon = PDMPSamplersRBridge._omrf_position_anchor_horizon!(
        manager, state, BouncyParticle(9, 0.0), :main)
    @test horizon == 1.0
    move_forward_time!(state, horizon, BouncyParticle(9, 0.0))
    PDMPSamplers.adapt!(Xoshiro(1), adapter, state,
        BouncyParticle(9, 0.0), cv, nothing; phase=:main, stats,
        event_type=:anchor_selection_boundary, anchor_boundary_won=true)
    @test PDMPSamplers.did_gradient_adapt(adapter)
    @test manager.active_idx == 2
    @test state.ξ.θ == before[2]
    @test state.free == before[3]
    @test state.stored_velocity == before[4]
    @test manager.preparations == preparations
    @test cv.anchor == manager.oracle.anchor ==
        PDMPSamplers._trajectory_scale_anchor(cv.envelope)

    # The same frozen global bank has exact harmonic Voronoi deadlines for
    # Boomerang.  A transverse crossing is owned at x=1, while a tangent
    # contact and a frozen coordinate do not create computational events.
    harmonic_flow = Boomerang(Diagonal(ones(9)), zeros(9), 0.0)
    harmonic_state = StickyPDMPState(0.0,
        SkeletonPoint(zeros(9), [2.0; zeros(8)]), trues(9))
    manager.active_idx = 1
    harmonic_horizon =
        PDMPSamplersRBridge._omrf_position_anchor_horizon!(
            manager, harmonic_state, harmonic_flow, :main)
    @test harmonic_horizon ≈ pi / 6 atol=1e-13
    @test @allocated(PDMPSamplersRBridge._omrf_position_anchor_horizon!(
        manager, harmonic_state, harmonic_flow, :main)) == 0
    constant, cosine, sine =
        PDMPSamplersRBridge._omrf_harmonic_anchor_difference_coefficients(
            manager, harmonic_state, harmonic_flow, 1, 2)
    for t in (0.0, 0.17, harmonic_horizon, 1.1)
        moved = move_forward_time(harmonic_state, t, harmonic_flow)
        explicit = PDMPSamplersRBridge._omrf_anchor_distance2(
            manager, moved.ξ.x, manager.entries[2].state) -
            PDMPSamplersRBridge._omrf_anchor_distance2(
                manager, moved.ξ.x, manager.entries[1].state)
        @test constant + cosine * cos(t) + sine * sin(t) ≈ explicit atol=2e-13
    end
    harmonic_before = (copy(harmonic_state.ξ.x), copy(harmonic_state.ξ.θ),
        copy(harmonic_state.free), copy(harmonic_state.stored_velocity))
    move_forward_time!(harmonic_state, harmonic_horizon, harmonic_flow)
    PDMPSamplers.adapt!(Xoshiro(2), adapter, harmonic_state,
        harmonic_flow, cv, nothing; phase=:main, stats,
        event_type=:anchor_selection_boundary, anchor_boundary_won=true)
    @test manager.active_idx == 2
    @test harmonic_state.ξ.θ[1] ≈ sqrt(3.0) atol=2e-13
    @test harmonic_state.free == harmonic_before[3]
    @test harmonic_state.stored_velocity == harmonic_before[4]
    @test PDMPSamplersRBridge._omrf_position_anchor_horizon!(
        manager, harmonic_state, harmonic_flow, :main) ≈ 2pi / 3 atol=2e-13

    tangent = StickyPDMPState(0.0,
        SkeletonPoint(zeros(9), [1.0; zeros(8)]), trues(9))
    manager.active_idx = 1
    PDMPSamplers.refresh_anchor_owned!(cv, manager.entries[1].state.anchor)
    @test isinf(PDMPSamplersRBridge._omrf_position_anchor_horizon!(
        manager, tangent, harmonic_flow, :main))
    frozen_harmonic = StickyPDMPState(0.0,
        SkeletonPoint(zeros(9), [0.0; zeros(8)]), trues(9))
    frozen_harmonic.free[1] = false
    frozen_harmonic.stored_velocity[1] = 3.0
    @test isinf(PDMPSamplersRBridge._omrf_position_anchor_horizon!(
        manager, frozen_harmonic, harmonic_flow, :main))
    reversed = deepcopy(tangent)
    reversed.ξ.θ .*= -1
    @test PDMPSamplersRBridge._omrf_position_anchor_index(
        manager, tangent.ξ.x) ==
        PDMPSamplersRBridge._omrf_position_anchor_index(
            manager, reversed.ξ.x)
    # Artificial grid subdivision does not enter the selection rule.
    subdivision_origin = StickyPDMPState(0.0,
        SkeletonPoint([0.2; zeros(8)], [0.7; zeros(8)]), trues(9))
    subdivided = deepcopy(subdivision_origin)
    for endpoint_time in (0.07, 0.19, 0.41, 0.83)
        direct = move_forward_time(subdivision_origin, endpoint_time,
            harmonic_flow)
        move_forward_time!(subdivided,
            endpoint_time - subdivided.t[], harmonic_flow)
        @test PDMPSamplersRBridge._omrf_position_anchor_index(manager,
            direct.ξ.x) == PDMPSamplersRBridge._omrf_position_anchor_index(
                manager, subdivided.ξ.x)
    end
    # Two competitors cross the current cell simultaneously and remain tied
    # on the outgoing diagonal; the established lowest-index convention owns
    # the boundary without physical displacement.
    original_third = manager.entries[3]
    manager.entries[3] = PDMPSamplersRBridge.OMRFAnchorEntry(
        PDMPSamplersRBridge._prepare_omrf_anchor(
            manager, [0.0, 2.0, zeros(7)...]), 0)
    simultaneous = StickyPDMPState(0.0,
        SkeletonPoint(zeros(9), [2.0, 2.0, zeros(7)...]), trues(9))
    manager.active_idx = 1
    @test PDMPSamplersRBridge._omrf_position_anchor_horizon!(manager,
        simultaneous, harmonic_flow, :main) ≈ pi / 6 atol=2e-13
    @test manager.pending_boundary_idx == 2
    manager.entries[3] = original_third
    linear_zz = StickyPDMPState(0.0,
        SkeletonPoint(zeros(9), [1.0; zeros(8)]), trues(9))
    manager.active_idx = 1
    @test PDMPSamplersRBridge._omrf_position_anchor_horizon!(
        manager, linear_zz, ZigZag(9), :main) == 1.0

    # Boomerang's Gaussian-reference correction is applied once to the
    # aggregate deterministic event gradient.  It is absent from the sampled
    # likelihood residual, whose value is dynamics-independent.
    correction_state = deepcopy(tangent)
    correction_state.ξ.x .= range(-0.3, 0.5; length=9)
    raw_deterministic = zeros(9)
    PDMPSamplers.deterministic_gradient!(raw_deterministic, cv,
        correction_state.ξ.x)
    gradient_cache = (∇ϕx=zeros(9), z=zeros(9), tmp=zeros(9))
    corrected_once = PDMPSamplers.compute_gradient!(correction_state, cv,
        harmonic_flow, gradient_cache)
    @test corrected_once ≈ raw_deterministic .-
        harmonic_flow.Γ * (correction_state.ξ.x .- harmonic_flow.μ)
    residual_bps = zeros(9); residual_boomerang = zeros(9)
    manager.oracle(residual_bps, correction_state.ξ.x, [2],
        manager.oracle.anchor)
    manager.oracle(residual_boomerang, correction_state.ξ.x, [2],
        manager.oracle.anchor)
    @test residual_boomerang == residual_bps

    # Production draw/refinement regression at a preserved complete live
    # candidate.  This invokes the specialized cell-marked draw_subset!, its
    # lazy selected-factor point bound, exact selected refinement, residual
    # oracle, and final thinning in their production order.
    function validate_live_context(cv, oracle, state, flow, D, M)
        env = cv.envelope
        scales = env.component_scales!
        cv.anchor == oracle.anchor == scales.base.anchor ||
            throw(ArgumentError("live anchors disagree"))
        counts = [count(==(component), vec(env.groups))
            for component in eachindex(env.totals)]
        counts == Int.(scales.component_totals) ||
            throw(ArgumentError("component membership counts disagree"))
        residual = zeros(length(state.ξ.x))
        oracle(residual, state.ξ.x, cv.subset, oracle.anchor)
        exact = 9 * PDMPSamplersRBridge._omrf_exact_selected_residual_abs_rate(
            state, residual, flow, oracle.residual_context.touched_positions)
        M + 1e-12 >= D + exact ||
            throw(ArgumentError("selected-factor bound is too small"))
        true
    end

    law_state = deepcopy(state)
    law_state.ξ.x[1] = 0.72
    law_state.ξ.θ .= range(-0.8, 0.9; length=9)
    PDMPSamplers.set_active_set!(cv, law_state.free)
    PDMPSamplers.component_cell_scales!(cv.envelope.cell_scales,
        cv.envelope, law_state, BouncyParticle(9, 0.0), 0.0, 0.2)
    B = PDMPSamplers.screening_residual_bound(
        cv.envelope, law_state, BouncyParticle(9, 0.0), 0.1)
    candidate = move_forward_time(law_state, 0.1, BouncyParticle(9, 0.0))
    deterministic = zeros(9)
    PDMPSamplers.deterministic_gradient!(deterministic, cv, candidate.ξ.x)
    D = PDMPSamplers.λ(candidate, deterministic, BouncyParticle(9, 0.0))
    exact = zeros(9); residual = zeros(9)
    expected = zeros(9)
    for factor in 1:9
        fill!(residual, 0.0)
        manager.oracle(residual, candidate.ξ.x, [factor], manager.oracle.anchor)
        @. exact = deterministic + 9residual
        expected[factor] = PDMPSamplers.λ(
            candidate, exact, BouncyParticle(9, 0.0)) / (9 * (D + B))
    end
    trials = 80_000
    accepted = zeros(Int, 9)
    early_rejections = 0
    residual_marks_before = cv.envelope.component_scales!.component_marks
    rng = Xoshiro(7401)
    for _ in 1:trials
        M = PDMPSamplers.draw_subset!(rng, cv, D, B)
        R = PDMPSamplers.subsampling_residual_subset_bound(
            manager.oracle, candidate, deterministic, BouncyParticle(9, 0.0),
            D, M, cv.subset, 9)
        if iszero(R) || (R < M && rand(rng) * M > R)
            early_rejections += 1
            continue
        end
        T = PDMPSamplers.subsampling_subset_bound(
            manager.oracle, candidate, deterministic, BouncyParticle(9, 0.0),
            D, R, cv.subset, 9)
        (iszero(T) || (T < R && rand(rng) * R > T)) && continue
        fill!(residual, 0.0)
        manager.oracle(residual, candidate.ξ.x, cv.subset,
            manager.oracle.anchor)
        @. exact = deterministic + 9residual
        A = PDMPSamplers.λ(candidate, exact, BouncyParticle(9, 0.0))
        rand(rng) * T <= A && (accepted[cv.subset[1]] += 1)
    end
    observed = accepted ./ trials
    standard_error = sqrt.(expected .* (1 .- expected) ./ trials)
    @test maximum(abs.(observed .- expected) ./
        max.(standard_error, 1 / trials)) < 5.0
    residual_draws = cv.envelope.component_scales!.component_marks -
        residual_marks_before
    @test 0 < residual_draws < trials
    @test early_rejections > 0
    println("live_selection_max_z=",
        maximum(abs.(observed .- expected) ./
            max.(standard_error, 1 / trials)),
        " residual_mixture_draws=", residual_draws,
        " early_refinement_rejections=", early_rejections)

    # The actual m=4 person-node design used by dependent Zig-Zag: enumerate
    # all 126 subsets independently, verify the estimator mean against the
    # native full gradient, and verify draw_subset!'s implemented size-biased
    # subset probabilities (including its deterministic/residual mixture).
    batch_ctx, batch_names =
        PDMPSamplersRBridge._new_stan_subsampling_context(
            lib, full, prior, 4)
    batch_spec = copy(spec)
    batch_spec["bound_type"] = "node_local_norm_low_rank"
    batch_weights = zeros(9, 3)
    for node in 1:3, person in 1:3
        radius = sqrt(sum(abs2, X[person, neighbour]
            for neighbour in 1:3 if neighbour != node))
        batch_weights[node, person] = 1.0
        batch_weights[3 + node, person] = radius
        batch_weights[6 + node, person] = radius * radius
    end
    batch_spec["norm_weights"] = batch_weights
    batch_sub = Dict{String,Any}("n_observations"=>9,
        "subsample_size"=>4, "anchor"=>zeros(9),
        "residual_envelope"=>batch_spec)
    batch_model, _, _, _, batch_manager =
        PDMPSamplersRBridge._build_stan_subsampling_model(
            batch_ctx, batch_names, batch_sub, zeros(9), nothing, falses(9);
            anchor_capacity=1)
    batch_cv = batch_model.grad
    zig_flow = ZigZag(9)
    batch_state = deepcopy(law_state)
    batch_state.ξ.θ .= [1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0]
    PDMPSamplers.set_active_set!(batch_cv, batch_state.free)
    # Use the production screening entry point: the lazy tree owns cached
    # state/flow used by subsequent selected-factor refinement.
    batch_B = PDMPSamplers.screening_residual_bound(
        batch_cv.envelope, batch_state, zig_flow, 0.0)
    batch_deterministic = zeros(9)
    PDMPSamplers.deterministic_gradient!(batch_deterministic, batch_cv,
        batch_state.ξ.x)
    batch_D = PDMPSamplers.λ(batch_state, batch_deterministic, zig_flow)
    subsets = NTuple{4,Int}[]
    for i in 1:6, j in (i + 1):7, k in (j + 1):8, l in (k + 1):9
        push!(subsets, (i, j, k, l))
    end
    @test length(subsets) == 126
    estimator_mean = zeros(9)
    enumerated_intensity = zeros(9)
    batch_residual = zeros(9)
    sampled_gradient = similar(batch_residual)
    for subset in subsets
        fill!(batch_residual, 0.0)
        batch_manager.oracle(batch_residual, batch_state.ξ.x,
            collect(subset), batch_manager.oracle.anchor)
        @. sampled_gradient = batch_deterministic + (9 / 4) * batch_residual
        estimator_mean .+= sampled_gradient
        enumerated_intensity .+= max.(0.0,
            batch_state.ξ.θ .* sampled_gradient)
    end
    estimator_mean ./= length(subsets)
    enumerated_intensity ./= length(subsets)
    native = PDMPSamplersRBridge.build_omrf_full_model(
        spec, batch_names, 9)
    native_gradient = zeros(9)
    native.grad.f(native_gradient, batch_state.ξ.x)
    @test estimator_mean ≈ native_gradient atol=2e-12
    # Exercise the actual lazy size-biased draw/refinement law.  Some residual
    # mixture paths are rejected within the component hierarchy before an
    # exact factor is evaluated, so the static generic subset-probability
    # formula is intentionally not used here.  The contribution of a sampled
    # candidate to coordinate i is (D+B) λ_i/M; its expectation must equal
    # the independently enumerated uniform-batch marked intensity.  This is
    # distinct from the full-gradient rate because positive-part formation is
    # nonlinear; unbiasedness of the underlying gradient is checked above.
    batch_trials = 80_000
    batch_rng = Xoshiro(7402)
    intensity_sum = zeros(9)
    intensity_square_sum = zeros(9)
    for _ in 1:batch_trials
        M = PDMPSamplers.draw_subset!(batch_rng, batch_cv, batch_D, batch_B)
        iszero(M) && continue
        R = PDMPSamplers.subsampling_residual_subset_bound(
            batch_manager.oracle, batch_state, batch_deterministic, zig_flow,
            batch_D, M, batch_cv.subset, 9 / 4)
        (iszero(R) || (R < M && rand(batch_rng) * M > R)) && continue
        T = PDMPSamplers.subsampling_subset_bound(
            batch_manager.oracle, batch_state, batch_deterministic, zig_flow,
            batch_D, R, batch_cv.subset, 9 / 4)
        (iszero(T) || (T < R && rand(batch_rng) * R > T)) && continue
        fill!(batch_residual, 0.0)
        batch_manager.oracle(batch_residual, batch_state.ξ.x,
            batch_cv.subset, batch_manager.oracle.anchor)
        @. sampled_gradient = batch_deterministic + (9 / 4) * batch_residual
        @inbounds for coordinate in eachindex(sampled_gradient)
            contribution = (batch_D + batch_B) * max(0.0,
                batch_state.ξ.θ[coordinate] * sampled_gradient[coordinate]) / T
            intensity_sum[coordinate] += contribution
            intensity_square_sum[coordinate] += contribution^2
        end
    end
    intensity_mean = intensity_sum ./ batch_trials
    intensity_variance = max.(0.0,
        intensity_square_sum ./ batch_trials .- intensity_mean .^ 2)
    intensity_se = sqrt.(intensity_variance ./ batch_trials)
    batch_z = maximum(abs.(intensity_mean .- enumerated_intensity) ./
        max.(intensity_se, 1e-12))
    println("batch4_intensity_mean=", intensity_mean,
        " enumerated=", enumerated_intensity, " se=", intensity_se)
    @test batch_z < 5.0
    println("batch4_accepted_intensity_max_z=", batch_z)

    # The test-only live-context oracle detects each requested deliberate
    # perturbation before it can be mistaken for a passing algebraic identity.
    M = PDMPSamplers.draw_subset!(rng, cv, D, B)
    R = PDMPSamplers.subsampling_residual_subset_bound(
        manager.oracle, candidate, deterministic, BouncyParticle(9, 0.0),
        D, M, cv.subset, 9)
    @test validate_live_context(cv, manager.oracle, candidate,
        BouncyParticle(9, 0.0), D, R)
    totals = copy(cv.envelope.component_scales!.component_totals)
    cv.envelope.component_scales!.component_totals[1] += 1
    @test_throws ArgumentError validate_live_context(cv, manager.oracle,
        candidate, BouncyParticle(9, 0.0), D, R)
    cv.envelope.component_scales!.component_totals .= totals
    @test_throws ArgumentError validate_live_context(cv, manager.oracle,
        candidate, BouncyParticle(9, 0.0), D, D)
    anchor0 = copy(cv.envelope.component_scales!.base.anchor)
    cv.envelope.component_scales!.base.anchor[1] += 0.1
    @test_throws ArgumentError validate_live_context(cv, manager.oracle,
        candidate, BouncyParticle(9, 0.0), D, R)
    cv.envelope.component_scales!.base.anchor .= anchor0

    # A prediction is not permission to commit.  The exact reviewer fixture
    # (anchors 0/2, x=0, v=+1) is interrupted at t=.25 by each relevant
    # physical or ordinary boundary class.  Each interruption invalidates the
    # target and selection is recomputed from the actual outgoing state.
    resize!(manager.entries, 2)
    manager.entries[1].state = PDMPSamplersRBridge._prepare_omrf_anchor(
        manager, [0.0; zeros(8)])
    manager.entries[2].state = PDMPSamplersRBridge._prepare_omrf_anchor(
        manager, [2.0; zeros(8)])
    for interrupted_by in (:reflect, :refresh, :sticky_freeze,
            :sticky_release, :horizon_hit, :terminal_checkpoint)
        manager.active_idx = 1
        manager.distance_scales .= 1.0
        interrupted = StickyPDMPState(0.0,
            SkeletonPoint(zeros(9), [1.0; zeros(8)]), trues(9))
        @test PDMPSamplersRBridge._omrf_position_anchor_horizon!(manager,
            interrupted, BouncyParticle(9, 0.0), :main) == 1.0
        move_forward_time!(interrupted, 0.25, BouncyParticle(9, 0.0))
        interrupted.ξ.θ[1] = -1.0
        @test !PDMPSamplersRBridge._select_omrf_anchor!(manager, cv,
            interrupted.ξ.x; phase=:main, state=interrupted,
            flow=BouncyParticle(9, 0.0), boundary_won=false,
            event_type=interrupted_by)
        @test manager.active_idx == 1
        @test manager.pending_boundary_idx == 0
        @test isinf(PDMPSamplersRBridge._omrf_position_anchor_horizon!(manager,
            interrupted, BouncyParticle(9, 0.0), :main))
    end

    # Adaptive dispatch must never enter the position branch.  A deliberately
    # staged position target is therefore left untouched.
    manager.selection_policy = :adaptive
    manager.main_refresh_distance = 1.0
    manager.pending_boundary_idx = 2
    manager.pending_boundary_time = 1.0
    generation = manager.boundary_generation
    PDMPSamplersRBridge._select_omrf_anchor!(manager, cv, zeros(9);
        phase=:main, state, flow=BouncyParticle(9, 0.0),
        boundary_won=true, event_type=:anchor_selection_boundary)
    @test manager.pending_boundary_idx == 2
    @test manager.boundary_generation == generation
    @test manager.active_idx == 1
    manager.selection_policy = :position

    # At the exact boundary, the outgoing anchor is installed at zero elapsed
    # time.  Recomputing from that outgoing cell cannot zero-loop.
    manager.active_idx = 2
    PDMPSamplers.refresh_anchor_owned!(cv, manager.entries[2].state.anchor)
    @test PDMPSamplersRBridge._omrf_position_anchor_horizon!(
        manager, state, BouncyParticle(9, 0.0), :main) == Inf
    tied = copy(state.ξ.x)
    @test PDMPSamplersRBridge._omrf_position_anchor_index(manager, tied) == 1
    coincident = PDMPSamplersRBridge._prepare_omrf_anchor(manager,
        manager.entries[2].state.anchor)
    push!(manager.entries, PDMPSamplersRBridge.OMRFAnchorEntry(coincident, 0))
    @test isinf(PDMPSamplersRBridge._omrf_position_anchor_horizon!(
        manager, state, BouncyParticle(9, 0.0), :main))

    # Real sampler-entry scheduler regression.  Closely spaced fixed anchors
    # force several exact computational boundaries before ordinary BPS events;
    # the run must progress without retained preparation or a zero-time loop.
    path_model, _, _, path_adapter, path_manager = built(5)
    path_manager.selection_policy = :position
    path_manager.main_refresh_distance = 0.03 # ignored by the fixed policy
    path_manager.distance_scales .= 1.0
    empty!(path_manager.entries)
    for coordinate in (-0.250, -0.246, -0.242, -0.238, -0.234)
        anchor = [coordinate; zeros(8)]
        PDMPSamplersRBridge._insert_omrf_anchor!(path_manager,
            PDMPSamplersRBridge._prepare_omrf_anchor(path_manager, anchor))
    end
    path_manager.active_idx = 1
    PDMPSamplers.refresh_anchor_owned!(path_model.grad,
        path_manager.entries[1].state.anchor)
    path_adapter.update_dt = Inf
    prepared_before_path = path_manager.preparations
    scheduler_records = NamedTuple[]
    PDMPSamplers.set_phase_step_observer!(
        (stage, phase, observed_state, criterion_horizon, adapter_horizon,
                event_type, observed_stats) -> push!(scheduler_records, (
            stage, phase, physical_time=Float64(observed_state.t[]),
            returned_horizon=Float64(adapter_horizon),
            active_anchor=path_manager.active_idx,
            pending_target=path_manager.pending_boundary_idx,
            event_type,
            rejection=Bool(observed_stats.last_rejected))))
    path_result = try
        pdmp_sample(
            SkeletonPoint([-0.251; zeros(8)], [2.0; fill(0.05, 8)]),
            BouncyParticle(9, 0.5), path_model,
            GridThinningStrategy(; N=4, t_max=0.01, lazy=false,
                bound_violation=:throw), 0.0, 0.04, 0.0;
            adapter=path_adapter, progress=false, seed=7511,
            statistic_counter=PDMPSamplers.DevelStatisticCounter)
    finally
        PDMPSamplers.set_phase_step_observer!(nothing)
    end
    path_trace = path_result.traces[1]
    @test path_manager.selection_boundaries >= 2
    @test path_manager.preparations == prepared_before_path
    @test path_manager.active_idx ==
        PDMPSamplersRBridge._omrf_position_anchor_index(path_manager,
            path_result.terminal_states[1].position)
    @test last(path_trace.times) > first(path_trace.times)
    @test all(diff(path_trace.times) .>= 0)
    @test count(iszero, diff(path_trace.times)) <= 1
    @test path_result.stats[1].sticky_events == 0
    boundary_records = filter(r -> r.stage === :after_adapt &&
        r.event_type === :anchor_selection_boundary, scheduler_records)
    @test length(boundary_records) == path_manager.selection_boundaries
    @test all(!r.rejection for r in boundary_records)
    @test all(r.pending_target == 0 for r in boundary_records)
    @test all(diff(getproperty.(boundary_records, :physical_time)) .> 0)

    # Real harmonic sampler-entry scheduling: a tightly spaced prepared bank
    # forces genuine source changes while ordinary proposal horizons can still
    # interrupt and invalidate predictions.
    boom_model, _, _, boom_adapter, boom_manager = built(5)
    boom_manager.selection_policy = :position
    boom_manager.distance_scales .= 1.0
    empty!(boom_manager.entries)
    for coordinate in (-0.010, -0.005, 0.0, 0.005, 0.010)
        anchor = [coordinate; zeros(8)]
        PDMPSamplersRBridge._insert_omrf_anchor!(boom_manager,
            PDMPSamplersRBridge._prepare_omrf_anchor(boom_manager, anchor))
    end
    boom_manager.active_idx = 1
    PDMPSamplers.refresh_anchor_owned!(boom_model.grad,
        boom_manager.entries[1].state.anchor)
    boom_adapter.update_dt = Inf
    boom_preparations = boom_manager.preparations
    boom_result = pdmp_sample(
        SkeletonPoint([-0.012; zeros(8)], [0.3; fill(0.01, 8)]),
        Boomerang(Diagonal(ones(9)), zeros(9), 0.0), boom_model,
        GridThinningStrategy(; N=4, t_max=0.004, lazy=false,
            bound_violation=:throw), 0.0, 0.08, 0.0;
        adapter=boom_adapter, progress=false, seed=7513,
        statistic_counter=PDMPSamplers.DevelStatisticCounter)
    @test boom_manager.selection_boundaries >= 2
    @test boom_manager.preparations == boom_preparations
    @test boom_result.terminal_states[1].t > 0.0
    @test boom_manager.active_idx ==
        PDMPSamplersRBridge._omrf_position_anchor_index(boom_manager,
            boom_result.terminal_states[1].position)

    # Logical checkpoint restoration: rebuild prepared quantities from values,
    # never from the restart endpoint or serialized callbacks/caches.
    restore_model, _, restore_names, restore_adapter, restore_manager = built(5)
    restore_manager.selection_policy = :position
    bank_spec = Dict{String,Any}(
        "version" => 1,
        "selection_policy" => "position",
        "tie_policy" => "lowest_index_outgoing",
        "parameter_order" => restore_names,
        "ordered_anchor_positions" => reduce(hcat,
            (entry.state.anchor for entry in path_manager.entries)),
        "distance_scales" => copy(path_manager.distance_scales))
    endpoint = path_result.terminal_states[1]
    @test PDMPSamplersRBridge._restore_omrf_position_bank!(
        restore_manager, restore_model.grad, bank_spec, endpoint.position)
    @test restore_manager.distance_scales == path_manager.distance_scales
    @test [entry.state.anchor for entry in restore_manager.entries] ==
        [entry.state.anchor for entry in path_manager.entries]
    @test restore_manager.active_idx ==
        PDMPSamplersRBridge._omrf_position_anchor_index(restore_manager,
            endpoint.position)
    restored_state = PDMPState(endpoint.t,
        SkeletonPoint(copy(endpoint.position), copy(endpoint.physical_velocity)))
    reference_cv = path_model.grad
    restored_cv = restore_model.grad
    PDMPSamplersRBridge._install_omrf_anchor!(path_manager,
        path_manager.entries[path_manager.active_idx].state)
    PDMPSamplers.refresh_anchor_owned!(reference_cv,
        path_manager.entries[path_manager.active_idx].state.anchor)
    g_reference = zeros(9); g_restored = zeros(9)
    PDMPSamplers.deterministic_gradient!(g_reference, reference_cv,
        endpoint.position)
    PDMPSamplers.deterministic_gradient!(g_restored, restored_cv,
        endpoint.position)
    @test g_restored == g_reference
    reference_residual = zeros(9); restored_residual = zeros(9)
    reference_cv.residual_oracle(reference_residual, endpoint.position, [1],
        reference_cv.anchor)
    restored_cv.residual_oracle(restored_residual, endpoint.position, [1],
        restored_cv.anchor)
    @test restored_residual == reference_residual
    reference_scales = similar(reference_cv.envelope.totals)
    restored_scales = similar(restored_cv.envelope.totals)
    PDMPSamplers.component_cell_scales!(reference_scales,
        reference_cv.envelope, restored_state, endpoint.flow, 0.0, 0.005)
    PDMPSamplers.component_cell_scales!(restored_scales,
        restored_cv.envelope, restored_state, endpoint.flow, 0.0, 0.005)
    @test restored_scales == reference_scales
    continuation_point = SkeletonPoint(copy(endpoint.position),
        copy(endpoint.physical_velocity))
    continuation_algorithm() = GridThinningStrategy(; N=4, t_max=0.01,
        lazy=false, bound_violation=:throw)
    reference_next = pdmp_sample(continuation_point, deepcopy(endpoint.flow),
        path_model, continuation_algorithm(), endpoint.t, endpoint.t + 0.005,
        0.0; adapter=path_adapter, progress=false, seed=7512,
        statistic_counter=PDMPSamplers.DevelStatisticCounter).traces[1]
    restored_next = pdmp_sample(continuation_point, deepcopy(endpoint.flow),
        restore_model, continuation_algorithm(), endpoint.t,
        endpoint.t + 0.005, 0.0; adapter=restore_adapter, progress=false,
        seed=7512, statistic_counter=PDMPSamplers.DevelStatisticCounter).traces[1]
    @test restored_next.times == reference_next.times
    @test restored_next.positions == reference_next.positions
    @test restored_next.velocities == reference_next.velocities

end

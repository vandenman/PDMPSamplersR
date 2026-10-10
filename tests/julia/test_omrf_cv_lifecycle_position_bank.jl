using Test, Random, LinearAlgebra, PDMPSamplers, BridgeStan
include(normpath(joinpath(@__DIR__, "..", "..", "inst", "julia",
    "main_interface_function.jl")))
using .PDMPSamplersRBridge

function direct_spec(X, seen)
    Dict{String,Any}(
        "type"=>"omrf", "backend"=>"analytic", "X"=>X, "seen"=>seen,
        "thresholds"=>"thresholds_0", "interactions"=>"interactions_0",
        "bound_type"=>"direct_factor_roof_person_hcv",
        "factorization"=>"person", "factor_sampling"=>"size_biased",
        "prior_backend"=>"analytic", "use_hcv"=>true, "hcv_damping"=>Inf,
        "hcv_scope"=>"global",
        "analytic_prior"=>Dict{String,Any}("type"=>"gaussian",
            "threshold_alpha"=>1.0, "threshold_beta"=>1.0,
            "interaction_scale"=>1.0),
        "exact_selected_refinement"=>true,
        "categorical_residual_cap"=>true)
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
    spec = direct_spec(X, [3, 3, 3])
    N = size(X, 1)
    sub = Dict{String,Any}("n_observations"=>N, "subsample_size"=>1,
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

    function entry_case(; warmup=0.04, updates=true)
        model, _, names, adapter, manager = built(1)
        adapter.update_dt = updates ? 0.01 : Inf
        flow = BouncyParticle(9, 0.5)
        alg = GridThinningStrategy(; N=4, t_max=0.02, lazy=false,
            bound_violation=:throw)
        observed = Ref(false)
        observer = PDMPSamplers.SamplerObserver(; phase_start=
            (phase, state, live_model, flow, alg, cache) -> begin
                phase === :main || return
                cv = live_cv(live_model.grad)
                @test cv !== nothing
                @test cv.anchor == manager.oracle.anchor
                @test PDMPSamplers.trajectory_scale_anchor(cv.envelope) ==
                    manager.oracle.anchor
                @test cv.residual_oracle.f.anchor == manager.oracle.anchor
                observed[] = true
            end)
        pdmp_sample(SkeletonPoint(fill(0.15, 9), fill(0.4, 9)), flow,
            model, alg, 0.0, warmup + 0.02, warmup;
            adapter, progress=false, seed=7301,
            statistic_counter=PDMPSamplers.DevelStatisticCounter, observer)
        @test observed[]
        @test manager.oracle.anchor != zeros(9)
    end

    # Actual sampler-entry regressions: subsampled warmup and a zero-warmup
    # restart.  Capacity one forces bank replacement without an index change.
    entry_case(warmup=0.04, updates=true)
    entry_case(warmup=0.0, updates=false)

    # Position-bank finalization is explicit and independent of the adaptive
    # route's historical `main_refresh_distance == Inf` sentinel.
    fixed_model, _, fixed_names, fixed_adapter, fixed_manager = built(4)
    fixed_manager.selection_policy = :position
    fixed_manager.main_refresh_distance = 0.2
    fixed_adapter.update_dt = Inf
    pdmp_sample(SkeletonPoint(fill(0.15, 9), fill(0.4, 9)),
        BouncyParticle(9, 0.5), fixed_model,
        GridThinningStrategy(; N=4, t_max=0.02, lazy=false,
            bound_violation=:throw), 0.0, 0.06, 0.04;
        adapter=fixed_adapter, progress=false, seed=7302,
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
    @test PDMPSamplers.trajectory_scale_anchor(wrapped.grad.envelope) != forced
    @test PDMPSamplersRBridge._finish_omrf_warmup!(
        manager, wrapped.grad, forced)
    @test wrapped.grad.anchor == forced
    @test PDMPSamplers.trajectory_scale_anchor(wrapped.grad.envelope) == forced

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
    # the pre-proposal hook and must attach the live context before event
    # search.
    preproposal = PDMPSamplers.with_stats(
        model, PDMPSamplers.DevelStatisticCounter())
    PDMPSamplers.adaptation_horizon(adapter, state,
        BouncyParticle(9, 0.0), preproposal.grad, :main, stats)
    @test preproposal.grad.envelope.context === manager.oracle.residual_context
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
        PDMPSamplers.trajectory_scale_anchor(cv.envelope)

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
    # candidate.  This invokes the direct envelope's draw_subset!, exact
    # selected refinement, residual oracle, and final thinning in their
    # production order: a person J must be accepted with probability
    # λ(D + N r_J) / (N (D + B)).
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
    expected = zeros(N)
    for person in 1:N
        fill!(residual, 0.0)
        manager.oracle(residual, candidate.ξ.x, [person], manager.oracle.anchor)
        @. exact = deterministic + N * residual
        expected[person] = PDMPSamplers.λ(
            candidate, exact, BouncyParticle(9, 0.0)) / (N * (D + B))
    end
    trials = 80_000
    accepted = zeros(Int, N)
    rng = Xoshiro(7401)
    for _ in 1:trials
        M = PDMPSamplers.draw_subset!(rng, cv, D, B)
        R = PDMPSamplers.subsampling_residual_subset_bound(
            manager.oracle, candidate, deterministic, BouncyParticle(9, 0.0),
            D, M, cv.subset, N)
        (iszero(R) || (R < M && rand(rng) * M > R)) && continue
        T = PDMPSamplers.subsampling_subset_bound(
            manager.oracle, candidate, deterministic, BouncyParticle(9, 0.0),
            D, R, cv.subset, N)
        (iszero(T) || (T < R && rand(rng) * R > T)) && continue
        fill!(residual, 0.0)
        manager.oracle(residual, candidate.ξ.x, cv.subset,
            manager.oracle.anchor)
        @. exact = deterministic + N * residual
        A = PDMPSamplers.λ(candidate, exact, BouncyParticle(9, 0.0))
        rand(rng) * T <= A && (accepted[cv.subset[1]] += 1)
    end
    observed = accepted ./ trials
    standard_error = sqrt.(expected .* (1 .- expected) ./ trials)
    @test maximum(abs.(observed .- expected) ./
        max.(standard_error, 1 / trials)) < 5.0
    println("live_selection_max_z=",
        maximum(abs.(observed .- expected) ./
            max.(standard_error, 1 / trials)))

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
    scheduler_observer = PDMPSamplers.SamplerObserver(; step=
        (stage, phase, observed_state, criterion_horizon, adapter_horizon,
                event_type, observed_stats) -> push!(scheduler_records, (
            stage, phase, physical_time=Float64(observed_state.t[]),
            returned_horizon=Float64(adapter_horizon),
            active_anchor=path_manager.active_idx,
            pending_target=path_manager.pending_boundary_idx,
            event_type,
            rejection=Bool(observed_stats.last_rejected))))
    path_result = pdmp_sample(
        SkeletonPoint([-0.251; zeros(8)], [2.0; fill(0.05, 8)]),
        BouncyParticle(9, 0.5), path_model,
        GridThinningStrategy(; N=4, t_max=0.01, lazy=false,
            bound_violation=:throw), 0.0, 0.04, 0.0;
        adapter=path_adapter, progress=false, seed=7511,
        statistic_counter=PDMPSamplers.DevelStatisticCounter,
        observer=scheduler_observer)
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

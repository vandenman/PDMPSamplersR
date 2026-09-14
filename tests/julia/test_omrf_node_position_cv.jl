using Test, LinearAlgebra, Random, PDMPSamplers, BridgeStan
include(normpath(joinpath(@__DIR__, "..", "..", "inst", "julia",
    "main_interface_function.jl")))
using .PDMPSamplersRBridge

function node_position_spec(X, seen; prior_type="gaussian")
    N, P = size(X)
    membership = zeros(Int, P, N); nodes = Int[]
    minima = Vector{Vector{Float64}}(); maxima = similar(minima)
    offset = 0
    for node in 1:P
        local_groups = [1 + mod(i - 1, 2) for i in 1:N]
        for group in 1:2
            people = findall(==(group), local_groups)
            lo = zeros(P); hi = zeros(P)
            for neighbour in 1:P
                neighbour == node && continue
                lo[neighbour] = minimum(X[people, neighbour])
                hi[neighbour] = maximum(X[people, neighbour])
            end
            push!(nodes, node); push!(minima, lo); push!(maxima, hi)
        end
        membership[node, :] .= offset .+ local_groups
        offset += 2
    end
    prior = prior_type == "gaussian" ? Dict{String,Any}(
        "type"=>"gaussian", "threshold_alpha"=>1.0,
        "threshold_beta"=>1.0, "interaction_scale"=>1.0) :
        Dict{String,Any}("type"=>"node_scale_gaussian",
            "threshold_alpha"=>1.0, "threshold_beta"=>1.0,
            "interaction_scale"=>0.7, "global_scale"=>"global_scale",
            "node_scale"=>"node_scale", "global_scale_sd"=>1.2,
            "node_scale_sd"=>0.8)
    Dict{String,Any}(
        "type"=>"omrf", "backend"=>"analytic", "X"=>X, "seen"=>seen,
        "thresholds"=>"thresholds_0", "interactions"=>"interactions_0",
        "bound_type"=>"cell_marked_node_box_first_order",
        "factorization"=>"person_node", "factor_sampling"=>"size_biased",
        "prior_backend"=>"analytic", "use_hcv"=>false,
        "analytic_prior"=>prior, "exact_selected_refinement"=>true,
        "node_box_membership"=>membership, "node_box_nodes"=>nodes,
        "node_box_minima"=>reduce(vcat, permutedims.(minima)),
        "node_box_maxima"=>reduce(vcat, permutedims.(maxima)))
end

function explicit_gram(context, node, metric_variances)
    N = size(context.X, 1); q = context.seen[node] - 1
    L = length(context.incident_edges[node]); s = q + L
    Q = zeros(s, s)
    for person in 1:N
        Z = zeros(q, s)
        for u in 1:q
            Z[u, u] = 1.0
            for edge in 1:L
                Z[u, q + edge] = u * context.X[person,
                    context.incident_neighbours[node][edge]]
            end
        end
        support = vcat(context.threshold_indices[
            context.threshold_starts[node]:(context.threshold_starts[node]+q-1)],
            [context.interaction_indices[e] for e in context.incident_edges[node]])
        C = Diagonal(metric_variances[support])
        weight = tr(Z * C * Z')
        Q .+= weight .* (Z' * Z) ./ N
    end
    Q
end

@testset "research node-position first-order CV" begin
    root = normpath(joinpath(@__DIR__, "..", ".."))
    lib = PDMPSamplersRBridge._compile_model_with_header(
        joinpath(root, "inst", "stan", "omrf", "omrf.stan"),
        joinpath(root, "inst", "stan", "pdmp_subsample.hpp"))
    full = joinpath(root, "tests", "stan", "omrf_subsampling_full.json")
    prior = joinpath(root, "tests", "stan", "omrf_subsampling_prior.json")
    ctx, names = PDMPSamplersRBridge._new_stan_subsampling_context(
        lib, full, prior, 1)
    X = [0 1 2; 2 0 1; 1 2 0]; seen = [3, 3, 3]
    spec = node_position_spec(X, seen)
    d = 9; K = 8
    candidates = reshape(collect(range(-0.35, 0.42; length=d*K)), d, K)
    metric_scale = exp.(range(-0.4, 0.5; length=d))
    flow = PDMPSamplers.PreconditionedBPS(d;
        refresh_rate=0.7, scale=metric_scale)
    counts = PDMPSamplersRBridge.StanSubsamplingCallCounts()
    nodecv = PDMPSamplersRBridge._prepare_omrf_node_position_cv(
        spec, names, candidates, flow, counts)
    @test all(bank -> !bank.use_hcv, nodecv.bank.contexts)
    @test all(bank -> bank.likelihood_hessian === nothing, nodecv.bank.contexts)
    for node in 1:3
        @test nodecv.grams[node] ≈ explicit_gram(
            nodecv.bank.contexts[1], node, metric_scale .^ 2) atol=2e-12
    end

    subsampling = Dict{String,Any}("n_observations"=>9,
        "subsample_size"=>1, "residual_envelope"=>spec)
    model, _, _, adapter, manager =
        PDMPSamplersRBridge._build_stan_subsampling_model(
            ctx, names, subsampling, zeros(d), nothing, falses(d);
            anchor_capacity=K)
    manager.selection_policy = :node_position
    empty!(manager.entries); manager.active_idx = 0
    for source in 1:K
        prepared = PDMPSamplersRBridge._prepare_omrf_anchor(
            manager, view(candidates, :, source))
        manager.active_idx = PDMPSamplersRBridge._insert_omrf_anchor!(
            manager, prepared)
    end
    manager.oracle.node_position_cv = nodecv
    manager.node_position_policy_version = 1
    x = collect(range(-0.25, 0.31; length=d))
    PDMPSamplersRBridge._omrf_node_position_sources!(
        nodecv.requested_sources, nodecv, x)
    PDMPSamplersRBridge._install_omrf_node_position_sources!(
        manager, model.grad, nodecv.requested_sources)
    @test length(unique(nodecv.sources)) > 1

    # Exhaustive factor enumeration is an implementation-independent CV
    # identity: the mixed baseline plus all mixed-source residuals equals the
    # same native complete likelihood gradient.
    deterministic = zeros(d); residual = zeros(d); mixed_sum = zeros(d)
    model.grad.deterministic_gradient!(deterministic, x)
    for factor in 1:9
        model.grad.residual_oracle(residual, x, [factor], model.grad.anchor)
        mixed_sum .+= residual
    end
    native_context = PDMPSamplersRBridge._omrf_residual_context(
        spec, names, zeros(d), PDMPSamplersRBridge.StanSubsamplingCallCounts())
    native_prior = PDMPSamplersRBridge._omrf_analytic_prior(spec, native_context)
    native = PDMPSamplersRBridge.OMRFFullGradientProvider(
        native_context, native_prior)
    exact = zeros(d); native(exact, x)
    @test deterministic + mixed_sum ≈ exact atol=2e-9
    for mask_code in 0:3
        x_mask = copy(x)
        mask_code & 1 == 1 && (x_mask[7] = 0.0)
        mask_code & 2 == 2 && (x_mask[8] = 0.0)
        PDMPSamplersRBridge._omrf_node_position_sources!(
            nodecv.requested_sources, nodecv, x_mask)
        PDMPSamplersRBridge._install_omrf_node_position_sources!(
            manager, model.grad, nodecv.requested_sources)
        model.grad.deterministic_gradient!(deterministic, x_mask)
        fill!(mixed_sum, 0.0)
        for factor in 1:9
            model.grad.residual_oracle(
                residual, x_mask, [factor], model.grad.anchor)
            mixed_sum .+= residual
        end
        native(exact, x_mask)
        @test deterministic + mixed_sum ≈ exact atol=2e-9
    end

    # A common source exactly reproduces the existing global-anchor CV.
    fill!(nodecv.requested_sources, 3)
    PDMPSamplersRBridge._install_omrf_node_position_sources!(
        manager, model.grad, nodecv.requested_sources)
    global_oracle = PDMPSamplersRBridge.OMRFAnalyticSubsamplingOracle(
        ctx, spec, names, collect(view(candidates, :, 3)))
    global_d = zeros(d); global_r = zeros(d)
    global_oracle(global_d, x)
    model.grad.deterministic_gradient!(deterministic, x)
    @test deterministic ≈ global_d atol=2e-9
    for factor in (1, 5, 9)
        model.grad.residual_oracle(residual, x, [factor], model.grad.anchor)
        global_oracle(global_r, x, [factor], global_oracle.anchor)
        @test residual ≈ global_r atol=2e-9
    end

    # Mixed references include one shared interaction with different sources.
    nodecv.requested_sources .= [1, 8, 4]
    PDMPSamplersRBridge._install_omrf_node_position_sources!(
        manager, model.grad, nodecv.requested_sources)
    scales = model.grad.envelope.component_scales!
    @test scales.base.node_sources == [1, 8, 4]
    scales.base.categorical_residual_cap = true
    state = StickyPDMPState(0.0,
        SkeletonPoint(copy(x), collect(range(0.7, -0.6; length=d))),
        trues(d))
    PDMPSamplers.component_cell_scales!(model.grad.envelope.cell_scales,
        model.grad.envelope, state, flow, 0.0, 0.4)
    for factor in 1:9
        component = model.grad.envelope.groups[1, factor]
        for t in (0.0, 0.2, 0.4)
            PDMPSamplers.screening_residual_bound(
                model.grad.envelope, state, flow, t)
            point = PDMPSamplersRBridge._omrf_cell_marked_factor_point_bound(
                scales, factor)
            candidate = move_forward_time(state, t, flow)
            model.grad.residual_oracle(
                residual, candidate.ξ.x, [factor], model.grad.anchor)
            exact_rate = PDMPSamplersRBridge._omrf_exact_selected_residual_abs_rate(
                candidate, residual, flow,
                manager.oracle.residual_context.touched_positions)
            @test exact_rate <= point * (1 + 2e-7) + 2e-8
            @test point <= model.grad.envelope.scales[component] *
                (1 + 64eps(Float64))
        end
    end
    scales.base.categorical_residual_cap = false

    # Crossing ownership: an interrupted prediction is not committed; exact
    # and simultaneous boundaries select the deterministic outgoing sources.
    PDMPSamplersRBridge._omrf_node_position_sources!(
        nodecv.requested_sources, nodecv, state.ξ.x)
    PDMPSamplersRBridge._install_omrf_node_position_sources!(
        manager, model.grad, nodecv.requested_sources)
    horizon = PDMPSamplersRBridge._omrf_node_position_anchor_horizon!(
        manager, state, flow, :main)
    @test horizon >= 0
    saved_sources = copy(nodecv.sources)
    early = move_forward_time(state, min(horizon / 2, 0.01), flow)
    adapter.select_fn!(model.grad, early, flow, :main,
        PDMPSamplers.DevelStatisticCounter(), :reflection, false)
    @test nodecv.sources == PDMPSamplersRBridge._omrf_node_position_sources!(
        similar(nodecv.sources), nodecv, early.ξ.x)
    horizon2 = PDMPSamplersRBridge._omrf_node_position_anchor_horizon!(
        manager, early, flow, :main)
    if isfinite(horizon2)
        boundary = move_forward_time(early, horizon2, flow)
        boundary.t[] = early.t[] + horizon2
        position_before = copy(boundary.ξ.x); velocity_before = copy(boundary.ξ.θ)
        adapter.select_fn!(model.grad, boundary, flow, :main,
            PDMPSamplers.DevelStatisticCounter(),
            :anchor_selection_boundary, true)
        @test boundary.ξ.x == position_before
        @test boundary.ξ.θ == velocity_before
    end

    # Grid subdivision is computational only: selections along a prescribed
    # path depend on position, not how often the rule is queried.
    path_times = range(0.0, 1.0; length=41)
    direct_sources = [PDMPSamplersRBridge._omrf_node_position_sources!(
        similar(nodecv.sources), nodecv,
        state.ξ.x .+ t .* state.ξ.θ) for t in path_times]
    subdivided_sources = Vector{Vector{Int}}()
    for (left, right) in zip(path_times[1:end-1], path_times[2:end])
        for t in range(left, right; length=7)[1:end-1]
            PDMPSamplersRBridge._omrf_node_position_sources!(
                nodecv.requested_sources, nodecv,
                state.ξ.x .+ t .* state.ξ.θ)
        end
        push!(subdivided_sources,
            PDMPSamplersRBridge._omrf_node_position_sources!(
                similar(nodecv.sources), nodecv,
                state.ξ.x .+ right .* state.ξ.θ))
    end
    @test subdivided_sources == direct_sources[2:end]

    # Node prior adds scale coordinates to the target, never to selector
    # supports or likelihood residuals.
    node_names = vcat(names, ["global_scale"],
        ["node_scale.$j" for j in 1:3])
    node_spec = node_position_spec(X, seen; prior_type="node")
    node_candidates = vcat(candidates, zeros(4, K))
    node_bank = PDMPSamplersRBridge._prepare_omrf_node_position_cv(
        node_spec, node_names, node_candidates, ones(13))
    node_residual = zeros(13)
    PDMPSamplersRBridge._omrf_node_local_residual!(node_bank.bank,
        node_residual, vcat(x, [0.1, -0.2, 0.3, 0.4]), [1, 5, 9],
        [1, 8, 4]; factorization=:person_node)
    @test node_residual[10:13] == zeros(4)
    @test all(s -> all(<(10), s), node_bank.supports)

    for mask_code in 0:3
        x_node = vcat(x, [0.1, -0.2, 0.3, 0.4])
        mask_code & 1 == 1 && (x_node[7] = 0.0)
        mask_code & 2 == 2 && (x_node[8] = 0.0)
        sources = [1, 8, 4]
        baseline = zeros(13); residual_total = zeros(13)
        PDMPSamplersRBridge._omrf_node_local_baseline!(
            baseline, node_bank.bank, sources)
        PDMPSamplersRBridge._omrf_node_local_residual!(node_bank.bank,
            residual_total, x_node, collect(1:9), sources;
            factorization=:person_node)
        node_context = node_bank.bank.contexts[1]
        node_prior = PDMPSamplersRBridge._omrf_analytic_prior(
            node_spec, node_context)
        prior_gradient = zeros(13)
        PDMPSamplersRBridge._omrf_prior_gradient!(
            prior_gradient, node_prior, x_node)
        exact_node = zeros(13)
        PDMPSamplersRBridge.OMRFFullGradientProvider(
            node_context, node_prior)(exact_node, x_node)
        @test prior_gradient + baseline + residual_total ≈ exact_node atol=2e-9
        @test residual_total[10:13] == zeros(4)
    end

    # Logical restart rebuilds prepared quantities from the mathematical bank
    # and restores the selector before the first proposal.
    restart_ctx, restart_names =
        PDMPSamplersRBridge._new_stan_subsampling_context(
            lib, full, prior, 1)
    restart_model, _, _, _, restart_manager =
        PDMPSamplersRBridge._build_stan_subsampling_model(
            restart_ctx, restart_names, subsampling, zeros(d), nothing,
            falses(d); anchor_capacity=K)
    restart_manager.selection_policy = :node_position
    restart_spec = Dict{String,Any}(
        "version"=>1, "selection_policy"=>"node_position",
        "tie_policy"=>"lowest_index_outgoing",
        "parameter_order"=>names,
        "ordered_anchor_positions"=>candidates,
        "distance_scales"=>ones(d),
        "physical_velocity_variances"=>metric_scale .^ 2,
        "node_support_order"=>nodecv.supports, "policy_version"=>1)
    PDMPSamplersRBridge._restore_omrf_position_bank!(restart_manager,
        restart_model.grad, restart_spec, state.ξ.x)
    restored = restart_manager.oracle.node_position_cv
    @test restored.sources == PDMPSamplersRBridge._omrf_node_position_sources!(
        similar(restored.sources), restored, state.ξ.x)
    @test restored.supports == nodecv.supports
    @test all(isapprox(restored.grams[j], nodecv.grams[j]; atol=2e-12)
        for j in 1:3)
    PDMPSamplersRBridge._omrf_node_position_sources!(
        nodecv.requested_sources, nodecv, state.ξ.x)
    PDMPSamplersRBridge._install_omrf_node_position_sources!(
        manager, model.grad, nodecv.requested_sources)
    restart_model.grad.deterministic_gradient!(global_d, state.ξ.x)
    model.grad.deterministic_gradient!(deterministic, state.ξ.x)
    @test global_d ≈ deterministic atol=2e-9
    restart_horizon = PDMPSamplersRBridge._omrf_node_position_anchor_horizon!(
        restart_manager, state, flow, :main)
    live_horizon = PDMPSamplersRBridge._omrf_node_position_anchor_horizon!(
        manager, state, flow, :main)
    @test restart_horizon ≈ live_horizon
    # The complete no-change adapter path is allocation-free after prewarming.
    PDMPSamplersRBridge._select_omrf_node_position_core!(manager, model.grad,
        state.ξ.x, :main, state, false)
    live_scales = model.grad.envelope.component_scales!
    PDMPSamplers.component_cell_scales!(model.grad.envelope.cell_scales,
        model.grad.envelope, state, flow, 0.0, 0.2)
    cached_cells = live_scales.n_cached_cells
    cached_total = copy(live_scales.cell_total_masses)
    @test !PDMPSamplersRBridge._install_omrf_node_position_sources!(
        manager, model.grad, nodecv.sources)
    @test live_scales.n_cached_cells == cached_cells
    @test live_scales.cell_total_masses == cached_total
    # The real adapter horizon runs before the first event search.  It must
    # attach the copied envelope's observational counters immediately (not
    # after the first source change) and preserve the installed cap.
    live_scales.base.categorical_residual_cap = true
    preproposal_model = PDMPSamplers.with_stats(
        model, PDMPSamplers.DevelStatisticCounter())
    preproposal_scales = preproposal_model.grad.envelope.component_scales!
    @test preproposal_scales.base.context !== manager.oracle.residual_context
    @test preproposal_scales.base.categorical_residual_cap
    PDMPSamplers.adaptation_horizon(adapter, state, flow,
        preproposal_model.grad, :main, PDMPSamplers.DevelStatisticCounter())
    @test preproposal_scales.base.context === manager.oracle.residual_context
    @test preproposal_scales.base.categorical_residual_cap
    # A subsequent source change must retain authoritative counter ownership.
    old_sources = copy(nodecv.sources)
    nodecv.sources .= mod1.(nodecv.sources .+ 1, K)
    adapter.select_fn!(preproposal_model.grad, state, flow, :main,
        PDMPSamplers.DevelStatisticCounter(), :reflection, false)
    @test preproposal_scales.base.context === manager.oracle.residual_context
    nodecv.sources .= old_sources
    live_scales.base.categorical_residual_cap = false
    # Statistics wrapping deep-copies the envelope.  Matching source indices
    # must therefore still force ownership/generation synchronization before
    # its first proposal.
    wrapped_model = PDMPSamplers.with_stats(
        model, PDMPSamplers.DevelStatisticCounter())
    wrapped_scales = wrapped_model.grad.envelope.component_scales!
    @test wrapped_scales.base.node_bank !== nodecv.bank
    @test PDMPSamplersRBridge._install_omrf_node_position_sources!(
        manager, wrapped_model.grad, nodecv.sources)
    @test wrapped_scales.base.node_bank === nodecv.bank
    @test wrapped_scales.base.node_installation_generation ==
        nodecv.installation_generation
    @test wrapped_scales.base.node_sources == nodecv.sources
    @test @allocated(PDMPSamplersRBridge._install_omrf_node_position_sources!(
        manager, model.grad, nodecv.requested_sources)) == 0
    @test @allocated(PDMPSamplersRBridge._select_omrf_node_position_core!(
        manager, model.grad, state.ξ.x, :main, state, false)) == 0

    # Exact two-node synthetic geometry: simultaneous crossing at t=1,
    # persistent ties, and interrupted reflection/refresh/freeze/release all
    # invalidate the old owned deadline without artificial travel.
    for node in 1:2
        nodecv.supports[node] = [node]
        nodecv.linear_terms[node] = zeros(1, K)
        nodecv.constants[node, :] .= 100.0
        nodecv.linear_terms[node][1, 1] = 0.0
        nodecv.linear_terms[node][1, 2] = 2.0
        nodecv.constants[node, 1] = 0.0
        nodecv.constants[node, 2] = 4.0
    end
    nodecv.supports[3] = Int[]
    nodecv.sources .= 1
    synthetic = StickyPDMPState(0.0,
        SkeletonPoint(zeros(d), vcat([1.0, 1.0], zeros(d-2))), trues(d))
    synthetic.stored_velocity .= collect(range(-0.4, 0.5; length=d))
    @test PDMPSamplersRBridge._omrf_node_position_horizon_impl!(
        manager, synthetic) ≈ 1.0
    interrupted = move_forward_time(synthetic, 0.25, flow)
    interrupted.t[] = 0.25
    interrupted.ξ.θ[1:2] .*= -1
    PDMPSamplersRBridge._select_omrf_node_position!(manager, model.grad,
        interrupted.ξ.x; phase=:main, state=interrupted,
        boundary_won=false)
    @test nodecv.sources[1:2] == [1, 1]
    @test isinf(PDMPSamplersRBridge._omrf_node_position_horizon_impl!(
        manager, interrupted))
    interrupted.free[1:2] .= false
    interrupted.stored_velocity[1:2] .= 1.0
    interrupted.ξ.θ[1:2] .= 0.0 # actual frozen physical coordinates
    adapter.select_fn!(model.grad, interrupted, flow, :main,
        PDMPSamplers.DevelStatisticCounter(), :sticky_freeze, false)
    @test isinf(PDMPSamplersRBridge._omrf_node_position_horizon_impl!(
        manager, interrupted))
    interrupted.free[1:2] .= true
    interrupted.ξ.θ[1:2] .= interrupted.stored_velocity[1:2]
    interrupted.stored_velocity[1:2] .= 0.0
    adapter.select_fn!(model.grad, interrupted, flow, :main,
        PDMPSamplers.DevelStatisticCounter(), :sticky_release, false)
    @test PDMPSamplersRBridge._omrf_node_position_horizon_impl!(
        manager, interrupted) ≈ 0.75
    tied = move_forward_time(interrupted, 0.75, flow); tied.t[] = 1.0
    tied_position = copy(tied.ξ.x); tied_velocity = copy(tied.ξ.θ)
    tied_stored = copy(tied.stored_velocity); tied_free = copy(tied.free)
    tied.ξ.θ[1:2] .= 0.0
    PDMPSamplersRBridge._omrf_node_position_sources!(
        nodecv.requested_sources, nodecv, tied.ξ.x)
    @test nodecv.requested_sources[1:2] == [1, 1]
    tied.ξ.θ[1:2] .= 1.0
    PDMPSamplersRBridge._omrf_node_position_outgoing_sources!(
        nodecv.requested_sources, nodecv, tied.ξ.x, tied.ξ.θ)
    @test nodecv.requested_sources[1:2] == [2, 2]
    @test tied.ξ.x == tied_position
    @test tied.ξ.θ == tied_velocity
    @test tied.stored_velocity == tied_stored
    @test tied.free == tied_free

    # Prewarmed selector, crossing schedule, and certificate callbacks allocate
    # no memory and perform no retained observation scan.
    PDMPSamplersRBridge._omrf_node_position_sources!(
        nodecv.requested_sources, nodecv, state.ξ.x)
    PDMPSamplersRBridge._omrf_node_position_anchor_horizon!(
        manager, state, flow, :main)
    PDMPSamplers.component_cell_scales!(model.grad.envelope.cell_scales,
        model.grad.envelope, state, flow, 0.0, 0.2)
    scans = counts.native_full_data_scans
    @test @allocated(PDMPSamplersRBridge._omrf_node_position_sources!(
        nodecv.requested_sources, nodecv, state.ξ.x)) == 0
    @test @allocated(PDMPSamplersRBridge._omrf_node_position_anchor_horizon!(
        manager, state, flow, :main)) == 0
    @test @allocated(PDMPSamplers.component_cell_scales!(
        model.grad.envelope.cell_scales, model.grad.envelope,
        state, flow, 0.0, 0.2)) == 0
    @test counts.native_full_data_scans == scans
end

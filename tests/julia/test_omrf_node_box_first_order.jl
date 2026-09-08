using Test, Random, LinearAlgebra, PDMPSamplers, BridgeStan
include(normpath(joinpath(@__DIR__, "..", "..", "inst", "julia",
    "main_interface_function.jl")))
using .PDMPSamplersRBridge

function node_box_spec(X, seen; groups=2,
        bound_type="node_box_first_order")
    N, P = size(X)
    membership = zeros(Int, P, N)
    nodes = Int[]
    minima = Vector{Vector{Float64}}()
    maxima = Vector{Vector{Float64}}()
    offset = 0
    for node in 1:P
        local_membership = [1 + mod(person - 1, min(groups, N))
            for person in 1:N]
        group_count = maximum(local_membership)
        membership[node, :] .= offset .+ local_membership
        neighbours = [k for k in 1:P if k != node]
        for group in 1:group_count
            people = findall(==(group), local_membership)
            lower = zeros(P)
            upper = zeros(P)
            for neighbour in neighbours
                lower[neighbour] = minimum(X[people, neighbour])
                upper[neighbour] = maximum(X[people, neighbour])
            end
            push!(minima, lower)
            push!(maxima, upper)
            push!(nodes, node)
        end
        offset += group_count
    end
    return Dict{String,Any}(
        "type" => "omrf", "backend" => "analytic", "X" => X,
        "seen" => seen, "thresholds" => "thresholds_0",
        "interactions" => "interactions_0",
        "bound_type" => bound_type,
        "factorization" => "person_node",
        "factor_sampling" => "size_biased",
        "prior_backend" => "analytic", "use_hcv" => false,
        "analytic_prior" => Dict{String,Any}(
            "type" => "gaussian", "threshold_alpha" => 1.0,
            "threshold_beta" => 1.0, "interaction_scale" => 1.0),
        "exact_selected_refinement" => true,
        "node_box_membership" => membership,
        "node_box_nodes" => nodes,
        "node_box_minima" => reduce(vcat, permutedims.(minima)),
        "node_box_maxima" => reduce(vcat, permutedims.(maxima)))
end

@testset "OMRF node-box first-order Boomerang envelope" begin
    root = normpath(joinpath(@__DIR__, "..", ".."))
    lib = PDMPSamplersRBridge._compile_model_with_header(
        joinpath(root, "inst", "stan", "omrf", "omrf.stan"),
        joinpath(root, "inst", "stan", "pdmp_subsample.hpp"))
    full = joinpath(root, "tests", "stan", "omrf_subsampling_full.json")
    prior = joinpath(root, "tests", "stan", "omrf_subsampling_prior.json")
    ctx, names = PDMPSamplersRBridge._new_stan_subsampling_context(
        lib, full, prior, 1)
    X = [0 1 2; 2 0 1; 1 2 0]
    seen = [3, 3, 3]
    N, P = size(X)
    spec = node_box_spec(X, seen)
    anchor = zeros(9)
    oracle = PDMPSamplersRBridge.OMRFAnalyticSubsamplingOracle(
        ctx, spec, names, anchor)
    envelope = PDMPSamplersRBridge._build_omrf_residual_envelope(
        spec, oracle.anchor, oracle.residual_context)
    envelope2 = PDMPSamplersRBridge._build_omrf_residual_envelope(
        spec, copy(anchor), oracle.residual_context)
    cell_marked_spec = node_box_spec(X, seen;
        bound_type="cell_marked_node_box_first_order")
    cell_marked = PDMPSamplersRBridge._build_omrf_residual_envelope(
        cell_marked_spec, copy(anchor), oracle.residual_context)
    @test envelope isa PDMPSamplers.GroupedResidualEnvelope
    @test PDMPSamplers.n_observations(envelope) == length(X)
    @test envelope.component_scales! isa
        PDMPSamplersRBridge.OMRFNodeBoxFirstOrderScales
    @test cell_marked.component_scales! isa
        PDMPSamplersRBridge.OMRFCellMarkedNodeBoxFirstOrderScales
    @test envelope.component_scales!.displacement_lower !==
        envelope2.component_scales!.displacement_lower
    @test envelope.component_scales!.velocity_upper !==
        envelope2.component_scales!.velocity_upper

    # A full anchor bank recycles the inactive node-box context and grouped
    # envelope rather than reconstructing their O(NP) storage at every
    # refresh. This is also the allocation contract exercised in trajectories.
    subsampling = Dict{String,Any}(
        "n_observations" => N * P,
        "subsample_size" => 1,
        "residual_envelope" => spec)
    manager_ctx, manager_names =
        PDMPSamplersRBridge._new_stan_subsampling_context(
            lib, full, prior, 1)
    managed_model, _, _, managed_adapter, manager =
        PDMPSamplersRBridge._build_stan_subsampling_model(
            manager_ctx, manager_names, subsampling, anchor, nothing,
            falses(9); anchor_capacity=1)
    for shift in (0.01, 0.02, 0.03)
        prepared = PDMPSamplersRBridge._prepare_omrf_anchor(
            manager, anchor .+ shift)
        PDMPSamplersRBridge._insert_omrf_anchor!(manager, prepared)
    end
    recycled = manager.recycled_state
    recycled_context = recycled.residual_context
    recycled_envelope = recycled.envelope
    target_anchor = anchor .+ 0.04
    recycled_before = manager.recycled_preparations
    refreshed = PDMPSamplersRBridge._prepare_omrf_anchor(
        manager, target_anchor)
    @test manager.recycled_preparations == recycled_before + 1
    @test refreshed.residual_context === recycled_context
    @test refreshed.envelope === recycled_envelope
    @test refreshed.envelope.component_scales!.anchor == target_anchor

    # Regression for the live statistics-wrapped event path. Anchor-bank
    # selection must update the retained control variate supplied to the
    # adapter, rather than the unwrapped provider captured at construction.
    # Otherwise oracle and factor-certificate anchors diverge after a switch.
    live_model, _, _, live_adapter, live_manager =
        PDMPSamplersRBridge._build_stan_subsampling_model(
            manager_ctx, manager_names, subsampling, anchor, nothing,
            falses(9); anchor_capacity=2)
    wrapped_stats = PDMPSamplers.DevelStatisticCounter()
    wrapped_model = PDMPSamplers.with_stats(live_model, wrapped_stats)
    second_anchor = anchor .+ 0.08
    live_manager.main_refresh_distance = 0.01
    switching_state = StickyPDMPState(0.0,
        SkeletonPoint(copy(second_anchor), ones(9)), trues(9))
    switching_flow = BouncyParticle(9, 0.0)
    @test live_adapter.select_fn!(wrapped_model.grad, switching_state,
        switching_flow, :main, wrapped_stats)
    @test wrapped_model.grad.anchor == live_manager.oracle.anchor
    @test wrapped_model.grad.envelope.component_scales!.anchor ==
        live_manager.oracle.anchor
    @test wrapped_model.grad.residual_oracle.f.anchor ==
        live_manager.oracle.anchor
    # Restore a recycled owner, then measure the same prepared branch without
    # assigning its return value (which would box a local in global test code).
    manager.recycled_state = recycled
    allocated_refresh = @allocated PDMPSamplersRBridge._prepare_omrf_anchor(
        manager, target_anchor)
    @test allocated_refresh == 0

    # The stored factor order is node-major, matching the analytic oracle.
    for node in 1:P, person in 1:N
        factor = (node - 1) * N + person
        component = envelope.groups[1, factor]
        @test spec["node_box_nodes"][component] == node
        @test factor in envelope.members[component]
    end

    flow = PDMPSamplersRBridge.build_flow("AdaptiveBoomerang",
        Matrix(Diagonal(range(0.6, 1.4; length=9))), zeros(9);
        refresh_rate=0.1)
    rng = Xoshiro(913)
    position = 0.15randn(rng, 9)
    velocity = PDMPSamplers.initialize_velocity(rng, flow, 9)
    state = PDMPState(0.0,
        SkeletonPoint(copy(position), copy(velocity)))
    PDMPSamplers.initialize_flow_state!(state, flow)

    # The G=2 shared incident-edge traversal is exactly equivalent to the
    # generic component-at-a-time implementation for point and cell bounds.
    specialized = similar(envelope.scales)
    generic = similar(envelope.scales)
    envelope.component_scales!(specialized, state, flow, 0.17)
    @inbounds for component in eachindex(generic)
        generic[component] =
            PDMPSamplersRBridge._omrf_node_box_first_order_bound(
                envelope.component_scales!, component,
                envelope.component_scales!.component_nodes[component])
    end
    @test specialized == generic
    envelope.component_scales!(specialized, state, flow, -0.3, 0.6)
    @inbounds for component in eachindex(generic)
        generic[component] =
            PDMPSamplersRBridge._omrf_node_box_first_order_bound(
                envelope.component_scales!, component,
                envelope.component_scales!.component_nodes[component])
    end
    @test specialized == generic

    # Ordinary first-order decomposition: anchor deterministic term plus all
    # person-node residuals equals the complete Stan gradient.
    x = 0.2randn(rng, 9)
    deterministic = zeros(9)
    residual = zeros(9)
    complete = zeros(9)
    oracle(deterministic, x)
    oracle(residual, x, collect(1:(N * P)), oracle.anchor)
    PDMPSamplersRBridge._clear_gradient!(ctx, :full, complete, x)
    @test deterministic + residual ≈ complete atol=2e-9 rtol=2e-10

    # Exact harmonic interval regression, including extrema and wide cells.
    for coefficients in ((0.0, 1.0, 0.0), (0.3, -0.8, 1.2),
                         (-1.0, 0.0, 0.0)),
        interval in ((-pi/2, pi/2), (0.0, 2pi), (-4pi, 3pi), (-0.2, 0.4))
        c, a, b = coefficients
        left, right = interval
        lower, upper = PDMPSamplersRBridge._omrf_harmonic_cell_interval(
            c, a, b, left, right)
        for t in range(left, right; length=257)
            value = c + a*cos(t) + b*sin(t)
            @test lower <= value <= upper
        end
    end

    # Dense random factor/cell checks cover anchor crossings and harmonic
    # extrema without evaluating people in either envelope callback.
    selected = zeros(9)
    maximum_violation = -Inf
    maximum_factor_point_violation = -Inf
    maximum_cell_mark_violation = -Inf
    for iteration in 1:20_001
        state.ξ.x .= 0.4randn(rng, 9)
        state.ξ.θ .= PDMPSamplers.initialize_velocity(rng, flow, 9)
        left = iteration % 7 == 0 ? -pi/2 : 2pi*rand(rng) - pi
        right = iteration % 11 == 0 ? left + 2pi : left + 0.01 + rand(rng)
        PDMPSamplers.component_cell_scales!(envelope.cell_scales,
            envelope, state, flow, left, right)
        cell_marked.component_scales!.n_cached_cells = 0
        PDMPSamplers.component_cell_scales!(cell_marked.cell_scales,
            cell_marked, state, flow, left, right)
        @test cell_marked.cell_scales == envelope.cell_scales
        copyto!(envelope.scales, envelope.cell_scales)
        factor = rand(rng, 1:(N * P))
        t = left + rand(rng) * (right - left)
        candidate = move_forward_time(state, t, flow)
        fill!(selected, 0.0)
        oracle(selected, candidate.ξ.x, [factor], oracle.anchor)
        exact = PDMPSamplersRBridge._omrf_exact_selected_residual_abs_rate(
            candidate, selected, flow,
            oracle.residual_context.touched_positions)
        bound = PDMPSamplers.observation_residual_bound(envelope, factor)
        maximum_violation = max(maximum_violation, exact - bound)
        @test exact <= bound * (1 + 2e-7) + 2e-8

        # The cell-marked route retrieves the already constructed cell and
        # evaluates only the selected node/person at the point. Its two
        # thinning inequalities are checked separately.
        B_cell = PDMPSamplers.screening_residual_bound(
            cell_marked, state, flow, t)
        @test B_cell == dot(cell_marked.cell_scales, cell_marked.totals)
        point_bound =
            PDMPSamplersRBridge._omrf_cell_marked_factor_point_bound(
                cell_marked.component_scales!, factor)
        component = cell_marked.groups[1, factor]
        component_cell_bound = cell_marked.scales[component]
        maximum_factor_point_violation = max(
            maximum_factor_point_violation, exact - point_bound)
        maximum_cell_mark_violation = max(
            maximum_cell_mark_violation, point_bound - component_cell_bound)
        @test exact <= point_bound * (1 + 2e-7) + 2e-8
        @test point_bound <= component_cell_bound * (1 + 64eps(Float64))
    end
    @test maximum_violation <= 2e-8
    @test maximum_factor_point_violation <= 2e-8
    @test maximum_cell_mark_violation <= 0.0

    # Inactive interactions and cells crossing the anchor remain certified.
    free = trues(9)
    free[oracle.residual_context.interaction_indices[1]] = false
    frozen = StickyPDMPState(0.0,
        SkeletonPoint(copy(anchor), copy(velocity)), free)
    frozen.ξ.θ[oracle.residual_context.interaction_indices[1]] = 0.0
    PDMPSamplers.component_cell_scales!(envelope.cell_scales,
        envelope, frozen, flow, -0.4, 0.4)
    @test all(isfinite, envelope.cell_scales)
    @inbounds for component in eachindex(generic)
        generic[component] =
            PDMPSamplersRBridge._omrf_node_box_first_order_bound(
                envelope.component_scales!, component,
                envelope.component_scales!.component_nodes[component])
    end
    @test envelope.cell_scales == generic

    # Component and within-box draws reproduce factor-level size-biased mass.
    B = PDMPSamplers.total_residual_bound(envelope, state, flow, 0.17)
    cv = SubsampledControlVariate((out, x) -> fill!(out, 0.0), oracle,
        envelope, anchor, 1)
    expected = [PDMPSamplers.observation_residual_bound(envelope, factor) / B
        for factor in 1:(N * P)]
    counts = zeros(Int, N * P)
    draw_rng = Xoshiro(77)
    for _ in 1:100_000
        PDMPSamplers.draw_subset!(draw_rng, cv, 0.0, B)
        counts[cv.subset[1]] += 1
    end
    empirical = counts ./ sum(counts)
    @test maximum(abs.(empirical .- expected)) < 0.007

    # Component-first selection uses the cached closed-cell mass. Component
    # sizes are unequal in this N=3/G=2 fixture, so this also guards against
    # accidentally sampling components uniformly.
    cell_cv = SubsampledControlVariate((out, x) -> fill!(out, 0.0), oracle,
        cell_marked, anchor, 1)
    cell_marked.component_scales!.n_cached_cells = 0
    PDMPSamplers.component_cell_scales!(cell_marked.cell_scales,
        cell_marked, state, flow, 0.0, 0.4)
    B_cell = PDMPSamplers.screening_residual_bound(
        cell_marked, state, flow, 0.17)
    cell_expected = [cell_marked.scales[cell_marked.groups[1, factor]] /
        B_cell for factor in 1:(N * P)]
    cell_counts = zeros(Int, N * P)
    component_counts = zeros(Int, length(cell_marked.totals))
    cell_rng = Xoshiro(78)
    for _ in 1:100_000
        PDMPSamplers.draw_subset!(cell_rng, cell_cv, 0.0, B_cell)
        factor = cell_cv.subset[1]
        cell_counts[factor] += 1
        component_counts[cell_marked.groups[1, factor]] += 1
    end
    cell_empirical = cell_counts ./ sum(cell_counts)
    @test maximum(abs.(cell_empirical .- cell_expected)) < 0.007
    component_expected = cell_marked.scales .* cell_marked.totals ./ B_cell
    @test maximum(abs.(component_counts ./ sum(component_counts) .-
        component_expected)) < 0.007

    # The two thinning factors reduce each component-marked cell intensity to
    # its factor-level point certificate exactly, including unequal component
    # sizes. Zero-mass components contribute exactly zero.
    for factor in 1:(N * P)
        component = cell_marked.groups[1, factor]
        nc = cell_marked.totals[component]
        cell_bound = cell_marked.scales[component]
        point_bound =
            PDMPSamplersRBridge._omrf_cell_marked_factor_point_bound(
                cell_marked.component_scales!, factor)
        reached = B_cell * (nc * cell_bound / B_cell) * (1 / nc) *
            (point_bound / cell_bound)
        @test reached ≈ point_bound atol=2e-14 rtol=2e-14
    end
    synthetic_roofs = [0.0, 0.7, 1.1]
    synthetic_sizes = [2.0, 3.0, 1.0]
    synthetic_mass = dot(synthetic_roofs, synthetic_sizes)
    @test synthetic_sizes[1] * synthetic_roofs[1] / synthetic_mass == 0.0

    # Both deterministic and residual mixture branches use the same selected
    # cell certificate. A very large deterministic mass exercises the former;
    # D=0 exercises the latter and increments the component-mark count.
    marks_before = cell_marked.component_scales!.component_marks
    deterministic_rng = Xoshiro(190)
    PDMPSamplers.draw_subset!(deterministic_rng, cell_cv, 1e300, B_cell)
    @test cell_marked.component_scales!.component_marks == marks_before
    PDMPSamplers.draw_subset!(deterministic_rng, cell_cv, 0.0, B_cell)
    @test cell_marked.component_scales!.component_marks == marks_before + 1

    # With batch one, proposal mass times selection and exact thinning equals
    # the intended uniform-mixture marked intensity exactly.
    marked_sum = 0.0
    recovered_sum = 0.0
    scale = N * P
    for factor in 1:(N * P)
        candidate = move_forward_time(state, 0.17, flow)
        fill!(selected, 0.0)
        oracle(selected, candidate.ξ.x, [factor], oracle.anchor)
        exact = PDMPSamplers.λ(candidate,
            deterministic .+ scale .* selected, flow)
        selected_bound = scale *
            PDMPSamplers.observation_residual_bound(envelope, factor)
        probability = PDMPSamplers._selected_subset_probability(
            cv, 0.0, B, selected_bound)
        marked_sum += exact / (N * P)
        recovered_sum += B * probability * exact / selected_bound
    end
    @test recovered_sum ≈ marked_sum atol=2e-12 rtol=2e-12

    # Prewarmed point, cell, selection, and exact-refinement paths allocate
    # no memory. The diagnostic arrays above are deliberately outside this.
    PDMPSamplers.total_residual_bound(envelope, state, flow, 0.1)
    PDMPSamplers.component_cell_scales!(envelope.cell_scales,
        envelope, state, flow, 0.0, 0.3)
    PDMPSamplers.draw_subset!(draw_rng, cv, 0.0,
        dot(envelope.scales, envelope.totals))
    deterministic_state = zeros(9)
    oracle(deterministic_state, state.ξ.x)
    D = PDMPSamplers.λ(state, deterministic_state, flow)
    selected_roof = D + (N * P) *
        PDMPSamplers.observation_residual_bound(envelope, cv.subset[1])
    PDMPSamplers.subsampling_residual_subset_bound(oracle, state,
        deterministic_state, flow, D, selected_roof, cv.subset, N * P)
    @test @allocated(PDMPSamplers.total_residual_bound(
        envelope, state, flow, 0.1)) == 0
    @test @allocated(PDMPSamplers.component_cell_scales!(envelope.cell_scales,
        envelope, state, flow, 0.0, 0.3)) == 0
    @test @allocated(PDMPSamplers.draw_subset!(draw_rng, cv, 0.0,
        dot(envelope.scales, envelope.totals))) == 0
    @test @allocated(PDMPSamplers.subsampling_residual_subset_bound(
        oracle, state, deterministic_state, flow, D, selected_roof,
        cv.subset, N * P)) == 0
    cell_marked.component_scales!.n_cached_cells = 0
    PDMPSamplers.component_cell_scales!(cell_marked.cell_scales,
        cell_marked, state, flow, 0.0, 0.3)
    cell_B = PDMPSamplers.screening_residual_bound(
        cell_marked, state, flow, 0.1)
    PDMPSamplers.draw_subset!(cell_rng, cell_cv, 0.0, cell_B)
    @test @allocated(PDMPSamplers.component_cell_scales!(
        cell_marked.cell_scales, cell_marked, state, flow, 0.3, 0.6)) == 0
    @test @allocated(PDMPSamplers.screening_residual_bound(
        cell_marked, state, flow, 0.4)) == 0
    @test @allocated(PDMPSamplers.draw_subset!(
        cell_rng, cell_cv, 0.0, cell_B)) == 0
    println("node_box_maximum_certificate_violation=", maximum_violation)
    println("cell_marked_factor_maximum_certificate_violation=",
        maximum_factor_point_violation)
    println("cell_marked_cell_maximum_certificate_violation=",
        maximum_cell_mark_violation)
end

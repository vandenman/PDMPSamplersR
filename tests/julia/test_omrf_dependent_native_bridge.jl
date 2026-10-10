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

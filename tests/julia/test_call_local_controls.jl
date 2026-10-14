using Test, PDMPSamplers, LinearAlgebra
include(joinpath(@__DIR__, "..", "..", "inst", "julia", "main_interface_function.jl"))
const RB = PDMPSamplersRBridge

@testset "explicit controls ignore ambient environment and do not persist" begin
    withenv("PDMP_ADAPTIVE_BOOMERANG_LAMBDA_REF"=>"invalid",
            "PDMP_GRID_N_MIN"=>"invalid", "PDMP_GRID_BOUND_VIOLATION"=>"invalid",
            "PDMP_ADAPTIVE_WARMUP_STOP"=>"invalid",
            "PDMPSAMPLERS_HCV_CACHE"=>"0") do
        covariance = Matrix{Float64}(I, 2, 2)
        control = Dict(:adaptive_boomerang_refresh_rate=>0.7,
            :grid_n_min=>3, :grid_bound_violation=>"throw",
            :warmup_stop=>Dict(:min_time=>1.0, :max_time=>4.0, :stable_time=>0.5,
                              :min_events=>12, :check_every=>3))
        custom = RB.build_flow("AdaptiveBoomerang", covariance, zeros(2); control)
        ordinary = RB.build_flow("AdaptiveBoomerang", covariance, zeros(2))
        @test PDMPSamplers.refresh_rate(custom) == 0.7
        @test PDMPSamplers.refresh_rate(ordinary) == 0.1
        custom_grid = RB.build_algorithm("GridThinningStrategy";
            c0=0.1, d=2, grid_n=10, grid_t_max=1.0, control)
        ordinary_grid = RB.build_algorithm("GridThinningStrategy";
            c0=0.1, d=2, grid_n=10, grid_t_max=1.0)
        @test custom_grid.N_min == 3
        @test ordinary_grid.N_min == 5
        stop = RB.build_warmup_stop(4.0, control)
        @test stop.min_time == 1.0
        @test stop.max_time == 4.0
        @test stop.stable_time == 0.5
        @test stop.min_events == 12
        @test stop.check_every == 3
        @test isnothing(RB.build_warmup_stop(4.0))
        @test RB.build_warmup_stop(4.0, Dict(:warmup_stop=>Dict())).min_time == 1.0
        cache = RB.OMRFHCVProductCache(2; enabled=false, research_diagnostics=true,
                                      check_fast_deterministic=true)
        ordinary_cache = RB.OMRFHCVProductCache(2)
        @test !cache.enabled && cache.research_diagnostics && cache.check_fast_deterministic
        @test ordinary_cache.enabled && !ordinary_cache.research_diagnostics &&
              !ordinary_cache.check_fast_deterministic
    end
end

@testset "compile cache paths use explicit directories and compiler arguments" begin
    mktempdir() do dir
        source = joinpath(dir, "test.stan")
        write(source, "parameters { real x; } model { x ~ normal(0,1); }")
        withenv("PDMPSAMPLERSR_BRIDGESTAN_CACHE_DIR"=>"ignored") do
            a = RB._cache_build_dir(source, String[], String[]; cache_root=joinpath(dir,"a"))
            b = RB._cache_build_dir(source, ["--O1"], String[]; cache_root=joinpath(dir,"a"))
            c = RB._cache_build_dir(source, String[], String[]; cache_root=joinpath(dir,"b"))
            @test dirname(a) == joinpath(dir,"a")
            @test a != b
            @test a != c
        end
    end
end

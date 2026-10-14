using Test, PDMPSamplers, BridgeStan
include(joinpath(@__DIR__, "..", "..", "inst", "julia", "main_interface_function.jl"))
const RB = PDMPSamplersRBridge

csv_chains(trace) = PDMPChains([trace], [nothing], [nothing], [nothing], [nothing])
csv_trace(times) = PDMPSamplers.PDMPTrace(times, reshape(copy(times), 1, :),
    ones(1, length(times)), BouncyParticle(1))

struct BrokenCSVTrace <: PDMPSamplers.AbstractPDMPTrace end
struct CSVDiscretizationError <: Exception end
PDMPSamplers.first_event_time(::BrokenCSVTrace) = 0.0
PDMPSamplers.last_event_time(::BrokenCSVTrace) = 1.0
Base.first(::BrokenCSVTrace) = (position=[0.0],)
PDMPSamplers.PDMPDiscretize(::BrokenCSVTrace, ::Float64) = throw(CSVDiscretizationError())
PDMPSamplers.adaptive_discretize(::BrokenCSVTrace; kwargs...) = throw(CSVDiscretizationError())

@testset "CSV discretization handles only zero-duration degeneracy" begin
    for times in ([0.0], [2.0], [0.0, 0.0]), dt in (0.01, 0.0)
        @test RB._draws_for_csv(csv_chains(csv_trace(times)), 1, dt) == reshape([first(times)], 1, 1)
    end
    # A positive-duration path without a bounce must still be discretized.
    for dt in (0.01, 0.0)
        @test size(RB._draws_for_csv(csv_chains(csv_trace([0.0, 0.1])), 1, dt), 1) > 1
        @test_throws BoundsError RB._draws_for_csv(csv_chains(csv_trace(Float64[])), 1, dt)
        @test_throws CSVDiscretizationError RB._draws_for_csv(csv_chains(BrokenCSVTrace()), 1, dt)
    end
end

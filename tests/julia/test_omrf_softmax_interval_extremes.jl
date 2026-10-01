using Test, Random, PDMPSamplers, BridgeStan
include(normpath(joinpath(@__DIR__, "..", "..", "inst", "julia",
    "main_interface_function.jl")))
using .PDMPSamplersRBridge

function exact_probability(logits, category)
    values = BigFloat[0; logits]
    m = maximum(values)
    probabilities = exp.(values .- m)
    Float64(probabilities[category + 1] / sum(probabilities))
end

@testset "OMRF softmax interval remains finite for extreme cells" begin
    ordinary_lo = [-2.0, 0.5, 3.0]
    ordinary_hi = [-1.0, 1.0, 4.0]
    extreme_lo = [-1000.0, -4.0, 900.0]
    extreme_hi = [-900.0, 3.0, 1000.0]
    for (lo, hi) in ((ordinary_lo, ordinary_hi),
                     (extreme_lo, extreme_hi))
        for category in 0:3
            lower, upper =
                PDMPSamplersRBridge._omrf_softmax_probability_interval(
                    lo, hi, category, 3)
            lower_logits = BigFloat.(hi)
            upper_logits = BigFloat.(lo)
            if category != 0
                lower_logits[category] = BigFloat(lo[category])
                upper_logits[category] = BigFloat(hi[category])
            end
            @test 0.0 <= lower <= upper <= 1.0 + 1e-14
            @test lower <= exact_probability(lower_logits, category)
            @test upper >= exact_probability(upper_logits, category)
        end
    end
    # Ordinary finite cells retain the pre-repair arithmetic exactly.
    for category in 0:3
        target_lo = category == 0 ? 0.0 : ordinary_lo[category]
        target_hi = category == 0 ? 0.0 : ordinary_hi[category]
        lower_den = 1.0; upper_den = 1.0
        for other in 0:3
            other == category && continue
            other_lo = other == 0 ? 0.0 : ordinary_lo[other]
            other_hi = other == 0 ? 0.0 : ordinary_hi[other]
            lower_den = PDMPSamplersRBridge._omrf_up(lower_den +
                exp(other_hi - target_lo))
            upper_den = PDMPSamplersRBridge._omrf_down(upper_den +
                exp(other_lo - target_hi))
        end
        legacy = (PDMPSamplersRBridge._omrf_down(inv(lower_den)),
            PDMPSamplersRBridge._omrf_up(inv(upper_den)))
        @test PDMPSamplersRBridge._omrf_softmax_probability_interval(
            ordinary_lo, ordinary_hi, category, 3) == legacy
    end
    rng = Xoshiro(9152026)
    for _ in 1:100
        centre = 1200 .* (rand(rng, 3) .- 0.5)
        width = 200 .* rand(rng, 3)
        lo = centre .- width
        hi = centre .+ width
        for category in 0:3
            lower, upper =
                PDMPSamplersRBridge._omrf_softmax_probability_interval(
                    lo, hi, category, 3)
            lower_logits = BigFloat.(hi)
            upper_logits = BigFloat.(lo)
            if category != 0
                lower_logits[category] = BigFloat(lo[category])
                upper_logits[category] = BigFloat(hi[category])
            end
            exact_lower = exact_probability(lower_logits, category)
            exact_upper = exact_probability(upper_logits, category)
            @test lower <= exact_lower || isapprox(lower, exact_lower;
                rtol=1e-12, atol=nextfloat(0.0))
            @test upper >= exact_upper || isapprox(upper, exact_upper;
                rtol=1e-12, atol=nextfloat(0.0))
        end
    end
end

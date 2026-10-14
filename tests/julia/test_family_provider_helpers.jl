using Test, PDMPSamplers, BridgeStan
include(joinpath(@__DIR__, "..", "..", "inst", "julia", "main_interface_function.jl"))
const RB = PDMPSamplersRBridge

@testset "reference-category softmax preserves input and supports aliasing" begin
    for logits in ([0.0], [30.0, -30.0], [-30.0, -29.0],
            [0.2, -1.3, 0.7], [1000.0, 999.0], [-1000.0, -999.0])
        original = copy(logits)
        out = similar(logits)
        @test RB._reference_category_softmax!(out, logits) === out
        @test logits == original
        high_precision = exp.(BigFloat.(logits))
        expected = Float64.(high_precision ./ (1 + sum(high_precision)))
        @test out ≈ expected rtol=2e-15 atol=0.0
        aliased = copy(logits)
        @test RB._reference_category_softmax!(aliased, aliased) === aliased
        @test isequal(out, aliased)
        @test all(isfinite, out)
    end
    # A zero reference logit is part of the denominator, even though omitted
    # from the output: three equally likely categories leave two probabilities.
    @test RB._reference_category_softmax!(zeros(2), zeros(2)) == fill(1/3, 2)
    matrix = [30.0 -30.0; 0.2 0.7]
    expected = RB._reference_category_softmax!(zeros(2), copy(matrix[1, :]))
    RB._reference_category_softmax!(@view(matrix[1, :]), @view(matrix[1, :]))
    @test matrix[1, :] == expected
    @test matrix[2, :] == [0.2, 0.7]
end

using Test
using LinearAlgebra

@testset "Boomerang bridge preserves diagonal covariance structure" begin
    covariance = Matrix(Diagonal([0.25, 1.5, 8.0, 0.04]))
    precision = PDMPSamplersRBridge._to_precision(covariance, 4)
    @test precision isa Diagonal
    @test precision.diag == 1.0 ./ diag(covariance)

    flow = PDMPSamplersRBridge.build_flow("Boomerang", precision, zeros(4))
    @test flow isa PDMPSamplers.Boomerang
    @test flow.Γ isa Diagonal
    @test flow.L isa Diagonal
    @test flow.ΣL isa Diagonal
    @test abs2.(flow.ΣL.diag) ≈ diag(covariance)

    bps = PDMPSamplersRBridge.build_flow(
        "PreconditionedBPS", precision, zeros(4); refresh_rate=0.7)
    zigzag = PDMPSamplersRBridge.build_flow(
        "PreconditionedZigZag", precision, zeros(4))
    @test bps.metric isa PDMPSamplers.DiagonalPreconditioner
    @test zigzag.metric isa PDMPSamplers.DiagonalPreconditioner
    @test bps.metric.scale ≈ sqrt.(diag(covariance))
    @test zigzag.metric.scale ≈ sqrt.(diag(covariance))
    @test !(bps.metric.scale == ones(4))

    dense_covariance = [1.0 0.2; 0.2 2.0]
    dense_precision = PDMPSamplersRBridge._to_precision(dense_covariance, 2)
    @test !(dense_precision isa Diagonal)
    @test Matrix(dense_precision) ≈ inv(dense_covariance)
end

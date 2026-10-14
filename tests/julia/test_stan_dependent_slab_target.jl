using Test, LinearAlgebra, BridgeStan, PDMPSamplers

@testset "Stan dependent slab matches manually constructed active faces" begin
    root = normpath(joinpath(@__DIR__, "..", ".."))
    stan = joinpath(root, "tests", "stan", "dependent_slab_face.stan")
    lib = BridgeStan.compile_model(stan)
    full_data = joinpath(root, "tests", "stan", "dependent_slab_full.json")
    full_model = PDMPModel(BridgeStan.StanModel(lib, full_data; warn=false))
    covariance = [1.1 0.35; 0.35 0.8]
    slab = DenseGaussianSlab([0.2, -0.1], covariance, [1, 2])
    prior = BernoulliModelPrior([0.5, 0.5])
    cases = (
        (BitVector([false, false, true]), "dependent_slab_empty.json"),
        (BitVector([true, false, true]), "dependent_slab_first.json"),
        (BitVector([false, true, true]), "dependent_slab_second.json"),
        (BitVector([true, true, true]), "dependent_slab_both.json"),
    )
    for (free, filename) in cases
        x = [0.45, -0.3, 1.4]
        free[1] || (x[1] = 0.0)
        free[2] || (x[2] = 0.0)
        target = DependentSlabTarget(3, full_model.grad, slab, prior;
            initial_free=free)
        corrected = zeros(3)
        target(corrected, x)
        manual_model = PDMPModel(BridgeStan.StanModel(
            lib, joinpath(root, "tests", "stan", filename); warn=false))
        manual = zeros(3)
        compute_gradient!(manual_model.grad, x, manual)
        @test corrected ≈ manual atol=2e-10 rtol=2e-10
        @test abs(corrected[3]) > 0.1
    end
end

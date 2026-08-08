using Test, LinearAlgebra, BridgeStan

json_vector(x) = "[" * join(string.(x), ",") * "]"
json_matrix(x) = "[" * join((json_vector(row) for row in eachrow(x)), ",") * "]"

@testset "RI-CLPM sufficient statistics equal the rowwise likelihood" begin
    root = normpath(joinpath(@__DIR__, "..", ".."))
    sufficient_stan = joinpath(root, "inst", "stan", "riclpm", "riclpm_sufficient.stan")
    direct_stan = joinpath(root, "tests", "stan", "riclpm_direct.stan")
    sufficient_lib = BridgeStan.compile_model(sufficient_stan)
    direct_lib = BridgeStan.compile_model(direct_stan)

    N, P, T = 7, 2, 3
    z = [sin(0.31 * n + 0.17 * j) + 0.08 * n * j / (N * P * T)
         for n in 1:N, j in 1:(P * T)]
    z_sum = vec(sum(z; dims=1))
    z_crossprod = z' * z
    common = "\"N\":$N,\"P\":$P,\"T\":$T," *
        "\"jitter\":1e-8,\"slab_mean\":[0.05,-0.03]," *
        "\"slab_cov\":[[0.20,0.04],[0.04,0.16]]," *
        "\"mean_prior_sd\":1.5,\"lag_prior_sd\":0.5," *
        "\"log_sd_prior_mean\":-0.4,\"log_sd_prior_sd\":0.6," *
        "\"rho_prior_sd\":0.7"
    sufficient_data = "{" * common *
        ",\"z_sum\":" * json_vector(z_sum) *
        ",\"z_crossprod\":" * json_matrix(z_crossprod) * "}"
    direct_data = "{" * common * ",\"z\":" * json_matrix(z) * "}"
    sufficient = BridgeStan.StanModel(sufficient_lib, sufficient_data; warn=false)
    direct = BridgeStan.StanModel(direct_lib, direct_data; warn=false)
    @test BridgeStan.param_unc_names(sufficient) == BridgeStan.param_unc_names(direct)
    @test BridgeStan.param_num(sufficient) == 19

    points = (
        zeros(19),
        collect(range(-0.18, 0.21; length=19)),
        [0.04 * sin(0.7 * j) - 0.03 * cos(0.2 * j) for j in 1:19],
    )
    for q in points
        sufficient_gradient = zeros(19)
        direct_gradient = zeros(19)
        sufficient_lp, _ = BridgeStan.log_density_gradient!(
            sufficient, q, sufficient_gradient; propto=false, jacobian=true)
        direct_lp, _ = BridgeStan.log_density_gradient!(
            direct, q, direct_gradient; propto=false, jacobian=true)
        @test sufficient_lp ≈ direct_lp atol=2e-9 rtol=2e-10
        @test sufficient_gradient ≈ direct_gradient atol=2e-8 rtol=2e-9
    end
end

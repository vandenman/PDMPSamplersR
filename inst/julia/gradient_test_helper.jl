using BridgeStan

function r_eval_stan_gradient(
        stan_file::String,
        data_json::String,
        theta::AbstractVector;
        hpp_path::String = "",
        indices_0based::AbstractVector = Int32[])
    theta = Vector{Float64}(theta)
    indices_0based = Vector{Int32}(indices_0based)

    lib = if !isempty(hpp_path)
        hpp_norm = replace(normpath(hpp_path), "\\" => "/")
        if endswith(stan_file, ".stan")
            BridgeStan.compile_model(stan_file;
                stanc_args=["--allow-undefined"],
                make_args=["USER_HEADER=$(hpp_norm)"])
        else
            stan_file
        end
    else
        if endswith(stan_file, ".stan")
            BridgeStan.compile_model(stan_file)
        else
            stan_file
        end
    end

    sm = BridgeStan.StanModel(lib, data_json; warn=false)

    if !isempty(indices_0based)
        dl = Libc.Libdl.dlopen(lib, Libc.Libdl.RTLD_NOLOAD | Libc.Libdl.RTLD_GLOBAL)
        set_fn = Libc.Libdl.dlsym(dl, :pdmp_set_subsample_indices)
        ccall(set_fn, Cvoid, (Ptr{Int32}, Int32), indices_0based, Int32(length(indices_0based)))
    end

    d = Int(BridgeStan.param_unc_num(sm))
    if isempty(theta)
        theta = ones(d) .* 0.1
    end

    grad = zeros(d)
    lp, _ = BridgeStan.log_density_gradient!(sm, theta, grad)

    Dict("lp" => lp, "grad" => grad, "d" => d, "theta" => theta)
end

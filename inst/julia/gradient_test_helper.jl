using BridgeStan

function r_eval_stan_gradient(stan_file::String, data_json::String, theta::Vector{Float64};
                              hpp_path::String="", indices_0based::Vector{<:Integer}=Int[])
    lib_path = if endswith(stan_file, ".stan")
        if !isempty(hpp_path)
            hpp_for_make = replace(normpath(hpp_path), "\\" => "/")
            BridgeStan.compile_model(stan_file;
                stanc_args=["--allow-undefined"],
                make_args=["USER_HEADER=$(hpp_for_make)"])
        else
            BridgeStan.compile_model(stan_file)
        end
    else
        stan_file
    end

    sm = BridgeStan.StanModel(lib_path, data_json; warn=false)
    d = Int(BridgeStan.param_unc_num(sm))

    if !isempty(indices_0based)
        set_fn = Libc.Libdl.dlsym(
            Libc.Libdl.dlopen(lib_path, Libc.Libdl.RTLD_NOLOAD | Libc.Libdl.RTLD_GLOBAL),
            :pdmp_set_subsample_indices
        )
        idx_buf = Int32.(indices_0based)
        ccall(set_fn, Cvoid, (Ptr{Int32}, Int32), idx_buf, Int32(length(idx_buf)))
    end

    if isempty(theta)
        theta = zeros(d)
    end

    grad = zeros(d)
    lp = BridgeStan.log_density_gradient!(sm, theta, grad)
    return Dict("d" => d, "lp" => lp, "grad" => grad, "theta" => copy(theta))
end

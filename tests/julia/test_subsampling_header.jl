using Test, Libdl

header = normpath(joinpath(@__DIR__, "..", "..", "inst", "stan",
    "pdmp_subsample.hpp"))

@testset "subset header setter/getter round trip" begin
    compilers = filter(!isnothing,
        [Sys.which("c++"), Sys.which("g++"), Sys.which("clang++")])
    cxx = isempty(compilers) ? nothing : first(compilers)
    if cxx === nothing
        @test_skip "no C++ compiler"
    else
        mktempdir() do directory
            source = joinpath(directory, "roundtrip.cpp")
            library = joinpath(directory, "roundtrip." * Libdl.dlext)
            open(source, "w") do io
                println(io, "#include \"", replace(header, "\\" => "\\\\"), "\"")
            end
            run(`$cxx -std=c++17 -shared -fPIC $source -o $library`)
            lib = Libdl.dlopen(library)
            set_fn = Libdl.dlsym(lib, :pdmp_set_subsample_indices)
            get_size = Libdl.dlsym(lib, :pdmp_get_subsample_size)
            get_index = Libdl.dlsym(lib, :pdmp_get_subsample_index)
            clear_fn = Libdl.dlsym(lib, :pdmp_clear_subsample_indices)
            indices = Int32[4, 1, 8]
            GC.@preserve indices begin
                ccall(set_fn, Cvoid, (Ptr{Int32}, Cint),
                    indices, length(indices))
            end
            @test ccall(get_size, Cint, ()) == 3
            @test [ccall(get_index, Cint, (Cint,), i) for i in 0:2] == indices
            ccall(clear_fn, Cvoid, ())
            @test ccall(get_size, Cint, ()) == 0
            Libdl.dlclose(lib)
        end
    end
end

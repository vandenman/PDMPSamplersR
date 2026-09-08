using Test, Random, LinearAlgebra, PDMPSamplers, BridgeStan
include(normpath(joinpath(@__DIR__, "..", "..", "inst", "julia", "main_interface_function.jl")))
using .PDMPSamplersRBridge

@testset "OMRF norm-local Boomerang callback" begin
    root = normpath(joinpath(@__DIR__, "..", ".."))
    lib = PDMPSamplersRBridge._compile_model_with_header(
        joinpath(root,"inst","stan","omrf","omrf.stan"),
        joinpath(root,"inst","stan","pdmp_subsample.hpp"))
    full = joinpath(root,"tests","stan","omrf_subsampling_full.json")
    prior = joinpath(root,"tests","stan","omrf_subsampling_prior.json")
    ctx,names = PDMPSamplersRBridge._new_stan_subsampling_context(lib,full,prior,1)
    X = [0 1 2; 2 0 1; 1 2 0]; P=3; N=3; anchor=zeros(9)
    weights=zeros(3P,N)
    for j=1:P,i=1:N
        r=sqrt(sum(abs2,X[i,k] for k=1:P if k!=j))
        weights[j,i]=1; weights[P+j,i]=r; weights[2P+j,i]=r*r
    end
    spec=Dict{String,Any}("type"=>"omrf","backend"=>"analytic","X"=>X,
        "seen"=>[3,3,3],"thresholds"=>"thresholds_0","interactions"=>"interactions_0",
        "weights"=>zeros(P,N*P),"node_weights"=>zeros(P,N*P),
        "norm_weights"=>weights,"bound_type"=>"node_local_norm_low_rank",
        "factorization"=>"person_node","factor_sampling"=>"size_biased",
        "prior_backend"=>"analytic","use_hcv"=>false)
    oracle=PDMPSamplersRBridge.OMRFAnalyticSubsamplingOracle(ctx,spec,names,anchor)
    legacy_env=PDMPSamplersRBridge._build_omrf_residual_envelope(spec,anchor,oracle.residual_context)
    compact_spec=copy(spec); delete!(compact_spec,"weights"); delete!(compact_spec,"node_weights")
    env=PDMPSamplersRBridge._build_omrf_residual_envelope(compact_spec,anchor,oracle.residual_context)
    env2=PDMPSamplersRBridge._build_omrf_residual_envelope(compact_spec,copy(anchor),oracle.residual_context)
    @test env.component_scales! isa PDMPSamplersRBridge.OMRFNormLocalScales
    @test env.component_scales!.threshold_coefficients !== env2.component_scales!.threshold_coefficients
    @test env.component_scales!.coordinate_displacements !==
        env2.component_scales!.coordinate_displacements
    @test env.component_scales!.coordinate_velocities !==
        env2.component_scales!.coordinate_velocities
    @test env.weights == legacy_env.weights
    @test env.component_blocks == legacy_env.component_blocks
    @test env.totals == legacy_env.totals
    @test env.alias_tables == legacy_env.alias_tables
    flow=PDMPSamplersRBridge.build_flow("AdaptiveBoomerang",
        Matrix(Diagonal(range(.6,1.4;length=9))),zeros(9);refresh_rate=.1)
    rng=Xoshiro(901)
    state=PDMPState(0.0,SkeletonPoint(.1randn(rng,9),PDMPSamplers.initialize_velocity(rng,flow,9)))
    PDMPSamplers.initialize_flow_state!(state,flow)
    scalar_scales=similar(env.scales)
    for t in (-0.7,0.0,0.23,1.4)
        for node in 1:P
            A,B,C,D=PDMPSamplersRBridge._omrf_norm_node_terms_at(
                env.component_scales!,state,flow,node,t)
            PDMPSamplersRBridge._fill_omrf_norm_scales!(
                scalar_scales,node,P,A,B,C,D)
        end
        PDMPSamplers.component_scales!(env.scales,env,state,flow,t)
        @test env.scales == scalar_scales
    end
    PDMPSamplers.component_cell_scales!(env.cell_scales,env,state,flow,0.0,.3)
    PDMPSamplers.component_scales!(env.scales,env,state,flow,.1)
    PDMPSamplers.component_cell_scales!(legacy_env.cell_scales,legacy_env,state,flow,0.0,.3)
    PDMPSamplers.component_scales!(legacy_env.scales,legacy_env,state,flow,.1)
    @test env.cell_scales == legacy_env.cell_scales
    @test env.scales == legacy_env.scales
    cv=SubsampledControlVariate((out,x)->fill!(out,0.0),oracle,env,anchor,1)
    legacy_cv=SubsampledControlVariate((out,x)->fill!(out,0.0),oracle,legacy_env,anchor,1)
    B=sum(env.scales .* env.totals); r1=Xoshiro(77); r2=Xoshiro(77)
    draws1=[(PDMPSamplers.draw_subset!(r1,cv,0.0,B);cv.subset[1]) for _=1:200]
    draws2=[(PDMPSamplers.draw_subset!(r2,legacy_cv,0.0,B);legacy_cv.subset[1]) for _=1:200]
    @test draws1 == draws2
    hcv_spec=copy(compact_spec); hcv_spec["bound_type"]="norm_local_hcv_compact"
    hcv_spec["use_hcv"]=true; hcv_spec["hcv_scope"]="global"
    hcv_spec["hcv_damping"]=Inf; hcv_spec["hcv_after_warmup"]=false
    hcv_weights=zeros(4P,N)
    for j=1:P,i=1:N
        r=sqrt(sum(abs2,X[i,k] for k=1:P if k!=j))
        for power=0:3; hcv_weights[power*P+j,i]=r^power; end
    end
    hcv_spec["norm_weights"]=hcv_weights
    hcv_oracle=PDMPSamplersRBridge.OMRFAnalyticSubsamplingOracle(ctx,hcv_spec,names,anchor)
    hcv_env=PDMPSamplersRBridge._build_omrf_residual_envelope(
        hcv_spec,anchor,hcv_oracle.residual_context)
    legacy_hcv=copy(hcv_spec);legacy_hcv["weights"]=zeros(P,N*P)
    legacy_hcv["node_weights"]=zeros(P,N*P)
    legacy_hcv_env=PDMPSamplersRBridge._build_omrf_residual_envelope(
        legacy_hcv,anchor,hcv_oracle.residual_context)
    @test size(hcv_env.weights)==(4P,N)
    @test hcv_env.weights==legacy_hcv_env.weights
    @test hcv_env.totals==legacy_hcv_env.totals
    @test hcv_env.alias_tables==legacy_hcv_env.alias_tables
    bpflow=PDMPSamplersRBridge.build_flow("BouncyParticle",
        Matrix(Diagonal(range(.6,1.4;length=9))),zeros(9);refresh_rate=.1)
    bpstate=PDMPState(0.0,SkeletonPoint(copy(state.ξ.x),copy(state.ξ.θ)))
    PDMPSamplers.initialize_flow_state!(bpstate,bpflow)
    PDMPSamplers.component_cell_scales!(hcv_env.cell_scales,hcv_env,bpstate,bpflow,0.0,.3)
    PDMPSamplers.component_scales!(hcv_env.scales,hcv_env,bpstate,bpflow,.1)
    PDMPSamplers.component_cell_scales!(legacy_hcv_env.cell_scales,legacy_hcv_env,bpstate,bpflow,0.0,.3)
    PDMPSamplers.component_scales!(legacy_hcv_env.scales,legacy_hcv_env,bpstate,bpflow,.1)
    @test hcv_env.cell_scales == legacy_hcv_env.cell_scales
    @test hcv_env.scales == legacy_hcv_env.scales
    hcv_cv=SubsampledControlVariate((out,x)->fill!(out,0.0),hcv_oracle,hcv_env,anchor,1)
    legacy_hcv_cv=SubsampledControlVariate((out,x)->fill!(out,0.0),hcv_oracle,legacy_hcv_env,anchor,1)
    hB=sum(hcv_env.scales .* hcv_env.totals); hr1=Xoshiro(79); hr2=Xoshiro(79)
    hdraws1=[(PDMPSamplers.draw_subset!(hr1,hcv_cv,0.0,hB);hcv_cv.subset[1]) for _=1:200]
    hdraws2=[(PDMPSamplers.draw_subset!(hr2,legacy_hcv_cv,0.0,hB);legacy_hcv_cv.subset[1]) for _=1:200]
    @test hdraws1 == hdraws2
    compact_sub=Dict{String,Any}("n_observations"=>N*P,
        "subsample_size"=>1,"residual_envelope"=>compact_spec)
    legacy_sub=Dict{String,Any}("n_observations"=>N*P,
        "subsample_size"=>1,"residual_envelope"=>spec)
    compact_model=PDMPSamplersRBridge._build_stan_subsampling_model(
        ctx,names,compact_sub,anchor,nothing,falses(9))[1]
    legacy_model=PDMPSamplersRBridge._build_stan_subsampling_model(
        ctx,names,legacy_sub,anchor,nothing,falses(9))[1]
    alg1=GridThinningStrategy(N=2,N_min=2,t_max=.25,lazy=false,
        bound_violation=:shrink);alg2=deepcopy(alg1)
    eflow1=deepcopy(flow); eflow2=deepcopy(flow)
    er1=Xoshiro(991);er2=Xoshiro(991)
    es1,em1,ea1,ec1,est1=PDMPSamplers.initialize_state(er1,eflow1,
        compact_model,alg1,0.0,SkeletonPoint(copy(state.ξ.x),copy(state.ξ.θ)))
    es2,em2,ea2,ec2,est2=PDMPSamplers.initialize_state(er2,eflow2,
        legacy_model,alg2,0.0,SkeletonPoint(copy(state.ξ.x),copy(state.ξ.θ)))
    @test es1.ξ.x == es2.ξ.x
    @test es1.ξ.θ ≈ es2.ξ.θ rtol=0 atol=eps(Float64)
    @test es1.t[] == es2.t[]
    for _=1:10
        out1=PDMPSamplers.next_event_time(er1,em1,eflow1,ea1,es1,ec1,est1,.25,false,:horizon_hit)
        out2=PDMPSamplers.next_event_time(er2,em2,eflow2,ea2,es2,ec2,est2,.25,false,:horizon_hit)
        @test out1[1:2]==out2[1:2]
        @test isequal(out1[3].∇ϕx,out2[3].∇ϕx)
        @test em1.grad.subset==em2.grad.subset
        PDMPSamplers.handle_event!(er1,out1[1],em1,eflow1,ea1,es1,ec1,
            out1[2],out1[3],est1)
        PDMPSamplers.handle_event!(er2,out2[1],em2,eflow2,ea2,es2,ec2,
            out2[2],out2[3],est2)
        @test es1.ξ.x ≈ es2.ξ.x rtol=0 atol=4eps(Float64)
        @test es1.ξ.θ ≈ es2.ξ.θ rtol=0 atol=4eps(Float64)
        @test es1.t[] == es2.t[]
    end
    @test @allocated(PDMPSamplers.component_cell_scales!(env.cell_scales,env,state,flow,0.0,.3)) == 0
    @test @allocated(PDMPSamplers.component_scales!(env.scales,env,state,flow,.1)) == 0
    # Closed cells spanning the anchor and harmonic extrema, plus inactive
    # interaction velocities, exercise the non-generic edge cases explicitly.
    for (left,right) in ((-pi/2,pi/2),(0.0,pi),(pi/2,3pi/2),(-2pi,2pi))
        PDMPSamplers.component_cell_scales!(env.cell_scales,env,state,flow,left,right)
        @test all(isfinite,env.cell_scales) && all(>=(0),env.cell_scales)
    end
    inactive_state=deepcopy(state); inactive_state.ξ.θ[7:9] .= 0
    PDMPSamplers.component_cell_scales!(env.cell_scales,env,inactive_state,flow,-.2,.4)
    @test all(isfinite,env.cell_scales) && all(>=(0),env.cell_scales)
    residual=zeros(9); max_violation=-Inf
    for iteration=1:20_001
        left=2pi*rand(rng)-pi; right=left + .001 + .7rand(rng)
        PDMPSamplers.component_cell_scales!(env.cell_scales,env,state,flow,left,right)
        copyto!(env.scales,env.cell_scales)
        factor=rand(rng,1:(N*P)); t=left+rand(rng)*(right-left)
        candidate=move_forward_time(state,t,flow)
        oracle(residual,candidate.ξ.x,[factor],anchor)
        exact=PDMPSamplers.λ(candidate,residual,flow)+PDMPSamplers.λ(candidate,-residual,flow)
        bound=PDMPSamplers.observation_residual_bound(env,factor)
        max_violation=max(max_violation,exact-bound)
        @test exact <= bound*(1+2e-7)+2e-8
    end
    println("maximum_certificate_violation=",max_violation)
    @test max_violation <= 2e-8
end

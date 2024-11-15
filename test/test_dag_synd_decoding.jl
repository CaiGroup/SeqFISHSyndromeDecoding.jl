
#include("simulation_fncs.jl")
#include("test_synd_decode_fncs.jl")

using SeqFISHSyndromeDecoding
using SeqFISHSyndromeDecoding: add_code_cols!, DotAdjacencyGraph, set_q, set_H


using DelimitedFiles
using Test
using DataFrames


q20_cb = readdlm("codes/Eng2019_647.csv", ',', UInt8)

cbs = [q20_cb]
pc_matrices = [[1 1 -1 -1;]]


@testset "all tests" for ii = 1:1

@testset "Test simulation set up" for (i, cb) in enumerate(cbs)
    ntargets = 10
    n = length(cb[1,:])
    q = length(unique(cb))
    #set_n(UInt8(n))
    set_q(UInt8(q))
    params = DecodeParams()
    set_zeros_probed(params, true)
    set_H(pc_matrices[i], params, cb)

    @test length(sim_true(ntargets).x) == ntargets
    @test length(encode(sim_true(ntargets), cb).x) == n*ntargets
    @test test_add_loc_errors_2(100, 0.1, cb)
    ndots = 100000
    rate = 0.5
    @test test_drop_random_dots(cb, ndots, rate) ≈ n*ndots*rate atol=10sqrt(ndots*rate*(1-rate))
    @inferred sim_true(100)
    @test length(draw_localization_errors(ntargets,0.1)[1]) == ntargets
    @test test_add_loc_errors(100, 0.1, cb)
    @test test_reconstruct_decode_message(1, cb)
    @test test_reconstruct_decode_message(2, cb)
    @test test_reconstruct_decode_message(100, cb)
end


@testset "Test DotAdjacencyGraph" for (i, cb) in enumerate(cbs)

    n = length(cb[1,:])
    q = length(unique(cb))
    #set_n(UInt8(n))
    set_q(UInt8(q))
    params = DecodeParams()
    set_zeros_probed(params, true)
    set_H(pc_matrices[i], params, cb)

    @test test_dag(300, cb, 0.05, 0.15, 0.1, 0)
    @test test_dag(300, cb, 0.05, 0.15, 0.1, 1)
    #@test test_dag_edges(300, cb)
    @test test_get_cw_pos_bnds(300, cb, 0.1, 0.3,1,ndrops=0)
    @test test_get_cw_pos_bnds(300, cb, 0.1, 0.3,1,ndrops=1)
end



@testset "test compute syndrome" begin
    for (i, cb) in enumerate(cbs)
        #n = length(cb[1,:])
        q = length(unique(cb))
        #set_n(UInt8(n))
        set_q(UInt8(q))
        params = DecodeParams()
        set_zeros_probed(params, true)
        set_H(pc_matrices[i], params, cb)
        ntargets = 50000
        @test test_get_cw_pos_bnds(ntargets, cb, 0.0, 0.0, 0.0)

        test_compute_syndromes(ntargets, cb, 0)
        #test_compute_syndromes(ntargets, cb, 1)

    end
end

println("full decode perfect")
@testset "full decode perfect" begin
    for (i, cb) in enumerate(cbs), ntargets in [1, 10, 100]
        params_ = DecodeParams()
        set_zeros_probed(params_, true)
        set_H(pc_matrices[i], params_, cb)
        H = pc_matrices[i]

        lat_thresh = 0.0
        z_thresh = 0.0
        ndrops = 0
        free_dot_energy = 5.0
        mip_sa_thresh = 80

        pnts, g = construct_test_dag(ntargets, 0, 0, 0, cb, ndrops)

        pnts.z = zeros(nrow(pnts))

        ndots = ntargets*length(cb[1,:])
        free_dot_cost = 5.0

        c_final = 2
        n_chains = 100
        l_chain = 300
        c₀ = free_dot_energy * 40
        lat_var_factor = 8000.0
        z_var_factor = 1.0
        lw_var_factor = 0.0
        s_var_factor = 0.0
        erasure_penalty = 4.0
        converge_thresh = 100 * ndots
        skip_thresh = 2000
        skip_density_thresh = 2000
        params = DecodeParams()
        set_xy_search_radius(params, lat_thresh)
        set_z_search_radius(params, z_thresh)
        set_lat_var_cost_coeff(params, lat_var_factor)
        set_z_var_cost_coeff(params, z_var_factor)
        set_lw_var_cost_coeff(params, lw_var_factor)
        set_s_var_cost_coeff(params, s_var_factor)
        set_n_allowed_drops(params, ndrops)
        set_zeros_probed(params, true)
        set_free_dot_cost(params, free_dot_energy)
        set_erasure_penalty(params, erasure_penalty)
        set_skip_thresh(params, skip_thresh)
        set_skip_density_thresh(params, skip_density_thresh)

        decode_syndromes!(pnts, cb, H, params)
        @test pnts.species == [Int(p) for p in pnts.decoded]
    end
end


println("full decode perfect search radius")
@testset "full decode perfect search radius" begin
    for (i, cb) in enumerate(cbs[1:end-1]), ntargets in [1, 10, 30, 100]
        n = length(cb[1,:])
        q = length(unique(cb))
        #set_n(UInt8(n))
        set_q(UInt8(q))
        
        H = pc_matrices[i]

        lat_thresh = 0.1
        z_thresh = 0.0
        ndrops = 0
        free_dot_energy = 1.0
        mip_sa_thresh = 800000000000

        

        ndots = ntargets*sum(cb[1,:] .!= 0)
        free_dot_cost = 1.0

        c_final = 2
        n_chains = 100
        l_chain = 300
        c₀ = free_dot_energy * 40
        lat_var_factor = 8000.0
        z_var_factor = 1.0
        lw_var_factor = 0.0
        s_var_factor = 0.0
        erasure_penalty = 0.0
        converge_thresh = 100 * ndots
        skip_thresh = 20000000000000
        skip_density_thresh = 20000000000009
        params = DecodeParams()
        set_lat_var_cost_coeff(params, lat_var_factor)
        set_xy_search_radius(params, lat_thresh)
        set_free_dot_cost(params, free_dot_cost)
        set_z_var_cost_coeff(params, z_var_factor)
        set_s_var_cost_coeff(params, s_var_factor)
        set_lw_var_cost_coeff(params, lw_var_factor)
        set_H(pc_matrices[i], params, cb)
        set_zeros_probed(params, true)
        set_skip_thresh(params, skip_thresh)
        set_skip_density_thresh(params, skip_density_thresh)
        set_H(pc_matrices[i], params, cb)
        pnts, g = construct_test_dag(ntargets, 0, 0, 0, cb, ndrops)

        pnts.z = zeros(nrow(pnts))
        decode_syndromes!(pnts, cb, H, params)
        @test pnts.species == [Int(p) for p in pnts.decoded]
    end
end

println("full decode drops")
@testset "full decode drops" begin
    for (i, cb) in enumerate(cbs), ntargets in [1, 2, 10, 100]
        H = pc_matrices[i]
        ndrops = 1

        encoded = construct_test_encoding(ntargets, cb)
        lat_thresh = 0.0
        z_thresh = 0
        n = length(cb[1,:])

        
        lat_thresh = 0.0
        z_thresh = 0.0
        free_dot_energy = 5.0
        mip_sa_thresh = 80

        ndots = ntargets*length(cb[1,:])
        free_dot_cost = 5.0

        c_final = 2
        n_chains = 100
        l_chain = 300
        c₀ = free_dot_energy * 40
        lat_var_factor = 8000.0
        z_var_factor = 1.0
        lw_var_factor = 0.0
        s_var_factor = 0.0
        erasure_penalty = 4.0
        converge_thresh = 100 * ndots
        skip_thresh = 2000
        skip_density_thresh=2000
        params = DecodeParams()
        set_lat_var_cost_coeff(params, lat_var_factor)
        set_xy_search_radius(params, lat_thresh)
        set_free_dot_cost(params, free_dot_cost)
        set_z_var_cost_coeff(params, z_var_factor)
        set_s_var_cost_coeff(params, s_var_factor)
        set_lw_var_cost_coeff(params, lw_var_factor)
        set_H(pc_matrices[i], params, cb)
        set_n_allowed_drops(params,1)
        set_zeros_probed(params, true)
        set_skip_thresh(params, skip_thresh)
        set_skip_density_thresh(params, skip_density_thresh)

        decode_syndromes!(encoded, cb, H, params)
        @test encoded.species == encoded.decoded
        

        select!(encoded, Not(["decoded", "mpath"]))
        sort!(encoded, :dot_ID)
        to_drop = rand(1:n , ntargets) + n*Array(0:(ntargets-1))
        deleteat!(encoded, to_drop)
        decode_syndromes!(encoded, cb, H, params)
        @test encoded.species == encoded.decoded
    end
end

end

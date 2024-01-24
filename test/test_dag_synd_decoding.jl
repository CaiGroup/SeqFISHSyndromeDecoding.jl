
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
    @test test_get_cw_pos_bnds(300, cb, 0.1, 0.3,1,0)
    @test test_get_cw_pos_bnds(300, cb, 0.1, 0.3,1,1)
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
        params = DecodeParams(
            lat_thresh,
            z_thresh,
            lat_var_factor,
            z_var_factor,
            lw_var_factor,
            s_var_factor,
            ndrops,
            true,
            mip_sa_thresh,
            free_dot_energy,
            n_chains,
            l_chain,
            c₀,
            (c₀/c_final-1)/log(n_chains),
            erasure_penalty,
            converge_thresh,
            skip_thresh,
            skip_density_thresh
        )
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
        params = DecodeParams(
            lat_thresh,
            z_thresh,
            lat_var_factor,
            z_var_factor,
            lw_var_factor,
            s_var_factor,
            ndrops,
            true,
            mip_sa_thresh,
            free_dot_energy,
            n_chains,
            l_chain,
            c₀,
            (c₀/c_final-1)/log(n_chains),
            erasure_penalty,
            converge_thresh,
            skip_thresh,
            skip_density_thresh
        )
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

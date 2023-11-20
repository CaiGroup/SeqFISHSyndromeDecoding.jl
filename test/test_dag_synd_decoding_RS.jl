
include("simulation_fncs.jl")

using SeqFISHSyndromeDecoding
include("test_synd_decode_fncs.jl")
using SeqFISHSyndromeDecoding: add_code_cols!, syndrome_find_barcodes!, DotAdjacencyGraph, set_q, set_H


using DelimitedFiles
using Test
using DataFrames
using CSV


#RS_q5_k2_cb = readdlm("RS_q5_k2_cb.csv", ',', UInt8)
#RS_q5_k2_H = readdlm("RS_q5_k2_H.csv", ',', UInt8)
RS_q7_k4_w4cb = readdlm("codes/RS_q7_k4_w4cb.csv", ',', UInt8)
RS_q7_k4_H = readdlm("codes/RS_q7_k4_H.csv", ',', UInt8)
RS_q11_k7_w4cb = readdlm("codes/RS_q11_k7_w4cb.csv", ',', UInt8)
RS_q11_k7_w5cb = readdlm("codes/RS_q11_k7_w5cb.csv", ',', UInt8)
RS_q11_k7_H = readdlm("codes/RS_q11_k7_H.csv", ',', UInt8)
RS_q11_k8_w4cb = readdlm("codes/RS_q11_k8_w4cb.csv", ',', UInt8)
RS_q11_k8_H = readdlm("codes/RS_q11_k8_H.csv", ',', UInt8)
hamming_merfish_cb = readdlm("codes/hamming_merfish_cb.csv", ',', UInt8)
hamming_merfish_H = readdlm("codes/hamming_merfish_H.csv", ',', UInt8)
RS_q8_n7_k4_H = readdlm("codes/RS_q8_n7_k4_H.csv", ',', String)
RS_q8_n7_k4_w4cb = readdlm("codes/RS_q8_n7_k4_w4cb.csv", ',', String)
RS_q9_n8_k5_H = readdlm("codes/RS_q9_n8_k5_H.csv", ',', String)
RS_q9_n8_k5_w4cb = readdlm("codes/RS_q9_n8_k5_w4cb.csv", ',', String)

RS_q11_n12_k9_H = readdlm("codes/RS_q11_n12_k9_H.csv", ',', UInt8)
RS_q11_n12_k9_w4cb = Matrix(select(DataFrame(CSV.File("codes/RS_q11_n12_k9_w4cb.csv")), Not(:gene_name)))

q9_n10_k7_w4_ncws1000_v1 = String.(Matrix(select(DataFrame(CSV.File("codes/q9_n10_k7_w4_ncws1000_v1.csv")), Not(:gene_name))))
RS_q9_n10_k7_H = readdlm("codes/RS_q9_n10_k7_H.csv", ',', String)

RS_q13_n12_k8_H = readdlm("codes/RS_q13_n12_k8_H.csv", ',', UInt8)
RS_q13_n12_k8_w5cb = Matrix(select(DataFrame(CSV.File("codes/RS_q13_n12_k8_w5cb.csv")), Not(:gene_name)))

#cbs = [RS_q7_k4_w4cb, RS_q11_k7_w4cb, RS_q11_k7_w5cb, RS_q11_k8_w4cb, hamming_merfish_cb, RS_q8_n7_k4_w4cb, RS_q9_n8_k5_w4cb, RS_q11_n12_k9_w4cb]
#pc_matrices = [RS_q7_k4_H, RS_q11_k7_H, RS_q11_k7_H, RS_q11_k8_H, hamming_merfish_H, RS_q8_n7_k4_H, RS_q9_n8_k5_H, RS_q11_n12_k9_H]

cbs = [RS_q7_k4_w4cb, RS_q11_k7_w4cb, RS_q11_k7_w5cb, RS_q11_k8_w4cb, RS_q8_n7_k4_w4cb, RS_q9_n8_k5_w4cb, RS_q11_n12_k9_w4cb, q9_n10_k7_w4_ncws1000_v1, RS_q13_n12_k8_w5cb]
pc_matrices = [RS_q7_k4_H, RS_q11_k7_H, RS_q11_k7_H, RS_q11_k8_H, RS_q8_n7_k4_H, RS_q9_n8_k5_H, RS_q11_n12_k9_H, RS_q9_n10_k7_H, RS_q13_n12_k8_H]

#cbs = [RS_q8_n7_k4_w4cb]
#pc_matrices = [RS_q8_n7_k4_H]

@testset "all tests" for ii = 1:1


@testset "Test simulation set up RS" for (i, cb) in enumerate(cbs)
    ntargets = 10
    n = length(cb[1,:])
    if typeof(cb[1,1]) == String
        w = sum(cb[1,:] .!= "0")
    else
        w = sum(cb[1,:] .!= 0)
    end
    q = length(unique(cb))
    len_cb = length(cb[:, 1])
    #set_n(UInt8(n))
    set_q(UInt8(q))
    params = DecodeParams()
    set_zeros_probed(params, false)
    set_H(pc_matrices[i], params, cb)

    @test length(sim_true(ntargets, len_cb).x) == ntargets
    @test length(encode(sim_true(ntargets, len_cb), cb).x) == w*ntargets
    @test test_add_loc_errors_2(100, 0.1, cb)
    ndots = 100000
    rate = 0.5
    @test test_drop_random_dots(cb, ndots, rate) ≈ w*ndots*rate atol=10sqrt(ndots*rate*(1-rate))
    @inferred sim_true(100, len_cb)
    @test length(draw_localization_errors(ntargets,0.1)[1]) == ntargets
    @test test_add_loc_errors(100, 0.1, cb)
    @test test_reconstruct_decode_message(1, cb)
    @test test_reconstruct_decode_message(2, cb)
    @test test_reconstruct_decode_message(100, cb)
end


@testset "Test DotAdjacencyGraph RS" for (i, cb) in enumerate(cbs)

    n = length(cb[1,:])
    q = length(unique(cb))
    #set_n(UInt8(n))
    set_q(UInt8(q))
    params = DecodeParams()
    set_zeros_probed(params, false)
    set_H(pc_matrices[i], params, cb)

    @test test_dag(300, cb, 0.05, 0.15, 0.1, 0)
    #@test test_dag(300, cb, 0.05, 0.15, 0.1, 1)
    #@test test_dag_edges(300, cb)
end


println("full decode perfect RS")
@testset "full decode perfect RS" begin
    for (i, cb) in enumerate(cbs), ntargets in [1, 10, 100, 1000]
        n = length(cb[1,:])
        q = length(unique(cb))
        #set_n(UInt8(n))
        set_q(UInt8(q))
        
        H = pc_matrices[i]

    
        lat_thresh = eps() #0.0
        z_thresh = eps() #0.0
        ndrops = 0
        free_dot_energy = 1.0
        mip_sa_thresh = 80        

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
        erasure_penalty = 1.0
        converge_thresh = 100 * ndots
        skip_thresh = 200000
        skip_density_thresh = 2000000
        params = DecodeParams()
        set_lat_var_cost_coeff(params, lat_var_factor)
        set_xy_search_radius(params, lat_thresh)
        set_free_dot_cost(params, free_dot_cost)
        set_z_var_cost_coeff(params, z_var_factor)
        set_s_var_cost_coeff(params, s_var_factor)
        set_lw_var_cost_coeff(params, lw_var_factor)
        set_H(pc_matrices[i], params, cb)
        set_zeros_probed(params, false)
        
        pnts, g = construct_test_dag(ntargets, 0, 0, 0, cb, ndrops)
        pnts.z = zeros(nrow(pnts))

        decode_syndromes!(pnts, cb, H, params)
        @test pnts.species == [Int(p) for p in pnts.decoded]
    end
end


println("full decode perfect RS search radius")
@testset "full decode perfect RS search radius" begin
    for (i, cb) in enumerate(cbs[1:end-1]), ntargets in [1, 10, 30] #, 100]
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
        set_zeros_probed(params, false)
        set_skip_thresh(params, skip_thresh)
        set_skip_density_thresh(params, skip_density_thresh)
        set_H(pc_matrices[i], params, cb)
        pnts, g = construct_test_dag(ntargets, 0, 0, 0, cb, ndrops)

        pnts.z = zeros(nrow(pnts))
        decode_syndromes!(pnts, cb, H, params)
        @test pnts.species == [Int(p) for p in pnts.decoded]
    end
end

println("full decode 1 drop RS")
@testset "full decode 1 drop RS" begin
    for (i, cb) in enumerate(cbs), ntargets in [1, 2, 10, 100]
        H = pc_matrices[i]
        if size(H)[1] > 1
            ndrops = 1

            encoded = construct_test_encoding(ntargets, cb)
            lat_thresh = 0.0
            z_thresh = 0
            if typeof(cb[1,1]) == String
                w = maximum(sum(cb .!= "0", dims=2))
            else
                w = maximum(sum(.~ iszero.(cb), dims=2))
            end

            to_drop = rand(1:w , ntargets) + w*Array(0:(ntargets-1))
            deleteat!(encoded, to_drop)
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
            skip_thresh = 200000000
            skip_density_thresh=200000000000
            params = DecodeParams(
                lat_thresh,
                z_thresh,
                lat_var_factor,
                z_var_factor,
                lw_var_factor,
                s_var_factor,
                ndrops,
                false,
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
            set_zeros_probed(params, false)
            set_skip_thresh(params, skip_thresh)
            decode_syndromes!(encoded, cb, H, params)
            @test encoded.species == encoded.decoded
        end
    end
end

println("full decode 2 drop RS")
@testset "full decode 2 drop RS" begin
    for (i, cb) in enumerate(cbs), ntargets in [1, 2, 10, 100]
        if size(H)[1] > 3
            ndrops = 2

            encoded = construct_test_encoding(ntargets, cb)
            lat_thresh = 0.0
            z_thresh = 0
            if typeof(cb[1,1]) == String
                w = maximum(sum(cb .!= "0", dims=2))
            else
                w = maximum(sum(.~ iszero.(cb), dims=2))
            end

            to_drop = rand(1:w , ntargets) + w*Array(0:(ntargets-1))
            deleteat!(encoded, to_drop)
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
            skip_thresh = 200000000
            skip_density_thresh=200000000000
            params = DecodeParams(
                lat_thresh,
                z_thresh,
                lat_var_factor,
                z_var_factor,
                lw_var_factor,
                s_var_factor,
                ndrops,
                false,
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
            set_zeros_probed(params, false)
            set_skip_thresh(params, skip_thresh)
            decode_syndromes!(encoded, cb, H, params)
            @test encoded.species == encoded.decoded
        end
    end
end



println("full decode drop RS search radius")
@testset "full decode drop RS search radius" begin
    for (i, cb) in enumerate(cbs), ntargets in [1, 10, 30] #, 100]
        H = pc_matrices[i]
        if size(H)[1] > 1
            n = length(cb[1,:])
            q = length(unique(cb))
            #set_n(UInt8(n))
            set_q(UInt8(q))
            
            lat_thresh = 0.1
            z_thresh = 0.0
            ndrops = 1
            free_dot_energy = 1.0
            mip_sa_thresh = 800000000000
            if typeof(cb[1,1]) == String
                w = maximum(sum(cb .!= "0", dims=2))
            else
                w = maximum(sum(.~ iszero.(cb), dims=2))
            end

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
            erasure_penalty = 1.0
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
            set_n_allowed_drops(params, ndrops)
            set_zeros_probed(params, false)
            set_skip_thresh(params, skip_thresh)
            set_skip_density_thresh(params, skip_density_thresh)
            
            set_H(pc_matrices[i], params, cb)

            #pnts, g = construct_test_dag(ntargets, 0, 0, 0, cb, ndrops)
            encoded = construct_test_encoding(ntargets, cb)

            to_drop = rand(1:w , ntargets) + w*Array(0:(ntargets-1))
            deleteat!(encoded, to_drop)        

            encoded.z = zeros(nrow(encoded))
            decode_syndromes!(encoded, cb, H, params)
            @test encoded.species == [Int(p) for p in encoded.decoded]
        end
    end
end

"""
println("full decode perfect RS barcode pairs w same coords")
@testset "full decode perfect RS barcode pairs w same coords" begin
    for ntargets in [1, 10, 100, 1000]
        n = length(RS_q11_k7_w4cb[1,:])
        q = length(unique(RS_q11_k7_w4cb))
        set_n(UInt8(n))
        set_q(UInt8(q))
        set_H(RS_q11_k7_H)
        H = RS_q11_k7_H

        lat_thresh = 0.0
        z_thresh = 0.0
        ndrops = 0
        free_dot_energy = 1.0
        mip_sa_thresh = 800000000

        pnts1, g = construct_test_dag(ntargets, 0, 0, 0, RS_q11_k7_w4cb, ndrops)
        pnts2, g2 = construct_test_dag(ntargets, 0, 0, 0, RS_q11_k7_w4cb, ndrops)
        pnts1.z = zeros(nrow(pnts1))
        pnts2.z = zeros(nrow(pnts2))

        sort!(pnts1, :dot_ID)
        sort!(pnts2, :dot_ID)

        pnts2.x .= pnts1.x
        pnts2.y .= pnts1.y
        pnts = vcat(pnts1, pnts2)
        #println("pnts:")
        #println(pnts)
        sort!(pnts, [:round, :pseudocolor])

        ndots = ntargets*sum(RS_q11_k7_w4cb[1,:] .!= 0)
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
            false,
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
        decode_syndromes!(pnts, RS_q11_k7_w4cb, H, params)
        println(sum(sort(pnts.decoded) .!= sort(pnts.decoded)))
        #println(hcat(dcd, sort(pnts.species))')
        sort!(pnts, [:x, :species, :decoded])
        println(pnts[!, [:x, :species, :decoded]])
        @test(all(sort(pnts.species) .== sort(pnts.decoded)))
       #@test all(sort(pnts.species) .== dcd)
    end
end
"""

end
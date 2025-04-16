using DataFrames
using CSV
using SeqFISHSyndromeDecoding
using GLPK
using Test
using DelimitedFiles

@testset "Test Real Reed-Solomon Encoded Data" begin
    cb = DataFrame(CSV.File("../example_data/full_RS_q11_k7_half_pool_cb.csv"))
    H = readdlm("../example_data/RS_q11_k7_H.csv", ',', UInt8)
    pnts = DataFrame(CSV.File("../example_data/example_RS_cell_points.csv"))

    filter!(pnt -> ~ismissing(pnt.pseudocolor), pnts)
    pnts.block = UInt8.(pnts.block)
    select!(pnts, Not([:ch,:hyb]))
    SeqFISHSyndromeDecoding.sort_readouts!(pnts)

    params = DecodeParams()
    set_zeros_probed(params, false)
    set_lat_var_cost_coeff(params, 7.0)
    set_z_var_cost_coeff(params, 0.0)
    set_lw_var_cost_coeff(params, 0.0)
    set_s_var_cost_coeff(params, 0.0)
    set_free_dot_cost(params, 1.0)
    set_n_allowed_drops(params, 0)
    
    set_xy_search_radius(params, 2)
    set_z_search_radius(params, 0.0);


    barcodes = decode_syndromes!(pnts, cb, H, params) #, GLPK.Optimizer)

    saved_results = DataFrame(CSV.File("../example_data/example_RS_results.csv"))
    sort!(saved_results, :gene_number)
    sort!(barcodes, :gene_number)
    @test saved_results[:,[:gene, :gene_number]] == barcodes[:,[:gene, :gene_number]]

end
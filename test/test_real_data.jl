using DataFrames
using CSV
using SeqFISHSyndromeDecoding
using GLPK


@testset "Test Real Data" begin
    cb = DataFrame(CSV.File("../example_data/codebook_ch_561.csv"))
    H = [1 1 -1 -1;]
    pnts = DataFrame(CSV.File("../example_data/example_cell_points.csv"))

    pnts.z = zeros(Float64, nrow(pnts))
    pnts.hyb = UInt8.(pnts.hyb)


    params = DecodeParams()

    set_lat_var_cost_coeff(params, 30.0)
    set_z_var_cost_coeff(params, 0.0)
    set_lw_var_cost_coeff(params, 8.0)
    set_s_var_cost_coeff(params, 0.0)
    set_free_dot_cost(params, 5.0)

    set_xy_search_radius(params, sqrt(5.0*size(H)[2]/30.0)*3)
    set_z_search_radius(params, 0.0)
    set_zeros_probed(params, true)


    barcodes = decode_syndromes!(pnts, cb, H, params) #, GLPK.Optimizer)

    saved_results = DataFrame(CSV.File("../example_data/example_results.csv"))
    sort!(barcodes, :gene_number)
    @test saved_results[:,[:gene, :gene_number]] == barcodes[:,[:gene, :gene_number]]

end
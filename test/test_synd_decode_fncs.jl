using SeqFISHSyndromeDecoding
using SeqFISHSyndromeDecoding: DotAdjacencyGraph, DotAdjacencyGraphBlankBlock,
make_cw_dict, get_cw_pos_inds, add_code_cols!, compute_syndromes, neighbors, get_number_of_dots

using DelimitedFiles
using Test
using Graphs: nv, ne
using DataFrames

#q20_cb = readdlm("Eng2019_647.csv", ',', UInt8)


function construct_test_encoding(n, cb)
    true_locs = sim_true(n, length(cb[:, 1]))
    encode(true_locs, cb)
end

"""
    test_drop_random_dots(cb, ntargets, drop_rate)ƒ
"""
function test_drop_random_dots(cb, ntargets, drop_rate)
    encoded = construct_test_encoding(ntargets, cb)
    n_initial = nrow(encoded)
    drop_random_dots!(encoded, drop_rate)
    n_dropped = n_initial - length(encoded.x)
    return n_dropped
end


function test_add_loc_errors(n,rstdv, cb)
    pnts = construct_test_encoding(n, cb)
    ndots = length(pnts.dot_ID)
    x_errors, y_errors = draw_localization_errors(ndots,0.1)
    pnts_1 = deepcopy(pnts)
    add_localization_errors!(pnts, x_errors, y_errors)
    pnts.x - pnts_1.x ≈ x_errors && pnts.y - pnts_1.y ≈ y_errors #approx to within floating point error
end


function test_add_loc_errors_2(n,rstdv,cb)
    pnts = construct_test_encoding(n, cb)
    pnts_1 = deepcopy(pnts)
    add_localization_errors!(pnts, rstdv)
    pnts.x ≠ pnts_1.x && pnts.y ≠ pnts_1.y
end


function test_dag_edges(ntargets, cb, ndrops, tforms=nothing)
    n, q, w = get_n_q_w(cb)
    pnts, g = construct_test_dag(ntargets, 0.0, 0.0, 0.0, cb, ndrops, tforms=tforms)
    sum_edges = sum([length(neighbors(g, d)) for d in 1:nrow(pnts)])
    #ne(g.g) == length(pnts.dot_ID)*(factorial(n)/(factorial(n-2)*2))/n
    sum_edges == length(pnts.dot_ID)*(factorial(n)/(factorial(n-2)*2))/n
end


"""
    construct_message(pnts :: DataFrame, mdot_inds, cb)

Converts a list of dots into a message vector.
"""
function construct_message(pnts :: DataFrame, mdot_inds, cb :: Array)
    n, q, w = get_n_q_w(cb)
    if typeof(cb[1,1]) == String
        mcw = fill("0", n)
    else
        mcw = fill(0x00, n)
    end
    for m_dot in mdot_inds
        mcw[pnts.pos[m_dot]] = pnts.coeff[m_dot]
    end
    return mcw
end

"""
function construct_message(pnts :: DataFrame, mdot_inds, n, q, w, symb_type :: DataType)
    if w == n
        mcw = fill(symb_type(q+1), n)
    else
        mcw = fill(symb_type(0), n)
    end
    #mdot_IDs = pnts.dot_ID[mdot_inds]
    for m_dot in mdot_inds
        mcw[pnts.pos[m_dot]] = pnts.coeff[m_dot]
    end
    
    mcw#, mdot_IDs
end
"""

function test_reconstruct_decode_message(ntargets, cb)
    true_locs = sim_true(ntargets, length(cb[:, 1]))
    encoded = encode(true_locs, cb)
    n, q, w = get_n_q_w(cb)
    cw_dict = make_cw_dict(cb)
    symb_type = typeof(cb[1,1])
    params = SeqFISHSyndromeDecoding.DecodeParams()
    set_zeros_probed(params, n == w)
    add_code_cols!(encoded, params)
    for i = 1:ntargets
        stop = w*i
        start = stop - w + 1
        mdots_inds = Array(start:stop)

        mcw = construct_message(encoded, mdots_inds, cb)
        for mdot in mdots_inds
            if mcw ∉ keys(cw_dict) || encoded.species[mdot] != cw_dict[mcw]
                println("mcw: ", mcw)
                mcw ∈ keys(cw_dict) ? println("species: ", cw_dict[mcw]) : println("mcw ∉ dict")
                println(encoded)
                println("mdot: ", mdot)
                return false
            end
        end
    end
    true
end

function encoded_2_dag!(pnts, cb, lat_thresh, z_thresh, ndrops, tforms=nothing)
    n, q, w = get_n_q_w(cb)
    params = SeqFISHSyndromeDecoding.DecodeParams()
    set_zeros_probed(params, n == w)

    if ("block" in names(pnts)) && ("pseudocolor" in names(pnts))
        sort!(pnts, [:block, :pseudocolor])
    else
        sort!(pnts, :hyb)
    end
    add_code_cols!(pnts, params)
    if isnothing(tforms)
        tforms_dict = nothing
    else
        tforms_dict = SeqFISHSyndromeDecoding.get_tform_dict(tforms)
    end
    if n > w
        DotAdjacencyGraphBlankBlock(pnts, lat_thresh, z_thresh, n, ndrops, w, tforms_dict)
    else
        return DotAdjacencyGraph(pnts, lat_thresh, z_thresh, n, ndrops, tforms=tforms_dict)
    end
end

function test_dag(ntargets, cb, rstdv, lat_thresh, z_thresh, ndrops, tforms=nothing)
    pnts, g = construct_test_dag(ntargets, rstdv, lat_thresh, z_thresh, cb, ndrops, tforms=tforms)
    nv(g.g) == length(pnts.dot_ID)
end

function construct_test_dag(ntargets, rstdv, lat_thresh, z_thresh, cb, ndrops; tforms=nothing)
    pnts = construct_test_encoding(ntargets, cb)
    n, q, w = get_n_q_w(cb)
    add_localization_errors!(pnts, rstdv)
    g = encoded_2_dag!(pnts, cb, lat_thresh, z_thresh, ndrops, tforms)
    [pnts, g]
end

function test_get_cw_pos_bnds(ntargets, cb, rstdv, lat_thresh, z_thresh; ndrops=0, tforms=nothing)
    pnts, g = construct_test_dag(ntargets, rstdv, lat_thresh, z_thresh, cb, ndrops, tforms=tforms)
    n = maximum(pnts.pos)
    for posᵢ in 1:n
        if pnts.pos[get_cw_pos_inds(g, posᵢ)] != fill(posᵢ, length(get_cw_pos_inds(g, posᵢ)))
            return false
        end
    end
    true
end

function test_compute_syndromes(ntargets :: Int64, cb, ndrops, tforms=nothing)
    n, q, w = get_n_q_w(cb)
    rstdv = 0.0
    lat_thresh = 0.0
    z_thresh = 0.0
    pnts, g = construct_test_dag(ntargets, rstdv, lat_thresh, z_thresh, cb, ndrops, tforms=tforms)
    syndromes, syndrome_coeff_positions = compute_syndromes(pnts, g)
    final_pos_dots = get_cw_pos_inds(g, n)
    fpd_sndrs = syndromes[final_pos_dots]
    ncws, npseudocolors = size(cb)
    fpd_pth_lns = [get_number_of_dots(p[1], npseudocolors) for p in syndrome_coeff_positions[final_pos_dots]]

    full_path_sums = []
    for i = eachindex(fpd_sndrs)
        full_len_paths = [fpd_sndrs[i][j] for j in 1:length(fpd_sndrs[i]) if fpd_pth_lns[i][j] == UInt8(n)]
        push!(full_path_sums, full_len_paths)
    end
    @test all(iszero.(full_path_sums))
end

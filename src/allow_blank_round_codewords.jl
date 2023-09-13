



function find_blank_round_codewords(
    pnts ::DataFrame,
    g :: DotAdjacencyGraph,
    cw_dict,
    w
)

    # initialize array for syndrome partial sums
    syndromes = fill(Vector{Vector{SyndromeComponent}}(), nrow(pnts)) #Vector{Vector{SyndromeComponent}}()
    syndrome_block_bounds = fill(Vector{Vector{Int64}}(), nrow(pnts))
    len_syndromes = fill(Vector{Int64}(), nrow(pnts))

    # get list of dots in rounds that could be the last of a codeword
    potential_barcode_final_dots = collect(g.cw_pos_bnds[w]:g.n)
    potential_barcode_final_dots_n_uncomputed_neighbors = fill(0, length(potential_barcode_final_dots))
    unprocessed_inrange_dots = fill(0, nrow(pnts))
    round_trees = [make_KDTree(pnts[get_cw_pos_inds(g, r), :]) for r in 1:g.n]

    # for each dot, count how many dots that could be the final dot of a barcode including it are in range
    # and how many dots that could be the final dot are in range of each other
    for dot in 1:nrow(pnts)
        r = pnts.pos[dot] 
        if r >= w
            unprocessed_inrange_dots[dot] += 1
            for ri in w:r
                ri_inrange = length(inrange(round_trees[ri], [pnts.x[dot], pnts.y[dot]], g.lat_thresh))
                potential_barcode_final_dots_n_uncomputed_neighbors[dot-g.cw_pos_bnds[w]-1] += ri_inrange
            end
        end

        ri = maximum([w, r+1])
        while ri < g.n
            ri_inrange = length(inrange(round_trees[ri], [pnts.x[dot], pnts.y[dot]], g.lat_thresh))
            unprocessed_inrange_dots[dot] += ri_inrange
            if r >= w
                potential_barcode_final_dots_n_uncomputed_neighbors[dot-g.cw_pos_bnds[w]-1] += ri_inrange
            end
            ri += 1
        end
    end


    while length(potential_barcode_final_dots) > 0
        dot_ind = argmin(potential_barcode_final_dots_n_uncomputed_neighbors[potential_barcode_final_dots .- (g.cw_pos_bnds[w]-1)])
        dot = potential_barcode_final_dots[dot_ind]    

    end

end



"""
    find_barcode_candidates()
"""
function find_barcode_candidates!(
    g :: DotAdjacencyGraph,
    pnts :: DataFrame,
    cw_dict :: Dict,
    w,
    syndromes,
    len_syndromes,
    unprocessed_inrange_dots,
    syndrome_block_bounds,
    dot
    )
    # compute dot syndromes
    # if allocated 
    if length(syndromes[dot]) > 0
        return
    # elif not allocated
    else
        # if all in range final round dots processed
        if unprocessed_inrange_dots[dot] == 0
            return 
        #elseif dot in first round
        
        else
            # get round
            r = pnts.pos[dot]
            #init nsynd counter
            nsyndc = fill(0, minimum([r,w-1]))
            if r < g.n - 3
                nsyndc[1] = 1
            end
            dot_synd_cs = Vector{{Vector{SyndromeComponent}}}()
            sizehint!(dot_synd_cs, length(nsyndc))
            dot_neighbors = neighbors(g, dot)

            start_ind = maximum([1, r - g.n + w])
            for dot in dot_neighbors
                find_barcode_candidates!(pnts, g, cw_dict, w, syndromes, len_syndromes, unprocessed_inrange_dots, syndrome_block_bounds, dot)
                for (i, lsa) in enumerate(len_syndromes[dot][start_ind:end])
                    nsyndc[i+start_ind] += lsa
                end
            end
            
            dot_synd_cs = fill.(pnts.sc[dot], nsyndc)
            @inbounds len_syndromes[dot] = nsyndc[1:(end-1)]
            
            # allocate syndome component array of length nsnd
            #@inbounds syndromes[dot] = fill(pnts.sc[dot], nsyndc)
            
            # for dot in neighbors
                
             

            nbr_block_bnds = fill(1, length(dot_neighbors)+1)

            # get neighbors of that dot
            for (i, neighbor) in enumerate(dot_neighbors)

                # add syndrome components to appropriate block
                @inbounds nbr_block_bnds[i+1] = nbr_block_bnds[i]+length(syndromes[neighbor])-1
                @inbounds syndromes[dot][nbr_block_bnds[i]:nbr_block_bnds[i+1]] += syndromes[neighbor]

                # keep track of indices of each neighbor's syndrome component block
                @inbounds nbr_block_bnds[i+1] += 1
            end
            @inbounds syndrome_block_bounds[dot] = nbr_block_bnds
        end
    end

    
    # trace barcodes that produced syndromes == 0
    trace_barcodes()

    return [barcode_candidates, gene_nums]
end


function trace_barcodes(

)
    barcode_candidates = []
    gene_nums = []
    for (i, s) in enumerate(syndromes[dot])
        if iszero(s)
            #println("found zero syndrome: ", dot, i)
            barcode_candidate = recursive_get_synd_neighbors_mem_eff(pnts, g, dot, i, syndromes, len_syndromes)
            #println("barcode_candidate: ", barcode_candidate)
            if ismissing(barcode_candidate)
                #println("missing barcode candidate")
                continue
            end
            cw = pnts.coeff[barcode_candidate]
            if cw in keys(cw_dict)
                push!(barcode_candidates, barcode_candidate)
                push!(gene_nums, cw_dict[cw])
            end
        end
    end

    return [barcode_candidates, gene_nums]
        
end
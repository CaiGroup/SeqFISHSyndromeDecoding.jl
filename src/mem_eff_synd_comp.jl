using Graphs

function find_barcodes_mem_eff(pnts ::DataFrame,
                               g :: DotAdjacencyGraph,
                               cw_dict
    )

    
    syndromes = fill(Vector{SyndromeComponent}, nrow(pnts)) #Vector{Vector{SyndromeComponent}}()
    syndrome_block_bounds = fill(Vector{Int64}, nrow(pnts))
    sizehint!(syndromes, nv(g.g))
    final_round_dots = get_cw_pos_inds(g.n)
    final_round_dots_start_m1 = minimum(final_round_dots) - 1


    unprocessed_inrange_dots = fill(0, nv(g.g))
    last_round_tree = make_KDTree(pnts[final_round_dots, :])

    # Count how many dots in other rounds are in range of each dot in the last round
    for dot in final_round_dots
        for round in 1:(g.n-1)
            unprocessed_inrange_dots[dot] += length(inrange(g.trees[round], [pnts.x[dot], pnts.y[dot]], g.lat_thresh))
        end 
    end

    # Count how many dots in the last round are within the search radius of each dot in previous round
    final_round_inrange_overlap_graph = SimpleWeightedGraph(nrow(pnts))
    final_round_inrange_computed_dots = unprocessed_inrange_dots[final_round_dots]
    for dot in 1:(g.cw_pos_bnd(g.n)-1)
        inrange_last_round_dots = inrange(last_round_tree, [pnts.x[dot], pnts.y[dot]], g.lat_thresh, true)
        unprocessed_inrange_dots[dot] = length(inrange_last_round_dots)
        for i in 1:length(inrange_last_round_dots)
            for j in (i+1):length(inrange_last_round_dots)
                doti = inrange_last_round_dots[i]
                dotj = inrange_last_round_dots[j]
                old_w = get_weight(final_round_inrange_overlap_graph, doti, dotj)
                add_edge!(final_round_inrange_overlap_graph, doti, dotj, old_w + 1)
            end
        end

    end

    unprocessed_last_round_dots = get_cw_pos_inds(g,g.n)
    barcode_candidates = []
    gene_nums = []
    # while there are dots in the last round that have not had their candidate barcodes searched for
    while len(unprocessed_last_round_dots) > 0
        # find the final roud dot that han the fewest number (or any tied for fewest) of dots in the previous rounds that have 
        # not yet contributed to a syndrome computation
        dot_ind = argmin(final_round_inrange_computed_dots) # unprocessed final round indexing
        dot = unprocessed_last_round_dots[dot_ind] # all points index

        # Compute syndromes for and find barcode candidates including that dot
        [dot_barcode_candidates, dot_gene_nums] = find_final_round_dot_barcode_candidates!(
            g :: DotAdjacencyGraph,
            pnts :: DataFrame,
            cw_dict :: Dict,
            syndromes,
            len_syndromes,
            unprocessed_inrange_dots,
            syndrome_block_bounds,
            dot
            )

        append!(barcode_candidates, dot_barcode_candidates)
        append!(gene_nums, dot_gene_nums)
    

        # update uncomputed inrange dots
        for neighbor in neighbors(final_round_inrange_overlap_graph) # all pnt indexing
            overlap = get_weight(final_round_inrange_overlap_graph, dot, neighbor)
            neighbor_ind = indexin(neighbor, final_round_inrange_computed_dots) # convert to unprocessed final round indexing
            final_round_inrange_computed_dots[neighbor_ind] -= overlap
        end

        # remove dot from list of unprocessed final round dots
        deleteat!(unprocessed_last_round_dots, dot_ind)

        #free space
        syndromes[dot] = []
        for round in 1:(g.n-1)
            for inrange_dot in inrange(g.trees[round], [pnts.x[dot], pnts.y[dot]], g.lat_thresh)
                @inbounds unprocessed_inrange_dots[inrange_dot] -= 1
                if unprocessed_inrange_dots[inrange_dot] == 0
                    syndromes[inrange_dot] = []
                    syndrome_block_bounds[inrange_dot] = []
                end
            end
        end
    end
    return [barcode_candidates, gene_nums]
end


"""
    find_final_round_dot_barcode_candidates()
"""
function find_final_round_dot_barcode_candidates!(
    g :: DotAdjacencyGraph,
    pnts :: DataFrame,
    cw_dict :: Dict,
    syndromes,
    len_syndromes,
    unprocessed_inrange_dots,
    syndrome_block_bounds,
    dot
    )
    # compute dot syndromes
    recursive_syndrome_computation!(pnts, g, syndromes, len_syndromes, unprocessed_inrange_dots,
        syndrome_block_bounds, dot)

    
    # trace barcodes that produced syndromes == 0
    barcode_candidates = []
    gene_nums = []
    for (i, s) in enumerate(syndromes[dot])
        if all(iszero.(s))
            barcode_candidate = recursive_get_synd_neighbors(pnts, g, dot, i, syndromes)
            cw = pnts.coeff[barcode_candidate]
            if cw in keys(cw_dict)
                push!(barcode_candidates, barcode_candidate)
                push!(gene_nums, cw_dict[barcode_candidate])
            end
        end
    end
    
    # decrease counters for number of unprocessed inrange final round dots for each previous round dot inrange 
    # of the dot inputed to this function
    for tree in g.trees
        early_dots = inrange(tree, [pnts.x[dot], pnts.y[dot]], g.lat_thresh)
        unprocessed_inrange_dots[early_dots] .-= 1
    end

    # clear allocated syndrome components for dots not in the last round that have had all final round dots in their range processed

    return [barcode_candidates, gene_nums]
end

"""
    recursive_syndrome_computation()
"""
function recursive_syndrome_computation!(
    pnts :: DataFrame,
    g :: DotAdjacencyGraph,
    syndromes,
    len_syndromes,
    unprocessed_inrange_dots,
    syndrome_block_bounds,
    dot
    )

    # if allocated 
    if len(syndromes[dot]) > 0
        return
    # elif not allocated
    else
        # if all in range final round dots processed
        if unprocessed_inrange_dots[dot] == 0
            return 
        #elseif dot in first round
        elseif dot < get_cw_pos_bounds
            #allocate and return
            @inbounds syndromes[dot] = [pnts.sc[dot]]
            return
        else
            #init nsynd counter
            nsyndc = 0
            dot_neighbors = neighbors(g, dot)
            for dot in dot_neighbors
                recursive_syndrome_computation!(pnts, g, syndromes, len_syndromes, unprocessed_inrange_dots, dot)
                @inbounds nsyndc += length(syndrome[dot])
            end
            # allocate syndome component array of length nsnd
            @inbounds syndromes[dot] = fill(pnts.sc[dot], nsyndc)
            @inbounds len_syndromes[dot] = nsyndc
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
end



function recursive_get_synd_neighbors(
    pnts :: DataFrame,
     g :: DotAdjacencyGraph,
     dot :: Int,
     synd_ind :: Int,
     syndromes
     )

    # if this is the last dot in the message, return number of dot in an array
    if dot < g.cw_pos_bnds[2+g.ndrops]    
        cpath = Int[dot]
        sizehint!(cpath, g.n)
        return cpath
    end

    #otherwise, get neighbor of the dot that produced the zero syndrome
    neighbor, neighbor_synd_ind = get_synd_neighbor(g, dot, synd_ind, syndromes)

    # if neighbor has been cleared, return missing

    # Add result to recursively defined array, and return
    push!(recursive_get_synd_neighbors(pnts, g, neighbor, neighbor_synd_ind, syndromes), dot)
end

"""
Helper Function used to trace back groups of dots that produce zero syndrome, and therefore are codewords
"""
function get_synd_neighbor_mem_eff(
    g :: DotAdjacencyGraph,
    dot :: Int,
    synd_ind :: Int,
    syndromes,
    len_syndromes
    )

    #@assert synd_ind != 1

    # keep track of number of syndromes from the dot that have been searched through
    #cum_n_syndromes = 1
    if dot < g.cw_pos_bnds[2]
        cum_n_syndromes = 1
    else
        cum_n_syndromes = 0
    end
    #cum_n_syndromes = dot < g.cw_pos_bnds[2+g.ndrops] ? 1 : 0

    #for each neighbor
    dot_neighbors = neighbors(g, dot)
    for neighbor in dot_neighbors
        # move the index tracker to the index just past the end of where syndrome
        # components from this neighbor are stored.
        @inbounds cum_n_syndromes += len_syndromes[neighbor] #length(syndromes[neighbor])

        #if the syndrome includes a component from this particular neighbor
        if synd_ind <= cum_n_syndromes
            # find the index of the neighbors syndrome component of interest
            @inbounds neighbor_synd_ind = synd_ind - cum_n_syndromes + length(syndromes[neighbor])
            return [neighbor, neighbor_synd_ind]
        end
    end
    error("We shouldn't have gotten here!")
end
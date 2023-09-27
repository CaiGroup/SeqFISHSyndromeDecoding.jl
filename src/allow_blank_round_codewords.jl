function find_blank_round_codewords(pnts ::DataFrame, g :: DotAdjacencyGraphBlankRound, cw_dict, w)
    # initialize array for syndrome partial sums
    syndromes = fill(Vector{Vector{SyndromeComponent}}(), nrow(pnts)) #Vector{Vector{SyndromeComponent}}()
    syndrome_block_sizes = fill(Vector{Vector{Int64}}(), nrow(pnts))
    barcode_candidates = []
    gene_nums = []

    # get list of dots in rounds that could be the last of a codeword
    potential_barcode_final_dot_range = g.first_potential_barcode_final_dot:nrow(pnts)
    potential_barcode_final_dots = collect(potential_barcode_final_dot_range)
    
    round_trees = []
    for r in 1:g.n
        if ismissing(g.cw_round_ranges[r])
            push!(round_trees, make_KDTree(pnts[1:0, :]))
        else
            push!(round_trees, make_KDTree(pnts[g.cw_round_ranges[r], :]))
        end
    end
    pot_final_dot_tree = make_KDTree(pnts[potential_barcode_final_dot_range,:])

    unprocessed_inrange_dots, potential_barcode_final_dots_n_uncomputed_neighbors = count_inrange_dots(pnts, g, round_trees, w)

    while length(potential_barcode_final_dots) > 0
        dot_ind = argmin(potential_barcode_final_dots_n_uncomputed_neighbors[potential_barcode_final_dots .- (g.first_potential_barcode_final_dot-1)])
        dot = potential_barcode_final_dots[dot_ind]
        find_barcode_candidates!(g, pnts, cw_dict, w, syndromes, unprocessed_inrange_dots, syndrome_block_sizes, barcode_candidates, gene_nums, dot)

        #update uncomputed neighbor counts for dots that may be the last in a barcode
        potential_barcode_final_dots_n_uncomputed_neighbors[inrange(pot_final_dot_tree, [pnts.x[dot], pnts.y[dot]], g.lat_thresh)] .-= 1

        # remove dot from list of unprocessed final round dots
        deleteat!(potential_barcode_final_dots, dot_ind)

        # delete pointers to allocated variables that are no longer needed
        #free_space!(pnts, g, unprocessed_inrange_dots, round_trees, syndromes, syndrome_block_sizes, dot)
    end
    return barcode_candidates, gene_nums
end


function count_inrange_dots(pnts, g, round_trees, w)
     # for each dot, count how many dots that could be the final dot of a barcode including it are in range
    # and how many dots that could be the final dot are in range of each other
    unprocessed_inrange_dots = fill(0, nrow(pnts))
    potential_barcode_final_dots_n_uncomputed_neighbors = fill(0, nrow(pnts)-g.first_potential_barcode_final_dot+1)
    for dot in 1:nrow(pnts)
        r = pnts.pos[dot] 
        if r >= w
            unprocessed_inrange_dots[dot] += 1
            for ri in w:r
                ri_inrange = length(inrange(round_trees[ri], [pnts.x[dot], pnts.y[dot]], g.lat_thresh))
                potential_barcode_final_dots_n_uncomputed_neighbors[dot-g.first_potential_barcode_final_dot+1] += ri_inrange
            end
        end

        ri = maximum([w, r+1])
        while ri <= g.n
            ri_inrange = length(inrange(round_trees[ri], [pnts.x[dot], pnts.y[dot]], g.lat_thresh))
            unprocessed_inrange_dots[dot] += ri_inrange
            if r >= w
                potential_barcode_final_dots_n_uncomputed_neighbors[dot-g.first_potential_barcode_final_dot+1] += ri_inrange
            end
            ri += 1
        end
    end
    return unprocessed_inrange_dots, potential_barcode_final_dots_n_uncomputed_neighbors
end

function free_space!(pnts, g, unprocessed_inrange_dots, round_trees, syndromes, syndrome_block_sizes, dot)
    for round in 1:(g.n-1)
        for inrange_dot in inrange(round_trees[round], [pnts.x[dot], pnts.y[dot]], g.lat_thresh)
            inrange_dot_ind = inrange_dot + g.cw_round_ranges[round][1] - 1
            unprocessed_inrange_dots[inrange_dot_ind] -= 1
            if unprocessed_inrange_dots[inrange_dot_ind] == 0
                syndromes[inrange_dot_ind] = []
                syndrome_block_sizes[inrange_dot_ind] = []
            end
        end
    end
end

"""
    find_barcode_candidates()
"""
function find_barcode_candidates!(
    g :: DotAdjacencyGraphBlankRound,
    pnts :: DataFrame,
    cw_dict :: Dict,
    w,
    syndromes,
    unprocessed_inrange_dots,
    syndrome_block_sizes,
    barcode_candidates,
    gene_nums,
    dot
    )
    r = pnts.pos[dot]
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
            nsyndc = fill(0, w)
            len_nbrs = [Int64[] for i in 1:(w-1)] # list comprehention ensures each element points to distinct array
            if r <= g.n - w + 1
                nsyndc[1] = 1
            end
            dot_synd_cs = Vector{Vector{SyndromeComponent}}[]
            sizehint!(dot_synd_cs, length(nsyndc))

            dot_neighbors = neighbors(g, dot)
            if length(dot_neighbors) > 0

                start_ind = maximum([1, r - (Int64(g.n) - w + 1)])
                for ndot in dot_neighbors
                    find_barcode_candidates!(g, pnts, cw_dict, w, syndromes, unprocessed_inrange_dots, syndrome_block_sizes, barcode_candidates, gene_nums, ndot)
                    append!.(len_nbrs[start_ind:end], length.(syndromes[ndot][start_ind:end]))
                end
                nsyndc[2:end] .+= sum.(len_nbrs)
            end
            # allocate syndome component array of length nsnd
            dot_synd_cs = [fill(pnts.sc[dot], nsyndci) for nsyndci in nsyndc]
            sum_dot_syndrome_components!(g, r, nsyndc, dot_neighbors, syndromes, dot_synd_cs, dot)
            # trace barcodes that produced syndromes == 0
            @inbounds syndrome_block_sizes[dot] = len_nbrs
            @inbounds syndromes[dot] = dot_synd_cs[1:(end-1)]
            map(is -> trace_barcode!(is, pnts, g, dot, cw_dict, w, syndromes, syndrome_block_sizes, barcode_candidates, gene_nums), enumerate(dot_synd_cs[end]))
        end
    end
end

function sum_dot_syndrome_components!(g, r, nsyndc, dot_neighbors, syndromes, dot_synd_cs, dot)
    synd_ind = fill(1, length(nsyndc))
    # for dot in neighbors
    start = maximum([2, g.w - (g.n - r)])
    for neighbor in dot_neighbors

        # add syndrome components to appropriate block
        @inbounds end_ind = synd_ind[2:end] .+ length.(syndromes[neighbor]) .- 1
        for i in start:length(dot_synd_cs)
            dot_synd_cs[i][synd_ind[i]:end_ind[i-1]] += syndromes[neighbor][i-1]
        end

        # keep track of indices of each neighbor's syndrome component block
        synd_ind[2:end] .= end_ind .+ 1
    end
    
end

function trace_barcode!(is, pnts, g, dot, cw_dict, w, syndromes, syndrome_block_sizes, barcode_candidates, gene_nums)
    (i, s) = is
    if iszero(s)
        barcode_candidate = recursive_get_synd_neighbors_blank_rounds(pnts, g, dot, i, 1, w, syndromes, syndrome_block_sizes)
        if ismissing(barcode_candidate)
            return
        end
        cw = zeros(UInt8, g.n)
        cw[pnts.round[barcode_candidate]] .= pnts.coeff[barcode_candidate]
        if cw in keys(cw_dict)
            push!(barcode_candidates, barcode_candidate)
            push!(gene_nums, cw_dict[cw])
        else
        end
    end
end

function recursive_get_synd_neighbors_blank_rounds(pnts, g, dot, synd_ind, recursion_depth, w, syndromes, syndrome_block_sizes)
    # if this is the last dot in the message, return number of dot in an array
    if recursion_depth == w 
        cpath = Int[dot]
        sizehint!(cpath, w)
        return cpath
    end

    #otherwise, get neighbor of the dot that produced the zero syndrome
    neighbor, neighbor_synd_ind = get_synd_neighbors_blank_round(g, dot, synd_ind, syndrome_block_sizes, recursion_depth, w)

    # if neighbor has been cleared, return missing
    if length(syndromes[neighbor]) == 0
        return missing
    end

    # Add result to recursively defined array, and return
    res = recursive_get_synd_neighbors_blank_rounds(pnts, g, neighbor, neighbor_synd_ind, recursion_depth+1, w, syndromes, syndrome_block_sizes)
    if ismissing(res)
        return missing
    else
        return push!(res, dot)
    end
end

"""
Helper Function used to trace back groups of dots that produce zero syndrome, and therefore are codewords
"""
function get_synd_neighbors_blank_round(g, dot, synd_ind, syndrome_block_sizes, recursion_depth, w)
    # keep track of number of syndromes from the dot that have been searched through
    cum_n_syndromes = 0

    #for each neighbor
    dot_neighbors = neighbors(g, dot)
    for (i, neighbor) in enumerate(dot_neighbors)
        # move the index tracker to the index just past the end of where syndrome
        # components from this neighbor are stored.
        
        cum_n_syndromes += syndrome_block_sizes[dot][w-recursion_depth][i]

        #if the syndrome includes a component from this particular neighbor
        if synd_ind <= cum_n_syndromes
            # find the index of the neighbors syndrome component of interest
            neighbor_synd_ind = synd_ind - cum_n_syndromes + syndrome_block_sizes[dot][w-recursion_depth][i]
            return [neighbor, neighbor_synd_ind]
        end
    end
    error("We shouldn't have gotten here!")
end
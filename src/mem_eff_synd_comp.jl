make_block_trees(pnts, g :: DotAdjacencyGraphRegistered2D) = [make_KDTree2D(pnts[get_cw_pos_inds(g, r), :]) for r in 1:g.n] 

make_block_trees(pnts, g :: DotAdjacencyGraphRegistered3D) = [make_KDTree3D(pnts[get_cw_pos_inds(g, r), :], g.lat_thresh, g.z_thresh) for r in 1:g.n]

function make_block_trees(pnts, g :: DotAdjacencyGraphPairwise2D)
    block_trees = []
    for block in 1:g.n
        block_pseudocolor_trees = []
        for pseudocolor in 1:q # q is global
            if sum(pnts.pos .== block .&& pnts.coeff .== (pseudocolor .% q)) != 0
            end
            push!(block_pseudocolor_trees, make_KDTree2D(pnts[pnts.pos .== block .&& pnts.coeff .== (pseudocolor .% q), :]))
        end
        push!(block_trees, block_pseudocolor_trees)
    end
    return block_trees
end

function make_block_trees(pnts, g :: DotAdjacencyGraphPairwise3D)
    block_trees = []
    for block in 1:g.n
        block_pseudocolor_trees = []
        for pseudocolor in 1:q # q is global
            push!(block_pseudocolor_trees, make_KDTree3D(pnts[pnts.pos .== block .&& pnts.coeff .== pseudocolor, :]))
        end
        push!(block_trees, block_pseudocolor_trees)
    end
    return block_trees
end

inrng(tree, dot, g :: DotAdjacencyGraphRegistered2D, tforms :: Nothing, block, adjust=true) = inrange(tree, [g.pnts.x[dot], g.pnts.y[dot]], g.lat_thresh, true)
inrng(tree, dot, g :: DotAdjacencyGraphRegistered3D, tforms :: Nothing, block, adjust=true) = inrange(tree, [g.pnts.x[dot], g.pnts.y[dot], g.pnts.z[dot]], g.lat_thresh, true)

"""
function inrng(tree, dot, g :: DotAdjacencyGraphRegistered2D, tforms, block_search)
    tform = tforms[g.pnts.block[dot], block_search]
    registered = Array(g.pnts[dot, [:x, :y]]) * tform[1:2, 1:2] .+ tform[1:2, 4]
    return inrange(tree, registered, g.lat_thresh, true)
end

function inrng(tree, dot, g :: DotAdjacencyGraphRegistered3D, tforms, block_search)
    tform = tforms[g.pnts.block[dot], block_search]
    registered = Array(g.pnts[dot, [:x, :y, :z]]) * tform[:, 1:3] .+ tform[:, 4]
    return inrange(tree, registered, g.lat_thresh, true)
end
"""

function inrng(block_trees, dot, g :: DotAdjacencyGraphPairwise2D, tforms :: Dict, nbr_block, adjust=true)
    inrange_dots = []
    searching_block = g.pnts.pos[dot]
    searching_pseudocolor = g.pnts.coeff[dot]
    search_pc_ind = searching_pseudocolor == 0 ? q : searching_pseudocolor
    for nbr_block_pseudocolor in 1:q #0:(q-1)
        nbr_rnd_pc_ind = nbr_block_pseudocolor == 0 ? q : nbr_block_pseudocolor
        tform = get_tform(tforms, searching_block, searching_pseudocolor, nbr_block, nbr_rnd_pc_ind % q)
        registered_dot_coordinates = Array(tform[1:2, 1:2] * Array(g.pnts[dot, [:x, :y]])) .+ tform[1:2, 4]
        
        block_pseudocolor_dots = inrange(block_trees[nbr_rnd_pc_ind], registered_dot_coordinates, g.lat_thresh, true)
        if length(block_pseudocolor_dots) > 0
            if adjust
                block_pseudocolor_dots .+= g.block_pc_block_starts[nbr_block, nbr_block_pseudocolor]
            end
        end
        inrange_dots = vcat(inrange_dots, block_pseudocolor_dots)
    end
    return inrange_dots
end

function find_barcodes_mem_eff(pnts ::DataFrame, g :: DotAdjacencyGraph, cw_dict, tforms)
    
    syndromes = fill(Vector{SyndromeComponent}(), nrow(pnts))
    syndrome_block_sizes = fill(Vector{Int64}(), nrow(pnts))
    sizehint!(syndromes, nv(g.g))
    final_block_dots = get_cw_pos_inds(g, g.n)

    unprocessed_inrange_dots = fill(0, nv(g.g))

    block_trees = make_block_trees(pnts, g)

    last_block_tree = block_trees[end]

    # Count how many dots in other blocks are in range of each dot in the last block
    for dot in final_block_dots
        for block in 1:(g.n-1)
            unprocessed_inrange_dots[dot] += length(inrng(block_trees[block], dot, g, tforms, block))
            #unprocessed_inrange_dots[dot] += length(inrange(block_trees[block], [pnts.x[dot], pnts.y[dot]], g.lat_thresh))
        end 
    end

    # Count how many dots in the last block are within the search radius of each dot in previous block
    for dot in 1:(g.cw_pos_bnds[g.n]-1)
        inrange_last_block_dots = inrng(last_block_tree, dot, g, tforms, g.n)
        #println("dot: $dot; inrange_last_block_dots: $inrange_last_block_dots")

        #inrange_last_block_dots = inrange(last_block_tree, [pnts.x[dot], pnts.y[dot]], g.lat_thresh, true)
        unprocessed_inrange_dots[dot] = length(inrange_last_block_dots)
    end
    #println("unprocessed_inrange_dots: $unprocessed_inrange_dots")
    # Count how many dots in the last block are in range of each other
    final_block_dot_n_uncomputed_neighbors = [length(inrng(last_block_tree, dot, g, tforms, g.n)) - 1 for dot in final_block_dots]
    #final_block_dot_n_uncomputed_neighbors = [length(inrange(last_block_tree, [pnts.x[dot], pnts.y[dot]], g.lat_thresh)) - 1 for dot in final_block_dots]
    #println("final_block_dot_n_uncomputed_neighbors")
    #println(final_block_dot_n_uncomputed_neighbors)

    unprocessed_last_block_dots = collect(get_cw_pos_inds(g,g.n))
    barcode_candidates = []
    gene_nums = []
    # while there are dots in the last block that have not had their candidate barcodes searched for
    while length(unprocessed_last_block_dots) > 0
        # find the final block dot that has the fewest number of other uncomputed final block dots in
        # its range
        dot_ind = argmin(final_block_dot_n_uncomputed_neighbors[unprocessed_last_block_dots .- (g.cw_pos_bnds[g.n]-1)])
        dot = unprocessed_last_block_dots[dot_ind] 

        # Compute syndromes for and find barcode candidates including that dot
        dot_barcode_candidates, dot_gene_nums = find_final_block_dot_barcode_candidates!(
            g :: DotAdjacencyGraph,
            pnts :: DataFrame,
            cw_dict :: Dict,
            syndromes,
            unprocessed_inrange_dots,
            syndrome_block_sizes,
            dot
            )

        append!(barcode_candidates, dot_barcode_candidates)
        append!(gene_nums, dot_gene_nums)
    
        #final_block_dot_n_uncomputed_neighbors[inrange(last_block_tree, [pnts.x[dot], pnts.y[dot]], g.lat_thresh)] .-= 1
        #println("inrng(last_block_tree, dot, g, tforms, g.n): ", inrng(last_block_tree, dot, g, tforms, g.n))
        #println("final_block_dot_n_uncomputed_neighbors")
        #println(final_block_dot_n_uncomputed_neighbors)
        final_block_dot_n_uncomputed_neighbors[inrng(last_block_tree, dot, g, tforms, g.n, false)] .-= 1

        # remove dot from list of unprocessed final block dots
        deleteat!(unprocessed_last_block_dots, dot_ind)        

        #free space
        #=
        for block in 1:(g.n-1)
            #for inrange_dot in inrange(block_trees[block], [pnts.x[dot], pnts.y[dot]], g.lat_thresh)
            for inrange_dot in inrng(block_trees[block], dot, g, tforms, block, true)
                if typeof(g) <: DotAdjacencyGraphRegistered2D || typeof(g) <: DotAdjacencyGraphRegistered3D
                    inrange_dot_ind = inrange_dot + g.cw_pos_bnds[block] - 1
                else
                    inrange_dot_ind = inrange_dot
                end
                unprocessed_inrange_dots[inrange_dot_ind] -= 1
                if unprocessed_inrange_dots[inrange_dot_ind] == 0
                    syndromes[inrange_dot_ind] = []
                    syndrome_block_sizes[inrange_dot_ind] = []
                end
            end
        end
        =#
    end
    return [barcode_candidates, gene_nums]
end


"""
    find_final_block_dot_barcode_candidates()
"""

function find_final_block_dot_barcode_candidates!(
    g :: DotAdjacencyGraph,
    pnts :: DataFrame,
    cw_dict :: Dict,
    syndromes,
    unprocessed_inrange_dots,
    syndrome_block_sizes,
    dot
    )
    # compute dot syndromes
    recursive_syndrome_computation!(pnts, g, syndromes, unprocessed_inrange_dots,
        syndrome_block_sizes, dot)

    # trace barcodes that produced syndromes == 0
    barcode_candidates = []
    gene_nums = []
    for (i, s) in enumerate(syndromes[dot])
        if iszero(s)
            barcode_candidate = recursive_get_synd_neighbors_mem_eff(pnts, g, dot, i, syndromes, syndrome_block_sizes)
            if ismissing(barcode_candidate)
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

"""
    recursive_syndrome_computation()
"""

function recursive_syndrome_computation!(
    pnts :: DataFrame,
    g :: DotAdjacencyGraph,
    syndromes,
    unprocessed_inrange_dots,
    syndrome_block_sizes,
    dot :: Integer
    )

    # if allocated 
    if length(syndromes[dot]) > 0
        return
    # elif not allocated
    else
        # if all in range final block dots processed
        if unprocessed_inrange_dots[dot] == 0
            return 
        #elseif dot in first block
        elseif dot < g.cw_pos_bnds[2]
            #allocate and return
            @inbounds syndromes[dot] = [pnts.sc[dot]]
            return
        else
            dot_neighbors = neighbors(g, dot)
            for dot in dot_neighbors
                recursive_syndrome_computation!(pnts, g, syndromes, unprocessed_inrange_dots, syndrome_block_sizes, dot)
            end

            len_nbrs = length.(syndromes[dot_neighbors])
            # allocate syndome component array of length nsnd
            @inbounds syndromes[dot] = fill(pnts.sc[dot], sum(len_nbrs)) #nsyndc)
                
            synd_ind = 1
            # for dot in neighbors
            for (i, neighbor) in enumerate(dot_neighbors)
                # add syndrome components to appropriate block
                @inbounds end_ind = synd_ind+length(syndromes[neighbor])-1
                @inbounds syndromes[dot][synd_ind:end_ind] += syndromes[neighbor]
                # keep track of indices of each neighbor's syndrome component block
                synd_ind = end_ind + 1
            end
            @inbounds syndrome_block_sizes[dot] = len_nbrs
        end
    end
end


function recursive_get_synd_neighbors_mem_eff(
    pnts :: DataFrame,
     g :: DotAdjacencyGraph,
     dot :: Int,
     synd_ind :: Int,
     syndromes,
     syndrome_block_sizes
     )

    # if this is the last dot in the message, return number of dot in an array
    if dot < g.cw_pos_bnds[2+g.ndrops]    
        cpath = Int[dot]
        sizehint!(cpath, g.n)
        return cpath
    end

    #otherwise, get neighbor of the dot that produced the zero syndrome
    neighbor, neighbor_synd_ind = get_synd_neighbors_mem_eff(g, dot, synd_ind, syndrome_block_sizes)

    # if neighbor has been cleared, return missing
    if length(syndromes[neighbor]) == 0
        return missing
    end

    # Add result to recursively defined array, and return
    res = recursive_get_synd_neighbors_mem_eff(pnts, g, neighbor, neighbor_synd_ind, syndromes, syndrome_block_sizes)
    if ismissing(res)
        return missing
    else
        return push!(res, dot)
    end
end


"""
Helper Function used to trace back groups of dots that produce zero syndrome, and therefore are codewords
"""

function get_synd_neighbors_mem_eff(
    g :: DotAdjacencyGraph,
    dot :: Int,
    synd_ind :: Int,
    syndrome_block_sizes
    )

    # keep track of number of syndromes from the dot that have been searched through
    cum_n_syndromes = 0

    #for each neighbor
    dot_neighbors = neighbors(g, dot)
    for (i, neighbor) in enumerate(dot_neighbors)
        # move the index tracker to the index just past the end of where syndrome
        # components from this neighbor are stored.
        @inbounds cum_n_syndromes += syndrome_block_sizes[dot][i]

        #if the syndrome includes a component from this particular neighbor
        if synd_ind <= cum_n_syndromes
            # find the index of the neighbors syndrome component of interest
            @inbounds neighbor_synd_ind = synd_ind - cum_n_syndromes + syndrome_block_sizes[dot][i]
            return [neighbor, neighbor_synd_ind]
        end
    end
    error("We shouldn't have gotten here!")
end
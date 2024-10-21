inrng(tree, dot, g :: DotAdjacencyGraphBlankBlock3D) = inrange(tree, [g.pnts.x[dot], g.pnts.y[dot], g.pnts.z[dot]], g.lat_thresh)
inrng(tree, dot, g :: DotAdjacencyGraphBlankBlock2D) = inrange(tree, [g.pnts.x[dot], g.pnts.y[dot]], g.lat_thresh)

function inrng(tree, dot, g :: DotAdjacencyGraphBlankBlock2D, tforms, block_search)
    if tforms == nothing
        return inrange(tree, [g.pnts.x[dot], g.pnts.y[dot]], g.lat_thresh, true)
    else
        tform = tforms[g.pnts.block[dot], block_search]
        registered = Array(g.pnts[dot, [:x, :y]]) * tform[1:2, 1:2] .+ tform[1:2, 4]
        return inrange(tree, registered, g.lat_thresh, true)
    end
end

function inrng(tree, dot, g :: DotAdjacencyGraphBlankBlock3D, tforms, block_search)
    if tforms == nothing
        return inrange(tree, [g.pnts.x[dot], g.pnts.y[dot], g.pnts.z[dot]], g.lat_thresh, true)
    else
        tform = tforms[g.pnts.block[dot], block_search]
        registered = Array(g.pnts[dot, [:x, :y, :z]]) * tform[:, 1:3] .+ tform[:, 4]
        return inrange(tree, registered, g.lat_thresh, true)
    end
end

"""
    find_blank_block_codewords(pnts ::DataFrame, g :: DotAdjacencyGraphBlankBlock, cw_dict, w)

To conserve memory when summing the parity check equation to identify codepath or paths correctable to codepaths, we want to have as 
few intermediate sums allocated at any time as is necessary. We accomplish this by keeping lists of the dots that may be the terminal
dot in a codepath and the dots that have not had their intermediate sums evaluated yet. At each iteration, we evaluate sums for the next
potentially terminal dot that has the fewest inrange potentially terminal dots that have not had their potential terminal codepaths evaluated.
After all potentially terminal dots within the search radius of a dot have had their codepaths identified, the dots intermediate sums are
no longer needed and are cleared.
"""
function find_blank_block_codewords(pnts ::DataFrame, g :: DotAdjacencyGraphBlankBlock, cw_dict, w, tforms)
    # initialize array for syndrome partial sums
    syndromes = fill(Vector{Vector{SyndromeComponent}}(), nrow(pnts)) #Vector{Vector{SyndromeComponent}}()
    syndrome_block_sizes = fill(Vector{Vector{Int64}}(), nrow(pnts))
    barcode_candidates = []
    gene_nums = []

    # get list of dots in blocks that could be the last of a codepath/barcode
    if g.first_potential_barcode_final_dot == nothing
        return barcode_candidates, gene_nums
    end
    potential_barcode_final_dot_range = g.first_potential_barcode_final_dot:nrow(pnts)
    potential_barcode_final_dots = collect(potential_barcode_final_dot_range)

    if typeof(g) ==  DotAdjacencyGraphBlankBlock2D
        make_KDTree = make_KDTree2D
    else
        make_KDTree(df) = make_KDTree3D(df, g.lat_thresh, g.z_thresh)
    end
    
    # make KDTree to search dots in each barcoding block
    block_trees = []
    for r in 1:g.n
        if ismissing(g.cw_block_ranges[r])
            push!(block_trees, make_KDTree(pnts[1:0, :]))
        else
            push!(block_trees, make_KDTree(pnts[g.cw_block_ranges[r], :]))
        end
    end
    pot_final_dot_tree = make_KDTree(pnts[potential_barcode_final_dot_range,:])

    unprocessed_inrange_pot_term_dots, potential_barcode_final_dots_n_uncomputed_neighbors = count_inrange_dots(pnts, g, block_trees, w, tforms)

    while length(potential_barcode_final_dots) > 0
        dot_ind = argmin(potential_barcode_final_dots_n_uncomputed_neighbors[potential_barcode_final_dots .- (g.first_potential_barcode_final_dot-1)])
        dot = potential_barcode_final_dots[dot_ind]
        find_barcode_candidates!(g, pnts, cw_dict, w, syndromes, unprocessed_inrange_pot_term_dots, syndrome_block_sizes, barcode_candidates, gene_nums, dot)

        #update uncomputed neighbor counts for dots that may be the last in a barcode
        potential_barcode_final_dots_n_uncomputed_neighbors[inrng(pot_final_dot_tree, dot, g)] .-= 1

        # remove dot from list of unprocessed final block dots
        deleteat!(potential_barcode_final_dots, dot_ind)

        # delete pointers to allocated variables that are no longer needed
        #free_space!(pnts, g, unprocessed_inrange_pot_term_dots, block_trees, syndromes, syndrome_block_sizes, dot, tforms)
    end
    return barcode_candidates, gene_nums
end


function count_inrange_dots(pnts, g, block_trees, w, tforms)
    # for each dot, count how many dots that could be the terminal dot of a codepath including it are in range
    # and how many dots that could be the final dot are in range of each other

    # this array gives the number of dots that could potentally terminate a codepath including 
    # it are inrange of each dot that have not been computed yet
    unprocessed_inrange_pot_term_dots = fill(0, nrow(pnts)) 
    
    # this array gives the number of potentially terminal dots that have not had their terminal paths summed that are 
    # neighbors to within a search radius of each potentially codepath terminating dot
    potential_barcode_final_dots_n_uncomputed_neighbors = fill(0, nrow(pnts)-g.first_potential_barcode_final_dot+1)

    for dot in 1:nrow(pnts)
        r = pnts.pos[dot] 
        if r >= w - g.ndrops
            unprocessed_inrange_pot_term_dots[dot] += 1 # this dot is a potential codepath terminus of its own
            for ri in (w-g.ndrops):r # dots in ri < (w-g.ndrops) cannot produce decodable paths when directly connected to the dot of concern
                # 
                ri_inrange = length(inrng(block_trees[ri], dot, g, tforms, ri))
                potential_barcode_final_dots_n_uncomputed_neighbors[dot-g.first_potential_barcode_final_dot+1] += ri_inrange
            end
        end

        ri = maximum([w - g.ndrops, r+1])
        while ri <= g.n
            ri_inrange = length(inrng(block_trees[ri], dot, g, tforms, ri))
            unprocessed_inrange_pot_term_dots[dot] += ri_inrange
            if r >= w -g.ndrops
                potential_barcode_final_dots_n_uncomputed_neighbors[dot-g.first_potential_barcode_final_dot+1] += ri_inrange
            end
            ri += 1
        end
    end
    return unprocessed_inrange_pot_term_dots, potential_barcode_final_dots_n_uncomputed_neighbors
end

function free_space!(pnts, g, unprocessed_inrange_pot_term_dots, block_trees, syndromes, syndrome_block_sizes, dot, tforms)

    for block in 1:(pnts.pos[dot]-1)
        for inrange_dot in  inrng(block_trees[block], dot, g, tforms, block)
            inrange_dot_ind = inrange_dot + g.cw_block_ranges[block][1] - 1
            unprocessed_inrange_pot_term_dots[inrange_dot_ind] -= 1
            if unprocessed_inrange_pot_term_dots[inrange_dot_ind] == 0
                syndromes[inrange_dot_ind] = []
                syndrome_block_sizes[inrange_dot_ind] = []
            end
        end
    end
end

"""
    find_barcode_candidates()

Each candidate barcode will be represented by "codepath" through the DAG of length at least w-ndrops. The pseudocolors/blocks in codepaths of length
w will satisfy the parity check equations. codepaths of length less than w will not satisfy the parity check equations, but missing symbol/dot may 
be found that match them to a closest codeword.
"""
function find_barcode_candidates!(
    g :: DotAdjacencyGraphBlankBlock,
    pnts :: DataFrame,
    cw_dict :: Dict,
    w,
    syndromes,
    unprocessed_inrange_pot_term_dots,
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
        # if all in range final block dots processed
        if unprocessed_inrange_pot_term_dots[dot] == 0
            return 
        #elseif dot unprocessed
        else
            # get block
            r = pnts.pos[dot]
            #init nsynd counter
            nsyndc = fill(0, w)
            len_nbrs = [Int64[] for i in 1:(w-1)] # list comprehention ensures each element points to distinct array
            if r <= g.n - w + 1 + g.ndrops
                nsyndc[1] = 1
            end
            dot_synd_cs = Vector{Vector{SyndromeComponent}}[]
            sizehint!(dot_synd_cs, length(nsyndc))

            dot_neighbors = neighbors(g, dot)
            if length(dot_neighbors) > 0

                start_ind = maximum([1, r - (Int64(g.n) - w  + g.ndrops + 1)])
                for ndot in dot_neighbors
                    if unprocessed_inrange_pot_term_dots[ndot] == 0
                        append!.(len_nbrs[start_ind:end], fill(0, length(len_nbrs[start_ind:end])))
                    else
                        find_barcode_candidates!(g, pnts, cw_dict, w, syndromes, unprocessed_inrange_pot_term_dots, syndrome_block_sizes, barcode_candidates, gene_nums, ndot)
                        append!.(len_nbrs[start_ind:end], length.(syndromes[ndot][start_ind:end]))
                    end
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
            for drops in 1:g.ndrops
                map(is -> trace_barcode_w_drops!(is, pnts, g, dot, cw_dict, w, syndromes, syndrome_block_sizes, barcode_candidates, gene_nums, drops), enumerate(dot_synd_cs[end-drops]))
            end
        end
    end
end

function sum_dot_syndrome_components!(g, r, nsyndc, dot_neighbors, syndromes, dot_synd_cs, dot)
    synd_ind = fill(1, length(nsyndc))
    # for dot in neighbors
    start = maximum([2, g.w - (g.n - r)])
    for neighbor in dot_neighbors
        # add syndrome components to appropriate block
        if length(syndromes[neighbor]) != 0
            @inbounds end_ind = synd_ind[2:end] .+ length.(syndromes[neighbor]) .- 1
            for i in start:length(dot_synd_cs)
                @inbounds dot_synd_cs[i][synd_ind[i]:end_ind[i-1]] += syndromes[neighbor][i-1]
            end

            # keep track of indices of each neighbor's syndrome component block
        @inbounds synd_ind[2:end] .= end_ind .+ 1
        end
    end
    
end

function trace_barcode!(is, pnts, g, dot, cw_dict, w, syndromes, syndrome_block_sizes, barcode_candidates, gene_nums)
    (i, s) = is
    if iszero(s)
        barcode_candidate = recursive_get_synd_neighbors_blank_blocks(pnts, g, dot, i, 1, w, syndromes, syndrome_block_sizes, 0)
        if ismissing(barcode_candidate)
            return
        end
        if typeof(pnts.coeff[1]) == String
            cw = fill("0", g.n)
        else
            cw = zeros(UInt8, g.n)
        end
        cw[pnts.block[barcode_candidate]] .= pnts.coeff[barcode_candidate]
        if cw in keys(cw_dict)
            push!(barcode_candidates, barcode_candidate)
            push!(gene_nums, cw_dict[cw])
        else
        end
    end
end

function trace_barcode_w_drops!(is, pnts, g, dot, cw_dict, w, syndromes, syndrome_block_sizes, barcode_candidates, gene_nums, drops)
    (i, s) = is
    
    barcode_candidate = recursive_get_synd_neighbors_blank_blocks(pnts, g, dot, i, 1, w, syndromes, syndrome_block_sizes, drops)
    if ismissing(barcode_candidate)
        return
    end
    if typeof(pnts.coeff[1]) <: AbstractString
        message = fill("0", g.n)
        message[pnts.block[barcode_candidate]] .= string.(pnts.coeff[barcode_candidate])
    else
        message = zeros(UInt8, g.n)
        message[pnts.block[barcode_candidate]] .= pnts.coeff[barcode_candidate]
    end
    r = find(BKTree_cb, message, drops)
    if length(r) == 1
        push!(barcode_candidates, barcode_candidate)
        push!(gene_nums, cw_dict[r[1][2]])
    end
end

function get_RS_error(bc :: Vector{SyndromeComponent}, ndrops)
    find(BKTree_cb, bc, ndrops)
end

function get_RS_error(H, S, ndrops)
    if ndrops == 1
        σ₀ = 1
        σ₁ = S[2]/S[1]
        σ(X) = σ₀ + σ₁*X
        dσ_dx(X) = σ₁
    end
    σS =  (σ₀ .* S) .+ (σ₁ .*S[[[length(S)]; collect(1:(length(S)-1))]])

    
end


function recursive_get_synd_neighbors_blank_blocks(pnts, g, dot, synd_ind, recursion_depth, w, syndromes, syndrome_block_sizes, ndrops)
    # if this is the last dot in the message, return number of dot in an array
    if recursion_depth == w - ndrops
        cpath = Int[dot]
        sizehint!(cpath, w - ndrops)
        return cpath
    end

    #otherwise, get neighbor of the dot that produced the zero syndrome
    neighbor, neighbor_synd_ind = get_synd_neighbors_blank_block(g, dot, synd_ind, syndrome_block_sizes, recursion_depth, w, ndrops)

    # if neighbor has been cleared, return missing
    if length(syndromes[neighbor]) == 0
        return missing
    end

    # Add result to recursively defined array, and return
    res = recursive_get_synd_neighbors_blank_blocks(pnts, g, neighbor, neighbor_synd_ind, recursion_depth+1, w, syndromes, syndrome_block_sizes, ndrops)
    if ismissing(res)
        return missing
    else
        return push!(res, dot)
    end
end

"""
Helper Function used to trace back groups of dots that produce zero syndrome, and therefore are codewords
"""
function get_synd_neighbors_blank_block(g, dot, synd_ind, syndrome_block_sizes, recursion_depth, w, ndrops)
    # keep track of number of syndromes from the dot that have been searched through
    cum_n_syndromes = 0

    #for each neighbor
    dot_neighbors = neighbors(g, dot)
    for (i, neighbor) in enumerate(dot_neighbors)
        # move the index tracker to the index just past the end of where syndrome
        # components from this neighbor are stored.
        
        cum_n_syndromes += syndrome_block_sizes[dot][w - ndrops - recursion_depth][i]

        #if the syndrome includes a component from this particular neighbor
        if synd_ind <= cum_n_syndromes
            # find the index of the neighbors syndrome component of interest
            neighbor_synd_ind = synd_ind - cum_n_syndromes + syndrome_block_sizes[dot][w - ndrops-recursion_depth][i]
            return [neighbor, neighbor_synd_ind]
        end
    end
    error("We shouldn't have gotten here!")
end
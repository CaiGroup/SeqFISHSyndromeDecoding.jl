using DelimitedFiles
using DataFrames
using LightGraphs
import LightGraphs.neighbors
using NearestNeighbors
using Distributions
using DataStructures
using Distributed

"""
    decode_syndromes!(
        pnts :: DataFrame,
        cb :: Matrix{UInt8},
        H :: Matrix,
        params :: DecodeParams
    )

Takes a DataFrame of aligned points with hyb, x, y, z, w, and s columns;
codebook matix; parity check matrix (H); and DecodeParams stucture with decoding parameters set.
Adds columns to input pnts matrix indicating the encoding round, syndrome component value,
and decoding result indicated as the row of the codebook matrix matched to.
Split up points into weakly connected components, then finds possible codeword messages
and runs simulated annealing to assign them. The pnts dataframe should have hybridization, x, y, and z columns
"""

function decode_syndromes!(pnts :: DataFrame, cb :: Matrix{UInt8}, H :: Matrix, params :: DecodeParams)
    println("start syndrome decoding")
    cpath_df = get_codepaths(pnts, cb, H, params)

    if nrow(cpath_df) == 0
        println("No viable code paths")
        return
    end

    return choose_optimal_codepaths(pnts, cb, H, params, cpath_df)
end


function obj_function(cpath, pnts, cw_n, params)
    lat_var_factor = params.lat_var_factor
    z_var_factor = params.z_var_factor
    lw_var_factor = params.lw_var_factor
    s_var_factor = params.s_var_factor
    dot_erasure_penalty = params.erasure_penalty

    lat_var_cost = (var(pnts.x[cpath]) + var(pnts.y[cpath])) * lat_var_factor
    z_var_cost = var(pnts.z[cpath]) * z_var_factor
    lw_var_cost = var(log2.(pnts.w[cpath])) * lw_var_factor
    s_var_cost = var(pnts.s[cpath]) * s_var_factor

    length(cpath) == cw_n ? erasure_cost = 0 : erasure_cost = length(cpath) * dot_erasure_penalty

    return lat_var_cost + z_var_cost + lw_var_cost + s_var_cost + erasure_cost
end

"""
"""
function check_inputs(pnts :: DataFrame, cb :: Matrix{UInt8}, H :: Matrix, params :: DecodeParams)
    alphabet = sort(unique(cb))
    q = UInt8(length(alphabet))
    n = length(cb[1,:])

    @assert alphabet[1] == 0x00
    @assert alphabet[q] < q
    @assert all(alphabet .>= 0)
    @assert maximum(pnts[!,"hyb"]) <= q*n
    @assert minimum(pnts[!,"hyb"]) > 0
    @assert params.ndrops <= size(H)[1]
    @assert params.ndrops >= 0
end

"""
get_codepaths(pnts :: DataFrame, cb :: Matrix{UInt8}, H :: Matrix, params :: DecodeParams)

Computes codepaths with syndrome decoding, removes codepaths that exceed the cost
of not decoding their component dots, and
and returns DataFrame of candidate codepaths.
"""
function get_codepaths(pnts :: DataFrame, cb :: Matrix{UInt8}, H :: Matrix, params :: DecodeParams)
    cb_dict = make_cw_dict(cb)
    alphabet = sort(unique(cb))
    q = UInt8(length(alphabet))
    n = length(cb[1,:])

    check_inputs(pnts, cb, H, params)

    #q = UInt8(length(unique(cb)))
    set_q(q)
    set_H(H)
    set_n(UInt8(n))
    if params.ndrops > 0
        get_decode_table()
    end

    add_code_cols!(pnts)
    g = DotAdjacencyGraph(pnts, params.lat_thresh, params.z_thresh, n, params.ndrops)

    cost(cpath) = obj_function(cpath, pnts, n, params)
    code_paths, values = syndrome_find_message_paths!(pnts, g, cb, params.ndrops)
    costs = cost.(code_paths)
    cpath_df = DataFrame(cpath = code_paths, cost = costs, value = values)
    sort!(cpath_df, :cost)
    cpath_df = remove_high_cost_cpaths(cpath_df, params.free_dot_cost, n, params.ndrops)
    cpath_df = threshold_cpaths(cpath_df, pnts, params.lat_thresh, params.z_thresh)

    return cpath_df
end

function choose_optimal_codepaths(pnts :: DataFrame, cb :: Matrix{UInt8}, H :: Matrix, params :: DecodeParams, cpath_df :: DataFrame)

    n = length(cb[1,:])
    ndots = nrow(pnts)
    cost(cpath) = obj_function(cpath, pnts, n, params)
    pnts[!,"decoded"] = fill(0, nrow(pnts))
    pnts[!, "mpath"] = [[] for i = 1:length(pnts.x)]

    cpath_df[!, "cost"] = cost.(cpath_df[!, "cpath"])
    sort!(cpath_df, :cost)
    cpath_df = remove_high_cost_cpaths(cpath_df, params.free_dot_cost, n, params.ndrops)
    cpath_df = threshold_cpaths(cpath_df, pnts, params.lat_thresh, params.z_thresh)

    #println("nrow(cpath_df): ", nrow(cpath_df))

    # build graph with by adding only edges in codepaths, and break into connected
    # components
    ccs = get_connected_components(cpath_df.cpath, nrow(pnts))
    cpath_df[!, "cc"] .= 0
    cpath_df[!, "cc_size"] .= 0

    n_ccs = length(ccs)
    println("n_ccs: ", n_ccs)

    mpaths = cpath_df[1:0, :]
    nmpaths = 0
    for (cc_i, cc) in enumerate(ccs)
        cpath_df[cc, "cc"] .= cc_i
        cpath_df[cc, "cc_size"] .= length(cc)
        cc_cpath_df = cpath_df[cc,:]
        println("cc: ", cc_i)
        println("npaths: ", nrow(cc_cpath_df))
        if nrow(cc_cpath_df) > params.skip_thresh
            println("skip")
            continue
        end

        ndots_cc = length(cc)
        if nrow(cc_cpath_df) == 1
            low_cost_state = [1]
        # ToDO: test 2 cpath cc case
        elseif nrow(cc_cpath_df) == 2
            costs = cc_cpath_df[!,"cost"]
            low_cost_state = (costs .== minimum(costs))
        else
            cpath_nbrs, cpath_partial_conflicts, cpath_partial_conflict_transitions = get_cpath_conflict_graph_remove_redundant_cpaths!(cc_cpath_df, ndots, n)
            if nrow(cc_cpath_df) < params.mip_sa_thresh
                low_cost_state = MIP_solve(cc_cpath_df, cpath_nbrs)
            else
                low_cost_state = simulated_annealing(cc_cpath_df, cpath_nbrs, cpath_partial_conflicts, cpath_partial_conflict_transitions, params, ndots_cc)
            end
        end

        mpath_df = cc_cpath_df[low_cost_state,:]
        nmpaths += nrow(mpath_df)
        append!(mpaths, mpath_df)

        for mpath_row in eachrow(mpath_df)
            pnts.decoded[mpath_row.cpath] .= mpath_row.value
            for dt in mpath_row.cpath
                pnts.mpath[dt] = mpath_row.cpath
            end
        end
    end
    println("found $nmpaths mpaths")
    mpaths
end

"""
    make_KDTree(pnts :: DataFrame)

Generate KDTree to aid in building adjacency graphs.
"""
make_KDTree(pnts :: DataFrame) = KDTree(Array([pnts.x pnts.y]'))

"""
    make_cw_dict(cb)

Generate dictionary of target sequences corresponding to each codeword.
"""
function make_cw_dict(cb)
    cw_dict = Dict()
    for i in 1:size(cb)[1]
        cw_dict[cb[i,:]] = i
    end
    cw_dict
end

"""
    add_code_cols!(pnts :: DataFrame, code :: String)

Add columns giving the coefficient and position the dot encodes in a codeword
and its associated syndrome component.
"""
function add_code_cols!(pnts :: DataFrame)
    pos = get_pos.(pnts.hyb)
    coeff = get_coeff.(pnts.hyb, pos)

    sc = SyndromeComponent.(coeff, pos)
    pnts.pos = pos
    pnts.coeff = coeff
    pnts.sc = sc

    pnts.decoded = zeros(length(pnts.x))
    pnts.mpath = [[] for i = 1:length(pnts.x)]
end


"""
    DotAdjacencyGraph(g :: SimpleDiGraph
                      cw_pos_bnds :: Tuple{Int64}
                      n :: Int8
                      q :: Int8)

Structure for storing the dot adjacency graph with some parameters
"""
struct DotAdjacencyGraph
    g :: SimpleDiGraph
    cw_pos_bnds :: Array{Int64}
    n :: Int8
    trees :: Vector{KDTree}
    lat_thresh :: Float64
    pnts :: DataFrame
    ndrops :: Int64
end

"""
    DotAdjacencyGraph(g :: SimpleDiGraph
                      cw_pos_bnds :: Tuple{Int64}
                      n :: Int8
                      q :: Int8)

Construct a dot adjacency graph where dots close to each other have directed
edges pointing towards the dot representing an earlier symbor
"""
function DotAdjacencyGraph(pnts :: DataFrame, lat_thresh :: Real, z_thresh :: Real, n, ndrops) #, q)
    g = SimpleDiGraph(nrow(pnts))
    #tree = make_KDTree(pnts)

    # Find the indices of dots representing each place, cᵢ, in a codeword start.
    cw_pos_bnds = get_cw_pos_bounds(pnts, n)
    trees = []
    for round in 1:n
        start_round = maximum([round-1-ndrops, 1])
        start_pnt = cw_pos_bnds[start_round]
        end_pnt = (cw_pos_bnds[round]-1)
        push!(trees, make_KDTree(pnts[start_pnt:end_pnt, :]))
    end
    #trees = [make_KDTree(pnts[(cw_pos_bnds[round-1-ndrops]):(cw_pos_bnds[round]-1),:]) for round in 1:n]

    # now add edges to the graph
    """
    for cᵢ in 2:n
        for pnt in cw_pos_bnds[cᵢ]:(cw_pos_bnds[cᵢ+1]-1)
            neighbors = inrange(tree, [pnts.x[pnt], pnts.y[pnt]], lat_thresh, true)
            for neighbor in neighbors
                if neighbor < cw_pos_bnds[cᵢ]
                    if "z" in names(pnts)
                        if abs(pnts.z[pnt] - pnts.z[neighbor] <= z_thresh)
                            add_edge!(g, pnt, neighbor)
                        end
                    else
                        add_edge!(g, pnt, neighbor)
                    end
                end
            end
        end
    end
    """

    DotAdjacencyGraph(g, cw_pos_bnds, n, trees, lat_thresh, pnts, ndrops)
end

"""
function DotAdjacencyGraph(pnts :: DataFrame, g :: SimpleDiGraph, n :: Int)
    cw_pos_bnds = get_cw_pos_bounds(pnts, n)
    DotAdjacencyGraph(g, cw_pos_bnds, n)
end
"""

"""
"""
function get_cw_pos_bounds(pnts, n)
    cw_pos_bnds = [1]
    for cᵢ = 1:(n-1)
        #if findfirst(x -> x>qm1*cᵢ, pnts.hyb) == nothing
        if findfirst(x -> x>cᵢ, pnts.pos) == nothing
            push!(cw_pos_bnds, length(pnts.hyb)+1)
        else
            push!(cw_pos_bnds, findfirst(x -> x>cᵢ, pnts.pos))
        end
    end
    push!(cw_pos_bnds,length(pnts.hyb)+1)

    return cw_pos_bnds
end

"""
    neighbors(g :: DotAdjacencyGraph, n)

Define SimpleDiGraph neighbors function for DotAdjacencyGraph
"""
#neighbors(g :: DotAdjacencyGraph, n) = neighbors(g.g, n)
function neighbors(g :: DotAdjacencyGraph, n)
    nbrs = inrange(g.trees[g.pnts.pos[n]], [g.pnts.x[n], g.pnts.y[n]], g.lat_thresh, true)
    pnts_prior_rnds = g.cw_pos_bnds[maximum([g.pnts.pos[n]-1-g.ndrops, 1])] - 1
    return nbrs .+ pnts_prior_rnds
end

"""
    get_cw_pos_inds(g :: DotAdjacencyGraph, pos :: Int)

Helper function to get the dots in the graph that represent a symbol in a given position of the codeword.
"""
function get_cw_pos_inds(g :: DotAdjacencyGraph, pos :: Int)
    return g.cw_pos_bnds[pos]:(g.cw_pos_bnds[pos+1]-1)
end

"""
"""
function syndrome_find_message_paths!(pnts ::DataFrame, g :: DotAdjacencyGraph, cb ::Matrix{UInt8}, ndrops)
    cw_dict = make_cw_dict(cb)
    cpaths, decode_cands = syndrome_find_message_paths!(pnts, g, cb, ndrops, cw_dict)
    return cpaths, decode_cands
end

"""
Computes syndrome of every path in the message graph, determines which ones represent valid
barcodes, and returns dataframe of message path candidates.
"""
function syndrome_find_message_paths!(pnts ::DataFrame,
                                      g :: DotAdjacencyGraph,
                                      cb ::Matrix{UInt8},
                                      ndrops,
                                      cw_dict
                                      )
    syndromes, syndrome_path_lengths, syndrome_coeff_positions = compute_syndromes(pnts, g)#, code)

    cpaths, decode_cands = find_code_paths!(
                    g,
                    pnts,
                    cw_dict,
                    syndromes,
                    syndrome_path_lengths,
                    syndrome_coeff_positions,
                    ndrops
                    )
    return cpaths, decode_cands
end

"""
    compute_syndromes(pnts :: DataFrame, g :: DotAdjacencyGraph)

Calculate the syndromes for paths in the message dag
"""
function compute_syndromes(pnts :: DataFrame, g :: DotAdjacencyGraph)
    syndromes, syndrome_path_length, syndrome_coeff_positions = init_syndromes(pnts, g)

    for cw_pos in 1:g.n
        cw_pos_inds = get_cw_pos_inds(g, cw_pos)

        # for each dot representing a candidate coefficient
        for dot in cw_pos_inds
            synd_ind = 2

            # get neighbors of that dot
            for neighbor in neighbors(g, dot)
                end_ind = synd_ind+length(syndromes[neighbor])-1
                syndromes[dot][synd_ind:end_ind] += syndromes[neighbor]
                syndrome_path_length[dot][synd_ind:end_ind] += syndrome_path_length[neighbor]
                syndrome_coeff_positions[dot][synd_ind:end_ind] += syndrome_coeff_positions[neighbor]
                synd_ind = end_ind + 1
            end
        end
    end
    return syndromes, syndrome_path_length, syndrome_coeff_positions
end

"""
find_nsnds(g :: DotAdjacencyGraph)

Counts the number of paths in the dot adjacency graph that run through each dot to
determine the size of syndrome component arrays to preallocate.
"""
function find_nsnds(g :: DotAdjacencyGraph)
    n = g.n
    n_pnts = nrow(g.pnts)#nv(g.g)
    nsnds = fill(1,n_pnts)

    for bc_round in 0x02:g.n
        for pnt in g.cw_pos_bnds[bc_round]:(g.cw_pos_bnds[bc_round+0x01]-1)
            nbrs = neighbors(g, pnt)
            for nbr in nbrs
                nsnds[pnt] += nsnds[nbr]
            end
        end
    end
    nsnds
end

"""
init_syndromes(pnts :: DataFrame, g :: DotAdjacencyGraph)

Pre-allocates arrays of syndrome partial sums.
"""
function init_syndromes(pnts :: DataFrame, g :: DotAdjacencyGraph)

    # initialize each syndrome partial sum with the contribution from each dot
    nsnds = find_nsnds(g)
    syndromes = Vector{Vector{SyndromeComponent}}()
    syndrome_path_length = Vector{Vector{UInt8}}()
    syndrome_coeff_positions = Vector{Vector{UInt8}}()
    sizehint!(syndromes, nv(g.g))
    for (pnt, sc) in enumerate(pnts.sc)
        push!(syndromes, fill(sc, nsnds[pnt]))
        push!(syndrome_path_length, fill(0x01, nsnds[pnt]))
        pos_pow = UInt8(0x02 ^ (pnts.pos[pnt] - 0x01))
        push!(syndrome_coeff_positions, fill(pos_pow, nsnds[pnt]))
    end
    (syndromes, syndrome_path_length, syndrome_coeff_positions)
end

"""
    find_code_paths!(
                     g :: DotAdjacencyGraph,
                     pnts :: DataFrame,
                     cw_dict :: Dict,
                     syndromes,
                     syndrome_path_lengths :: Vector{Vector{UInt8}},
                     syndrome_coeff_positions :: Vector{Vector{UInt8}}
                    )

Uses syndrome decoding to identify codepaths in the message dag: both perfect and
correctable. Adds columns to the pnts DataFrame listing each dot's perfect
codepaths and correctable codepaths.
"""
function find_code_paths!(
                         g :: DotAdjacencyGraph,
                         pnts :: DataFrame,
                         cw_dict :: Dict,
                         syndromes,
                         syndrome_path_lengths :: Vector{Vector{UInt8}},
                         syndrome_coeff_positions :: Vector{Vector{UInt8}},
                         ndrops
                         )

    cpaths = Vector{Int}[]
    decode_cands = Int[]

    cw_n_symbols = length(collect(keys(cw_dict))[1])


    full_bin_pos_indicator = 0x00
    for i = 0x00:UInt8(cw_n_symbols-1)
        full_bin_pos_indicator += 0x02^i
    end

    for dot_ind in 1:nrow(pnts)
        for (synd_ind, path_length) in enumerate(syndrome_path_lengths[dot_ind])
            if path_length == cw_n_symbols && iszero(syndromes[dot_ind][synd_ind])
                code_path = recursive_get_synd_neighbors(pnts, g, dot_ind, synd_ind, syndromes)
                message = pnts.coeff[code_path]
                if message in keys(cw_dict)
                    push!(cpaths, code_path)
                    push!(decode_cands, cw_dict[message])
                end
            elseif path_length >= cw_n_symbols - ndrops && path_length < cw_n_symbols
                s = syndromes[dot_ind][synd_ind]
                dot_pos_sum = syndrome_coeff_positions[dot_ind][synd_ind]

                # use bitwise masking to get the position of the missing dot
                # positions are bitwise "one-hot" encoded from sums of powers of 2
                drop_pos_pow = dot_pos_sum ⊻ full_bin_pos_indicator
                if drop_pos_pow > 32
                    println("dot_pos_sum: $dot_pos_sum")
                end

                result = check_mpath_decodable(drop_pos_pow, s)

                if result.decodable
                    code_path = recursive_get_synd_neighbors(pnts, g, dot_ind, synd_ind, syndromes)
                    message = pnts.coeff[code_path]
                    #convert drop_pos to standard integer format
                    drop_pos = UInt8(log2(drop_pos_pow)) + 0x01

                    insert!(message, drop_pos, result.coeff)
                    if message in keys(cw_dict)
                    # add imperfect codeword
                        push!(cpaths, code_path)
                        push!(decode_cands, cw_dict[message])
                    end
                end
            end
        end
    end
    return cpaths, decode_cands
end

"""
Helper Function used to trace back groups of dots that produce zero syndrome, and therefore are codewords
"""
function recursive_get_synd_neighbors(
    pnts :: DataFrame,
     g :: DotAdjacencyGraph,
     dot :: Int,
     synd_ind :: Int,
     syndromes
     )

    # if this is the last dot in the message, return number of dot in an array
    if synd_ind == 1
        return Int[dot]
    end

    #otherwise, get neighbor of the dot that produced the zero syndrome
    neighbor, neighbor_synd_ind = get_synd_neighbor(pnts, g, dot, synd_ind, syndromes)

    # Add result to recursively defined array, and return
    push!(recursive_get_synd_neighbors(pnts, g, neighbor, neighbor_synd_ind, syndromes), dot)
end

"""
Helper Function used to trace back groups of dots that produce zero syndrome, and therefore are codewords
"""
function get_synd_neighbor(
    pnts :: DataFrame,
    g :: DotAdjacencyGraph,
    dot :: Int,
    synd_ind :: Int,
    syndromes
    )

    @assert synd_ind != 1

    # keep track of number of syndromes from the dot that have been searched through
    cum_n_syndromes = 1

    #for each neighbor
    dot_neighbors = neighbors(g, dot)
    for neighbor in dot_neighbors
        # move the index tracker to the index just past the end of where syndrome
        # components from this neighbor are stored.
        cum_n_syndromes += length(syndromes[neighbor])

        #if the syndrome includes a component from this particular neighbor
        if synd_ind <= cum_n_syndromes

            # find the index of the neighbors syndrome component of interest
            neighbor_synd_ind = synd_ind - cum_n_syndromes + length(syndromes[neighbor])
            return [neighbor, neighbor_synd_ind]
        end
    end
end

"""
remove_high_cost_cpaths(cpath_df :: DataFrame, free_dot_cost, n :: Int, ndrops :: Int)

remove candidate codepaths that have a higher cost than the free dot cost of their component
dots.
"""
function remove_high_cost_cpaths(cpath_df :: DataFrame, free_dot_cost, n :: Int, ndrops :: Int)
    ncpath_b4 = nrow(cpath_df)
    cutoff = searchsorted(cpath_df.cost, n*free_dot_cost).stop

    # first remove paths that are more costly than a not decoding the number of dots
    # in a full length code path
    if cutoff < nrow(cpath_df)
        cpath_df = cpath_df[1:cutoff, :]
    end

    if cutoff == 0
        return cpath_df
    end

    # now remove code paths with erasures that are more costly than not score_decoding
    # their component dots.
    #i = searchsorted(cpath_df.cost, (n-ndrops)*free_dot_cost).stop
    i = searchsorted(cpath_df.cost, (n-ndrops)*free_dot_cost).start
    n_erasure_paths_removed = 0

    while i <= nrow(cpath_df)
        if cpath_df.cost[i] >= length(cpath_df.cpath[i])*free_dot_cost
            delete!(cpath_df, i)
            n_erasure_paths_removed+= 1
        else
            i += 1
        end
    end
    return cpath_df
end

"""
threshold_cpaths(cpaths_df, pnts, lat_thresh, z_thresh)

remove cpaths from the cpaths_df where the distance any pair of its dots exceeds lat_thresh on the imaging plane
or z_thresh along the z axis
"""
function threshold_cpaths(cpaths_df, pnts, lat_thresh, z_thresh)
    row = 1
    while row <= nrow(cpaths_df)
        xs = pnts.x[cpaths_df.cpath[row]]
        ys = pnts.y[cpaths_df.cpath[row]]
        zs = pnts.z[cpaths_df.cpath[row]]
        exceeds_threshold = false
        len_cp = length(xs)
        for i = 1:(len_cp-1), j = (i+1):len_cp
            z_diff = abs(zs[i] - zs[j])
            lat_diff = sqrt((xs[i]-xs[j])^2 + (ys[i]-ys[j])^2)
            if lat_diff > lat_thresh || z_diff > z_thresh
                exceeds_threshold = true
                break
            end
        end
        if exceeds_threshold
            delete!(cpaths_df, row)
        else
            row += 1
        end
    end
    cpaths_df
end

"""
get_connected_components(cpaths, ndots)

Finds connected components of candidate codepaths that are connected by their
conflicting dots.
"""
function get_connected_components(cpaths, ndots)
    g = SimpleGraph(ndots)
    for cpath in cpaths
        for dt in 1:(length(cpath)-1)
            add_edge!(g, cpath[dt], cpath[dt+1])
        end
    end

    dot_ccs = connected_components(g)
    for dcc in dot_ccs
        @assert issorted(dcc)
    end

    cpath_ccs = [Vector{Int}() for i = 1:length(dot_ccs)]
    for (cpath_i, cpath) in enumerate(cpaths)
        for (cc_i, dt_cc) in enumerate(dot_ccs)
            result = searchsorted(dt_cc, cpath[1])
            if ~isempty(result)
                push!(cpath_ccs[cc_i], cpath_i)
                break
            end
        end
    end
    ccs = filter(cc->~isempty(cc), cpath_ccs)
end

"""
get_cpath_conflict_graph_remove_redundant_cpaths!(cpaths_df, ndots, n)

redundant paths are subpaths with drops of longer paths. If such a subpath has
higher cost than its super path, and the dot(s) by which the paths differ is not
shared with any other codepaths that are not subpaths of the same super path,
than the subpath is strictly worse than the super path and we can remove it.
Removing these subpaths will reduce size of the connected components and easier to
find optimal solutions.

Once redundant subpaths are removed, return list of connecte

cpath_nbr_cpath_indices is a vector of vectors. The nested vectors contains the indices
of cpaths from the input cpaths_df that conflict with the cpath that has the
same index in cpaths_df as the nested vector in cpath_nbr_cpath_indices
"""
function get_cpath_conflict_graph_remove_redundant_cpaths!(cpaths_df, ndots, n)
    println("Building Adjacency Graph")
    cpaths = cpaths_df.cpath
    costs = cpaths_df.cost
    n_cpaths = length(cpaths)

    # get the number of code paths that include each dot
    dot_ncpaths = zeros(Int, ndots)
    for cpath in cpaths
        for dot in cpath
            dot_ncpaths[dot] += 1
        end
    end

    # preallocate list of the code paths that include each dot
    dot_path_list = Vector{Int}[]
    for i in 1:ndots
        push!(dot_path_list, Int[])
        sizehint!(dot_path_list[i], dot_ncpaths[i])
    end

    # get list of the code paths that include each dot
    for (i, cpath) in enumerate(cpaths)
        for dt in cpath
            push!(dot_path_list[dt], i)
        end
    end

    # preallocate list of dictionaries of neighbors for each cpath
    cpath_neighbors_dict = Dict{Int, Bool}[]
    for (cpath_i, cpath) in enumerate(cpaths)
        cpath_i_dict = Dict{Int, Bool}()
        sizehint!(cpath_i_dict, sum(dot_ncpaths[cpath] .- 1))
        push!(cpath_neighbors_dict, cpath_i_dict)
    end

    # fill dictionaries of neighbors for each cpath
    nedges = 0
    for dt in 1:ndots
        for i in 1:(dot_ncpaths[dt]-1), j in (i+1):dot_ncpaths[dt]
            cpath_i = dot_path_list[dt][i]
            cpath_j = dot_path_list[dt][j]
            if ~haskey(cpath_neighbors_dict[cpath_i], cpath_j)
                cpath_neighbors_dict[cpath_i][cpath_j] = true
                cpath_neighbors_dict[cpath_j][cpath_i] = false
                nedges += 1
            end
        end
    end

    # make neighbors Array
    cpath_nbr_cpath_indices = Vector{Int}[]
    cpath_sub_cpath_indices = Vector{Int}[]
    for i = 1:n_cpaths
        nbr_cpath_indices = collect(keys(cpath_neighbors_dict[i]))
        push!(cpath_sub_cpath_indices, [])
        cpath = cpaths[i]
        if length(cpath) == n
            for nbr_cpath_ind in nbr_cpath_indices
                nbr_cpath = cpaths[nbr_cpath_ind]
                if length(nbr_cpath) < n && nbr_cpath ⊆ cpath
                    push!(cpath_sub_cpath_indices[i], nbr_cpath_ind)
                end
            end
        end
        sort!(cpath_sub_cpath_indices[i])
        sort!(nbr_cpath_indices)
        push!(cpath_nbr_cpath_indices, nbr_cpath_indices)
    end
    #return cpath_nbr_cpath_indices, Vector{Int}[], Vector{Int}[]

    #check if cpath_sub_cpath_indices are redundant
    #println("check if subpaths are redundant")
    redundant_subcodepath_indices = []
    partial_conflicts = [Int[] for i = 1:n_cpaths]
    partial_conflict_transitions = [Int[] for i = 1:n_cpaths]
    for (i, cpath) in enumerate(cpaths)
        for j in cpath_sub_cpath_indices[i]
            difference = setdiff(cpath, cpaths[j]) # get the indices of the missing dots from the subcpath
            if costs[j] > costs[i] && all([dot_path_list[missing_dot] ⊆ cpath_nbr_cpath_indices[j] for missing_dot in difference])
                # the jth sub codepath is redundant, we can remove it.
                push!(redundant_subcodepath_indices, j)
                #println("j: $j, ", cpaths[j])
                #println("i: $i, ", cpaths[i])
            else # find the dots that conflict with the whole path, but not some subpath, and the subpaths they do not conflict with
                for missing_dot in difference
                    for k in dot_path_list[missing_dot] # for each index, k, of a cpath that passes through the missing dot
                        if isempty(cpaths[j] ∩ cpaths[k])
                            #ToDO: before supporting codes that allow for two drops, need to handle more cases here
                            #println("Partial Conflict i: ", cpaths_df.cpath[i], ", j: ", cpaths_df.cpath[j])
                            push!(partial_conflicts[k], i)
                            #println("Transition k: ", cpaths_df.cpath[k], ", j: ", cpaths_df.cpath[j])
                            push!(partial_conflict_transitions[k], j)
                            #filter!(x -> x != i, cpath_nbr_cpath_indices[k])
                        end
                    end
                end
            end
        end
    end

    #remove redundant codepaths, and rename indices of the codepaths with larger indices
    sort!(redundant_subcodepath_indices, rev = true)
    redundant_subcodepath_indices = unique(redundant_subcodepath_indices) # this is necessary if two dots in the same hyb are close enough to compete for a space in a codepath.
    for rsbcp_ind in redundant_subcodepath_indices
        remove_redundant_cpath_update_nbr_inds!(cpath_nbr_cpath_indices, rsbcp_ind)
        remove_redundant_cpath_update_nbr_inds!(partial_conflicts, rsbcp_ind)
        remove_redundant_cpath_update_nbr_inds!(partial_conflict_transitions, rsbcp_ind)
        delete!(cpaths_df, rsbcp_ind)
    end

    return cpath_nbr_cpath_indices, partial_conflicts, partial_conflict_transitions#, cpath_sub_cpath_indices




    #println("Neighbor Dictionary: ")
    #println(cpath_neighbors_dict)
    #println()
    #=
    #preallocate list of edges
    edge_list = Edge[]
    sizehint!(edge_list, nedges)

    # get list of edges
    for (cpath, nbr_dict) in enumerate(cpath_neighbors_dict)
        for nbr in keys(nbr_dict)
            if nbr_dict[nbr]
                #println("found conflit between edge $cpath: ", cpaths[cpath], " and $nbr: ", cpaths[nbr])
                push!(edge_list, Edge(cpath, nbr))
            end
        end
    end=#


    #edgelist = Edge[]
    #sizehint!(edgelist, sum((dot_ncpaths.^2)))
    #=
    println("Building Adjacency Graph")

    for i = 1:n_cpaths
        for dt in cpaths[i]
            for j in dot_path_list[dt]
                #if (j, i) ∉ edgelist #edges(A)
                #if (i, j) ∉ edges(A)
                add_edge!(A, i, j)
                #push!(edgelist, Edge(i, j))
                #end
            end
            push!(dot_path_list[dt], i)
        end
    end
    =#

    #println("build graph from iterator")
    #println(edge_list)
    #A = SimpleGraphFromIterator(edge_list)


    #return A, ccs
end


function get_cpath_conflict_graph_remove_redundant_cpaths_rA!(cpaths_df, ndots, n)# :: Vector{Vector{Int}})
    cpath_nbr_cpath_indices, partial_conflicts, partial_conflict_transitions = get_cpath_conflict_graph_remove_redundant_cpaths!(cpaths_df, ndots, n)
    ncpaths= length(cpath_nbr_cpath_indices)

    A = falses(ncpaths, ncpaths)
    for i in 1:ncpaths, j in cpath_nbr_cpath_indices[i]
        A[i,j] = 1
    end
    return A
end

function remove_redundant_cpath_update_nbr_inds!(nbr_array, redundant_codepath_ind)
    deleteat!(nbr_array, redundant_codepath_ind)
    for i in 1:length(nbr_array)
        j = 1
        while j <= length(nbr_array[i])
            if nbr_array[i][j] == redundant_codepath_ind
                deleteat!(nbr_array[i], j)
            elseif nbr_array[i][j] > redundant_codepath_ind
                nbr_array[i][j] -= 1
                j += 1
            else
                j += 1
            end
        end
    end
end

function get_cpath_conflict_graph2(cpaths_in :: DataFrame, n, ndrops)# :: Vector{Vector{Int}})
    println("Building Adjacency Graph")
    #println("cpaths:")
    #println(cpaths)
    cpaths = deepcopy(cpaths_in)
    n_cpaths = nrow(cpaths)

    println("n_cpaths: $n_cpaths")
    A = SimpleGraph(n_cpaths)

    cpaths["i"]= Array(1:n_cpaths)

    start_dot_ind = 0
    last_dot = 0
    for i in 1:(n-ndrops)
        sort!(cpaths, :cpath)
        for (j, r) in enumerate(eachrow(cpaths))
            #println("j: $j, r: $r")
            #println(r.cpath)
            dotj = popfirst!(r.cpath)
            if dotj == last_dot
                for k = start_dot_ind:j
                    if (k,j) ∉ edges(A)
                        add_edge!(A, k, j)
                    end
                end
            else
                last_dot = dotj
                start_dot_ind = j
            end
        end
    end

end

function get_cpath_conflict_graph_pairwise(cpaths)
    println("Building Adjacency Graph")
    #println("cpaths:")
    #println(cpaths)
    n_cpaths = length(cpaths)
    println("n_cpaths: $n_cpaths")
    A = SimpleGraph(n_cpaths)

    for i = 1:n_cpaths, j = (i+1):n_cpaths # ToDo: use smarter more efficient algorithm later
        for dt in cpaths[i]
            if dt in cpaths[j]
                add_edge!(A, i, j)
                break
            end
        end
    end

    println("Done building adjacency Graph")
    return A
end

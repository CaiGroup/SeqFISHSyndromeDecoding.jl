using DelimitedFiles
using DataFrames
using Graphs
import Graphs.neighbors
using NearestNeighbors
using DataStructures
using Statistics
using Clustering
using JuMP
using GLPK
#using Revise

"""
    decode_syndromes!(
        pnts :: DataFrame,
        cb,
        H :: Matrix,
        params :: DecodeParams
        optimizer = GLPK.Optimizer
    )

Arguments
- `pnts`: DataFrame of seqFISH psfs. Must include columns:
	- `x` : x spatial coordinate of psfs
	- `y` : y spatial coordinate of psfs
    - `z` : z spatial coordinate of psfs
	- `s` : the sigma width parameter of the psfs
    - `w` : the weight (or brightness) of the psfs
    Additionally, there the data frame must either have columns
    - `round` : the barcoding round in which the psf was found
    - `pseudocolor` : the pseudocolor of the barcoding round in which the psf was found
    or the round and pseudocolor can be computed from the hybridization
    - `hyb` : the hybridization in which the dot was found
    where round = ceil(hyb / q), pseudocolor = (hyb - (round - 1) * q) % q, and q is the number of pseudocolors.

- `cb` : The codebook.
- `H` : The parity check Matrix
- `params` : DecodeParams object holding the parameters for decoding
- `optimizer` : solver for integer linear programming optimizations. Uses the open source GLPK optimizer by default, but allows faster commercial optimizers to be used if necessary.
A list of supported solvers is available [here](https://jump.dev/JuMP.jl/stable/installation/#Supported-solvers)


Takes a DataFrame of aligned points with hyb, x, y, z, w, and s columns;
codebook matix; parity check matrix (H); and DecodeParams stucture with decoding parameters set.
Adds columns to input pnts matrix indicating the encoding round, syndrome component value,
and decoding result indicated as the row of the codebook matrix matched to.
Split up points into weakly connected components, then finds possible codeword messages
and runs simulated annealing to assign them. The pnts dataframe should have hybridization, x, y, and z columns
"""
function decode_syndromes!(pnts :: DataFrame, cb, H :: Matrix, params :: DecodeParams; optimizer = GLPK.Optimizer, tforms = nothing, obj_func = obj_function)
    #println("start syndrome decoding")
    if tforms == nothing # pnts preregistered
        cpath_df = get_codepaths(pnts, cb, H, params)
    else
        cpath_df = get_codepaths(pnts, cb, H, params, tforms)
    end

    if typeof(cpath_df) != DataFrame || nrow(cpath_df) == 0
        println("No viable barcodes")
        return
    end

    return choose_optimal_codepaths(pnts, cb, H, params, cpath_df, optimizer, tforms=tforms, obj_func=obj_function)
end

"""
function decode_syndromes!(pnts :: DataFrame, cb_df :: DataFrame, H :: Matrix, params :: DecodeParams)

    cb = Matrix(UInt8.(cb_df[!, 2:end]))
    gene = cb_df[!, "Gene"]
    value = Array(1:length(gene))
    decoded = decode_syndromes!(pnts :: DataFrame, cb :: Matrix, H :: Matrix, params :: DecodeParams)
    gene_value_df = DataFrame("Gene" => gene, "value" => value)
    decoded_joined = rightjoin(gene_value_df, decoded, on=:value)

end
"""

function obj_function(cpath, pnts, cw_w, params, tforms :: Nothing)
    lat_var_factor = params.lat_var_factor
    z_var_factor = params.z_var_factor
    lw_var_factor = params.lw_var_factor
    s_var_factor = params.s_var_factor
    dot_erasure_penalty = params.erasure_penalty

    lat_var_cost = (var(pnts.x[cpath]) + var(pnts.y[cpath])) * lat_var_factor
    z_var_cost = var(pnts.z[cpath]) * z_var_factor
    lw_var_cost = var(log2.(pnts.w[cpath])) * lw_var_factor
    s_var_cost = var(pnts.s[cpath]) * s_var_factor

    erasure_cost = (cw_w - length(cpath)) * dot_erasure_penalty

    return lat_var_cost + z_var_cost + lw_var_cost + s_var_cost + erasure_cost
end

function obj_function(cpath, pnts, cw_w, params, tforms :: Dict)
    lat_var_factor = params.lat_var_factor
    z_var_factor = params.z_var_factor
    lw_var_factor = params.lw_var_factor
    s_var_factor = params.s_var_factor
    dot_erasure_penalty = params.erasure_penalty

    cpath_pnts = pnts[cpath, :]
    lat_dists_sq = []
    zdists_sq = []
    for i in 1:(nrow(cpath_pnts) - 1), j in (i+1):nrow(cpath_pnts)
        tform = get_tform(tforms, cpath_pnts.pos[i], cpath_pnts.coeff[i], cpath_pnts.pos[j], cpath_pnts.coeff[j])
        dist = (Array(cpath_pnts[j, [:x, :y, :z]]) .- (tform[:,1:3] * Array(cpath_pnts[i, [:x, :y, :z]]) .+ tform[:,4])).^2
        push!(lat_dists_sq, sum(dist[1:2]))
        push!(zdists_sq, dist[3])
    end

    lat_var_cost = mean(lat_dists_sq) * lat_var_factor
    z_var_cost = mean(zdists_sq) * z_var_factor
    lw_var_cost = var(log2.(pnts.w[cpath])) * lw_var_factor
    s_var_cost = var(pnts.s[cpath]) * s_var_factor

    erasure_cost = (cw_w - length(cpath)) * dot_erasure_penalty

    return lat_var_cost + z_var_cost + lw_var_cost + s_var_cost + erasure_cost
end


"""
"""
function check_inputs(pnts :: DataFrame, cb :: Matrix, H :: Matrix, params :: DecodeParams)
    alphabet = sort(unique(cb))
    q = UInt8(length(alphabet))
    n = length(cb[1,:])
    if ~params.zeros_probed && ~("round" in names(pnts)) && ~("pseudocolor" in names(pnts))
        println("names(pnts): ", names(pnts))
        error("'round' and 'pseudocolor' columns must be included to decode experiments where zeros are not probed.")
    end
    @assert alphabet[1] == 0x00 || alphabet[1] == "0"
    if ~(typeof(alphabet[1]) <: AbstractString)
        @assert alphabet[q] < q
        @assert all(alphabet .>= 0)
    end
    if "hyb" in names(pnts)
        @assert maximum(pnts[!,"hyb"]) <= q*n
        @assert minimum(pnts[!,"hyb"]) > 0
    else
        @assert "round" in names(pnts)
        @assert "pseudocolor" in names(pnts)
    end
    if ~params.zeros_probed
        @assert typeof(cb) == typeof(H)
        if size(H)[1] < 2*params.ndrops
            error("Reed-Solomon Codes require 2 parity check symbols for every error corrected. Your code has ", size(H)[1], " parity check symbols, but you requested correction of up to ", params.ndrops, " errors.")
        end
    end
    if "round" in names(pnts)
        @assert maximum(pnts.round) <= n
    end
    if "pseudocolor" in names(pnts)
        if params.zeros_probed
            @assert maximum(pnts.pseudocolor) <= q
        else
            @assert maximum(pnts.pseudocolor) < q
        end
    end

    @assert params.ndrops <= size(H)[1]
    @assert params.ndrops >= 0
end

function sort_readouts!(pnts :: DataFrame)
    if "round" in names(pnts) && "pseudocolor" in names(pnts)
        sort!(pnts, [:round, :pseudocolor])
    elseif "hyb" in names(pnts)
        sort!(pnts, :hyb)
    else
        error("no way to sort points")
    end
end

"""
    get_codepaths(pnts :: DataFrame, cb_df :: DataFrame, H :: Matrix, params :: DecodeParams)

Arguments
- `pnts`: DataFrame of seqFISH psfs. Must include columns:
	- `x` : x spatial coordinate of psfs
	- `y` : y spatial coordinate of psfs
    - `z` : z spatial coordinate of psfs
	- `s` : the sigma width parameter of the psfs
    - `w` : the weight (or brightness) of the psfs
    Additionally, there the data frame must either have columns
    - `round` : the barcoding round in which the psf was found
    - `pseudocolor` : the pseudocolor of the barcoding round in which the psf was found
    or the round and pseudocolor can be computed from the hybridization
    - `hyb` : the hybridization in which the dot was found
    where round = ceil(hyb / q), pseudocolor = (hyb - (round - 1) * q) % q, and q is the number of pseudocolors.

- `cb` : The codebook.
- `H` : The parity check Matrix
- `params` : DecodeParams object holding the parameters for decoding

Computes codepaths with syndrome decoding, removes codepaths that exceed the cost
of not decoding their component dots, and
and returns DataFrame of candidate codepaths.
"""
function get_codepaths(pnts :: DataFrame, cb_df :: DataFrame, H :: Matrix, params :: DecodeParams, tforms=nothing)
    if typeof(cb_df[2, 2]) <: AbstractString
        cb = Matrix(string.(cb_df[!, 2:end]))
    else
        cb = Matrix(UInt8.(cb_df[!, 2:end]))
    end
    cpaths = get_codepaths(pnts :: DataFrame, cb :: Matrix, H :: Matrix, params :: DecodeParams, tforms)
    if typeof(cpaths) == DataFrame
        gene = cb_df[!, 1]
        gene_number = Array(1:length(gene))
        gene_df = DataFrame("gene_name" => gene, "gene_number" => gene_number)
        println(first(cpaths, 5))
        decoded_joined = rightjoin(gene_df, cpaths, on=:gene_number)
        return decoded_joined
    else
        return cpaths
    end
end

"""
    get_codepaths(pnts :: DataFrame, cb :: Matrix, H :: Matrix, params :: DecodeParams)

Computes codepaths with syndrome decoding, removes codepaths that exceed the cost
of not decoding their component dots, and returns DataFrame of candidate codepaths.

Arguments
- `pnts`: DataFrame of seqFISH psfs. Must include columns:
	- `x` : x spatial coordinate of psfs
	- `y` : y spatial coordinate of psfs
    - `z` : z spatial coordinate of psfs
	- `s` : the sigma width parameter of the psfs
    - `w` : the weight (or brightness) of the psfs
- `cb` : The codebook.
- `H` : The parity check Matrix
- `params` : DecodeParams object holding the parameters for decoding
"""
function get_codepaths(pnts :: DataFrame, cb :: Matrix, H :: Matrix, params :: DecodeParams, tforms=nothing)

    cb_dict = make_cw_dict(cb)
    alphabet = sort(unique(cb))
    q = UInt8(length(alphabet))
    #n = length(cb[1,:])
    if params.zeros_probed
        w = length(cb[1,:])
    else
        if typeof(cb[1,1]) <: AbstractString
            w = minimum(sum(cb .!= "0", dims=2))
        else
            w = minimum(sum(cb .!= 0, dims=2))
        end
    end

    if typeof(tforms) == DataFrame
        tforms_dict = get_tform_dict(tforms)
    else
        tforms_dict = nothing
    end

    if typeof(H[1,1]) <: AbstractString
        H = string.(H)
        cb = string.(cb)
    end

    check_inputs(pnts, cb, H, params)
    sort_readouts!(pnts)

    set_q(q)
    set_H(H, params, cb)
    #set_n(UInt8(n))
    #set_k(n-size(H)[1])
    if params.ndrops > 0 && params.zeros_probed
        get_decode_table()
    end

    #break into clusters
    if nrow(pnts) > 3
        pnts_mat = Array([pnts.x pnts.y pnts.z]')
        dbr = dbscan(pnts_mat, params.lat_thresh, min_neighbors=1, min_cluster_size=w-params.ndrops)
        if typeof(dbr) <: AbstractArray
            dbscan_clusters = [sort(vcat(dbc.core_indices,  dbc.boundary_indices)) for dbc in dbr]
        else
            dbscan_clusters = [sort(vcat(dbc.core_indices,  dbc.boundary_indices)) for dbc in dbr.clusters]
        end
    elseif nrow(pnts) == 3
        dbscan_clusters = [[1,2,3]]
    else
        dbscan_clusters = []
    end

    find_cluster_cpaths = function(dbscan_cluster)
        clust_pnts = pnts[dbscan_cluster, :]
        #println("cluster size: ", nrow(clust_pnts))
        sort_readouts!(clust_pnts)
        add_code_cols!(clust_pnts)
        g = DotAdjacencyGraph(clust_pnts, params, n, w, tforms_dict)
        #g = DotAdjacencyGraph(clust_pnts, params.lat_thresh, params.z_thresh, n, params.ndrops)

        cost(cpath) = obj_function(cpath, clust_pnts, w, params, tforms_dict)

        code_paths, gene_nums = syndrome_find_barcodes!(clust_pnts, g, cb, params.ndrops, w, tforms_dict)
        costs = cost.(code_paths)

        cpath_df = DataFrame(cpath = code_paths, cost = costs, gene_number = gene_nums)
        sort!(cpath_df, :cost)
        cpath_df = remove_high_cost_cpaths(cpath_df, params.free_dot_cost, w, params.ndrops)
        cpath_df = threshold_cpaths(cpath_df, clust_pnts, params.lat_thresh, params.z_thresh, tforms_dict)
        cpath_df[!,"x"] = mean.([clust_pnts.x[cpath] for cpath in cpath_df.cpath])
        cpath_df[!,"y"] = mean.([clust_pnts.y[cpath] for cpath in cpath_df.cpath])
        cpath_df[!,"z"] = mean.([clust_pnts.z[cpath] for cpath in cpath_df.cpath])
        replace!.(i->dbscan_cluster[i], cpath_df.cpath)

        return cpath_df
    end
    cpath_df = vcat(map(find_cluster_cpaths, dbscan_clusters)...)
    return cpath_df
end



"""
    choose_optimal_codepaths(pnts :: DataFrame, cb_df :: DataFrame, H :: Matrix, params :: DecodeParams, cpath_df :: DataFrame, optimzer, ret_discarded :: Bool=false)

Arguments
- `pnts`: DataFrame of seqFISH psfs. Must include columns:
	- `x` : x spatial coordinate of psfs
	- `y` : y spatial coordinate of psfs
    - `z` : z spatial coordinate of psfs
	- `s` : the sigma width parameter of the psfs
    - `w` : the weight (or brightness) of the psfs
    Additionally, there the data frame must either have columns
    - `round` : the barcoding round in which the psf was found
    - `pseudocolor` : the pseudocolor of the barcoding round in which the psf was found
    or the round and pseudocolor can be computed from the hybridization
    - `hyb` : the hybridization in which the dot was found
    where round = ceil(hyb / q), pseudocolor = (hyb - (round - 1) * q) % q, and q is the number of pseudocolors.

- `cb` : The codebook.
- `H` : The parity check Matrix
- `params` : DecodeParams object holding the parameters for decoding
- `cpath_df` : A DataFrame output from the `get_codepaths` function, defined above.
- `optimizer` : solver for integer linear programming optimizations. A list of supported solvers is available [here](https://jump.dev/JuMP.jl/stable/installation/#Supported-solvers)
- `ret_discarded` : if true, return data frame of candidate codepaths that were discarded for being in too large/too dense of a conflict network

Choose best codepaths from previouly found candidates that may have been found with less strict parameters. Reevaluates the costs for each candidate and trims according
to the passed parameters.

"""
function choose_optimal_codepaths(pnts :: DataFrame, cb_df :: DataFrame, H :: Matrix, params :: DecodeParams, cpath_df :: DataFrame, optimizer; ret_discarded :: Bool=false, tforms=nothing, obj_func=obj_function)
    if any(typeof(cb_df[2:end, 2:end]) .<: AbstractString) #typeof(cb_df[2, 2]) <: AbstractString
        cb = Matrix(string.(cb_df[!, 2:end]))
    else
        cb = Matrix(UInt8.(cb_df[!, 2:end]))
    end
    #gene = cb_df[!, 1]
    #gene_number = Array(1:length(gene))
    decoded, discarded_cpaths = choose_optimal_codepaths(pnts, cb, H, params, cpath_df, optimizer, tforms=tforms)
    #gene_df = DataFrame("gene_name" => gene, "gene_number" => gene_number)
    #decoded_joined = rightjoin(gene_df, decoded, on=:gene_number)
    if ret_discarded
        #return decoded_joined, discarded_cpaths
        return decoded, discarded_cpaths
    else
        #return decoded_joined
        return decoded
    end
end

function choose_optimal_codepaths(pnts :: DataFrame, cb :: Matrix, H :: Matrix, params :: DecodeParams, cpath_df :: DataFrame, optimizer; tforms=nothing, obj_func=obj_function)
    alphabet = sort(unique(cb))
    q = UInt8(length(alphabet))
    set_q(q)
    set_H(H, params, cb)
    sort_readouts!(pnts)
    #n = length(cb[1,:])
    ndots = nrow(pnts)
    if params.zeros_probed
        w = length(cb[1,:])
        filter!(:cpath => cpath -> length(cpath) >= w - params.ndrops, cpath_df)
    else
        w = minimum(sum(cb .!= "0" .&& cb .!= 0, dims=2))
        filter!(:cpath => cpath -> length(cpath) >= w - params.ndrops, cpath_df)
    end

    if typeof(tforms) == DataFrame
        tforms_dict = get_tform_dict(tforms)
    else
        tforms_dict = nothing
    end

    cost(cpath) = obj_func(cpath, pnts, w, params, tforms_dict)

    pnts[!,"decoded"] = fill(0, nrow(pnts))
    pnts[!, "mpath"] = [[] for i = 1:length(pnts.x)]
    add_code_cols!(pnts)
    cpath_df[!, "cost"] = cost.(cpath_df[!, "cpath"])
    sort!(cpath_df, :cost)
    cpath_df = remove_high_cost_cpaths(cpath_df, params.free_dot_cost, w, params.ndrops)
    cpath_df = threshold_cpaths(cpath_df, pnts, params.lat_thresh, params.z_thresh, tforms_dict)

    # build graph with by adding only edges in codepaths, and break into connected
    # components
    ccs = get_connected_components(cpath_df.cpath, nrow(pnts))
    cpath_df[!, "cc"] .= 0
    cpath_df[!, "cc_size"] .= 0

    n_ccs = length(ccs)
    #println("n_ccs: ", n_ccs)

    mpaths = cpath_df[1:0, :]
    dense_cpaths = cpath_df[1:0, :]
    nmpaths = 0
    for (cc_i, cc) in enumerate(ccs)
        cpath_df[cc, "cc"] .= cc_i
        cpath_df[cc, "cc_size"] .= length(cc)
        cc_cpath_df = cpath_df[cc,:]
        area = (maximum(cc_cpath_df[:,"x"]) - minimum(cc_cpath_df.x))*(maximum(cc_cpath_df.y) - minimum(cc_cpath_df.y))

        ndots_cc = length(cc)
        if nrow(cc_cpath_df) == 1
            low_cost_state = [1]
        # ToDO: test 2 cpath cc case
        elseif nrow(cc_cpath_df) == 2
            costs = cc_cpath_df[!,"cost"]
            low_cost_state = (costs .== minimum(costs))
        elseif nrow(cc_cpath_df)/area > params.skip_thresh && nrow(cc_cpath_df) > params.skip_density_thresh
            println("skip ", cc_i, " size ", length(cc), " area: $area")
            append!(dense_cpaths, cc_cpath_df)
            continue
        else
            cpath_nbrs, cpath_partial_conflicts, cpath_partial_conflict_transitions = get_cpath_conflict_graph_remove_redundant_cpaths!(cc_cpath_df, ndots, n)

            # get heuristic start point?
            if nrow(cc_cpath_df) < params.mip_sa_thresh
                low_cost_state = MIP_solve(cc_cpath_df, cpath_nbrs, params.free_dot_cost, optimizer)
            else
                low_cost_state = simulated_annealing(cc_cpath_df, cpath_nbrs, cpath_partial_conflicts, cpath_partial_conflict_transitions, params, ndots_cc)
            end
        end

        mpath_df = cc_cpath_df[low_cost_state,:]
        nmpaths += nrow(mpath_df)
        append!(mpaths, mpath_df)

        for mpath_row in eachrow(mpath_df)
            pnts.decoded[mpath_row.cpath] .= mpath_row.gene_number
            for dt in mpath_row.cpath
                pnts.mpath[dt] = mpath_row.cpath
            end
        end
    end
    #println("found $nmpaths mpaths")
    return mpaths, dense_cpaths
end



"""
    make_KDTree(pnts :: DataFrame)

Generate KDTree to aid in building adjacency graphs.
"""
make_KDTree2D(pnts :: DataFrame) = KDTree(Float64[pnts.x pnts.y]')

function make_KDTree2D(pnts :: Matrix)
    unregistered_pnts = Array([pnts.x pnts.y])
    registered_pnts = unregistered_pnts * tform[1:2, 1:2] .+ tform[1:2,4]
    return KDTree(Array([registered_pnts]'))
end

"""
    make_KDTree3D(pnts :: DataFrame)

Generate KDTree to aid in building adjacency graphs.
"""
function make_KDTree3D(pnts :: DataFrame, lat_thresh, z_thresh)
    z_scaled = pnts.z .* lat_thresh ./ z_thresh
    return KDTree(Array(Float64[pnts.x pnts.y z_scaled]'))
end

function make_KDTree3D(pnts :: Matrix, lat_thresh, z_thresh)
    z_scaled = pnts[:,3] .* lat_thresh ./ z_thresh
    pnts[:,3] .= z_scaled
    return KDTree(Float64.(Array(pnts')))
end

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
    if "round" in names(pnts)
        pnts.pos = UInt8.(pnts.round)
    else
        pnts.pos = get_pos.(pnts.hyb)
    end

    if "pseudocolor" in names(pnts)
        if q == 8 || q == 9
            pnts.coeff = map(c -> pseudocolor_2_savestring[c], pnts.pseudocolor)
        else
            pnts.coeff = UInt8.(pnts.pseudocolor)
        end
    else
        pnts.coeff = get_coeff.(pnts.hyb, pnts.pos)
    end

    pnts.sc = SyndromeComponent.(pnts.coeff, pnts.pos)
    pnts.decoded = zeros(length(pnts.x))
    pnts.mpath = [[] for i = 1:length(pnts.x)]
end

abstract type abstractDotAdjacencyGraph end

abstract type DotAdjacencyGraph <: abstractDotAdjacencyGraph end

abstract type DotAdjacencyGraph2D <: DotAdjacencyGraph end

abstract type DotAdjacencyGraph3D <: DotAdjacencyGraph end

"""
    DotAdjacencyGraph(g :: SimpleDiGraph
                      cw_pos_bnds :: Tuple{Int64}
                      n :: Int8
                      q :: Int8)

Structure for storing the dot adjacency graph with some parameters
"""
struct DotAdjacencyGraphRegistered2D <: DotAdjacencyGraph2D
    g :: SimpleDiGraph
    cw_pos_bnds :: Array{Int64}
    n :: Int8
    trees :: Vector{KDTree}
    lat_thresh :: Float64
    pnts :: DataFrame
    ndrops :: Int64
end

struct DotAdjacencyGraphPairwise2D <: DotAdjacencyGraph2D
    g :: SimpleDiGraph
    cw_pos_bnds :: Array{Int64}
    n :: Int8
    trees :: Matrix{KDTree}
    lat_thresh :: Float64
    pnts :: DataFrame
    ndrops :: Int64
    round_pc_block_starts :: Matrix
end

struct DotAdjacencyGraphRegistered3D <: DotAdjacencyGraph3D
    g :: SimpleDiGraph
    cw_pos_bnds :: Array{Int64}
    n :: Int8
    trees :: Vector{KDTree}
    lat_thresh :: Float64
    z_thresh :: Float64
    pnts :: DataFrame
    ndrops :: Int64
end

struct DotAdjacencyGraphPairwise3D <: DotAdjacencyGraph3D
    g :: SimpleDiGraph
    cw_pos_bnds :: Array{Int64}
    n :: Int8
    trees :: Matrix{KDTree}
    lat_thresh :: Float64
    z_thresh :: Float64
    pnts :: DataFrame
    ndrops :: Int64
    round_pc_block_starts :: Matrix
end

function DotAdjacencyGraph(pnts :: DataFrame, params :: DecodeParams, n, w, tforms=nothing)
    if params.zeros_probed
        return DotAdjacencyGraph(pnts, params.lat_thresh, params.z_thresh, n, params.ndrops, tforms=tforms)
    else
        return DotAdjacencyGraphBlankRound(pnts, params.lat_thresh, params.z_thresh, n, params.ndrops, w, tforms)
    end
end

"""
    DotAdjacencyGraph(g :: SimpleDiGraph
                      cw_pos_bnds :: Tuple{Int64}
                      n :: Int8
                      q :: Int8)

Construct a dot adjacency graph where dots close to each other have directed
edges pointing towards the dot representing an earlier symbor
"""
function DotAdjacencyGraph(pnts :: DataFrame, lat_thresh :: Real, z_thresh :: Real, n, ndrops; tforms=nothing)
    g = SimpleDiGraph(nrow(pnts))

    data_2d = length(unique(pnts.z)) == 1
    # Find the indices of dots representing each place, cᵢ, in a codeword start.
    cw_pos_bnds = get_cw_pos_bounds(pnts, n)
    if isnothing(tforms)
        trees = []
    else
        trees = Matrix{KDTree}(undef, n, q)
        round_pseudocolor_start_positions = Matrix{Union{Int64, Nothing}}(undef, n, q)
    end

    for round in 1:n
        start_round = maximum([round-1-ndrops, 1])
        start_pnt = cw_pos_bnds[start_round]
        end_pnt = (cw_pos_bnds[round]-1)
        if isnothing(tforms) 
            if data_2d
                push!(trees, make_KDTree2D(pnts[start_pnt:end_pnt, :]))
            else
                push!(trees, make_KDTree3D(pnts[start_pnt:end_pnt, :], lat_thresh, z_thresh))
            end
        else
            for searching_pseudocolor in 0:(q-1)
                registered_pnts = copy(pnts[start_pnt:end_pnt, :])
                for nbr_round in start_round:(round-1)
                    for nbr_round_pseudocolor in 0:(q-1)
                        tform = get_tform(tforms, nbr_round, nbr_round_pseudocolor, round, searching_pseudocolor)
                        rows_of_interest = registered_pnts.pos .== nbr_round .&& registered_pnts.coeff .== nbr_round_pseudocolor
                        nbr_rnd_pc_pnts = registered_pnts[rows_of_interest, :]
                        registered_pnts[rows_of_interest, [:x, :y, :z]] .= (Array(tform[:, 1:3] * Array(nbr_rnd_pc_pnts[:, [:x, :y, :z]])') .+ tform[:, 4])'
                    end
                end
                spc_ind = searching_pseudocolor == 0 ? q : searching_pseudocolor
                if data_2d
                    #push!(trees, KDTree(registered_pnts'))
                    trees[round, spc_ind] = KDTree(Matrix(registered_pnts[:,[:x, :y]])')
                    #trees[round, spc_ind] = KDTree(transformed_coords')
                else
                    #push!(trees, MakeKDTree3D(registered_pnts, lat_thresh, z_thresh))
                    trees[round, spc_ind] = MakeKDTree3D(registered_pnts, lat_thresh, z_thresh)
                end
                #round_pseudocolor_start_positions[round, spc_ind] = start_pnt - 1 #findfirst(d -> d.pos==round && d.coeff ==searching_pseudocolor, eachrow(pnts))

                #get the adjustment factor (global start index minus 1) for dots in each pseudocolor-round.
                round_pseudocolor_start_positions[round, spc_ind] = findfirst(d -> d.pos==round && d.coeff ==searching_pseudocolor, eachrow(pnts))
                if ~isnothing(round_pseudocolor_start_positions[round, spc_ind])
                    round_pseudocolor_start_positions[round, spc_ind] -= 1
                end
            end
        end
    end

    if isnothing(tforms)
        if data_2d
            return DotAdjacencyGraphRegistered2D(g, cw_pos_bnds, n, trees, lat_thresh, pnts, ndrops)
        else
            return DotAdjacencyGraphRegistered3D(g, cw_pos_bnds, n, trees, lat_thresh, z_thresh, pnts, ndrops)
        end
    else
        if data_2d
            return DotAdjacencyGraphPairwise2D(g, cw_pos_bnds, n, trees, lat_thresh, pnts, ndrops, round_pseudocolor_start_positions)
        else
            return DotAdjacencyGraphPairwise3D(g, cw_pos_bnds, n, trees, lat_thresh, z_thresh, pnts, ndrops, round_pseudocolor_start_positions)
        end
    end
end

function get_tform_dict(tforms :: DataFrame)
    return Dict([((rsrc=row.r_src, pcsrc=row.pc_src, rdst=row.r_dst, pcdst=row.pc_dst), row.tform) for row in eachrow(tforms)])
end

function get_tform(tforms :: Dict, src_round, src_pseudocolor, dst_round, dst_pseudocolor)
    #return tforms[(tforms.r_src .== nbr_round) .&& (tforms.pc_src .== (nbr_round_pseudocolor .% q)) .&& (tforms.r_dst .== searching_round) .&& (tforms.pc_dst .== (searching_pseudocolor .% q)), "tform"][1]
    return tforms[(rsrc=src_round, pcsrc = src_pseudocolor, rdst = dst_round, pcdst = dst_pseudocolor)]
end

abstract type DotAdjacencyGraphBlankRound <: abstractDotAdjacencyGraph end


struct DotAdjacencyGraphBlankRound2D <: DotAdjacencyGraphBlankRound
    g :: SimpleDiGraph
    cw_round_ranges
    n :: Int8
    trees :: Vector{KDTree}
    lat_thresh :: Float64
    pnts :: DataFrame
    ndrops :: Int64
    w :: Int64
    first_potential_barcode_final_dot
end

struct DotAdjacencyGraphBlankRound3D <: DotAdjacencyGraphBlankRound
    g :: SimpleDiGraph
    cw_round_ranges
    n :: Int8
    trees :: Vector{KDTree}
    lat_thresh :: Float64
    z_thresh :: Float64
    pnts :: DataFrame
    ndrops :: Int64
    w :: Int64
    first_potential_barcode_final_dot
end

function DotAdjacencyGraphBlankRound(pnts :: DataFrame, lat_thresh :: Real, z_thresh :: Real, n, ndrops, w, tforms=nothing)

    g = SimpleDiGraph(nrow(pnts))
    data_2d = length(unique(pnts.z)) == 1
    if data_2d
        make_KDTree = make_KDTree2D
    else
        make_KDTree(df) = make_KDTree3D(df, lat_thresh, z_thresh)
    end

    # Find the indices of dots representing each place, cᵢ, in a codeword start.
    cw_round_ranges = get_cw_round_ranges(pnts, n)

    #make KDTrees to search for dots that may be neigbors to dots each barcoding round that 
    trees = []

    # Dots in any previous round may be neighbors with dots in rounds up to the n-w+1st round in paths that may represent a decodable barcode
    for round in 1:(n-w+1) 
        if ismissing(cw_round_ranges[round]) # if there are no dots in the barcoding round
            push!(trees, make_KDTree(pnts[1:0, :])) #make empty KBTree
        else
            # make KDTree that searches dots in all previous rounds
            end_pnt = (cw_round_ranges[round][1]-1)
            push!(trees, make_KDTree(pnts[1:end_pnt, :])) 
        end
    end

    # dots in the n-w+2st or greater barcoding round cannot form paths in the DAG representing valid barcodes with early dots 
    # since no path containing directed edges will be long enough (paths must have length of at least w-ndrops to be decodable)
    for (i, round) in enumerate((n-w+2):n) 
        if ismissing(cw_round_ranges[round]) # if no dots in round
            push!(trees, make_KDTree(pnts[1:0, :])) # make empty KDTree
        else
            # make KDTree searching previous rounds starting at maximum([w-(n-round)-1-ndrops,1]), which may
            # be included in paths of decodable length when directed edges connect dots in both rounds
            start_pnt = find_previous_round_start(cw_round_ranges, maximum([w-(n-round)-1-ndrops,1]))
            end_pnt = (cw_round_ranges[round][1]-1)
            push!(trees, make_KDTree(pnts[start_pnt:end_pnt, :]))
        end
    end

    # note: row indices of dots in the pnts DataFrame are sorted by barcoding round
    # find the index of the first dot that may be the dot of highest barcoding round in a decodable path through the DAG.
    # all dots of higher index also may be the dot of highest barcoding round in a decodable path through the DAG.
    first_potential_barcode_final_dot=find_previous_round_start(cw_round_ranges,w-ndrops)
    if data_2d
        DotAdjacencyGraphBlankRound2D(g, cw_round_ranges, n, trees, lat_thresh, pnts, ndrops, w, first_potential_barcode_final_dot)
    else
        DotAdjacencyGraphBlankRound3D(g, cw_round_ranges, n, trees, lat_thresh, z_thresh, pnts, ndrops, w, first_potential_barcode_final_dot)
    end

end

function find_previous_round_start(cw_round_ranges, r0)
    for r in r0:length(cw_round_ranges)
        start = cw_round_ranges[r]
        ismissing(start) ? r += 1 : return start[1]
    end
end

"""
"""
function get_cw_pos_bounds(pnts, n)
    cw_pos_bnds = [1]
    for cᵢ = 1:(n-1)
        if findfirst(x -> x>cᵢ, pnts.pos) == nothing
            push!(cw_pos_bnds, nrow(pnts)+1)
        else
            push!(cw_pos_bnds, findfirst(x -> x>cᵢ, pnts.pos))
        end
    end
    push!(cw_pos_bnds,nrow(pnts)+1)
    return cw_pos_bnds
end

"""
"""
function get_cw_round_ranges(pnts, n)
    cw_round_ranges = []
    sizehint!(cw_round_ranges, n)

    for cᵢ = 1:n
        if findfirst(x -> x==cᵢ, pnts.round) == nothing
            push!(cw_round_ranges, missing)
        else
            push!(cw_round_ranges, findfirst(x -> x==cᵢ, pnts.round):findlast(x -> x==cᵢ, pnts.round))
        end
    end
    return cw_round_ranges
end

"""
    neighbors(g :: abstractDotAdjacencyGraph2D, n)

Define SimpleDiGraph neighbors function for DotAdjacencyGraph
"""

function neighbors(g :: DotAdjacencyGraphRegistered2D, n)
    nbrs = inrange(g.trees[g.pnts.pos[n]], [g.pnts.x[n], g.pnts.y[n]], g.lat_thresh, true)
    pnts_prior_rnds = g.cw_pos_bnds[maximum([g.pnts.pos[n]-1-g.ndrops, 1])] - 1
    return nbrs .+ pnts_prior_rnds
end

function neighbors(g :: DotAdjacencyGraphPairwise2D, dot)
    round = g.pnts.pos[dot]
    pseudocolor = g.pnts.coeff[dot] == 0 ? q : g.pnts.coeff[dot]
    nbrs = inrange(g.trees[round, pseudocolor], [g.pnts.x[dot], g.pnts.y[dot]], g.lat_thresh, true)
    #nbrs = inrange(g.trees[g.pnts.pos[n]], [g.pnts.x[n], g.pnts.y[n]], g.lat_thresh, true)
    pnts_prior_rnds = g.cw_pos_bnds[maximum([round-1-g.ndrops, 1])] - 1
    return nbrs .+ pnts_prior_rnds
end

"""
    neighbors(g :: abstractDotAdjacencyGraph3D, n)

Define SimpleDiGraph neighbors function for DotAdjacencyGraph
"""

function neighbors(g :: DotAdjacencyGraph3D, n)
    nbrs = inrange(g.trees[g.pnts.pos[n]], [g.pnts.x[n], g.pnts.y[n], g.pnts.z[n]*g.lat_thresh/g.z_thresh], g.lat_thresh, true)
    pnts_prior_rnds = g.cw_pos_bnds[maximum([g.pnts.pos[n]-1-g.ndrops, 1])] - 1
    return nbrs .+ pnts_prior_rnds
end

function neighbors(g :: DotAdjacencyGraphPairwise3D, dot)
    round = g.pnts.round[dot]
    pseudocolor = g.pnts.pseudocolor[dot]
    nbrs = inrange(g.trees[round, pseudocolor], [g.pnts.x[dot], g.pnts.y[dot], g.pnts.z[dot]*g.lat_thresh/g.z_thresh], g.lat_thresh, true)
    #nbrs = inrange(g.trees[g.pnts.pos[n]], [g.pnts.x[n], g.pnts.y[n]], g.lat_thresh, true)
    pnts_prior_rnds = g.cw_pos_bnds[maximum([round-1-g.ndrops, 1])] - 1
    return nbrs .+ pnts_prior_rnds
end

"""
    neighbors(g :: DotAdjacencyGraphBlankRound2D, n)

Define SimpleDiGraph neighbors function for DotAdjacencyGraph
"""
function neighbors(g :: DotAdjacencyGraphBlankRound2D, dot)
    nbrs = inrange(g.trees[g.pnts.pos[dot]], [g.pnts.x[dot], g.pnts.y[dot]], g.lat_thresh, true)
    return nbr_index_round_to_global!(nbrs, g, dot)
end

"""
    neighbors(g :: DotAdjacencyGraphBlankRound, n)

Define SimpleDiGraph neighbors function for DotAdjacencyGraph
"""
function neighbors(g :: DotAdjacencyGraphBlankRound3D, dot)
    nbrs = inrange(g.trees[g.pnts.pos[dot]], [g.pnts.x[dot], g.pnts.y[dot], g.pnts.z[dot]], g.lat_thresh, true)
    return nbr_index_round_to_global!(nbrs, g, dot)
end

function nbr_index_round_to_global!(nbrs, g, dot)
    if g.pnts.pos[dot] >  g.n - g.w + 2 + g.ndrops
        pnts_prior_rnds = find_previous_round_start(g.cw_round_ranges, g.w - g.ndrops - (g.n - g.pnts.round[dot]+1))
        nbrs .+= pnts_prior_rnds - 1
    end
    return nbrs
end

"""
    get_cw_pos_inds(g :: DotAdjacencyGraph, pos :: Int)

Helper function to get the dots in the graph that represent a symbol in a given position of the codeword.
"""
function get_cw_pos_inds(g :: DotAdjacencyGraph, pos :: Integer)
    return g.cw_pos_bnds[pos]:(g.cw_pos_bnds[pos+1]-1)
end

"""
"""
function syndrome_find_barcodes!(pnts ::DataFrame, g :: abstractDotAdjacencyGraph, cb ::Matrix, ndrops, w, tforms=nothing)
    cw_dict = make_cw_dict(cb)
    if typeof(g) <: DotAdjacencyGraphBlankRound
        cpaths, decode_cands = find_blank_round_codewords(pnts ::DataFrame, g :: DotAdjacencyGraphBlankRound, cw_dict, w, tforms)
    elseif ndrops == 0
        cpaths, decode_cands = find_barcodes_mem_eff(pnts, g, cw_dict, tforms)
    else
        cpaths, decode_cands = syndrome_find_barcodes!(pnts, g, ndrops, cw_dict, tforms)
    end
    return cpaths, decode_cands
end

"""
Computes syndrome of every path in the message graph, determines which ones represent valid
barcodes, and returns dataframe of message path candidates.
"""
function syndrome_find_barcodes!(pnts ::DataFrame,
                                      g :: DotAdjacencyGraph,
                                      ndrops :: Int,
                                      cw_dict :: Dict,
                                      tforms
                                      )
    syndromes, syndrome_coeff_positions = compute_syndromes(pnts, g)

    cpaths, decode_cands = find_code_paths!(
                    g,
                    pnts,
                    cw_dict,
                    syndromes,
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
    syndromes, syndrome_coeff_positions = init_syndromes(pnts, g)

    for cw_pos in 1:g.n
        cw_pos_inds = get_cw_pos_inds(g, cw_pos)

        # for each dot representing a candidate coefficient
        for dot in cw_pos_inds
            synd_ind = cw_pos <= (1 + g.ndrops) ? 2 : 1

            # get neighbors of that dot
            for neighbor in neighbors(g, dot)
                @inbounds end_ind = synd_ind+length(syndromes[neighbor])-1
                @inbounds syndromes[dot][synd_ind:end_ind] += syndromes[neighbor]
                @inbounds syndrome_coeff_positions[dot][synd_ind:end_ind] += syndrome_coeff_positions[neighbor]
                synd_ind = end_ind + 1
            end
        end
    end
    return syndromes, syndrome_coeff_positions
end

"""
    find_nsnds(g :: DotAdjacencyGraph)

Counts the number of paths in the dot adjacency graph that run through each dot to
determine the size of syndrome component arrays to preallocate.
"""
function find_nsnds(g :: DotAdjacencyGraph)
    n = g.n
    n_pnts = nrow(g.pnts)
    #nsnds = fill(1,n_pnts)
    nsnds = fill(0,n_pnts)
    nsnds[1:(g.cw_pos_bnds[2+g.ndrops]-1)] .= 1

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
    syndrome_coeff_positions = Vector{Vector{UInt8}}()
    sizehint!(syndromes, nv(g.g))
    for (pnt, sc) in enumerate(pnts.sc)
        push!(syndromes, fill(sc, nsnds[pnt]))
        pos_pow = UInt8(0x02 ^ (pnts.pos[pnt] - 0x01))
        push!(syndrome_coeff_positions, fill(pos_pow, nsnds[pnt]))
    end
    (syndromes, syndrome_coeff_positions)
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
        #for (synd_ind, path_length) in enumerate(syndrome_path_lengths[dot_ind])
        for (synd_ind, path_barcoding_rounds) in enumerate(syndrome_coeff_positions[dot_ind])
            if path_barcoding_rounds == full_bin_pos_indicator && iszero(syndromes[dot_ind][synd_ind])
                code_path = recursive_get_synd_neighbors(pnts, g, dot_ind, synd_ind, syndromes)
                @assert length(code_path) == (length(g.cw_pos_bnds) -1)
                message = pnts.coeff[code_path]
                if message in keys(cw_dict)
                    push!(cpaths, code_path)
                    push!(decode_cands, cw_dict[message])
                end
            else
                ndots = get_number_of_dots(path_barcoding_rounds, cw_n_symbols)
                if ndots >= cw_n_symbols - ndrops && ndots < cw_n_symbols
                    s = syndromes[dot_ind][synd_ind]
                    #dot_pos_sum = syndrome_coeff_positions[dot_ind][synd_ind]

                    # use bitwise masking to get the position of the missing dot
                    # positions are bitwise "one-hot" encoded from sums of powers of 2
                    drop_pos_pow = path_barcoding_rounds ⊻ full_bin_pos_indicator
                    if drop_pos_pow > 32
                        println("dot_pos_sum: $dot_pos_sum")
                    end

                    result = check_mpath_decodable(drop_pos_pow, s)

                    if result.decodable
                        code_path = recursive_get_synd_neighbors(pnts, g, dot_ind, synd_ind, syndromes)
                        @assert length(code_path) >= (length(g.cw_pos_bnds) -1 - g.ndrops)
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
    end
    return cpaths, decode_cands
end

"""
Helper function to get number of dots in path using bitwise operations
"""
function get_number_of_dots(pos_indicator, n_barcoding_rounds)
    ndots = 0
    for r in 1:n_barcoding_rounds
        bc_round_has_dot = (((2^(r-1)) & pos_indicator) > 0)
        ndots += bc_round_has_dot
    end
    return ndots
end

"""
Helper Function used to trace back groups of dots that produce zero syndrome, and therefore are codewords
"""
function recursive_get_synd_neighbors(
    pnts :: DataFrame,
     g :: DotAdjacencyGraph,
     dot :: Int,
     synd_ind :: Int,
     syndromes,
     recursion_depth  = 1:: Int
     )

    # if this is the last dot in the message, return number of dot in an array
    if synd_ind == 1 && dot < g.cw_pos_bnds[2+g.ndrops] && recursion_depth >= (length(g.cw_pos_bnds)-1-g.ndrops)
        cpath = Int[dot]
        sizehint!(cpath, g.n)
        return cpath
    end

    #otherwise, get neighbor of the dot that produced the zero syndrome
    neighbor, neighbor_synd_ind = get_synd_neighbor(g, dot, synd_ind, syndromes)

    # Add result to recursively defined array, and return
    push!(recursive_get_synd_neighbors(pnts, g, neighbor, neighbor_synd_ind, syndromes, recursion_depth+1), dot)
end

"""
Helper Function used to trace back groups of dots that produce zero syndrome, and therefore are codewords
"""
function get_synd_neighbor(
    g :: DotAdjacencyGraph,
    dot :: Int,
    synd_ind :: Int,
    syndromes
    )

    #@assert synd_ind != 1

    # keep track of number of syndromes from the dot that have been searched through
    #cum_n_syndromes = 1
    if dot < g.cw_pos_bnds[2+g.ndrops]
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
        cum_n_syndromes += length(syndromes[neighbor])

        #if the syndrome includes a component from this particular neighbor
        if synd_ind <= cum_n_syndromes

            # find the index of the neighbors syndrome component of interest
            neighbor_synd_ind = synd_ind - cum_n_syndromes + length(syndromes[neighbor])
            return [neighbor, neighbor_synd_ind]
        end
    end
    error("We shouldn't have gotten here!")
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
    i = searchsorted(cpath_df.cost, (n-ndrops)*free_dot_cost).start
    n_erasure_paths_removed = 0

    while i <= nrow(cpath_df)
        if cpath_df.cost[i] >= length(cpath_df.cpath[i])*free_dot_cost
            deleteat!(cpath_df, i)
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
function threshold_cpaths(cpaths_df, pnts, lat_thresh, z_thresh, tforms :: Nothing)
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
            deleteat!(cpaths_df, row)
        else
            row += 1
        end
    end
    cpaths_df
end

function threshold_cpaths(cpaths_df, pnts, lat_thresh, z_thresh, tforms :: Dict)
    row = 1
    while row <= nrow(cpaths_df)
        xs = pnts.x[cpaths_df.cpath[row]]
        ys = pnts.y[cpaths_df.cpath[row]]
        zs = pnts.z[cpaths_df.cpath[row]]
        coords = pnts[cpaths_df.cpath[row], [:x, :y, :z]]
        rounds = pnts.pos[cpaths_df.cpath[row]]
        pseudocolors = pnts.coeff[cpaths_df.cpath[row]]
        exceeds_threshold = false
        len_cp = length(xs)
        for i = 1:(len_cp-1), j = (i+1):len_cp
            tform = get_tform(tforms, rounds[i], pseudocolors[i], rounds[j], pseudocolors[j])
            registered_coords =  tform[:, 1:3] * Array(coords[i, :]) + tform[:,4]
            z_diff = abs(zs[j] - registered_coords[3])
            lat_diff = sqrt((xs[j]-registered_coords[1])^2 + (ys[j]-registered_coords[2])^2)
            if lat_diff > lat_thresh || z_diff > z_thresh
                exceeds_threshold = true
                break
            end
        end
        if exceeds_threshold
            deleteat!(cpaths_df, row)
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

    #check if cpath_sub_cpath_indices are redundant
    redundant_subcodepath_indices = []
    partial_conflicts = [Int[] for i = 1:n_cpaths]
    partial_conflict_transitions = [Int[] for i = 1:n_cpaths]
    for (i, cpath) in enumerate(cpaths)
        for j in cpath_sub_cpath_indices[i]
            difference = setdiff(cpath, cpaths[j]) # get the indices of the missing dots from the subcpath
            if costs[j] > costs[i] && all([dot_path_list[missing_dot] ⊆ cpath_nbr_cpath_indices[j] for missing_dot in difference])
                # the jth sub codepath is redundant, we can remove it.
                push!(redundant_subcodepath_indices, j)
            else # find the dots that conflict with the whole path, but not some subpath, and the subpaths they do not conflict with
                for missing_dot in difference
                    for k in dot_path_list[missing_dot] # for each index, k, of a cpath that passes through the missing dot
                        if isempty(cpaths[j] ∩ cpaths[k])
                            #ToDO: before supporting codes that allow for two drops, need to handle more cases here
                            push!(partial_conflicts[k], i)
                            push!(partial_conflict_transitions[k], j)
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
        deleteat!(cpaths_df, rsbcp_ind)
    end

    return cpath_nbr_cpath_indices, partial_conflicts, partial_conflict_transitions
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

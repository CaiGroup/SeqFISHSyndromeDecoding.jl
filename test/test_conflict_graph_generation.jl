using Random

Random.seed!(1)

#push!(LOAD_PATH, pwd())
#using syndrome_decoding
using SeqFISHPointDecoding
#using syndrome_decoding: add_code_cols!, DotAdjacencyGraph, objective_function,
#remove_high_cost_cpaths, codebooks, make_cw_dict, syndrome_find_message_paths!,
#get_cpath_conflict_graph_remove_redundant_cpaths!, get_cpath_conflict_graph_pairwise, get_cpath_conflict_graph2

using SeqFISHPointDecoding: add_code_cols!, DotAdjacencyGraph, objective_function,
remove_high_cost_cpaths, codebooks, make_cw_dict, syndrome_find_message_paths!,
get_cpath_conflict_graph_remove_redundant_cpaths!, get_cpath_conflict_graph_pairwise, get_cpath_conflict_graph2


using DelimitedFiles
using simulation_fncs
using DataFrames
using Profile
using LightGraphs

rstdv = 0.01
lat_thresh = 0.04
z_thresh = 1.0
drop_rate = 0.0
ndrops = 1
free_dot_energy = 5.0

code = "E2019"
ntargets = 10

ndots = 4 * ntargets

q20_cb = readdlm("Eng2019_647.csv", ',', UInt8)

true_locs = sim_true(ntargets)
pnts = encode(true_locs, q20_cb)
pnts.z = zeros(Float64, nrow(pnts))

drop_random_dots!(pnts, drop_rate)
add_localization_errors!(pnts, rstdv)
sort!(pnts, :hyb)

c_final = 2
l_chain = round(Int, ndots * log(ndots) * 100)
c₀ = free_dot_energy * 20

n_chains = 50
sa_params = (
    n_chains = n_chains,
    l_chain = l_chain,
    free_dot_cost = free_dot_energy,
    cooling_factor = c₀,#ndots/10,
    cooling_timescale = (c₀/c_final-1)/log(n_chains),#1 / (100*(1 + log(ndots))),#1/(10*ndots),
    lat_var_factor = 300.0, #150.0,
    z_var_factor = 1.0,
    erasure_penalty = 3.0,
    converge_thresh = 100 * ndots
)

println("start syndrome decoding")
#sa_params :: NamedTuple{(:n_iterations, :free_dot_cost, :cooling_factor, :cooling_timescale, :var_factor, :error_penalty, :converge_thresh), Tuple{Int64,Float64,Float64,Float64,Float64,Float64,Int64}}
#)
add_code_cols!(pnts, code)
cb = codebooks[code]
cb_dict = make_cw_dict(cb)
q = length(unique(cb))
w = sum(cb[1,:] .!= 0 )
n = length(cb[1,:])
g = DotAdjacencyGraph(pnts, lat_thresh, z_thresh, n)

cost(cpath) = objective_function(cpath, pnts, n, sa_params)

# Split graph into connected components
#wccs = weakly_connected_components(g.g)
code_paths, values = syndrome_find_message_paths!(pnts, g, code, ndrops)

costs = cost.(code_paths)

cpath_df = DataFrame(cpath = code_paths, cost = costs, value = values)

sort!(cpath_df, :cost)

cpath_df = remove_high_cost_cpaths(cpath_df, sa_params.free_dot_cost, n, ndrops)

if nrow(cpath_df) == 0
    println("No viable code paths")
    return
end

function time_get_cpath_conflict_graph_remove_redundant_cpaths(cpaths, ndots, n)
    @time cpath_nbrs, partial_conflicts, partial_conflict_transitions = get_cpath_conflict_graph_remove_redundant_cpaths!(cpaths, ndots, n)
    return (cpath_nbrs, partial_conflicts, partial_conflict_transitions)
end

function profile_get_cpath_conflict_graph(cpaths, ndots)
    @profile nbr_dicts = get_cpath_conflict_graph(cpaths.cpath, nrow(pnts))
    return nbr_dicts
end

function time_get_cpath_conflict_graph_pairwise(cpaths)#, ndots)
    #@time As, ccs = get_cpath_conflict_graph_pairwise(cpaths.cpath)#, nrow(pnts))
    @time As = get_cpath_conflict_graph_pairwise(cpaths.cpath)#, nrow(pnts))
    return As
end

function time_get_cpath_conflict_graph2(cpaths, ndots, ndrops)
    @time As, ccs = get_cpath_conflict_graph2(cpaths, nrow(pnts), ndrops)
end


println("pairwise")
A = time_get_cpath_conflict_graph_pairwise(cpath_df)

println("remove redundant cpaths")
cpath_nbrs, partial_conflicts, partial_conflict_transitions = time_get_cpath_conflict_graph_remove_redundant_cpaths(cpath_df, nrow(pnts), n)
#profile_get_cpath_conflict_graph(cpath_df, nrow(pnts))
#Profile.print()
#time_get_cpath_conflict_graph(cpath_df, nrow(pnts))


for i = 1:length(cpath_nbrs)
    @assert cpath_nbrs[i] == neighbors(A, i)
end

println("passed asserts")

#println()
#, nrow(pnts))
#time_get_cpath_conflict_graph_new(cpath_df, nrow(pnts))
#=

println()
println("3")
time_get_cpath_conflict_graph2(cpath_df, nrow(pnts), ndrops)
time_get_cpath_conflict_graph2(cpath_df, nrow(pnts), ndrops)
#@time As, ccs = get_cpath_conflict_graph(cpath_df.cpath, nrow(pnts))
=#

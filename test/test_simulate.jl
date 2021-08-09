include("simulation_fncs.jl")

using DelimitedFiles
using SeqFISHPointDecoding: DotAdjacencyGraph, syndrome_find_message_paths!, add_code_cols!
using LightGraphs
using SeqFISHPointDecoding

q20_cb = readdlm("Eng2019_647.csv", ',', UInt8)

ntargets = 10
cb = q20_cb

rstdv = 0.01
thresh = 0.2
drop_rate = 0.1
ndrops = 1
free_dot_energy = 5.0

n, q, w = get_n_q_w(cb)
SeqFISHPointDecoding.set_n(UInt8(n))
SeqFISHPointDecoding.set_q(UInt8(q))
SeqFISHPointDecoding.set_H([1 1 -1 -1;])
SeqFISHPointDecoding.get_decode_table()
#for i = 1:100
true_locs = sim_true(ntargets)
pnts = encode(true_locs, cb)
drop_random_dots!(pnts, drop_rate)
add_localization_errors!(pnts, rstdv)
sort!(pnts, :hyb)
add_code_cols!(pnts :: DataFrame)
g = DotAdjacencyGraph(pnts, thresh, 0.0, n)

code_paths, values = syndrome_find_message_paths!(pnts, g, cb, 1)

ndots_correct = sum(pnts.species .== pnts.decoded)
total_dots = nv(g.g)

println("$ndots_correct correct of $total_dots total")

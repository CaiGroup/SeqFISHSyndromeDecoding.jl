
#using syndrome_decoding: syndrome_find_message_paths!
using SeqFISHPointDecoding: syndrome_find_message_paths!
using test_synd_decode_fncs
using rule_out_cpaths: sort_cpaths_by_cost!, find_min_var_consensus_cpaths!,
rule_out_high_cost_cpaths!, init_min_cost_consensus_state, get_cpath_list
using simulation_fncs

using DataFrames
using DelimitedFiles


var_factor = 200.0
error_penalty = 3.0
free_dot_energy = 5.0

ntargets = 20
rstdv = 0.1
thresh = 0.4
code = "T2020"
drop_rate = 0.0

ndots = ntargets*5

cb = readdlm("q16_n5_w5_md2.txt", UInt8)


sa_params = (
    n_iterations = ndots*1000000,
    free_dot_energy = free_dot_energy,
    cooling_factor = free_dot_energy*8,
    cooling_timescale = 1/(10*ndots),
    var_factor = var_factor,
    error_penalty = error_penalty,
    converge_thresh = 100*ndots
)


pnts, g = construct_test_dag(ntargets, rstdv, thresh, code)

#=
true_locs = sim_true(ntargets)
pnts = encode(true_locs, cb)
drop_random_dots!(pnts, drop_rate)
add_localization_errors!(pnts, rstdv)
sort!(pnts, :hyb)
=#

syndrome_find_message_paths!(pnts, g, code)

sort_cpaths_by_cost!(pnts, var_factor, error_penalty)

all_cpaths_b4 = get_cpath_list(pnts)

println("b4: ", sum(nrow.(pnts.cpath_costs)), "; ", length(all_cpaths_b4))
println(sum(length.(pnts.perfect_cpaths)))
println(sum(length.(pnts.error_cpaths)))

rule_out_high_cost_cpaths!(pnts, free_dot_energy)

all_cpaths_after = get_cpath_list(pnts)

println("after: ", sum(nrow.(pnts.cpath_costs)), "; ", length(all_cpaths_after))
println(sum(length.(pnts.perfect_cpaths)))
println(sum(length.(pnts.error_cpaths)))

#mvar_consensus_cpaths = find_min_var_consensus_cpaths!(pnts)
#println(score_decoding(true_locs, pnts))
#println("sum matched dots: ", sum(length.(mvar_consensus_cpaths)))

#state = init_min_cost_consensus_state(pnts, mvar_consensus_cpaths, free_dot_energy, sa_params)

println("Done")

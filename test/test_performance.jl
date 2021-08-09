using Random

Random.seed!(1)
include("simulation_fncs.jl")

using SeqFISHPointDecoding

using DelimitedFiles
using Profile
using DataFrames

#q11_cb = readdlm("trimmed_RS_q11_n8_w5_md4_codewords.txt", UInt8)
#q12_cb = readdlm("ST2018_cb.csv", ',', UInt8)[1:10421,:]
#q16_cb = readdlm("q16_n5_w5_md2.txt", UInt8)
q20_cb = readdlm("Eng2019_647.csv", ',', UInt8)



cb = q20_cb
#cb = q16_cb
#cb = q12_cb
#q = 11
#n = 8
rstdv = 0.03
#rstdv = 0.0
#rstdv = 0.025
#thresh = 0.4
lat_thresh = 0.13#0.8
#lat_thresh = 0.0

z_thresh = 1.0
#drop_rate = 0.0#5
drop_rate = 0.05
ndrops = 1
free_dot_energy = 5.0
mip_sa_thresh = 80

ntargets = 200



#ndots = 5 * ntargets
ndots = 4 * ntargets

true_locs = sim_true(ntargets)#, 10421)
pnts = encode(true_locs, cb)
pnts[!, "z"] .= 0.0 #zeros(Float64, nrow(pnts))
pnts[!,"w"] .= 1.0
pnts[!, "s"] .= 0.0

c_final = 2
n_chains = 100
l_chain = 300#100#
c₀ = free_dot_energy * 40
lat_var_factor = 8000.0
z_var_factor = 1.0
lw_var_factor = 0.0
s_var_factor = 0.0
erasure_penalty = 4.0
converge_thresh = 100 * ndots
skip_thresh = 2000

params = SeqFISHPointDecoding.DecodeParams(
    lat_thresh,
    z_thresh,
    lat_var_factor,
    z_var_factor,
    lw_var_factor,
    s_var_factor,
    ndrops,
    mip_sa_thresh,
    free_dot_energy,
    n_chains,
    l_chain,
    c₀,
    (c₀/c_final-1)/log(n_chains),
    erasure_penalty,
    converge_thresh,
    skip_thresh
)

drop_random_dots!(pnts, drop_rate)
add_localization_errors!(pnts, rstdv)
sort!(pnts, :hyb)
println()
pnts0 = deepcopy(pnts)

println("Decode Syndrome")

println("No negative control")
@time mpaths = decode_syndromes!(pnts, q20_cb, [1 1 -1 -1;], params)



println(score_decoding(true_locs, pnts))

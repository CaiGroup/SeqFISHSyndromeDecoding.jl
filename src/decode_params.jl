export set_xy_search_radius, set_z_search_radius, set_n_allowed_drops, set_lat_var_cost_coeff,
set_z_var_cost_coeff, set_lw_var_cost_coeff, set_s_var_cost_coeff,
set_mip_sa_thresh, set_free_dot_cost, set_n_chains, set_l_chains, set_cooling_factor,
set_cooling_timescale, set_erasure_penalty, set_converge_multiplier, set_skip_thresh

mutable struct DecodeParams
    lat_thresh :: Float64
    z_thresh :: Float64
    lat_var_factor :: Float64
    z_var_factor :: Float64
    lw_var_factor :: Float64
    s_var_factor :: Float64
    ndrops :: Int64
    mip_sa_thresh :: Int64
    free_dot_cost :: Float64
    n_chains :: Int64
    l_chain :: Int64
    cooling_factor :: Float64
    cooling_timescale :: Float64
    erasure_penalty :: Float64
    converge_multiplier :: Int64
    skip_thresh :: Int64
end


"""
    DecodeParams()

Constructor for DecodeParams object that sets default values for every parameter. Dost function parameters must be set
by the user. If they are not, the decoding functions will throw an error.
"""
function DecodeParams()
    prms = DecodeParams(fill(0.0, 16)...)
    c_final = 1
    # set simmulated annealing parameters default values
    set_mip_sa_thresh(prms, 80)
    set_free_dot_cost(prms, 5.0)
    set_n_chains(prms, 100)
    set_l_chains(prms, 300)
    set_cooling_factor(prms, 200.0)
    ct = (prms.cooling_factor/c_final -1)/log(prms.n_chains)
    set_cooling_timescale(prms, ct :: Real)
    set_erasure_penalty(prms, 4.0)
    set_converge_multiplier(prms, 10000)
    set_skip_thresh(prms, 2000)
    return prms
end


"""
    set_xy_search_radius(prms :: DecodeParams, r :: Real)

Arguments:
    -`prms`: DecodeParams Object
    -`r`: The lateral KD Search radius in pixels for aligned dots in previous barcoding rounds.
"""
set_xy_search_radius(prms :: DecodeParams, r :: Real) = (prms.lat_thresh = Float64(r))

"""
    set_z_search_radius(prms :: DecodeParams, r :: Real)

Arguments:
    -`prms`: DecodeParams Object
    -`r`: The z KD Search radius in slices for aligned dots in previous barcoding rounds.
"""
set_z_search_radius(prms :: DecodeParams, r :: Real) = (prms.z_thresh = Float64(r))

"""
    set_n_allowed_drops(prms :: DecodeParams, d :: Integer)

Arguments:
    -`prms`: DecodeParams Object
    -`r`: The number of allowed drops in a barcode. Only 0 and 1 are currently supported. The parity check matrix must have redundancy to support correcting for one drop, otherwise errors will occur.
"""
set_n_allowed_drops(prms :: DecodeParams, d :: Integer) = (prms.ndrops = d)

"""
    set_lat_var_cost_coeff(prms :: DecodeParams, lvf :: Real)

Arguments:
    -`prms`: DecodeParams Object
    -`lvf`: The lateral variance penalty coefficient in the cost function.
"""
set_lat_var_cost_coeff(prms :: DecodeParams, lvf :: Real) = (prms.lat_var_factor = lvf)

"""
    set_z_var_cost_coeff(prms :: DecodeParams, zvf :: Real)

Arguments:
    -`prms`: DecodeParams Object
    -`r`: The z variance penalty coefficient in the cost function.
"""
set_z_var_cost_coeff(prms :: DecodeParams, zvf :: Real) = (prms.z_var_factor = Float64(zvf))

"""
    set_lw_var_cost_coeff(prms :: DecodeParams, lwvf :: Real)

Arguments:
    -`prms`: DecodeParams Object
    -`lwvf`: The log weight variance penalty coefficient in the cost function.
"""
set_lw_var_cost_coeff(prms :: DecodeParams, lwvf :: Real) = (prms.lw_var_factor = Float64(lwvf))

"""
    set_s_var_cost_coeff(prms :: DecodeParams, sf :: Real)

Arguments:
    -`prms`: DecodeParams Object
    -`svf`: The σ variance penalty coefficient in the cost function.
"""
set_s_var_cost_coeff(prms :: DecodeParams, svf :: Real) = (prms.s_var_factor = Float64(svf))

"""
    set_erasure_penalty(prms :: DecodeParams, ep :: Real)

The penalty that each dot dropped from a barcode adds to the cost function of the corrected barcode is ep. For example a barcode with one drop adds 1*ep to the cost of the barcode.
"""
set_erasure_penalty(prms :: DecodeParams, ep :: Real) = (prms.erasure_penalty = Float64(ep))


"""
    set_free_dot_cost(prms :: DecodeParams, fdc :: Real)

Arguments:
    -`prms`: DecodeParams Object
    -`fdc`: The cost of not decoding a dot in the barcode candidate network
"""
set_free_dot_cost(prms :: DecodeParams, fdc :: Real) = (prms.free_dot_cost = Float64(fdc))

"""
    set_mip_sa_thresh(prms :: DecodeParams, mst :: Real)

Arguments:
    -`prms`: DecodeParams Object
    -`mst`: Connected networks of more barcode candidates than this threshold will use simulated annealing to approximate the optimal barcodes.
"""
set_mip_sa_thresh(prms :: DecodeParams, mst :: Integer) = (prms.mip_sa_thresh = mst)

"""
    set_skip_thresh(prms :: DecodeParams, st :: Integer)

Arguments:
    -`prms`: DecodeParams Object
    -`st`: conflicting networks of barcode candidates containing more barcode candidates than this threshold will be discarded.
"""
set_skip_thresh(prms :: DecodeParams, st :: Integer) = (prms.skip_thresh = st)

# Simulated Annealing Parameters

"""
    set_n_chains(prms :: DecodeParams, nc :: Integer)

Sets the number of markov chains used int the simulated annealing algorithm. Each
Markov chain operates at the same "temperature" or control parameter
"""
set_n_chains(prms :: DecodeParams, nc :: Integer) = (prms.n_chains = nc)

"""
    set_l_chains(prms :: DecodeParams, lc :: Integer)

Sets the number of estimated Metropolis steps that the rejectionless markov chains take.
"""
set_l_chains(prms :: DecodeParams, lc :: Integer) = (prms.l_chain = lc)

"""
    set_cooling_factor(prms :: DecodeParams, cf :: Real)

Sets the start "temperature" or control parameter of the simulated annealing system
"""
set_cooling_factor(prms :: DecodeParams, cf :: Real) = (prms.cooling_factor = Float64(cf))

"""
    set_cooling_timescale(prms :: DecodeParams, ct :: Real)
    
Sets the decay rate of the simulated annealing temperature between markov chains
"""
set_cooling_timescale(prms :: DecodeParams, ct :: Real) = (prms.cooling_timescale = Float64(ct))


"""
    set_converge_multiplier(prms :: DecodeParams, ct :: Integer)

If the rejectionless simulated annealing algorithm favors the current state by
more than probablity (ct/(ct+1)) over all transitioning to any other state, converge
"""
set_converge_multiplier(prms :: DecodeParams, ct :: Integer) = (prms.converge_multiplier = ct)

export set_xy_search_radius, set_z_search_radius, set_n_allowed_drops, set_lat_var_cost_coeff,
set_z_var_cost_coeff, set_lw_var_cost_coeff, set_s_var_cost_coeff, set_zeros_probed, set_free_dot_cost, set_erasure_penalty, set_skip_thresh, set_skip_density_thresh

mutable struct DecodeParams
    lat_thresh :: Float64
    z_thresh :: Float64
    lat_var_factor :: Float64
    z_var_factor :: Float64
    lw_var_factor :: Float64
    s_var_factor :: Float64
    ndrops :: Int64
    zeros_probed :: Bool
    free_dot_cost :: Float64
    erasure_penalty :: Float64
    skip_thresh :: Int64
    skip_density_thresh :: Float64
end


"""
    DecodeParams()

Constructor for DecodeParams object that sets default values for every parameter. Cost function parameters must be set
by the user. If they are not, the decoding functions will throw an error.
"""
function DecodeParams()
    prms = DecodeParams(fill(0.0, 12)...)
    set_zeros_probed(prms, true)
    set_free_dot_cost(prms, 1.0)
    set_erasure_penalty(prms, 1.0)
    set_skip_thresh(prms, 2000)
    set_skip_density_thresh(prms, 50)
    return prms
end


"""
    set_xy_search_radius(prms :: DecodeParams, r :: Real)

# Arguments
- `prms`: DecodeParams Object
- `r`: The lateral KDTree Search radius in pixels for aligned dots in previous barcoding blocks.
"""
set_xy_search_radius(prms :: DecodeParams, r :: Real) = (prms.lat_thresh = Float64(r))

"""
    set_z_search_radius(prms :: DecodeParams, r :: Real)

# Arguments
- `prms`: DecodeParams Object
- `r`: The z KDTree Search radius in slices for aligned dots in previous barcoding blocks.
"""
set_z_search_radius(prms :: DecodeParams, r :: Real) = (prms.z_thresh = Float64(r))

"""
    set_n_allowed_drops(prms :: DecodeParams, d :: Integer)

# Arguments
- `prms`: DecodeParams Object
- `r`: The number of allowed drops in a barcode. Only 0 and 1 are currently supported. The parity check matrix must have redundancy to support correcting for one drop, otherwise errors will occur.
"""
set_n_allowed_drops(prms :: DecodeParams, d :: Integer) = (prms.ndrops = d)

"""
    set_zeros_probed(prms :: DecodeParams, zp :: Bool)

# Arguments
- `prms`: DecodeParams Object
- `zp`: Whether or not zeros in codewords are probed in the experiment. Drops are not supported when false.
"""
set_zeros_probed(prms :: DecodeParams, zp :: Integer) = (prms.zeros_probed = zp)
set_zeros_probed(prms :: DecodeParams, zp :: Bool) = (prms.zeros_probed = zp)


"""
    set_lat_var_cost_coeff(prms :: DecodeParams, lvf :: Real)

# Arguments
- `prms`: DecodeParams Object
- `lvf`: The lateral variance penalty coefficient in the cost function.
"""
set_lat_var_cost_coeff(prms :: DecodeParams, lvf :: Real) = (prms.lat_var_factor = lvf)

"""
    set_z_var_cost_coeff(prms :: DecodeParams, zvf :: Real)

# Arguments
- `prms`: DecodeParams Object
- `r`: The z variance penalty coefficient in the cost function.
"""
set_z_var_cost_coeff(prms :: DecodeParams, zvf :: Real) = (prms.z_var_factor = Float64(zvf))

"""
    set_lw_var_cost_coeff(prms :: DecodeParams, lwvf :: Real)

# Arguments
- `prms`: DecodeParams Object
- `lwvf`: The log weight variance penalty coefficient in the cost function.
"""
set_lw_var_cost_coeff(prms :: DecodeParams, lwvf :: Real) = (prms.lw_var_factor = Float64(lwvf))

"""
    set_s_var_cost_coeff(prms :: DecodeParams, sf :: Real)

# Arguments
- `prms`: DecodeParams Object
- `svf`: The σ variance penalty coefficient in the cost function.
"""
set_s_var_cost_coeff(prms :: DecodeParams, svf :: Real) = (prms.s_var_factor = Float64(svf))

"""
    set_erasure_penalty(prms :: DecodeParams, ep :: Real)

The penalty that each dot dropped from a barcode adds to the cost function of the corrected barcode is ep. For example a barcode with one drop adds 1*ep to the cost of the barcode.

# Arguments
- `prms`: DecodeParams Object
- `ep`: penalty for decoding a barcode with a dropped dot
"""
set_erasure_penalty(prms :: DecodeParams, ep :: Real) = (prms.erasure_penalty = Float64(ep))


"""
    set_free_dot_cost(prms :: DecodeParams, fdc :: Real)

# Arguments
- `prms`: DecodeParams Object
- `fdc`: The cost of not decoding a dot in the barcode candidate network
"""
set_free_dot_cost(prms :: DecodeParams, fdc :: Real) = (prms.free_dot_cost = Float64(fdc))


"""
    set_skip_thresh(prms :: DecodeParams, st :: Integer)

# Arguments
- `prms`: DecodeParams Object
- `st`: conflicting networks of barcode candidates containing more barcode candidates than this threshold will be discarded.
"""
set_skip_thresh(prms :: DecodeParams, st :: Integer) = (prms.skip_thresh = st)

"""
    set_skip_density_thresh(prms :: DecodeParams, sdt :: Real)

# Arguments:
- `prms`: DecodeParams Object
- `sdt`: conflicting networks of barcode candidates that have a higher ratio of candidate barcodes to area of their bounding box in pixels than this threshold will be discarded.
"""
set_skip_density_thresh(prms :: DecodeParams, sdt :: Real) = (prms.skip_density_thresh = Float64(sdt))
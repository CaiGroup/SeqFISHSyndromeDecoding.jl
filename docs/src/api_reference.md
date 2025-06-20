
# API Reference

## Setting Parameters

To set parameters, use a DecodeParams object

```@docs
DecodeParams
```

## Required parameters

Users are required to set the some parameters for decoding using the following methods:

```@docs
set_xy_search_radius
set_z_search_radius
set_n_allowed_drops
set_lat_var_cost_coeff
set_z_var_cost_coeff
set_lw_var_cost_coeff
set_s_var_cost_coeff
```

### Parameters with default values 

The rest of the parameters are given a default value by the [`DecodeParams`](@ref) Constructor, but may be set using the following methods:

```@docs
set_zeros_probed
set_erasure_penalty
set_free_dot_cost
set_skip_thresh
set_skip_density_thresh
```

## Running Decoding

```@docs
get_codepaths
choose_optimal_codepaths
decode_syndromes!
```


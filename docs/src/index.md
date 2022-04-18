# SeqFISHSyndromeDecoding.jl


## Introduction
This package provides functions to decode seqFISH points using [Syndrome Decoding](https://en.wikipedia.org/wiki/Decoding_methods#Syndrome_decoding).

In SeqFISH Experiments, target sequences are probed in 


## Installation

To install, from the julia REPL run:

```
julia> using Pkg
julia> Pkg.add(url="https://github.com/CaiGroup/SeqFISHSyndromeDecoding")
```

Alternatively, open the package manager by typing
```
julia>] 
```
then
```
Pkg> add "https://github.com/CaiGroup/SeqFISHSyndromeDecoding"
```

## Example



## API Reference

### Setting Parameters

To set parameters, use a Parameters Class

```@docs
DecodeParams
```

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

The rest of the parameters are given a default value by the `DecodeParams` Constructor, but may be set using the following methods:

```@docs
set_erasure_penalty
set_free_dot_cost
set_mip_sa_thresh
set_skip_thresh
set_skip_density_thresh
```

Simulated Annealing Parameters:

```@docs
set_n_chains
set_l_chains
set_cooling_factor
set_cooling_timescale
set_converge_multiplier
```


### Running Decoding

```@docs
get_codepaths
choose_optimal_codepaths
```

```@docs
decode_syndromes!
```


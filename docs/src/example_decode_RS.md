# Decoding Example: Reed-Solomon Codes

*note: this example was generated using data and a jupyter notebook that are freely available at the [SeqFISHSyndromeDecoding github repository](https://github.com/CaiGroup/SeqFISHSyndromeDecoding) and on [Google Colab](https://colab.research.google.com/github/CaiGroup/SeqFISHSyndromeDecoding.jl/blob/master/example_notebook/colab/example_decode_RS_colab.jl.ipynb).*

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
using DataFrames
using CSV
using SeqFISHSyndromeDecoding
using GLPK
using DelimitedFiles
```    

This notebook shows demonstrates how to use SeqFISHSyndromeDecoding.jl. The smallest cell with the fewest dots in our Reed-Solomon encoded experiment, chosen for computational convienience. We also reduce computation time by using the highest lateral position variance reported in our manuscript with half the search radius used in the manuscript computations. The larger positional variance penalty would prohibit most additional candidate barcodes found with larger search radius.

First, load the codebook that we will use to decode our sample data.


```julia
cb = DataFrame(CSV.File("../example_data/full_RS_q11_k7_half_pool_cb.csv"))
println(first(cb, 5))
```

    DataFrame
    Row  │ gene      block1  block2  block3  block4  block5  block6  block7  block8  block9  block10 
         │ String31  Int64   Int64   Int64   Int64   Int64   Int64   Int64   Int64   Int64   Int64   
    ─────┼───────────────────────────────────────────────────────────────────────────────────────────
       1 │ Aars           2       8       5       0       0       0       0       0       2        0
       2 │ Aco2          10       7       0       0       0       0       8       0       0        6
       3 │ Actn4          9       7       0       0       9       2       0       0       0        0
       4 │ Aebp1          2       0       1       3       0       0       0       2       0        0
       5 │ Aqp1          10       0       3       0       0       0       9       1       0        0
    

Define the [parity check matrix](https://en.wikipedia.org/wiki/Parity-check_matrix) for the codebook


```julia
H = readdlm("../example_data/RS_q11_k7_H.csv", ',', UInt8)
```




    3×10 Matrix{UInt8}:
     0x02  0x04  0x08  0x05  0x0a  0x09  0x07  0x03  0x06  0x01
     0x04  0x05  0x09  0x03  0x01  0x04  0x05  0x09  0x03  0x01
     0x08  0x09  0x06  0x04  0x0a  0x03  0x02  0x05  0x07  0x01



We can verify that H is actually the parity check matrix of the codebook.


```julia
all(H * Matrix(cb[:,2:end])' .% 11 .== 0)
```




    true



Next we can load the aligned points from each hybridization for our example cell.


```julia
pnts = DataFrame(CSV.File("../example_data/example_RS_cell_points.csv"))
filter!(pnt -> ~ismissing(pnt.pseudocolor), pnts)
pnts.block = UInt8.(pnts.block)
select!(pnts, Not([:ch,:hyb]))
SeqFISHSyndromeDecoding.sort_readouts!(pnts)
println(first(pnts, 5))
```

    DataFrame
    Row  │ Column1  x        y        s        w        z        pos    pseudocolor  block  cellid 
         │ Int64    Float64  Float64  Float64  Float64  Float64  Int64  Float64?     UInt8  Int64  
    ─────┼─────────────────────────────────────────────────────────────────────────────────────────
       1 │       0   94.938  1884.65  1.21555  267.268      0.0      6          1.0      1       6
       2 │    3815  105.651  1873.63  1.32669  253.041      0.0      6          1.0      1       6
       3 │       8  109.834  1836.6   2.0      304.417      0.0      6          1.0      1       6
       4 │      26  110.489  1873.2   1.47603  884.79       0.0      6          1.0      1       6
       5 │      33  110.8    1885.17  1.50077  571.21       0.0      6          1.0      1       6
    

Next we initialize a ```DecodeParams``` object, and set the parameters


```julia
params = DecodeParams()

set_zeros_probed(params, false)
set_lat_var_cost_coeff(params, 7.0)
set_z_var_cost_coeff(params, 0.0)
set_lw_var_cost_coeff(params, 0.0)
set_s_var_cost_coeff(params, 0.0)
set_free_dot_cost(params, 1.0)
set_n_allowed_drops(params, 0)

set_xy_search_radius(params, 2)
set_z_search_radius(params, 0.0);
```

We can then decode


```julia
barcodes = decode_syndromes!(pnts, cb, H, params);
println(first(barcodes, 5))
```


     DataFrame
     Row │ gene              gene_number  cpath                         cost      x        y        z    cc     cc_size 
         │ String31?         Any         Any                          Float64  Any     Any     Any Int64 Int64   
    ─────┼──────────────────────────────────────────────────────────────────────────────────────────────────────────────
       1 │ St13              220          [7127, 9160, 11942, 13432]    0.209086  114.12   1890.67  0.0      1       11
       2 │ Ybx1              263          [13625, 16733, 21893, 24052]  0.621507  113.441  1891.12  0.0      1       11
       3 │ negative_control  299          [6, 7923, 10905, 16128]       0.790891  110.843  1889.99  0.0      1       11
       4 │ Msn               112          [5598, 7267, 19024, 23060]    3.25967   111.708  1890.39  0.0      1       11
       5 │ Rarg              179          [11, 4496, 16538, 22858]      1.42675   114.981  1918.55  0.0      2        1
    

Alternatively, if we aren't sure what parameters we want to use, we can save time by splitting decode_syndromes! into its two steps. First we can identify barcode candidates with the ```get_codepaths``` (named for the paths that candidate barcodes take the the decoding graph in figure 1a) function using the least strict parameter set that we are interested in.


```julia
candidates = get_codepaths(pnts, cb, H, params);
println(first(candidates, 5))

```

    1m5×7 DataFrame
    Row  │ gene              gene_number  cpath                        cost      x        y        z  
         │ String31?         Any          Any                          Any       Any      Any      Any
    ─────┼─────────────────────────────────────────────────────────────────────────────────────────────
       1 │ Ercc1             40           [4965, 12136, 16528, 19370]  0.286652  106.317  1843.72  0.0
       2 │ Vim               261          [9474, 14973, 16125, 21132]  0.651009  108.019  1846.97  0.0
       3 │ Ercc1             40           [4971, 12138, 16529, 19371]  1.27286   108.758  1849.4   0.0
       4 │ negative_control  971          [763, 4117, 6469, 19371]     3.28321   108.719  1849.74  0.0
       5 │ Slc38a2           209          [1571, 7619, 12630, 14130]   1.16154   108.261  1854.62  0.0
    

We can then use the ```choose_optimal_codepaths``` function to find the same barcodew that we found earlier


```julia
barcodes_again = choose_optimal_codepaths(pnts, cb, H, params, candidates, GLPK.Optimizer)
barcodes == barcodes_again
```




    true



We can now also try choosing candidates using stricter parameters. This saves computation time by reducing the number of times that we have to run ```get_codepaths```.


```julia
strict_params = DecodeParams()

set_zeros_probed(strict_params, false)
set_lat_var_cost_coeff(strict_params, 10.0)
set_z_var_cost_coeff(strict_params, 0.0)
set_lw_var_cost_coeff(strict_params, 0.0)
set_s_var_cost_coeff(strict_params, 0.0)
set_free_dot_cost(strict_params, 1.0)
set_n_allowed_drops(strict_params, 0)

set_xy_search_radius(strict_params, 2)
set_z_search_radius(strict_params, 0.0);

stricter_barcodes = choose_optimal_codepaths(pnts, cb, H, strict_params, candidates, GLPK.Optimizer)
println(first(stricter_barcodes, 5))
```

    DataFrame
     Row │ gene              gene_number  cpath                      cost      x        y        z    cc     cc_size
         │ String31?         Any          Any                        Float64   Any      Any      Any  Int64  Int64   
    ─────┼───────────────────────────────────────────────────────────────────────────────────────────────────────────
       1 │ negative_control  299          [6, 7923, 10905, 16128]    1.12984   110.843  1889.99  0.0      1        2
       2 │ Rarg              179          [11, 4496, 16538, 22858]   2.03821   114.981  1918.55  0.0      2        1
       3 │ Khdrbs1           73           [13, 15343, 19604, 23071]  0.231704  117.603  1892.83  0.0      3        1
       4 │ Spp1              214          [18, 10922, 11844, 14141]  2.79514   119.057  1867.67  0.0      4        1
       5 │ Spp1              214          [19, 10921, 11846, 14142]  2.42644   119.15   1862.08  0.0      5        2
    

We can compare the decoding results using the two different sets of parameters. For brevity, we use gene encoding barcodes found in decoding runs that include searches for negative control barcodes, which differs from the procedure described in our manuscript in which datasets are also decoded with the negative control codewords ommitted from the codebook.


```julia
println("Number of gene encoding barcodes: ", sum(barcodes.gene .!= "negative_control"))
estimated_false_discovery_rate = sum(barcodes.gene .== "negative_control")*sum(cb.gene .!= "negative_control")/sum(cb.gene .== "negative_control")/sum(barcodes.gene .!= "negative_control")
println("Estimated False Discovery rate: ", estimated_false_discovery_rate)
```

    Number of gene encoding barcodes: 1737
    Estimated False Discovery rate: 0.01654845665984571
    


```julia
println("Number of gene encoding barcodes: ", sum(stricter_barcodes.gene .!= "negative_control"))
estimated_false_discovery_rate = sum(stricter_barcodes.gene .== "negative_control")*sum(cb.gene .!= "negative_control")/sum(cb.gene .== "negative_control")/sum(stricter_barcodes.gene .!= "negative_control")
println("Estimated False Discovery rate: ", estimated_false_discovery_rate)
```

    Number of gene encoding barcodes: 1303
    Estimated False Discovery rate: 0.011712467380864363
    

The less strict parameter set decodes about 40% more gene encoding barcodes at a cost of having twice the estimated false discovery rate. Since the estimated false positive rate is still small, it is probably an acceptable trade off.

To save your results, use the ```CSV.write``` command.


```julia
CSV.write("example_RS_results.csv", barcodes)
```




    "example_RS_results.csv"




# Decoding Example: SeqFISH+

*note: this example was generated using data and a jupyter notebook that are freely available at the [SeqFISHSyndromeDecoding github repository](https://github.com/CaiGroup/SeqFISHSyndromeDecoding) and can be run on [google colab](https://colab.research.google.com/github/CaiGroup/SeqFISHSyndromeDecoding.jl/blob/master/example_notebook/colab/example_decode_colab.jl.ipynb)*

```julia
using DataFrames
using CSV
using SeqFISHSyndromeDecoding
using GLPK
```

This notebook shows demonstrates how to use SeqFISHSyndromeDecoding. The example data was taken from the 561 channel of cell number 8 in position 4 of replicate 2 of the [2019 SeqFISH+](https://doi.org/10.1038/s41586-019-1049-y) NIH3T3 cell experiment. This particular subset of the data was chosen for its small size.

First load the codebook that we will use to decode our sample data.


```julia
cb = DataFrame(CSV.File("../example_data/codebook_ch_561.csv"))
println(first(cb, 5))
```
```
5×5 DataFrame
 Row │ gene  block1  block2  block3  block4 
     │ String31   Int64   Int64   Int64   Int64  
─────┼───────────────────────────────────────────
   1 │ Atp6v0e2       17      16      13       0
   2 │ Pclo           14      14       4       4
   3 │ Higd1a          0       1       1       0
   4 │ Srrm1           6      14       4      16
   5 │ Mapk8ip3       14      19      13       0
```


Define the [parity check matrix](https://en.wikipedia.org/wiki/Parity-check_matrix) for the codebook


```julia
H = [1 1 -1 -1;]
```




    1×4 Matrix{Int64}:
     1  1  -1  -1



We can verify that H is actually the parity check matrix of the codebook.


```julia
all(H * Matrix(cb[:,2:end])' .% 20 .== 0)
```




    true



Next we can load the aligned points from each hybridization for our example cell.


```julia
pnts = DataFrame(CSV.File("../example_data/example_cell_points.csv"))
println(first(pnts, 5))
```
```
5×5 DataFrame
 Row │ hyb    x        y        s        w        
     │ Int64  Float64  Float64  Float64  Float64  
─────┼────────────────────────────────────────────
   1 │     1  767.664  1463.64     1.22   633.14
   2 │     1  759.413  1534.17     1.22  1118.99
   3 │     1  757.458  1501.22     1.22   866.506
   4 │     1  808.817  1400.84     1.22  1292.37
   5 │     1  804.688  1448.16     1.22  1734.99
```


The SeqFISHSyndromeDecoding package requires that they hybridization column be UInt8s (to increase efficiency), and that
there be a z column (for generality to 3d data)


```julia
pnts.z = zeros(Float64, nrow(pnts))
pnts.hyb = UInt8.(pnts.hyb);
```

Next we initialize a [`DecodeParams`](@ref) object, and set the parameters


```julia
params = DecodeParams()

set_lat_var_cost_coeff(params, 8.0)
set_z_var_cost_coeff(params, 0.0)
set_lw_var_cost_coeff(params, 3.2)
set_s_var_cost_coeff(params, 0.0)
set_free_dot_cost(params, 1.0)

set_xy_search_radius(params, sqrt(size(H)[2]/6.0)*3)
set_z_search_radius(params, 0.0);
```

We can then decode


```julia
barcodes = decode_syndromes!(pnts, cb, H, params)
first(barcodes, 5)
```
```
5×9 DataFrame
 Row │ gene  gene_number  cpath                      cost      x        y        z    cc     cc_size 
     │ String31?  Int64        Array…                     Float64   Any      Any      Any  Int64  Int64   
─────┼────────────────────────────────────────────────────────────────────────────────────────────────────
   1 │ Higd1a               3  [4000, 4203, 8632, 16480]  1.06754   990.275  1517.47  0.0   1406        1
   2 │ Srrm1                4  [1133, 7101, 9168, 15521]  2.14305   764.521  1474.04  0.0     95        2
   3 │ Srrm1                4  [1020, 7057, 9256, 15609]  1.17522   873.084  1535.57  0.0    390        1
   4 │ Srrm1                4  [1177, 7120, 9365, 15689]  2.8198    882.568  1479.06  0.0    431        1
   5 │ Srrm1                4  [1180, 7119, 9368, 15691]  0.241011  869.379  1540.38  0.0    432        1
```




Alternatively, if we aren't sure what parameters we want to use, we can save time by splitting [`decode_syndromes!`](@ref) into its two steps. First we can identify barcode candidates with the [`get_codepaths`](@ref) (named for the paths that candidate barcodes take the the decoding graph in figure 1a) function using the least strict parameter set that we are interested in.


```julia
candidates = get_codepaths(pnts, cb, H, params)
println(first(candidates, 5))
```
```
5×6 DataFrame
 Row │ cpath                       cost       gene_number  x        y        z   
     │ Array…                      Any        Int64        Any      Any      Any 
─────┼───────────────────────────────────────────────────────────────────────────
   1 │ [4034, 7497, 9354, 14834]   0.0322611          738  792.754  1543.43  0.0
   2 │ [2479, 6290, 9655, 15972]   0.0770821         1300  876.973  1460.59  0.0
   3 │ [1484, 8356, 10071, 16418]  0.0919919          271  884.483  1527.45  0.0
   4 │ [1148, 7491, 11855, 13651]  0.139732          2287  833.953  1495.81  0.0
   5 │ [2959, 6294, 8792, 13517]   0.168859          1755  845.245  1458.23  0.0
```


We can then use the [`choose_optimal_codepaths`](@ref) function to find the same barcodew that we found earlier


```julia
barcodes_again = choose_optimal_codepaths(pnts, cb, H, params, candidates, GLPK.Optimizer)
barcodes == barcodes_again
```




    true



We can now also try choosing candidates using stricter parameters. This saves computation time by reducing the number of times that we have to run [`get_codepaths`](@ref).


```julia
strict_params = DecodeParams()
set_lat_var_cost_coeff(strict_params, 12.0)
set_z_var_cost_coeff(strict_params, 0.0)
set_lw_var_cost_coeff(strict_params, 4.8)
set_s_var_cost_coeff(strict_params, 0.0)
set_free_dot_cost(strict_params, 1.0)

set_xy_search_radius(strict_params, sqrt(size(H)[2]/6.0)*3)
set_z_search_radius(strict_params, 0.0);


stricter_barcodes = choose_optimal_codepaths(pnts, cb, H, strict_params, candidates, GLPK.Optimizer)
println(first(stricter_barcodes, 5))
```
```
5×9 DataFrame
 Row │ gene  gene_number  cpath                      cost      x        y        z    cc     cc_size 
     │ String31?  Int64        Array…                     Float64   Any      Any      Any  Int64  Int64   
─────┼────────────────────────────────────────────────────────────────────────────────────────────────────
   1 │ Higd1a               3  [4000, 4203, 8632, 16480]  1.60131   990.275  1517.47  0.0   1091        1
   2 │ Srrm1                4  [1020, 7057, 9256, 15609]  1.76283   873.084  1535.57  0.0    291        1
   3 │ Srrm1                4  [1133, 7101, 9168, 15521]  3.21457   764.521  1474.04  0.0    312        1
   4 │ Srrm1                4  [1180, 7119, 9368, 15691]  0.361517  869.379  1540.38  0.0    326        1
   5 │ Wbp11                7  [3862, 7010, 9228, 14463]  1.26753   777.584  1569.5   0.0   1049        1
```


We can compare the decoding results using the two different sets of parameters.


```julia
println("Number of gene encoding barcodes: ", sum(barcodes.gene .!= "negative_control"))
estimated_false_discovery_rate = sum(barcodes.gene .== "negative_control")*sum(cb.gene .!= "negative_control")/sum(cb.gene .== "negative_control")/sum(barcodes.gene .!= "negative_control")
println("Estimated False Discovery rate: ", estimated_false_discovery_rate)
```

    Number of gene encoding barcodes: 1474
    Estimated False Discovery rate: 0.024709855479406056



```julia
println("Number of gene encoding barcodes: ", sum(stricter_barcodes.gene .!= "negative_control"))
estimated_false_discovery_rate = sum(stricter_barcodes.gene .== "negative_control")*sum(cb.gene .!= "negative_control")/sum(cb.gene .== "negative_control")/sum(stricter_barcodes.gene .!= "negative_control")
println("Estimated False Discovery rate: ", estimated_false_discovery_rate)
```

    Number of gene encoding barcodes: 1118
    Estimated False Discovery rate: 0.015330875292705262


The less strict parameter set decodes about 40% more gene encoding barcodes at a cost of having twice the estimated false discovery rate. Since the estimated false discovery rate is still small, it is probably an acceptable trade off.

To save your results, use the [`CSV.write`](https://csv.juliadata.org/stable/writing.html) command.


```julia
CSV.write("example_results.csv", barcodes)
```




    "example_results.csv"



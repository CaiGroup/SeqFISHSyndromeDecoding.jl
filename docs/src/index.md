# SeqFISHSyndromeDecoding.jl


## Introduction
This package provides functions to decode seqFISH points using [syndrome decoding](https://en.wikipedia.org/wiki/Decoding_methods#Syndrome_decoding).

The main functions are [`get_codepaths`](@ref) which uses dynamic programming to find barcode candidates and [`choose_optimal_codepaths`](@ref) which runs an integer programming optimization to choose the best non-conflicting candidates. These functions break the algorithm into two steps because we recommend that users try multiple values for the integer programming objective function parameters, as different values of the parameters will work best for different datasets.

## Contents
```@contents
Pages = ["installation.md", "example_decode.md", "api_reference.md"]
Depth = 3
```
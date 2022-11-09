# SeqFISHSyndromeDecoding.jl


## Introduction
This package provides functions to decode seqFISH points using [syndrome decoding](https://en.wikipedia.org/wiki/Decoding_methods#Syndrome_decoding).

In addition to the dynamic programming and integer programming algorithms described in our preprint, this package implements a simulated annealing solver
to resolve conflicting barcode candidates following [*Simulated Annealing: Theory and Applications, 1987*](https://books.google.com/books/about/Simulated_Annealing_Theory_and_Applicati.html?id=-IgUab6Dp_IC) by P.J. van Laarhoven, E.H. Aarts. For
the work described in our manuscript, we opted to use a commercial solver for all optimizations since it is faster. We still include our simulated annealing
solution here since it may be useful to some users. 


## Contents
```@contents
Pages = ["installation.md", "example_decode.md", "api_reference.md"]
Depth = 3
```
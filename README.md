

# SeqFISHSyndromeDecoding

This julia package implements the decoding method described in our preprint, "[Untangling overlapping barcodes in spatial genomics data](https://doi.org/10.1101/2025.06.10.658913)." The main repository for the preprint is available [here](https://github.com/CaiGroup/UntanglingBarcodes) and includes jupyter notebooks for [generating codes](https://github.com/CaiGroup/UntanglingBarcodes/tree/main/codebook_generation).

## Documentation
Documenatation is available [here](https://caigroup.github.io/SeqFISHSyndromeDecoding.jl/).

## Installation
Julia package for decoding SeqFISH experiments using [syndrome decoding](https://en.wikipedia.org/wiki/Decoding_methods#Syndrome_decoding). To install, navigate to the package directory in a terminal window, enter a Julia session, then type:

```
using Pkg
Pkg.add(".")
```

Alternativesly, open the Julia package manager by typing:

```
]
```

then type:

<pre> <code> add . </code> </pre>

and press enter. Exit the Julia package manager by pressing the backspace key.

A third installation is to type from a regular julia session:
```
using Pkg
Pkg.add("https://github.com/CaiGroup/SeqFISHSyndromeDecoding")
```

# Example notebooks

Example notebooks for [decoding data](https://colab.research.google.com/github/CaiGroup/SeqFISHSyndromeDecoding.jl/blob/master/example_notebook/colab/example_decode_colab.jl.ipynb) from the original [RNA SeqFISH+ paper](https://doi.org/10.1038/s41586-019-1049-y) and for [decoding Reed-Solomon encoded data](https://colab.research.google.com/github/CaiGroup/SeqFISHSyndromeDecoding.jl/blob/master/example_notebook/colab/example_decode_RS_colab.jl.ipynb) from our [preprint](https://doi.org/10.1101/2025.06.10.658913).

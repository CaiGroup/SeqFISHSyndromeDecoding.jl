push!(LOAD_PATH,"../src/")
using Documenter
using SeqFISHSyndromeDecoding

_PAGES = [
    "index.md",
    "installation.md",
    "example_decode.md",
    "example_decode_RS.md",
    "api_reference.md"
]

makedocs(
    sitename = "SeqFISHSyndromeDecoding",
    format = Documenter.HTML(prettyurls=false),
    modules = [SeqFISHSyndromeDecoding],
    pages = _PAGES
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#

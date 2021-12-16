push!(LOAD_PATH,"../src/")
using Documenter
using SeqFISHSyndromeDecoding

makedocs(
    sitename = "SeqFISHSyndromeDecoding",
    format = Documenter.HTML(),
    modules = [SeqFISHSyndromeDecoding]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#

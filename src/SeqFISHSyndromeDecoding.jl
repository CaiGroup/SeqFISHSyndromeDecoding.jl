module SeqFISHSyndromeDecoding

export decode_syndromes!, DecodeParams, get_codepaths, choose_optimal_codepaths
include("decode_params.jl")
include("coding_types.jl")
include("syndrome_decoding.jl")
include("01_MIP_solve.jl")
include("mem_eff_synd_comp.jl")
include("allow_blank_block_codewords.jl")
include("PrimePowerSyndromeComponent.jl")

end # module

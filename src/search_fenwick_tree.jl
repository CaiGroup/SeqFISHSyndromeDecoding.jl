using DataStructures

"""
    search(ft :: FenwickTree{Float64}, draw_val :: Float64)

Do a binary search for a bin in a binary index (fenwick) tree.
"""

function search(ft :: FenwickTree, val)
    @assert typeof(val) == eltype(ft)

    search_bit= 2^floor(Int, log2(ft.n))
    sig_bits = 0
    sum = 0

    while search_bit > 0
        search_ind = search_bit + sig_bits
        if search_ind < length(ft)
            search_val = sum + ft.bi_tree[search_ind]
            if val > search_val
                sig_bits = search_ind
                sum = search_val
            end
        end

        search_bit = search_bit >> 1
    end
    sig_bits + 1
end

using Test

#include("search_fenwick_tree.jl")
using SeqFISHSyndromeDecoding: search
using DataStructures

ft = FenwickTree(ones(12))

res = search(ft, 1.5)

println("res: $res")

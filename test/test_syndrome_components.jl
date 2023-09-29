

using SeqFISHSyndromeDecoding: ℤnRingElem, set_q, set_n, set_H, SyndromeComponent, check_mpath_decodable, get_decode_table
using SeqFISHSyndromeDecoding
using Test
using DelimitedFiles

Eng2019_ontarget = readdlm("Eng2019_647.csv", ',', UInt8)
RS_q5_k2_cb = readdlm("RS_q5_k2_cb.csv", ',', UInt8)
RS_q5_k2_H = readdlm("RS_q5_k2_H.csv", ',', UInt8)

RS_q7_k4_H = readdlm("RS_q7_k4_H.csv", ',', UInt8)
RS_q7_k4_w3cb = readdlm("RS_q7_k4_w3cb.csv", ',', UInt8)
RS_q7_k4_w4cb = readdlm("RS_q7_k4_w4cb.csv", ',', UInt8)
RS_q7_k4_w5cb = readdlm("RS_q7_k4_w5cb.csv", ',', UInt8)

hamming_merfish_cb = readdlm("hamming_merfish_cb.csv", ',', UInt8)
hamming_merfish_H = readdlm("hamming_merfish_H.csv", ',', UInt8)

cbs = [Eng2019_ontarget]
Hs = [[1 1 -1 -1;]]

cbs_zeros_unprobed = [RS_q5_k2_cb, RS_q7_k4_w3cb, RS_q7_k4_w4cb, RS_q7_k4_w5cb, hamming_merfish_cb]
Hs_zeros_unprobed = [RS_q5_k2_H, RS_q7_k4_H, RS_q7_k4_H, RS_q7_k4_H, hamming_merfish_H]

@testset "all syndrome type tests" begin

@testset "ℤnRingElem Arithmetic" begin
    for q = 0x01:0x14
        @testset "q = $q" begin
            for i = 0x00:q, j = 0x00:q
                set_q(q)
                a = ℤnRingElem(i)
                b = ℤnRingElem(j)
                c = a + b
                d = a * b
                @test c.v == (i + j) % q
                @test d.v == (i * j) % q
            end
        end
    end
end

@testset "Test parity check verification" begin
    for (i, cb) in enumerate(cbs)
        n = UInt8(length(cb[1,:]))
        pos = Array(0x01:n)
        q = UInt8(maximum(cb) + 1)
        set_q(q)
        set_n(n)
        set_H(Hs[i])
        for coeffs = eachrow(cb)
            syndrome = sum(SyndromeComponent.(coeffs, pos))
            @test iszero(syndrome)
        end
    end
end

@testset "test decode table drops true positive" begin
    for (i, cb) in enumerate(cbs)
        get_decode_table()
        n = UInt8(length(cb[1,:]))
        pos = Array(0x01:n)
        q = UInt8(maximum(cb) + 1)
        set_q(q)
        set_n(n)
        set_H(Hs[i])
        nrows = length(cb[:,1])
        for row = 1:nrows, drop_pos = pos
            coeffs = deepcopy(cb[row, :])
            dropped_coeff = coeffs[drop_pos]
            coeffs[drop_pos] = 0x00
            s = sum(SyndromeComponent.(coeffs, pos))
            result = check_mpath_decodable(0x02^(drop_pos-0x01), s)
            @test result.decodable && result.coeff == dropped_coeff
        end
    end
end


@testset "test decode table test no false positive" begin
    for (i, cb) in enumerate(cbs)
        ndecodable = 0
        n = UInt8(length(cb[1,:]))
        pos = Array(0x01:n)
        q = UInt8(maximum(cb) + 1)
        set_q(q)
        set_n(n)
        set_H(Hs[i])
        qm1 = q - 0x01
        nrows = length(cb[:,1])
        R = CartesianIndices(Tuple([0x00:qm1 for i = 1:(n-1)]))
        for drop_pos = pos
            for I in R
                message_tup = Tuple(I)
                message = [UInt8(symb) for symb in message_tup]
                insert!(message, drop_pos, 0x00)
                s = sum(SyndromeComponent.(message, pos))
                result = check_mpath_decodable(0x02^(drop_pos-0x01), s)
                if result.decodable
                    message[drop_pos] = result.coeff
                    ndecodable += 1
                end
            end
        end
        @test ndecodable == Int(n)*Int64(q)^(Int64(n)-1)
    end
end


@testset "Test parity check verification blank" begin
    for (i, cb) in enumerate(cbs_zeros_unprobed)
        n = UInt8(length(cb[1,:]))
        pos = Array(0x01:n)
        q = UInt8(maximum(cb) + 1)
        set_q(q)
        set_n(n)
        set_H(Hs_zeros_unprobed[i])
        for coeffs = eachrow(cb)
            syndrome = sum(SyndromeComponent.(coeffs, pos))
            @test iszero(syndrome)
        end
    end
end


"""
@testset "Q11N8SC decode table drops true positive" begin
    pos = collect(0x01: 0x08)
    for row = 1:1267, drop_pos = 0x01:0x08
        if q11_cb[row, drop_pos] != 0
            coeffs = deepcopy(q11_cb[row, :])
            dropped_coeff = coeffs[drop_pos]
            coeffs[drop_pos] = 0x00
            s = sum(Q11N8SC.(coeffs, pos))
            result = Q11N8SC_Mod.check_mpath_decodable(0x02^(drop_pos-0x01), s)
            @test result.decodable && result.coeff == dropped_coeff
        end
    end
end

zero_combinations = []
for i=1:6,j=(i+1):7,k=(j+1):8
    push!(zero_combinations, (i, j, k))
end
"""


end

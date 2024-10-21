
import Base: +, -, *, ^, ==, !=, /, iszero
using BKTrees

"""
    ℤnRingElem

Abstract type of a residue ring where elements are integers [0, n-1] and
the addition and multiplication operations are modulo n.
"""
struct ℤnRingElem
    v :: UInt8
end

# global variables are specific to the code. These are set by the function syndrome_decoding
q = 0x00
H = zeros(Int64, (2,2))
n = 0x00
k = 0x00


function set_q(_q :: UInt8)
    global q = _q
    global qm1 = q - 0x01
    init_ones_and_twos_masks()
end


function set_H(_H :: Matrix, params :: DecodeParams, cb)
    if ~params.zeros_probed
        if (q == 8 || q == 9)
            global H = FFExtElemExpForm.(_H)
            if params.ndrops > 0
                cws = [String.(collect(cw)) for cw in eachrow(cb)]
                global BKTree_cb = BKTree((x, y) -> sum(x .!= y), cws)
            end
        else
            global H = _H
            if params.ndrops > 0
                cws = [collect(cw) for cw in eachrow(cb)]
                global BKTree_cb = BKTree((x, y) -> sum(x .!= y), cws)
            end
        end
    else
        global H = _H
    end
    global n = UInt8(size(H)[2])
    global k = UInt8(n-size(H)[1])
    init_ones_and_twos_masks()
end

function init_ones_and_twos_masks()
    if q == 9 && n != 0x00
        #initials ones place and twos place masks for bitwise operations
        global ones_mask = vcat(fill(BitArray([1,0]),2*(n-k))...)
        global twos_mask = .~ones_mask
    end
end

add_mod(a :: ℤnRingElem, b :: ℤnRingElem, mod :: Unsigned) = (a.v + b.v) % mod

add_mod(a :: Unsigned, b :: ℤnRingElem, mod :: Unsigned) = (a + b.v) % mod

mult_mod(a :: ℤnRingElem, b :: ℤnRingElem, mod :: Unsigned) = (a.v * b.v) % mod

mult_mod(a :: Unsigned, b :: ℤnRingElem, mod :: Unsigned) = (a * b.v) % mod

inv₊(x :: ℤnRingElem) =  ℤnRingElem((q - x.v) % q)
inv₊(x :: UInt8) =  ℤnRingElem((q - x) % q)


+(a :: Unsigned, b :: ℤnRingElem) = typeof(b)(add_mod(a, b, q))

+(a :: ℤnRingElem, b :: Unsigned) = typeof(a)(add_mod(b, a, q))

+(a :: ℤnRingElem, b :: ℤnRingElem) = a + b.v

-(a :: ℤnRingElem, b :: ℤnRingElem) = a + inv₊(b)

*(a :: Unsigned, b :: ℤnRingElem) = typeof(b)(mult_mod(a, b, q))

*(a :: ℤnRingElem, b :: Unsigned) = typeof(a)(mult_mod(b, a, q))

*(a :: ℤnRingElem, b :: ℤnRingElem) = a * b.v

/(a :: ℤnRingElem, b :: ℤnRingElem) = ℤnRing_div_table[a, b]

/(a :: Unsigned, b :: ℤnRingElem) = ℤnRing_div_table[typeof(b)(a), b]

/(a :: ℤnRingElem, b :: Unsigned) = ℤnRing_div_table[a, typeof(a)(b)]

#function *(a :: ℤnRingElem, b :: ℤnRingElem)
    #@assert typeof(a) == typeof(b)
    #a*b.v
#end

function exp_mod(a :: ℤnRingElem, order :: Unsigned, q :: Unsigned)
    r = 0x01
    for i in 0x01:order
        r = (r*a.v)%q
    end
    typeof(a)(r)
end

^(a :: ℤnRingElem, order :: Unsigned) = exp_mod(a, order, q)

iszero(x :: ℤnRingElem) = iszero(x.v)

abstract type SyndromeComponent end

==(a :: SyndromeComponent, b :: SyndromeComponent) = a.s == b.s

!=(a :: SyndromeComponent, b :: SyndromeComponent) = a.s != b.s

struct SyndromeComponentℤnRing <: SyndromeComponent
    s :: Tuple{Vararg{ℤnRingElem}}
end

function SyndromeComponent(coeff :: UInt8, pos :: UInt8)
    pos_fncs = func_from_H_val.(H)
    res = [pos_fncs[i,pos](ℤnRingElem(coeff)) for i = 1:(Base.size(H)[1])]
    return SyndromeComponentℤnRing(Tuple(res))
end

function func_from_H_val(h_val)
    if h_val == 0
        return (x) -> x-x
    elseif h_val == 1
        return identity
    elseif h_val == -1
        return inv₊
    elseif typeof(h_val) == FFExtElemExpForm 
        return (x) -> (h_val) * x
    elseif h_val > 1 
        return (x) -> ℤnRingElem(UInt8(h_val)) * x
    else
        raise(error("parity check function value $h_val not supported."))
    end
end

+(a :: SyndromeComponentℤnRing, b :: SyndromeComponentℤnRing) = typeof(a)(a.s .+ b.s)

iszero(x :: SyndromeComponent) = all(iszero.(x.s))

inv₊(s :: SyndromeComponentℤnRing) = typeof(s)(inv₊.(s.s))

get_pos(hyb :: Integer) = ceil(UInt8, hyb / q)

get_coeff(hyb :: Integer, pos :: Integer) = UInt8((hyb - (pos - 0x01) * q) % q)


function get_decode_table()
    global decode_table = Dict{Tuple{UInt8, SyndromeComponentℤnRing},  UInt8}()
    #global decode_table = Dict{Tuple{UInt8, UInt8},  UInt8}()
    #for coeff = 0x01:0x13, pos = 0x01:0x04
    for coeff = 0x01:(q - 0x01), pos = 0x01:n #0x04
        s = SyndromeComponent(coeff, pos)

        #@assert (pos, inv₊(s)) ∉ keys(E2019_decode_table)
        pos_indicator = 0x02^(pos-0x01)
        inv_s = inv₊(s)
        #svs = ([re.v for re in s.s]...)

        @assert (pos_indicator, inv_s) ∉ keys(decode_table)
        #@assert (pos_indicator, svs) ∉ keys(decode_table)

        #E2019_decode_table[(pos, inv(s))] = coeff
        # use "one-hot" binary encoding for position in the key/
        # this makes syndrome decoding of message dag more efficient
        #decode_table[(0x02^(pos-0x01), inv₊(s))] = coeff
        decode_table[(pos_indicator, inv_s)] = coeff

        #decode_table[(pos_indicator, svs)] = coeff
    end
    decode_table
end


"""
    check_mpath_decodable(drop_pos :: UInt8, s :: SyndromeComponentℤnRing)

Checks to see if a message path with a drop can be decoded.
Attempts to fill in the drop with syndrome decoding and decoding table.

If the drop is correctable, returns the pseudocolor/coeff of the dropped encoding block

otherwise, returns false
"""
function check_mpath_decodable(drop_pos :: UInt8, s :: SyndromeComponentℤnRing)
    if iszero(s)
        return (decodable = true, coeff = 0x00)
    end
    missing_coeff = get(decode_table, (drop_pos, s), "Undecodable")
    #missing_coeff = get(decode_table, (drop_pos, s.v), "Undecodable")
    if missing_coeff != "Undecodable"
        return (decodable = true, coeff = missing_coeff)
    else
        return (decodable = false, coeff = nothing)
    end
end

function make_div_table_Prime_Field(q)
    global ℤnRing_div_table = ℤnRingElem.(zeros(UInt8,q,q))
    for i = 0x00:UInt8(q-1), j = i:UInt8(q-1)
        res = ℤnRingElem(i)*ℤnRingElem(j)
        ℤnRing_div_table[res,i+1] = ℤnRingElem(j)
        ℤnRing_div_table[res,j+1] = ℤnRingElem(i)
    end
end



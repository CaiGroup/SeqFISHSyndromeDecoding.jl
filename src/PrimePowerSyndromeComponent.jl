
const savestring_2_exp = Dict(
    [
        ("0",Int8(-1)),
        ("1",Int8(0)),
        ("α1",Int8(1)),
        ("α",Int8(1)),
        ("α2",Int8(2)),
        ("α²",Int8(2)),
        ("α3",Int8(3)),
        ("α³",Int8(3)),
        ("α4",Int8(4)),
        ("α⁴",Int8(4)),
        ("α5",Int8(5)),
        ("α⁵",Int8(5)),
        ("α6",Int8(6)),
        ("α⁶",Int8(6)),
        ("α7",Int8(7)),
        ("α⁷",Int8(7))
    ]
)

const pseudocolor_2_savestring = Dict(
    [
        (0, "0"),
        (1, "1"),
        (2, "α1"),
        (2, "α"),
        (3, "α2"),
        (3, "α²"),
        (4, "α3"),
        (4, "α³"),
        (5, "α4"),
        (5, "α⁴"),
        (6, "α5"),
        (6, "α⁵"),
        (7, "α6"),
        (7, "α⁶"),
        (8, "α7"),
        (8, "α⁷")
    ]
)


if q == 9
    global p = 3
elseif q == 8
    global p = 2
end

function SyndromeComponent(coeff :: String, pos :: UInt8)
    if typeof(H[1,pos]) <: AbstractString
        res = [FFExtElemExpForm(H[i,pos])*FFExtElemExpForm(coeff) for i = 1:(size(H)[1])]
    else
        res = [H[i,pos]*FFExtElemExpForm(coeff) for i = 1:(size(H)[1])]
    end
    if q == 8
        return SyndromeComponentExtensionField2(res)
    elseif q == 9
        return SyndromeComponentExtensionField3(res)
    end
end

struct FFExtElemExpForm
    pow :: Int8
end

FFExtElemExpForm(elem :: String) = FFExtElemExpForm(savestring_2_exp[elem])


function *(a :: FFExtElemExpForm, b :: FFExtElemExpForm)
    if a.pow == Int8(-1) || b.pow == Int8(-1)
        return FFExtElemExpForm(Int8(-1))
    else
        return FFExtElemExpForm( (a.pow + b.pow) %  qm1)
    end
end

struct SyndromeComponentExtensionField2 <: SyndromeComponent
    s :: BitArray #Tuple{Vararg{Tuple{Vararg{ℤnRingElem}}}}
end

function SyndromeComponentExtensionField2(input :: Vector{FFExtElemExpForm})
    SyndromeComponentExtensionField2(vcat(to_ff2_3_poly_form.(input)...))
end

+(a :: SyndromeComponentExtensionField2, b :: SyndromeComponentExtensionField2) = SyndromeComponentExtensionField2(a.s .⊻ b.s)

struct FF2ExtensionElemPolyForm 
    v :: BitArray
end

struct SyndromeComponentExtensionField3 <: SyndromeComponent
    s :: BitArray #Tuple{Vararg{Tuple{Vararg{ℤnRingElem}}}}
end

function SyndromeComponentExtensionField3(input :: Vector{FFExtElemExpForm})
    SyndromeComponentExtensionField3(vcat(to_ff3_2_poly_form.(input)...))
end

function +(a :: SyndromeComponentExtensionField3, b :: SyndromeComponentExtensionField3)
    and_res = a.s .& b.s
    #and_swap = ((and_res << 1) .& twos_mask) .| ((and_res >> 1) .& ones_mask)
    and_swap = ((and_res >> 1) .& twos_mask) .| ((and_res << 1) .& ones_mask)
    xor_res = a.s .⊻ b.s
    not_or_res = .~(a.s .| b.s)
    #ones_or_res = ((xor_res .& (not_or_res >> 1)) .& ones_mask)
    #twos_or_res = ((xor_res .& (not_or_res << 1)) .& twos_mask)
    ones_or_res = ((xor_res .& (not_or_res << 1)) .& ones_mask)
    twos_or_res = ((xor_res .& (not_or_res >> 1)) .& twos_mask)
    or_res = ones_or_res .| twos_or_res
    SyndromeComponentExtensionField3(or_res .| and_swap)
end

FF2_ext3_poly_form_dict = Dict([
        (Int8(-1), BitArray([0,0,0])),
        (Int8(0), BitArray([1,0,0])),
        (Int8(1),BitArray([0,1,0])),
        (Int8(2), BitArray([0,0,1])),
        (Int8(3), BitArray([1,1,0])),
        (Int8(4), BitArray([0,1,1])),
        (Int8(5), BitArray([1,1,1])),
        (Int8(6), BitArray([1,0,1]))
    ]
)

FF3_ext2_poly_form_dict = Dict([
    (Int8(-1), BitArray([0,0,0,0])), # 0 + 0*α
    (Int8(0), BitArray([1,0,0,0])),  # 1 + 0*α
    (Int8(1),BitArray([0,0,1,0])),   # 0 + 1*α
    (Int8(2), BitArray([1,0,1,0])),  # 1 + 1*α
    (Int8(3), BitArray([1,0,0,1])),  # 1 + 2*α
    (Int8(4), BitArray([0,1,0,0])),  # 2 + 0*α
    (Int8(5), BitArray([0,0,0,1])),  # 0 + 2*α
    (Int8(6), BitArray([0,1,0,1])),  # 2 + 2*α
    (Int8(7), BitArray([0,1,1,0])),  # 2 + 1*α
    ]
)

to_ff3_2_poly_form(elem ::FFExtElemExpForm) = FF3_ext2_poly_form_dict[elem.pow]
to_ff2_3_poly_form(elem ::FFExtElemExpForm) = FF2_ext3_poly_form_dict[elem.pow]
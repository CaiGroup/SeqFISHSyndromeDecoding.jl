#module simulation_fncs

#export sim_true, encode, drop_random_dots!, draw_off_target_hyb_dots!, draw_rand_dots,
#draw_localization_errors, add_localization_errors!, score_decoding, score_dir

using DataFrames
using Distributions
using CSV

"""
    sim_true(n)

Simulate n objects in a 1 x 1 field of view. The species of each object is a
random integer from 1 to 1267. It's the location is drawn from the uniform
distribution in the field of view. Returns a data frame of each object's species
x, and y position.
"""
sim_true(n) = DataFrame(species = rand(1:1267, n), x = rand(n), y = rand(n))

sim_true(n, max_species) = DataFrame(species = rand(1:max_species, n), x = rand(n), y = rand(n))

function get_n_q_w(cb)
    ncws, n = size(cb)
    q = length(unique(cb))
    cw_nonzeros = [sum(cb[i, :] .!= -) for i in ncws]
    w =  maximum(cw_nonzeros)
    #w = sum(cb[1,:] .!= 0)
    [n, q, w]
end

"""
    encode(true_locs :: DataFrame, cb :: Array{Int8,2})

Given a data frame of true locations, and a codebook, return a dataframe giving
the hybridization of the encoded dots, their species, their true locations, and
coefficient.
"""
function encode(true_locs :: DataFrame, cb)
    n, q, wt = get_n_q_w(cb)
    n = UInt8(n)
    q = UInt8(q)
    wt = UInt8(wt)
    #n = UInt8(length(cb[1,:]))
    #wt = sum(cb[1,:] .!= 0)
    #q = UInt8(length(unique(cb)))

    qm1 = UInt8(q-1)

    npts = nrow(true_locs) #length(true_locs.species)
    dot = Array(1:(npts*n))
    species = Int64[]
    x = Float64[]
    y = Float64[]
    #hyb = Int8[]
    hyb = UInt8[]
    for i = 1:npts
        target = true_locs[i, :]
        codeword = cb[target.species,:]
        #for pos = 1:n
        for pos = 0x01:n
            if codeword[pos] != 0 && n > wt
                push!(species, target.species)
                push!(x, target.x)
                push!(y, target.y)
                dot_hyb = qm1*(pos-1)+codeword[pos]
                push!(hyb, dot_hyb)
                #push!(hyb, UInt8(dot_hyb))
            elseif n == wt
                push!(species, target.species)
                push!(x, target.x)
                push!(y, target.y)
                dot_hyb = q*(pos-1)+codeword[pos]
                #codeword[pos] == 0x00 ? dot_hyb += 0x10 : nothing
                codeword[pos] == 0x00 ? dot_hyb += UInt8(q) : nothing
                push!(hyb, dot_hyb)
                #push!(hyb, UInt8(dot_hyb))
            end
        end
    end
    pnts = DataFrame(dot_ID=dot,hyb=hyb,species=species,x=x,y=y)
    pnts[!, "z"] .= 1.0
    pnts[!, "w"] .= 1.0
    pnts[!, "s"] .= 1.0
    pnts
end


"""
    drop_random_dots!(pnts :: DataFrame, drop_rate)

Adds encoding errors to a simulated dataset by dropping dots at rate drop_rate
"""
function drop_random_dots!(pnts :: DataFrame, drop_rate)
    @assert drop_rate ≥ 0
    @assert drop_rate < 1

    to_drop = rand(length(pnts.x)) .≤ drop_rate
    delete!(pnts, to_drop)
end

"""
    draw_off_target_hyb_dots!(true_locs :: DataFrame, off_target_rate)

Draw the objects that should appear in off target hybs and in which off target
hybs they should appear in.
"""
function draw_off_target_hyb_dots!(true_locs :: DataFrame, off_target_rate)
    @assert off_target_rate ≥ 0
    n_objects = length(true_locs)

    objects = rand(1:n_objects, round(off_target_rate*n_objects))
    hybs = rand(1:80, length(off_target_objects))
    x = true_locs.x[off_target_objects]
    y = true_locs.y[off_target_objects]

    return DataFrame(species = -objects, hyb = hybs, x = x, y = y)
end

"""
    draw_rand_dots(n)

Draw random off target dots that simulate probes sticking to some unintended
object that is not of interest.
"""
draw_rand_dots(n) = DataFrame(species=0, hyb=rand(1:80), x=rand(n), y=rand(n))


"""
    draw_localization_errors(n, rstdv)

Draw vectors of x and y localization errors for n points.
"""
function draw_localization_errors(n, rstdv)
    r_error = rand(Normal(0,rstdv), n)
    θ = 2π*rand(Float64, n)

    x_error = r_error.*sin.(θ)
    y_error = r_error.*cos.(θ)
    return x_error, y_error
end

"""
    add_localization_errors(pnts :: DataFrame, x_errors, y_errors)
"""
function add_localization_errors!(pnts :: DataFrame, x_errors, y_errors)
    pnts.x .+= x_errors
    pnts.y .+= y_errors
end

"""
    add_localization_errors(pnts :: DataFrame)
"""
function add_localization_errors!(pnts :: DataFrame, rstdv)
    x_errors, y_errors = draw_localization_errors(length(pnts.x), rstdv)
    add_localization_errors!(pnts, x_errors, y_errors)
end


"""
    score_decoding(true_locs :: DataFrame, decoded :: DataFrame)

Dates in a data frame of true locations that generated a simulation, and the
decoded results. scores the decoding based on the number target species counts it
correctly and incorrectly identified.

Returns:
    n_correct: number of decoded gene counts that agree with the simulated true number of counts
    n_error: number of decoded gene counts that disagree with the simulated true number of counts
"""
function score_decoding(true_locs :: DataFrame, decoded :: DataFrame)

    max_species = Int(maximum([maximum(true_locs.species), maximum(decoded.decoded)]))
    # get true count vector
    #true_counts = fill(0, 1267)
    true_counts = fill(0, max_species)
    for species_ID in true_locs.species
        true_counts[species_ID] += 1
    end

    # get decoded count vector
    last_mpath = Int[]
    sort!(decoded, :mpath)
    #decode_counts = fill(0, 1267)
    decode_counts = fill(0, max_species)
    for i in 1:length(decoded.mpath)
        if decoded.mpath[i] != last_mpath && decoded.decoded[i] != 0
            decode_counts[Int(decoded.decoded[i])] += 1
            last_mpath = decoded.mpath[i]
        end
    end

    #find difference between true and decoded count vectors
    diff = true_counts .- decode_counts

    # negative components of the difference vector represent incorrect counts
    # the number of errors is the absolute value of the sum the negative components
    # of the difference vector
    #n_error = abs(sum(diff[diff .< 0]))

    n_false_positive = abs(sum(diff[diff .< 0]))
    n_false_negative = abs(sum(diff[diff .> 0]))

    # the number of correct counts is the total sum of decoded counts minus
    # the number of  incorrectly decoded counts found above.
    #n_correct = sum(decode_counts) - n_error

    n_targets = length(true_locs.species)

    n_correct = n_targets - n_false_negative

    return n_targets, n_correct, n_false_positive#n_error
end

"""
    score_decoding(true_locs :: String, decoded :: String)

Wrapper for above definition of score decoding that takes in file path for input files,
loads them, and then inputs them into the above definition.
"""
function score_decoding(true_locs_path :: String, decoded_path :: String)
    true_locs_df = DataFrame(CSV.File(true_locs_path))
    decoded_df = DataFrame(CSV.File(decoded_path))
    score_decoding(true_locs_df, decoded_df)
end

"""

"""
function score_dir(true_locs_dir_path :: String, decoded_dir_path :: String)

    # get lists of files to score
    orig_dir = pwd()
    cd(true_locs_dir_path)
    true_locs_paths = true_locs_dir_path .* readdir()
    cd(orig_dir)
    cd(decoded_dir_path)
    decoded_paths = decoded_dir_path .* readdir()
    cd(orig_dir)

    # initialize arrays to store scores
    n_files = length(true_locs_paths)
    file_n_correct = zeros(Int, n_files)
    file_n_error = zeros(Int, n_files)
    file_n_dots = zeros(Int, n_files)
    file_n_targets = zeros(Int, n_files)

    # get scores
    for i = 1:n_files
        n_targets, n_correct, n_error = score_decoding(true_locs_paths[i], decoded_paths[i])
        file_n_correct[i] = n_correct
        file_n_error[i] = n_error
        file_n_targets[i] = n_targets
    end

    scores = DataFrame(n_targets = file_n_targets,
                        n_correct = file_n_correct,
                        n_error = file_n_error)

    return scores
end

#end

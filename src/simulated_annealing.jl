using DataFrames
using DataStructures
#using Distributions
using Statistics: var
using SparseArrays
using LinearAlgebra
using LightGraphs

include("search_fenwick_tree.jl")



"""
    simulated_annealing!(pnts :: DataFrame, sa_params)

Partitions connected message dag into message paths using simulaterd annealing.
"""
#=
function simulated_annealing!(pnts :: DataFrame, sa_params)
    cpaths, cpath_costs = get_cpath_list(pnts)
    simulated_annealing!(cpaths, cpath_costs, pnts, sa_params)
end
=#

#=
function simulated_annealing!(cpaths,# :: Vector{Vector{Int}},
                              cpath_costs,# :: Vector{Float64},
                              A :: SimpleGraph,
                              pnts :: DataFrame,
                              sa_params)

=#

function simulated_annealing(cpaths_df,
                              A :: SimpleGraph,
                              sa_params,
                              ndots)

    simulated_annealing_rejectionless(cpaths_df, A, sa_params, ndots)
    #simulated_annealing_metropolis(cpaths_df, A, sa_params)
end

function simulated_annealing(cpaths_df :: DataFrame,
    cpath_conflict_cpath_indices :: Vector{Vector{Int64}},
    cpath_partial_conflict_indices :: Vector{Vector{Int64}},
    cpath_partial_conflict_transition_indices :: Vector{Vector{Int64}},
                             sa_params,
                             ndots :: Int)
    println("Begin Simulated Annealing")
    println("n cpaths: ", nrow(cpaths_df))
    simulated_annealing_rejectionless_disrupt_subcpath(cpaths_df, cpath_conflict_cpath_indices, cpath_partial_conflict_indices,
    cpath_partial_conflict_transition_indices, sa_params, ndots)
    #simulated_annealing_metropolis(cpaths_df, A, sa_params)
end

function simulated_annealing_rejectionless_disrupt_subcpath(cpaths_df :: DataFrame,
                                                            cpath_conflict_cpath_indices :: Vector{Vector{Int64}},
                                                            cpath_partial_conflict_indices :: Vector{Vector{Int64}},
                                                            cpath_partial_conflict_transition_indices :: Vector{Vector{Int64}},
                                                            sa_params,
                                                            ndots :: Int)
    # number each code path by increasing cost
    ncpaths = nrow(cpaths_df)
    #if nrow(cpaths_df) == 0
    if ncpaths == 0
        return
    end

    state = state_rejectionless_smtrns(cpaths_df, cpath_conflict_cpath_indices, cpath_partial_conflict_indices,
    cpath_partial_conflict_transition_indices, sa_params)

    lchain = round(Int, sa_params.l_chain * ndots * log(ndots+1))

    #no_change_streak = 0
    max_iters = true
    for mc in 1:sa_params.n_chains
        # set the control parameter
        c = sa_params.cooling_factor/(1.0 + sa_params.cooling_timescale*log(mc))

        #println("Starting chain $mc")
        reset_chain!(state, c)
        expected_metropolis_itertions = 0

        #for i = 1:(sa_params.l_chain)
        i = 0
        while expected_metropolis_itertions < lchain# sa_params.l_chain

            draw = rand()*state.sum_weights
            tcp = search(state.wts_ft, draw)

            # remove transition costs of disrupting the disrupted cpaths
            transition_update!(state, tcp, c, sa_params)

            check_lowest_cost!(state)

            expected_metropolis_itertions += expected_metropolis_repetitions(state)

            if termination_condition(state, sa_params)
                max_iters = false
                return state.state
            end
            i += 1
        end
    end
    println("maximum iterations: $lchain")
    #end
    state.state
end


mutable struct state_rejectionless_smtrns
    state :: BitArray{1}
    weights :: Vector{Float64}
    wts_ft :: FenwickTree
    sum_weights :: Float64
    current_cost :: Float64
    lowest_cost :: Float64
    lowest_cost_state :: BitArray{1}
    ΔCost_transition_cpath :: Vector{Float64}
    cpath_conflict_cpath_indices :: Vector{Vector{Int64}}
    cpath_partial_conflict_indices :: Vector{Vector{Int64}}
    cpath_partial_conflict_transition_indices :: Vector{Vector{Int64}}
    ΔCost_freedots_to_cpath :: Vector{Float64}
    cpath_lengths :: Vector{Int}
    ndots :: Int
end

function state_rejectionless_smtrns(cpaths_df :: DataFrame,
                                    cpath_conflict_cpath_indices :: Vector{Vector{Int64}},
                                    cpath_partial_conflict_indices :: Vector{Vector{Int64}},
                                    cpath_partial_conflict_transition_indices :: Vector{Vector{Int64}},
                                    sa_params)
    # number each code path by increasing cost
    n_cpaths = nrow(cpaths_df)

    if n_cpaths == 0
        return
    end

    cpath_lengths = length.(cpaths_df.cpath)
    dots = unique(vcat(cpaths_df.cpath...))
    ndots = length(dots)#nrow(pnts)

    ΔCost_freedots_to_cpath = cpaths_df.cost .- cpath_lengths.*sa_params.free_dot_cost


    # initialize state vector to all zeros
    state = falses(n_cpaths)

    #current_cost = sa_params.free_dot_cost*nrow(pnts)
    current_cost = sa_params.free_dot_cost*ndots
    lowest_cost = current_cost
    lowest_cost_state = falses(n_cpaths)

    #prealloacate, defined in reset_chain!
    ΔCost_transition_cpath = zeros(Float64, n_cpaths)#deepcopy(ΔCost_freedots_to_cpath)
    weights = zeros(Float64, n_cpaths) #exp.(-ΔCost_transition_cpath./c)
    wts_ft = FenwickTree(weights)
    sum_weights = 0.0

    state_rejectionless_smtrns(
        state,
        weights,
        wts_ft,
        sum_weights,
        current_cost,
        lowest_cost,
        lowest_cost_state,
        ΔCost_transition_cpath,
        cpath_conflict_cpath_indices,
        cpath_partial_conflict_indices,
        cpath_partial_conflict_transition_indices,
        ΔCost_freedots_to_cpath,
        cpath_lengths,
        ndots
    )
end

function transition_update!(s :: state_rejectionless_smtrns, tcp :: Int, c :: Float64, sa_params)
    update_weights_costs!(s, s.cpath_conflict_cpath_indices, tcp, c, sa_params)
    #transition_partial_conflict(s, tcp, c, sa_params)

end


function transition_partial_conflict(s :: state_rejectionless_smtrns, tcp, c, sa_params)
    #update_weights_costs!(s, s.cpath_partial_conflict_indices, tcp, c, sa_params)

    # ToDO: Something is wrone with this transition subroutine

    for (i, subcpath_ind) in enumerate(s.cpath_partial_conflict_indices[tcp])
        # add not conflicting subpath to the state
        if tcp == subcpath_ind
            update_weights_costs!(s, s.cpath_conflict_cpath_indices, s.cpath_partial_conflict_transition_indices[i], c, sa_params)
        end
    end


end


function update_weights_costs!(s :: state_rejectionless_smtrns, nbrs ::Vector{Vector{Int}}, tcp :: Int, c :: Float64, sa_params)
    update_neighbor_weights!(s, nbrs, tcp, c)

    #s.current_cost += s.ΔCost_freedots_to_cpath[tcp]
    s.state[tcp] = true

    s.ΔCost_transition_cpath[tcp] = 0
    s.current_cost = sa_params.free_dot_cost*s.ndots + s.state⋅s.ΔCost_freedots_to_cpath

    new_weight = 0
    diff_weight = 0 - s.weights[tcp]
    inc!(s.wts_ft, tcp, diff_weight)
    s.weights[tcp] = new_weight
    s.sum_weights += diff_weight
end

function update_neighbor_weights!(s :: state_rejectionless_smtrns, nbrs :: Vector{Vector{Int}}, tcp :: Int, c :: Float64)

    # remove transition costs of disrupting the disrupted cpaths
    for d_ind in nbrs[tcp]#neighbors(s.A, tcp)
        if s.state[d_ind] # if the conflicting cpath is in the state
            # remove the cost of disrupting the disrupted cpath for each of its neighbors
            for nbr in nbrs[d_ind]#neighbors(s.A, d_ind)
                s.ΔCost_transition_cpath[nbr] += s.ΔCost_freedots_to_cpath[d_ind]
                #current_cost += ΔCost_freedots_to_cpath[d_ind]
                #println("update current cost remove disruption cost: $current_cost")
                new_weight = exp(-s.ΔCost_transition_cpath[nbr]/c)
                diff_weight = new_weight - s.weights[nbr]
                #println("neighbor of disrupted diff_weight: $diff_weight")

                inc!(s.wts_ft, nbr, diff_weight)

                s.weights[nbr] = new_weight
                s.sum_weights += diff_weight
            end
            #update the state and transition cost of the disrupted cpath
            s.state[d_ind] = 0
            s.ΔCost_transition_cpath[d_ind] += s.ΔCost_freedots_to_cpath[d_ind] - s.ΔCost_freedots_to_cpath[tcp]
            #s.current_cost -= s.ΔCost_freedots_to_cpath[d_ind]
            #println("increased current cost: $current_cost")
            #println("update current cost disruption: $current_cost")

            # update sum weights and fenwick tree weights
            new_weight = exp(-s.ΔCost_transition_cpath[d_ind]/c)
            diff_weight = new_weight - s.weights[d_ind]
            inc!(s.wts_ft, d_ind, diff_weight)
            s.weights[d_ind] = new_weight
            s.sum_weights += diff_weight
            #println("disrupted diff_weight: $diff_weight")

        # if the neighbor is not in the the state, update the cost of transitioning to it to include the cost of disrupting the new codepath in the state
        else
            s.ΔCost_transition_cpath[d_ind] -= s.ΔCost_freedots_to_cpath[tcp]

            new_weight = exp(-s.ΔCost_transition_cpath[d_ind]/c)
            diff_weight = new_weight - s.weights[d_ind]
            inc!(s.wts_ft, d_ind, diff_weight)
            s.weights[d_ind] = new_weight
            s.sum_weights += diff_weight
            #println("disrupted diff_weight: $diff_weight")
        end
    end
end

function reset_chain!(s :: state_rejectionless_smtrns, c :: Float64)
    #reset ΔCost_transition_cpath
    s.ΔCost_transition_cpath .= s.ΔCost_freedots_to_cpath
    for (cpath_ind, cpath_bool) in enumerate(s.state)
        if cpath_bool
            for d_ind in s.cpath_conflict_cpath_indices[cpath_ind]#neighbors(s.A, cpath_ind)
                s.ΔCost_transition_cpath[d_ind] -= s.ΔCost_freedots_to_cpath[cpath_ind]
            end
            s.ΔCost_transition_cpath[cpath_ind] = 0
        end
    end

    s.weights = exp.(-s.ΔCost_transition_cpath./c)
    s.weights[s.state] .= 0
    s.wts_ft = FenwickTree(s.weights)
    s.sum_weights = sum(s.weights)
end

function check_lowest_cost!(s :: state_rejectionless_smtrns)
    if s.current_cost < s.lowest_cost
        s.lowest_cost = s.current_cost
        s.lowest_cost_state .= s.state
    end
end

function expected_metropolis_repetitions(s :: state_rejectionless_smtrns)
    maximum([1, length(s.state)/s.sum_weights])
end


# Free all dots from disrupted codepaths

function simulated_annealing_rejectionless(cpaths_df,
                              A :: SimpleGraph,
                              sa_params,
                              ndots)
    # number each code path by increasing cost
    ncpaths = nrow(cpaths_df)
    #if nrow(cpaths_df) == 0
    if ncpaths == 0
        return
    end

    state = state_rejectionless(cpaths_df, A, sa_params)

    lchain = round(Int, sa_params.l_chain * ndots * log(ndots+1))

    #no_change_streak = 0
    max_iters = true
    for mc in 1:sa_params.n_chains
        # set the control parameter
        c = sa_params.cooling_factor/(1.0 + sa_params.cooling_timescale*log(mc))

        #println("Starting chain $mc")
        reset_chain!(state, c)
        expected_metropolis_itertions = 0

        #for i = 1:(sa_params.l_chain)
        i = 0
        while expected_metropolis_itertions < lchain# sa_params.l_chain

            draw = rand()*state.sum_weights
            tcp = search(state.wts_ft, draw)

            # remove transition costs of disrupting the disrupted cpaths
            transition_update!(state, tcp, c, sa_params)

            check_lowest_cost!(state)

            expected_metropolis_itertions += expected_metropolis_repetitions(state)

            if termination_condition(state, sa_params)
                max_iters = false
                return state.state
            end
            i += 1
        end
        #println("expected metroplis iterations: ", expected_metropolis_itertions)
        #println("expected proportion rejected: ", (expected_metropolis_itertions-i)/expected_metropolis_itertions)
        #println("expected prop rejections:", expected_rejections/(sa_params.l_chain+expected_rejections))
    end
    println("maximum iterations: $lchain")
    #end
    state.state
end


mutable struct state_rejectionless
    state :: BitArray{1}
    weights :: Vector{Float64}
    wts_ft :: FenwickTree
    sum_weights :: Float64
    current_cost :: Float64
    lowest_cost :: Float64
    lowest_cost_state :: BitArray{1}
    ΔCost_transition_cpath :: Vector{Float64}
    A :: SimpleGraph
    ΔCost_freedots_to_cpath :: Vector{Float64}
    cpath_lengths :: Vector{Int}
    ndots :: Int
end

function state_rejectionless(cpaths_df :: DataFrame,
                              A :: SimpleGraph,
                              sa_params)
    # number each code path by increasing cost
    n_cpaths = nrow(cpaths_df)

    if n_cpaths == 0
        return
    end

    cpath_lengths = length.(cpaths_df.cpath)
    dots = unique(vcat(cpaths_df.cpath...))
    ndots = length(dots)#nrow(pnts)

    ΔCost_freedots_to_cpath = cpaths_df.cost .- cpath_lengths.*sa_params.free_dot_cost


    # initialize state vector to all zeros
    state = falses(n_cpaths)

    #current_cost = sa_params.free_dot_cost*nrow(pnts)
    current_cost = sa_params.free_dot_cost*ndots
    lowest_cost = current_cost
    lowest_cost_state = falses(n_cpaths)

    #prealloacate, defined in reset_chain!
    ΔCost_transition_cpath = zeros(Float64, n_cpaths)#deepcopy(ΔCost_freedots_to_cpath)
    weights = zeros(Float64, n_cpaths) #exp.(-ΔCost_transition_cpath./c)
    wts_ft = FenwickTree(weights)
    sum_weights = 0.0

    state_rejectionless(
        state,
        weights,
        wts_ft,
        sum_weights,
        current_cost,
        lowest_cost,
        lowest_cost_state,
        ΔCost_transition_cpath,
        A,
        ΔCost_freedots_to_cpath,
        cpath_lengths,
        ndots
    )
end

function transition_update!(s :: state_rejectionless, tcp :: Int, c :: Float64, sa_params)
    update_weights_costs!(s, tcp, c, sa_params)
end

function update_weights_costs!(s :: state_rejectionless, tcp :: Int, c :: Float64, sa_params)
    update_neighbor_weights!(s, tcp, c)

    #s.current_cost += s.ΔCost_freedots_to_cpath[tcp]
    s.state[tcp] = true

    s.ΔCost_transition_cpath[tcp] = 0
    s.current_cost = sa_params.free_dot_cost * s.ndots + s.state ⋅ s.ΔCost_freedots_to_cpath

    new_weight = 0
    diff_weight = 0 - s.weights[tcp]
    inc!(s.wts_ft, tcp, diff_weight)
    s.weights[tcp] = new_weight
    s.sum_weights += diff_weight
end

function update_neighbor_weights!(s :: state_rejectionless, tcp :: Int, c :: Float64)
    # remove transition costs of disrupting the disrupted cpaths
    for d_ind in neighbors(s.A, tcp)
        if s.state[d_ind] # if the conflicting cpath is in the state
            # remove the cost of disrupting the disrupted cpath for each of its neighbors
            for nbr in neighbors(s.A, d_ind)
                s.ΔCost_transition_cpath[nbr] += s.ΔCost_freedots_to_cpath[d_ind]
                #current_cost += ΔCost_freedots_to_cpath[d_ind]
                #println("update current cost remove disruption cost: $current_cost")
                new_weight = exp(-s.ΔCost_transition_cpath[nbr]/c)
                diff_weight = new_weight - s.weights[nbr]
                #println("neighbor of disrupted diff_weight: $diff_weight")

                inc!(s.wts_ft, nbr, diff_weight)

                s.weights[nbr] = new_weight
                s.sum_weights += diff_weight
            end
            #update the state and transition cost of the disrupted cpath
            s.state[d_ind] = 0
            s.ΔCost_transition_cpath[d_ind] += s.ΔCost_freedots_to_cpath[d_ind] - s.ΔCost_freedots_to_cpath[tcp]
            #println("increased current cost: $current_cost")
            #println("update current cost disruption: $current_cost")

            # update sum weights and fenwick tree weights
            new_weight = exp(-s.ΔCost_transition_cpath[d_ind]/c)
            diff_weight = new_weight - s.weights[d_ind]
            inc!(s.wts_ft, d_ind, diff_weight)
            s.weights[d_ind] = new_weight
            s.sum_weights += diff_weight
            #println("disrupted diff_weight: $diff_weight")

        else
            s.ΔCost_transition_cpath[d_ind] -= s.ΔCost_freedots_to_cpath[tcp]

            new_weight = exp(-s.ΔCost_transition_cpath[d_ind]/c)
            diff_weight = new_weight - s.weights[d_ind]
            inc!(s.wts_ft, d_ind, diff_weight)
            s.weights[d_ind] = new_weight
            s.sum_weights += diff_weight
            #println("disrupted diff_weight: $diff_weight")
        end
    end
end

function reset_chain!(s :: state_rejectionless, c :: Float64)
    #reset ΔCost_transition_cpath
    s.ΔCost_transition_cpath .= s.ΔCost_freedots_to_cpath
    for (cpath_ind, cpath_bool) in enumerate(s.state)
        if cpath_bool
            for d_ind in neighbors(s.A, cpath_ind)
                s.ΔCost_transition_cpath[d_ind] -= s.ΔCost_freedots_to_cpath[cpath_ind]
            end
            s.ΔCost_transition_cpath[cpath_ind] = 0
        end
    end

    s.weights = exp.(-s.ΔCost_transition_cpath./c)
    s.weights[s.state] .= 0
    s.wts_ft = FenwickTree(s.weights)
    s.sum_weights = sum(s.weights)
end

function check_lowest_cost!(s :: state_rejectionless)
    if s.current_cost < s.lowest_cost
        s.lowest_cost = s.current_cost
        s.lowest_cost_state .= s.state
    end
end

function expected_metropolis_repetitions(s :: state_rejectionless)
    #maximum([1, length(s.state)/s.sum_weights])
    exp_reps_from_weight = length(s.state)/s.sum_weights
    exp_reps_from_weight > 1 ? exp_reps_from_weight : 1
end

function simulated_annealing_metropolis(cpaths_df :: DataFrame,
                              A :: SimpleGraph,
                              sa_params)

    # number each code path by increasing cost
    n_cpaths = nrow(cpaths_df)

    if n_cpaths == 0
        return
    end

    cpath_lengths = length.(cpaths_df.cpath)
    dots = unique(vcat(cpaths_df.cpath...))
    ndots = length(dots)#nrow(pnts)

    uniform_cpath_draw = DiscreteUniform(1,n_cpaths)

    ΔCost_freedots_to_cpath = cpaths_df.cost .- cpath_lengths.*sa_params.free_dot_cost


    # initialize state vector to all zeros
    state = falses(n_cpaths)

    #current_cost = sa_params.free_dot_cost*nrow(pnts)
    current_cost = sa_params.free_dot_cost*ndots
    lowest_cost = current_cost
    lowest_cost_state = falses(n_cpaths)

    ΔCost_transition_cpath = deepcopy(ΔCost_freedots_to_cpath)


    max_iters = true
    for mc in 1:sa_params.n_chains
        # set the control parameter
        c = sa_params.cooling_factor/(1 + log(mc*sa_params.l_chain))

        println("Starting chain $mc")
        no_change_streak = 0
        n_reject = 0

        weights = exp.(-ΔCost_transition_cpath./c)
        #wts_ft = FenwickTree(weights)
        sum_weights = sum(weights)
        for i = 1:(sa_params.l_chain)

            # draw candidate codepath to transition to
            tcp = rand(uniform_cpath_draw)
            while state[tcp] == 1
                tcp = rand(uniform_cpath_draw)
            end

            # calculate disruption cost or transitioning to candidate code path
            cost_of_disruption = 0

            for d_ind in neighbors(A, tcp)
                if state[d_ind]
                    cost_of_disruption -= ΔCost_freedots_to_cpath[d_ind]
                end
            end

            if ΔCost_freedots_to_cpath[tcp] < cost_of_disruption
                metropolis_accept!(tcp, state, A)
                current_cost += ΔCost_freedots_to_cpath[tcp]

            elseif rand()*exp(-cost_of_disruption/c) < exp(-ΔCost_freedots_to_cpath[tcp]/c)
                metropolis_accept!(tcp, state, A)
                current_cost += ΔCost_freedots_to_cpath[tcp]+cost_of_disruption
            else
                n_reject += 1
            end

            #current_cost += ΔCost_freedots_to_cpath[tcp]
            #ΔCost_transition_cpath[tcp] = 0


            #new_weight = 1
            #diff_weight = 1 - weights[tcp]
            #inc!(wts_ft, tcp, diff_weight)
            #weights[tcp] = new_weight
            #sum_weights += diff_weight

            if current_cost < lowest_cost
                lowest_cost = current_cost
                lowest_cost_state = deepcopy(state)
            end

            if termination_condition(state,
                                     cpath_lengths,
                                     ndots,
                                     current_cost,
                                     sa_params,
                                     no_change_streak
                                     )#,
                max_iters = false
                #record_decode_results!(pnts, cpaths, state)
                return state
                #return cpaths[state]
            end
        end
        println("prop rejected: ", n_reject/sa_params.l_chain)
    end
    #if max_iters
    println("maximum iterations")
    #end
    state
    #record_decode_results!(pnts, cpaths, lowest_cost_state)
    #return cpaths[state]
    #println("final:")
    #println(state.lowest_energy_cws)
end

function metropolis_accept!(tcp, state, A)
    # add transition codeword to the state
    state[tcp] = true

    #remove conflicts from the state
    for d_ind in neighbors(A, tcp)
        if state[d_ind]
            state[d_ind] = 0
        end
    end
end

"""
function objective_fnc(n_free_dots, sum_cpath_vars :: Float64, n_error_dots, prms)
    energy = (prms.free_dot_energy * n_free_dots +
            prms.var_factor * sum_cpath_vars +
             prms.error_penalty * n_error_dots)
    return energy
end
"""


"""
    record_decode_results!(pnts :: DataFrame, cpaths :: Vector{Int})

Mark the final result at the termination of simulated annealing
"""
function record_decode_results!(pnts :: DataFrame, cpaths, state)
    mpaths = cpaths[state]
    #println("recording mpaths:")
    #println(mpaths)
    #println()
    decoded = zeros(length(mpaths))

    for mpath in mpaths
        if length(mpath) == 4
            mpath_vec = pnts.error_cpaths[mpath[1]]
            cand_vec = pnts.error_decode_cands[mpath[1]]
        else
            mpath_vec = pnts.perfect_cpaths[mpath[1]]
            cand_vec = pnts.perfect_decode_cands[mpath[1]]
        end

        for (i, dt_cpath) in enumerate(mpath_vec)
            if mpath_vec[i] == mpath
                pnts.decoded[mpath] .= cand_vec[i]
                decoded .= cand_vec[i]
                #println("decoded as: ", cand_vec[j])
                for dt in mpath
                    mpath_dots = Int[]
                    for dt2 in mpath
                        push!(pnts.mpath[dt], pnts.dot_ID[dt2])
                    end
                end
                break
            end
        end
    end
    return decoded, mpaths
end

#=
function termination_condition(state,
                               weights,
                               sum_weights,
                               path_lengths :: Vector{Int},
                               ndots :: Int,
                               current_cost :: Float64,
                               sa_params,
                               no_change_streak)#,
                               #cpaths)
=#
function termination_condition(s :: state_rejectionless, sa_params)

    #if #path_lengths ⋅ state == ndots
        #println("all m: ", state)
        #println("all dots matched")
        #=
        println("state: ", state)
        println("mpaths: ", cpaths[state])
        =#

        #return true
    if s.sum_weights < 0.000001
        println("converged")
        return true
    elseif s.current_cost < sa_params.free_dot_cost
        # ToDo: define minial energy parameter
        println("reached minimal cost: ")
        println(s.current_cost)
        return true

    else
        return false
    end
end

function termination_condition(s :: state_rejectionless_smtrns, sa_params)

    if s.sum_weights < 0.000001
        println("converged")
        return true
    elseif s.current_cost < sa_params.free_dot_cost
        # ToDo: define minial energy parameter
        println("reached minimal cost: ")
        println(s.current_cost)
        return true

    else
        return false
    end
end

function termination_condition(state,
                               cpath_lengths :: Vector{Int},
                               ndots :: Int,
                               current_cost :: Float64,
                               sa_params,
                               no_change_streak)#,
    if false #cpath_lengths ⋅ state == ndots
        #println("all m: ", state)
        println("all dots matched")
        #=
        println("state: ", state)
        println("mpaths: ", cpaths[state])
        =#
        return true

    elseif false#current_cost < sa_params.free_dot_cost
        # ToDo: define minial energy parameter
        println("reached minimal cost: ")
        println(current_cost)
        return true
    elseif no_change_streak > sa_params.converge_multiplier * ndots
        println("converged")
        #println("p stay in state: ", state ⋅ p)

        return true
    else
        return false
    end
end

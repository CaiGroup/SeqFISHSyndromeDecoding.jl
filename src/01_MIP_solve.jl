using JuMP
using GLPK
using Gurobi

function MIP_solve(cpath_df :: DataFrame, cpath_nbr_cpath_indices :: Vector{Vector{Int}}, optimizer)

    costs = cpath_df.cost .- 20.0 #5.0 .* length.(cpath_df.cpath) #.- maximum(cpath_df.cost)
    n = length(costs)

    model = Model(optimizer)
    if optimizer == Gurobi.Optimizer
        MOI.set(model, MOI.Silent(), true)
        MOI.set(model, MOI.NumberOfThreads(), 1)
    end
    MOI.set(model, MOI.TimeLimitSec(), 1200)
    @variable(model, x[1:n], Bin)
    for i = 1:n, j in cpath_nbr_cpath_indices[i]
        if i < j
            @constraint(model, x[i] + x[j] <= 1)
        end
    end
    @objective(model, Min, sum(x .* costs))

    optimize!(model)

    return round.(Bool, value.(x))
end


using JuMP
using GLPK

function MIP_solve(cpath_df :: DataFrame, cpath_nbr_cpath_indices :: Vector{Vector{Int}}; start_values=[false])

    costs = cpath_df.cost .- 5.0 .* length.(cpath_df.cpath)
    n = length(costs)

    model = Model(GLPK.Optimizer)
    #model = Model(Ipopt.Optimizer)
    @variable(model, x[1:n], Bin)

    if ~all(start_values .== false)
        set_start_value.(x, start_values)
    end
    
    #@constraint(model, A * x .<= b.*(b.-x))
    #@constraint(model, [i = 1:n, j = 1:n; A[i,j] == 1], x[i] + x[j] <= 1)
    for i = 1:n, j in cpath_nbr_cpath_indices[i]
        if i < j
            @constraint(model, x[i] + x[j] <= 1)
        end
    end
    #@constraint(model, x in MOI.ZeroOne())
    @objective(model, Min, sum(x .* costs))

    optimize!(model)

    return Bool.(value.(x))
end


function MIP_solve2(cpath_df :: DataFrame, cpath_nbr_cpath_indices :: Vector{Vector{Int}})

    costs = cpath_df.cost .- 5.0 .* length.(cpath_df.cpath)
    n = length(costs)

    model = Model(GLPK.Optimizer)
    #model = Model(Ipopt.Optimizer)
    @variable(model, x[1:n], Bin)
    #@constraint(model, A * x .<= b.*(b.-x))
    #@constraint(model, [i = 1:n, j = 1:n; A[i,j] == 1], x[i] + x[j] <= 1)
    for i = 1:n, j in cpath_nbr_cpath_indices[i]
        if i < j
            @constraint(model, x[i] + x[j] <= 1)
        end
    end
    #@constraint(model, x in MOI.ZeroOne())
    @objective(model, Min, sum(x .* costs))

    optimize!(model)

    return Bool.(value.(x))
end

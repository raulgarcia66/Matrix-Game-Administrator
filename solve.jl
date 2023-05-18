using Gurobi
using JuMP

function solve_game(A::Matrix{T}, c::Vector{U}, B::V; relax::Bool=false, TimeLimit::Int=300, MIPGap::W=0.05) where {T,U,V,W <: Real}
    
    num_rows = size(A, 1)

    #### Initialize model
    # model = Model(Gurobi.Optimizer)
    model = Model(optimizer_with_attributes(Gurobi.Optimizer, "MIPGap"=>MIPGap, "TimeLimit"=>TimeLimit))

    @variable(model, x[1:num_rows] >= 0)   # Don't need upper bound since sum(x) == 1 will enforce it
    @variable(model, r[1:num_rows], Bin)   # 1 means play is available
    @variable(model , z)

    @objective(model, Max, z)

    @constraint(model, [j = 1:size(A,2)], z - A[:,j]' * x <= 0 )
    @constraint(model, sum(c' * r) <= B)
    @constraint(model, sum(x) == 1)
    @constraint(model, x .<= r)

    #### Solve
    if relax
        relax_integrality(model)
    end
    optimize!(model)

    if termination_status(model) == INFEASIBLE_OR_UNBOUNDED || termination_status(model) == INFEASIBLE
        error("Model is infeasible.")
    elseif termination_status(model) == TIME_LIMIT && !has_values(model)
        error("No primal solution obtained within the time limit.")
    end

    return value.(x), value.(r), objective_value(model), termination_status(model), solve_time(model), relative_gap(model), node_count(model)
end


function solve_game(A::Matrix{T}, c::Vector{U}, B::V, r_fix::Vector{W}; relax::Bool=false, TimeLimit::Int=300, MIPGap::X=0.05) where {T,U,V,W,X <: Real}
    
    num_rows = size(A, 1)

    #### Initialize model
    # model = Model(Gurobi.Optimizer)
    model = Model(optimizer_with_attributes(Gurobi.Optimizer, "MIPGap"=>MIPGap, "TimeLimit"=>TimeLimit))

    @variable(model, x[1:num_rows] >= 0)   # Don't need upper bound since sum(x) == 1 will enforce it
    @variable(model, r[1:num_rows], Bin)   # 1 means play is available
    @variable(model , z)

    @objective(model, Max, z)

    @constraint(model, [j = 1:size(A,2)], z - A[:,j]' * x <= 0 )
    @constraint(model, sum(c' * r) <= B)
    @constraint(model, sum(x) == 1)
    @constraint(model, x .<= r)

    #### Fix values of r (1 or 0 will be fixed, all other values in r_fix are ignored)
    if r_fix != -ones(Int, length(c))
        for i in eachindex(r_fix)
            if r_fix[i] == 1
                @constraint(model, r[i] == 1)
            elseif r_fix[i] == 0
                @constraint(model, r[i] == 0)
            end
        end
    end

    #### Solve
    if relax
        relax_integrality(model)
    end
    optimize!(model)

    if termination_status(model) == INFEASIBLE_OR_UNBOUNDED || termination_status(model) == INFEASIBLE
        error("Model is infeasible.")
    elseif termination_status(model) == TIME_LIMIT && !has_values(model)
        error("No primal solution obtained within the time limit.")
    end

    return value.(x), value.(r), objective_value(model), termination_status(model), solve_time(model), relative_gap(model), node_count(model)
end
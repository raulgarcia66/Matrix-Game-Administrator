using Gurobi
using JuMP

function solve_game(A, c, B; relax=false, r_fix::Vector{T}=-ones(Int,length(c))) where {T <: Real}
    
    num_rows = size(A, 1)

    #### Initialize model
    model = Model(Gurobi.Optimizer)

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

    return value.(x), value.(r), objective_value(model), termination_status(model), solve_time(model), relative_gap(model), node_count(model) # model
end


function solve_game(A, c, B, r_fix::Vector{T}; relax=false) where {T <: Real}
    
    num_rows = size(A, 1)

    #### Initialize model
    model = Model(Gurobi.Optimizer)

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

    return value.(x), value.(r), objective_value(model), termination_status(model), solve_time(model), relative_gap(model), node_count(model) # model
end
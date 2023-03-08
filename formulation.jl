using Gurobi
using JuMP

function solve_game(A, c, B; relax=false)
    
    num_row_plays = size(A, 1)

    #### Initialize model
    model = Model(Gurobi.Optimizer)

    @variable(model, x[1:num_row_plays] >= 0)
    @variable(model, r[1:num_row_plays], Bin)   # 1 means play is available
    @variable(model , z)

    @objective(model, Max, z)

    @constraint(model, [j = 1:size(A,2)], z - A[:,j]' * x <= 0 )
    @constraint(model, sum(c' * r) <= B)
    @constraint(model, sum(x) == 1)
    @constraint(model, x .<= r)

    #### Regulatory constraints
    # E.g., Bound on the number of rows that can be purchased


    if relax
        relax_integrality(model)
    end
    optimize!(model)

    return model, value.(x), value.(r), value.(z)
end

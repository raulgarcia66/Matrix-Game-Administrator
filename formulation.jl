
using Gurobi
using JuMP

function compute_big_M(A)
    # Want z - dot(A[:,j], x) <= M[j]
    # z - dot(A[:,j], x) <= Largest {z} - Smallest {dot(A[:,j], x)}

    # Largest value z can be is max_ij A_[i,j]. This needs to be over all (i,j), correct?
    u = maximum( vec(A) )
    # Smallest value dot(A[:,j], x) can be is min_i A[i,j]. Can do this for each j.
    ℓ = map(j -> minimum(A[:,j]), 1:size(A,2))

    return u .- ℓ
end

function solve_game(A, Γ, cost_r, cost_s; relax=false)

    # Create some parameters
    num_row_plays = length(cost_r)
    num_col_plays = length(cost_s)
    
    #### Initialize model

    model = Model(Gurobi.Optimizer)

    @variable(model, x[1:num_row_plays] >= 0)
    @variable(model, r[1:num_row_plays], Bin)   # 1 means play is available
    @variable(model, s[1:num_col_plays], Bin)   # 0 means play is eliminated
    @variable(model , z)

    @objective(model, Max, Γ * z - cost_r' * (1 .- r) - cost_s' * (1 .- s) )

    #### Upper level play selection
    # E.g., Bound on the number of rows and/or columns that can be eliminated
    # E.g., Cannot eliminate all plays of column player

    #### Lower level strategy
    M = compute_big_M(A)
    @constraint(model, [j = 1:length(M)], z - A[:,j]' * x <= M[j] * (1 - s[j]) )
    @constraint(model, sum(x) == 1)
    @constraint(model, x .<= r)

    if relax
        relax_integrality(model)
    end
    optimize!(model)

    return model, value.(x), value.(r), value.(s), value.(z)
end

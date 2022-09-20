
using LinearAlgebra

# Parameters
Γ = 0.2
cost_r = fill(0.1, 5)
cost_r = rand(1:5, 5) ./ 10
cost_s = fill(0.1, 5)
cost_s = rand(1:5, 5) ./ 10
A = rand(1:9, 5, 5)
for i = 1:size(A,1)
    for j = i+1:size(A,2)
        A[i,j] = -A[j,i]
    end
end

# Solve
model, x, r, s, z = solve_game(A, Γ, cost_r, cost_s, relax=true)
model, x, r, s, z = solve_game(A, Γ, cost_r, cost_s)

termination_status(model)
objective_value(model)
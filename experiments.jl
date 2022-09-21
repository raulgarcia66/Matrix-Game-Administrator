
using LinearAlgebra
using Dates
include("formulation.jl")

# Parameters
Γ = 0.2
cost_r = fill(0.1, 5)
cost_r = rand(1:5, 5) ./ 10
cost_s = fill(0.1, 5)
cost_s = rand(1:5, 5) ./ 10
A = rand(-10:10, 5, 5)
for i = 1:size(A,1)
    for j = i+1:size(A,2)
        A[i,j] = -A[j,i]
    end
end

# Solve
model, x, r, s, z = solve_game(A, Γ, cost_r, cost_s)
# model, x, r, s, z = solve_game(A, Γ, cost_r, cost_s, relax=true)

term_status = termination_status(model)
obj_val = objective_value(model)
solution_time = solve_time(model)


#### Run experiments
# Make sure the A looks balanced
# A = rand(-10:10, 5, 5)
# for i = 1:size(A,1)
#     for j = i+1:size(A,2)
#         A[i,j] = -A[j,i]
#     end
# end
# A_ref = copy(A)

A_ref = [0 -5 -2 -4 -10; 5 -6 4 7 2; 2 -4 5 0 -8; 4 -7 0 -9 7; 10 -2 8 -7 0]
# for i = 1:size(A,1)
#     A_ref[i,i] = 0
# end

# Scale matrix by a range of values to find balance with other parameters
scales = [i for i = 1:25]

file_name = "./Experiments/Set $(Dates.today()) $(Dates.format(now(), "HH MM SS")).txt"
# file_name = "./Experiments/Relaxed Set $(Dates.today()) $(Dates.format(now(), "HH MM SS")).txt"
f = open(file_name, "w")
write(f, "Run\tGamma\tCost_r\tCost_s\tA\tTermination_status\tSolve_time\tObjective_val\tz\tx\tr\ts\n")

for i = 1:length(scales)
    Γ = 0.2
    # cost_r = fill(0.1, 5)
    # cost_r = rand(1:5, 5) ./ 10
    # Integral
    cost_r = fill(1, 5) * 5
    # cost_r = rand(1:5, 5)

    # cost_s = fill(0.1, 5)
    # cost_s = rand(1:5, 5) ./ 10
    # Integral
    cost_s = fill(1, 5) * 5
    # cost_s = rand(1:5, 5)

    A = copy(A_ref) * scales[i]

    model, x, r, s, z = solve_game(A, Γ, cost_r, cost_s)
    # model, x, r, s, z = solve_game(A, Γ, cost_r, cost_s, relax=true)
    term_status = termination_status(model)
    obj_val = objective_value(model)
    solution_time = solve_time(model)

    # Write to files
    # f = open(file_name, "a")
    write(f, "$i\t$Γ\t$cost_r\t$cost_s\t$A\t$term_status\t$solution_time\t$obj_val\t$z\t$x\t$r\t$s\n")
    flush(f)
end
close(f)


# Interesting case. s has two 1's and x his a singleton.
A_ref_12 = A_ref * 12

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

term_status = termination_status(model)
obj_val = objective_value(model)
solution_time = solve_time(model)

###################################################################################
###################################################################################

function run_experiments(A_ref, Γ, cost_r, cost_s;
        file_name="./Run of Experiments", scales=[1], mode="w", relax=false)

    f = open(file_name, mode)
    write(f, "Scale\tGamma\tCost_r\tCost_s\tA\tTermination_status\tSolve_time\tM\tObjective_val\tz\tx\tr\ts\n")

    for scale in scales
        A = copy(A_ref) * scale
    
        model, x, r, s, z, M = solve_game(A, Γ, cost_r, cost_s, relax=relax)
        term_status = termination_status(model)
        obj_val = objective_value(model)
        solution_time = solve_time(model)
    
        write(f, "$scale\t$Γ\t$cost_r\t$cost_s\t$A\t$term_status\t$solution_time\t$M\t$obj_val\t$z\t$x\t$r\t$s\n")
        flush(f)
    end
    close(f)
end

function create_matrix(scores, num_rows, num_cols; zero_diag::Bool=true)
    A = rand(scores, num_rows, num_cols)

    for i = 1:num_rows
        # Make symmetric
        for j = i+1:num_cols
            A[i,j] = -A[j,i]
        end

        if zero_diag == true && num_rows ==  num_cols
            A[i,i] = 0
        end
    end

    return A
end


###################################################################################
################################ Interesting Cases ################################
#### Sep 20
A = [0 -5 -2 -4 -10; 5 -6 4 7 2; 2 -4 5 0 -8; 4 -7 0 -9 7; 10 -2 8 -7 0]
# Interesting case. s has two 1's and x his a singleton.
A_12 = A * 12

#### Sep 27
A = [4 -9 9 -5 -3; 9 -3 3 0 2; -9 -3 3 -6 2; 5 0 6 6 8; 3 -2 -2 -8 -5]
# With M = 1e6, get round off errors
# M = 1000 and my indiviual M values yielded the same results
# LP solutions first decrease at the first s[j] value that becomes 0 in the MIP

A = [0 -4 -9 0 6; 4 0 7 -3 4; 9 -7 0 4 -8; 0 3 -4 0 10; -6 -4 8 -10 0]
# Diagonal is 0s.
# OBSERVATION: A fair game occurs when the matrix is negative symmetric (diagonals do not have to be 0)
# and all plays are available to the players.
A = [0 0 0 3 -1; 0 0 -6 -3 0; 0 6 0 8 -9; -3 3 -8 0 -6; 1 0 9 6 0]
A = [0 1 -1 1 -1; -1 0 1 -1 1; 1 -1 0 1 -1; -1 1 -1 0 1; 1 -1 1 -1 0]

#### Oct 7
# Must play equal number of rows and columns
A = [4 -9 9 -5 -3; 9 -3 3 0 2; -9 -3 3 -6 2; 5 0 6 6 8; 3 -2 -2 -8 -5]
A = [0 -5 -2 -4 -10; 5 -6 4 7 2; 2 -4 5 0 -8; 4 -7 0 -9 7; 10 -2 8 -7 0]
A = [0 0 0 3 -1; 0 0 -6 -3 0; 0 6 0 8 -9; -3 3 -8 0 -6; 1 0 9 6 0]
A = [0 1 -1 1 -1; -1 0 1 -1 1; 1 -1 0 1 -1; -1 1 -1 0 1; 1 -1 1 -1 0]
# Interesting. Optimal solution x is always play all plays an equal amount (ie, x = ones(m) / m ).
# Value of the game is always 0.

# Try different cost vectors
# Try constrains on plays. Rows eliminated >= or == or <= columns eliminated. (Seemingly fair)
# Look at nonsymmetric matrices and nonsquare matrices

###################################################################################
###################################################################################
# Scale matrix by a range of values to find balance with other parameters
scales = [i for i = 1:50]

Γ = 0.2

# cost_r = fill(0.1, 5)
cost_r = fill(1, 5) * 5
# cost_r = rand(1:5, 5)
# cost_r = rand(1:5, 5) ./ 10

# cost_s = fill(0.1, 5)
cost_s = fill(1, 5) * 5
# cost_s = rand(1:5, 5)
# cost_s = rand(1:5, 5) ./ 10

file_name = "./Experiments/Set $(Dates.today()) $(Dates.format(now(), "HH MM SS")).txt"
# file_name = "./Experiments/Relaxed Set $(Dates.today()) $(Dates.format(now(), "HH MM SS")).txt"

# A = create_matrix(-10:10, 5, 5, zero_diag=false)
# println("$A")

run_experiments(A, Γ, cost_r, cost_s, file_name=file_name, scales=scales, relax=false)

###################################################################################
#################################### Load Data ####################################

using DataFrames

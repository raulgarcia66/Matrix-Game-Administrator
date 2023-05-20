include("parameters.jl")
include("solve.jl")

##### Parameters
A = create_matrix(-10:10, 6, 7)
A = create_matrix(-2:10, 100)
A = create_matrix_symmetric(-10:10, 90)

num_rows = size(A,1)

# c = fill(2, num_rows)
c = rand(2:5, num_rows)

B = sum(c) / 5
B = 10

#### Solve
x, r, obj_val, term_status, soln_time, rel_gap, nodes = solve_game(A, c, B, TimeLimit=100)
B_used = r' * c
opt_val = copy(obj_val)

####
r_fix = ones(num_rows)   # values of 1 and 0 are fixed, all else ignored
x, r, obj_val, term_status, soln_time, rel_gap, nodes = solve_game(A, c, B, r_fix)
B_used = r' * c

####
x, r, obj_val, soln_time, soln_attempts = solve_game_naive(A, c, B, TimeLimit=10)
B_used = r' * c

####
x, r, obj_val, term_status, soln_time, gap, soln_attempts = solve_game_naive_test(A, c, B, opt_val, TimeLimit=10, seed = 2)
B_used = r' * c

# Use dual_objective_bound(model) to get dual bound. Can use this to compare with the lower bound of the naive method

###################################################################################
############################### Logging Experiments ###############################

num_rows_vec = [1000] # [100], [1000] [10000]
num_cols_vec = [10,100,1000] # [10000]
total_experiments_per_matrix_size = 1
budget_denominators = [1, 4/3, 2, 3, 4, 10]
matrix_entry_range = -100:100
costs_entry_range = 1:10

set_num = 3
mkpath("./Experiments/Set $set_num")
subpath = "./Experiments/Set $set_num/"

filenames = String[]
push!(filenames, "Set $set_num summary.txt")
push!(filenames, "Sets summary.txt")
for file in filenames
    if file == "Sets summary.txt"
        filename = "./Experiments/" * file
        f = open(filename, "a")
    else
        filename = subpath * file
        f = open(filename, "w")
    end
    write(f, "Set $set_num\n")
    write(f, "num_rows_vec = $num_rows_vec\n")
    write(f, "num_cols_vec = $num_cols_vec\n")
    write(f, "total_experiments_per_matrix_size = $total_experiments_per_matrix_size\n")
    write(f, "budget_fraction = $(1 ./ budget_denominators)\n")
    write(f, "matrix entry range = $matrix_entry_range\n")
    write(f, "costs_entry_range = $costs_entry_range\n")
    write(f, "\n")
    close(f)
end

for num_rows in num_rows_vec, num_cols in num_cols_vec
    A_vec = map(seed -> create_matrix(matrix_entry_range, num_rows, num_cols), Base.OneTo(total_experiments_per_matrix_size))
    c_vec = map(seed -> rand(costs_entry_range, num_rows), Base.OneTo(total_experiments_per_matrix_size))

    filename = subpath * "Matrices $num_rows by $num_cols.txt"
    f = open(filename, "a")
    write(f, "num_rows = $num_rows\n")
    write(f, "num_cols = $num_cols\n")
    write(f, "Matrix seed\tCosts seed\tB\tBudget fraction\tBudget spent\tObj val\tTermination status\tSolve time\tRelative gap\tNode count\n")
    for i in Base.OneTo(total_experiments_per_matrix_size)
    # for (i,A) in enumerate(A_vec)
    #     for (j,c) in enumerate(c_vec)
            
            B_vec = [sum(c_vec[i]) / d for d = budget_denominators]
            for (k,B) in enumerate(B_vec)
                x, r, obj_val, term_status, soln_time, rel_gap, nodes = solve_game(A_vec[i], c_vec[i], B)
                B_used = r' * c_vec[i]

                write(f, "$i\t$i\t$B\t$(1 / budget_denominators[k])\t$B_used\t$obj_val\t$term_status\t$soln_time\t$rel_gap\t$nodes\n")
                flush(f)
            end

        # end
    end
    close(f)
end


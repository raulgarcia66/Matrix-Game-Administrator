work_dir = pwd()
include(work_dir * "/src/parameters.jl")
include(joinpath(pwd(), "src/solve.jl"))

# TODO: Need to correct directories in experiments section

##### Parameters
# A = create_matrix(-10:10, 6, 7, seed=1)
# A = create_matrix(-2:10, 7, seed=5);
# A = create_matrix_symmetric(-10:10, 90);

entry_range = [i/20 for i = -7:10]
seed = 2
A = create_matrix(entry_range, 6, seed=seed)

num_rows = size(A,1)

# c = fill(2, num_rows)
c = create_cost_vector(2:5, num_rows, seed=seed);

B = sum(c) * 0.5
# B = sum(c)

########################################################################
x, r, obj_val, obj_bound, dual_obj, term_status, soln_time, rel_gap, nodes = solve_game(A, c, B, TimeLimit=30)
B_used = r' * c
R = filter(i -> r[i] == 1, eachindex(r))

opt_val, x_opt, r_opt, R_top, B_used_opt = copy(obj_val), copy(x), copy(r), copy(R), B_used

# Containment
c_patho = copy(c)
# c_patho = [4, 2, 2, 4, 3, 3, 5]

B_1 = 6
x, r, obj_val, obj_bound, dual_obj, term_status, soln_time, rel_gap, nodes = solve_game(A, c_patho, B_1, TimeLimit=30)
val_1 = copy(obj_val)
B_used = r' * c_patho
R_1 = filter(i -> r[i] == 1, eachindex(r))

B_2 = 11
x, r, obj_val, obj_bound, dual_obj, term_status, soln_time, rel_gap, nodes = solve_game(A, c_patho, B_2, TimeLimit=30)
val_2 = copy(obj_val)
B_used = r' * c_patho
R_2 = filter(i -> r[i] == 1, eachindex(r))

########################################################################
r_fix = ones(num_rows)   # values of 1 and 0 are fixed, all else ignored
x, r, obj_val, term_status, soln_time, rel_gap, nodes = solve_game(A, c, B, r_fix)
B_used = r' * c

# Submodularity testing
B = sum(c)
# R = {2,3}
r_fix = [0, 1, 1, 0, 0, 0]
x, r, obj_val, term_status, soln_time, rel_gap, nodes = solve_game(A, c, B, r_fix)
val_R = copy(obj_val)
# R = {1,2,3} (k=1)
r_fix = [1, 1, 1, 0, 0, 0]
x, r, obj_val, term_status, soln_time, rel_gap, nodes = solve_game(A, c, B, r_fix)
val_R_k = copy(obj_val)

val_R_k - val_R

# S = {2,3,5}
r_fix = [0, 1, 1, 1, 0, 0]
x, r, obj_val, term_status, soln_time, rel_gap, nodes = solve_game(A, c, B, r_fix)
val_S = copy(obj_val)
# S = {1,2,3,5} (k=1)
r_fix = [1, 1, 1, 1, 0, 0]
x, r, obj_val, term_status, soln_time, rel_gap, nodes = solve_game(A, c, B, r_fix)
val_S_k = copy(obj_val)

val_S_k - val_S

########################################################################
# x, r, obj_val, soln_time, soln_attempts = solve_game_naive(A, c, B, TimeLimit=10)
# B_used = r' * c

########################################################################
# x, r, obj_val, term_status, soln_time, gap, soln_attempts = solve_game_naive_compare(A, c, B, opt_val, TimeLimit=10, seed = 2)
# B_used = r' * c

########################################################################
# Greedy
B = 4
c = [1,1,1,1,1,1]
x, r, obj_val, term_status, time_elapsed, num_purchases, nodes, R, obj_val_vec, B_used_vec = solve_game_greedy(A, c, B, TimeLimit=30)
gains = compute_gains(obj_val_vec)

x
r
B_used = r' * c
B
obj_val
term_status
time_elapsed
num_purchases
nodes
R
obj_val_vec
B_used_vec

########################
# Greedy LP
x_lp, r_lp, obj_val_lp, term_status_lp, time_elapsed_lp, num_purchases_lp, r_vec_lp, obj_val_vec_lp, B_used_vec_lp = solve_game_greedy_LP(A, c, B, TimeLimit=30)
x_lp
r_lp
B_used_lp = r' * c
B
obj_val_lp
term_status_lp
time_elapsed_lp
num_purchases_lp
r_vec_lp
obj_val_vec_lp
B_used_vec_lp


###################################################################################
####################### Comparison with Greedy Approaches #######################

# On: [1000] × [100,1000];
# Done: [10] × [10,100,1000]; [100] × [10,100,1000]; [1000] × [10]
num_rows_vec = [1000]
num_cols_vec = [100, 1000]
MIPGap = 1E-2
TimeLimit = 1200
total_experiments_per_matrix_size = 5
budget_denominators = [1, 4/3, 2, 3, 4, 10]
matrix_entry_range = -10:100
costs_entry_range = 1:10

set_num = 5
subpath = "./Experiments/Row Variant/Set $set_num/"
# mkpath(subpath)

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
    write(f, "TimeLimit = $TimeLimit\n")
    write(f, "MIPGap = $MIPGap\n")
    write(f, "total_experiments_per_matrix_size = $total_experiments_per_matrix_size\n")
    write(f, "budget_fraction = $(1 ./ budget_denominators)\n")
    write(f, "matrix entry range = $matrix_entry_range\n")
    write(f, "costs_entry_range = $costs_entry_range\n")
    write(f, "\n")
    close(f)
end

test_MIP = false; test_naive = false; test_greedy = false; test_greedy_LP = true;
for num_rows in num_rows_vec, num_cols in num_cols_vec
    println("Num rows = $num_rows, Num cols = $num_cols")
    A_vec = map(seed -> create_matrix(matrix_entry_range, num_rows, num_cols, seed=seed), Base.OneTo(total_experiments_per_matrix_size))
    c_vec = map(seed -> create_cost_vector(costs_entry_range, num_rows, seed=seed), Base.OneTo(total_experiments_per_matrix_size))

    #### MIP
    if test_MIP
        filename = subpath * "Matrices $num_rows by $num_cols MIP.txt"
        f = open(filename, "a")
        write(f, "num_rows = $num_rows\n")
        write(f, "num_cols = $num_cols\n")
        write(f, "TimeLimit = $TimeLimit\n")
        write(f, "MIPGap = $MIPGap\n")
        write(f, "Matrix seed\tCosts seed\tB\tBudget fraction\tBudget spent\tObj val\tObj bound\tDual obj\tTermination status\tSolve time\tRelative gap\tNode count\n")
        # for (i,A) in enumerate(A_vec)
        #     for (j,c) in enumerate(c_vec)
        for i in Base.OneTo(total_experiments_per_matrix_size)
                
            B_vec = [sum(c_vec[i]) / d for d = budget_denominators]
            for (k,B) in enumerate(B_vec)
                x, r, obj_val, obj_bound, dual_obj, term_status, soln_time, rel_gap, nodes = solve_game(A_vec[i], c_vec[i], B, MIPGap=MIPGap, TimeLimit=TimeLimit)
                B_used = r' * c_vec[i]

                write(f, "$i\t$i\t$B\t$(1 / budget_denominators[k])\t$B_used\t$obj_val\t$obj_bound\t$dual_obj\t$term_status\t$soln_time\t$rel_gap\t$nodes\n")
                flush(f)
            end

            # end
        end
        close(f)
    end

    #### Naive approach
    if test_naive
        filename = subpath * "Matrices $num_rows by $num_cols naive.txt"
        f = open(filename, "a")
        write(f, "num_rows = $num_rows\n")
        write(f, "num_cols = $num_cols\n")
        write(f, "TimeLimit = $TimeLimit\n")
        write(f, "MIPGap = NA\n")
        write(f, "Matrix seed\tCosts seed\tB\tBudget fraction\tBudget spent\tBest obj val\tTime achieved\tSolve time\tSolution attempts\n")
        flush(f)
        for i in Base.OneTo(total_experiments_per_matrix_size)
            # println("Experiment $i of $total_experiments_per_matrix_size")
                
            B_vec = [sum(c_vec[i]) / d for d = budget_denominators]
            for (k,B) in enumerate(B_vec)
                # println("Budget fraction: $(1 / budget_denominators[k])")
                x, r, obj_val, time_achieved, time_elapsed, soln_attempts = solve_game_naive(A_vec[i], c_vec[i], B)
                B_used = r' * c_vec[i]

                write(f, "$i\t$i\t$B\t$(1 / budget_denominators[k])\t$B_used\t$obj_val\t$time_achieved\t$time_elapsed\t$soln_attempts\n")
                flush(f)
            end
        end
        close(f)
    end

    #### Greedy Algorithm
    if test_greedy
        filename = subpath * "Matrices $num_rows by $num_cols greedy.txt"
        f = open(filename, "a")
        write(f, "num_rows = $num_rows\n")
        write(f, "num_cols = $num_cols\n")
        write(f, "TimeLimit = $TimeLimit\n")
        write(f, "MIPGap = $MIPGap\n")
        write(f, "Matrix seed\tCosts seed\tB\tBudget fraction\tBudget spent\tObj val\tTermination status\tSolve time\tNum purchases\tNodes\n")
        flush(f)
        for i in Base.OneTo(total_experiments_per_matrix_size)
            # println("Experiment $i of $total_experiments_per_matrix_size")
                
            B_vec = [sum(c_vec[i]) / d for d = budget_denominators]
            for (k,B) in enumerate(B_vec)
                # println("Budget fraction: $(1 / budget_denominators[k])")
                x, r, obj_val, term_status, time_elapsed, num_purchases, nodes, _, _, _ = solve_game_greedy(A_vec[i], c_vec[i], B, MIPGap=MIPGap, TimeLimit=TimeLimit)
                B_used = r' * c_vec[i]

                write(f, "$i\t$i\t$B\t$(1 / budget_denominators[k])\t$B_used\t$obj_val\t$term_status\t$time_elapsed\t$num_purchases\t$nodes\n")
                flush(f)
            end
        end
        close(f)
    end

    #### Greedy Algorithm LP
    if test_greedy_LP
        filename = subpath * "Matrices $num_rows by $num_cols greedy LP.txt"
        f = open(filename, "a")
        write(f, "num_rows = $num_rows\n")
        write(f, "num_cols = $num_cols\n")
        write(f, "TimeLimit = $TimeLimit\n")
        write(f, "MIPGap = $MIPGap\n")
        write(f, "Matrix seed\tCosts seed\tB\tBudget fraction\tBudget spent\tObj val\tTermination status\tSolve time\tNum purchases\n")
        flush(f)
        for i in Base.OneTo(total_experiments_per_matrix_size)
            # println("Experiment $i of $total_experiments_per_matrix_size")
                
            B_vec = [sum(c_vec[i]) / d for d = budget_denominators]
            for (k,B) in enumerate(B_vec)
                # println("Budget fraction: $(1 / budget_denominators[k])")
                x, r, obj_val, term_status, time_elapsed, num_purchases, _, _, _ = solve_game_greedy_LP(A_vec[i], c_vec[i], B, MIPGap=MIPGap, TimeLimit=TimeLimit)
                B_used = r' * c_vec[i]

                write(f, "$i\t$i\t$B\t$(1 / budget_denominators[k])\t$B_used\t$obj_val\t$term_status\t$time_elapsed\t$num_purchases\n")
                flush(f)
            end
        end
        close(f)
    end
end


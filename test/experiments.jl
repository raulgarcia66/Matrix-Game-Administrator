include(".././src/parameters.jl")
include(".././src/solve.jl")

##### Parameters
A = create_matrix(-10:10, 6, 7)
A = create_matrix(-2:10, 100)
A = create_matrix_symmetric(-10:10, 90)

num_rows = size(A,1)

# c = fill(2, num_rows)
c = rand(2:5, num_rows)

B = sum(c) / 5
B = sum(c)

#### Solve
x, r, obj_val, obj_bound, dual_obj, term_status, soln_time, rel_gap, nodes = solve_game(A, c, B, TimeLimit=300)
B_used = r' * c
opt_val = copy(obj_val)

x_opt = copy(x)
r_opt = copy(r)
B_used_opt = B_used

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

####
x, r, obj_val, term_status, time_elapsed, num_purchases, r_vec, obj_val_vec, B_used_vec = solve_game_greedy_LP(A, c, B)
x
r
B_used = r' * c
B
obj_val
term_status
time_elapsed
num_purchases
r_vec
obj_val_vec
B_used_vec

x_gn, r_gn, obj_val_gn, term_status_gn, time_elapsed_gn, num_purchases_gn, nodes_gn, r_vec_gn, obj_val_vec_gn, B_used_vec_gn = solve_game_greedy_naive(A, c, B, TimeLimit=100)
x_gn
r_gn
B_used_gn = r' * c
B
obj_val_gn
term_status_gn
time_elapsed_gn
num_purchases_gn
r_vec_gn
obj_val_vec_gn
B_used_vec_gn


#### Submodularity testing
A = create_matrix(-10:10, 6, 7, seed=1)
num_rows = size(A,1)
c = rand(2:5, num_rows)  # c = [5 3 3 2 2 5]
B = sum(c)
r_fix = zeros(Int,num_rows)

r_vec = []
x_vec = []
obj_val_vec = []
B_used_vec = []
desc_vec = []

r_fix[1] = 1; r_fix[2] = 1; r_fix[3] = 0; r_fix[4] = 0; r_fix[6] = 0
r_fix[5] = 0; r_fix
x, r, obj_val, term_status, soln_time, rel_gap, nodes = solve_game(A,c,B,r_fix)
B_used = r' * c
push!(r_vec, r)
push!(x_vec, x)
push!(obj_val_vec, obj_val)
push!(B_used_vec, B_used)
# push!(desc_vec, "S ∪ {4}")

t = 4
r_vec[t] = r
obj_val_vec[t] = obj_val
B_used_vec[t] = B_used

desc_vec[3] = "R ∪ {5}"
desc_vec[4] = "S ∪ {5}"

obj_val_vec[2] - obj_val_vec[1]
obj_val_vec[4] - obj_val_vec[3]

#### Subset inclusivity testing
A = create_matrix(-10:10, 6, 7, seed=1)
num_rows = size(A,1)
c = rand(2:5, num_rows)  # c = [5 3 3 2 2 5]
B = sum(c) / 5  # 4
B_2 = sum(c) / 4  # 5

x, r, obj_val, term_status, soln_time, rel_gap, nodes, = solve_game(A, c, B)
B_used = r' * c

x_2, r_2, obj_val_2, term_status_2, soln_time_2, rel_gap_2, nodes_2 = solve_game(A, c, B_2)
B_used_2 = r_2' * c

##
B = sum(c)

r_fix = zeros(Int,num_rows)
r_fix[1] = 1; r_fix[2] = 0; r_fix[3] = 1; r_fix[4] = 0; r_fix[6] = 0
r_fix[5] = 0; r_fix
x, r, obj_val, term_status, soln_time, rel_gap, nodes = solve_game(A, c, B, r_fix)
B_used = r' * c

r_fix_2 = zeros(Int,num_rows)
r_fix_2[1] = 1; r_fix_2[2] = 0; r_fix_2[3] = 1; r_fix_2[4] = 0; r_fix_2[6] = 0
r_fix_2[5] = 1; r_fix_2
x_2, r_2, obj_val_2, term_status_2, soln_time_2, rel_gap_2, nodes_2 = solve_game(A, c, B, r_fix_2)
B_used_2 = r_2' * c

obj_val_2 - obj_val

r_fix_3 = zeros(Int,num_rows)
r_fix_3[1] = 0; r_fix_3[2] = 1; r_fix_3[3] = 1; r_fix_3[4] = 0; r_fix_3[6] = 0
r_fix_3[5] = 0; r_fix_3
x_3, r_3, obj_val_3, term_status_3, soln_time_3, rel_gap_3, nodes_3 = solve_game(A, c, B, r_fix_3)
B_used_3 = r_3' * c

r_fix_4 = zeros(Int,num_rows)
r_fix_4[1] = 0; r_fix_4[2] = 1; r_fix_4[3] = 1; r_fix_4[4] = 0; r_fix_4[6] = 0
r_fix_4[5] = 1; r_fix_4
x_4, r_4, obj_val_4, term_status_4, soln_time_4, rel_gap_4, nodes_4 = solve_game(A, c, B, r_fix_4)
B_used_4 = r_4' * c

obj_val_4 - obj_val_3


###################################################################################
############################### General Experiments ###############################

num_rows_vec = [10,100,1000] # [100], [1000] [10000]
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
    # TODO: If gonna use these, need to rerun them since initial batch didn't incorporate the seed
    A_vec = map(seed -> create_matrix(matrix_entry_range, num_rows, num_cols, seed=seed), Base.OneTo(total_experiments_per_matrix_size))
    c_vec = map(seed -> create_cost_vector(costs_entry_range, num_rows, seed=seed), Base.OneTo(total_experiments_per_matrix_size))

    filename = subpath * "Matrices $num_rows by $num_cols.txt"
    f = open(filename, "a")
    write(f, "num_rows = $num_rows\n")
    write(f, "num_cols = $num_cols\n")
    # write(f, "Matrix seed\tCosts seed\tB\tBudget fraction\tBudget spent\tObj val\tTermination status\tSolve time\tRelative gap\tNode count\n")
    write(f, "Matrix seed\tCosts seed\tB\tBudget fraction\tBudget spent\tObj val\tObj bound\tDual obj\tTermination status\tSolve time\tRelative gap\tNode count\n")
    for i in Base.OneTo(total_experiments_per_matrix_size)
    # for (i,A) in enumerate(A_vec)
    #     for (j,c) in enumerate(c_vec)
            
            B_vec = [sum(c_vec[i]) / d for d = budget_denominators]
            for (k,B) in enumerate(B_vec)
                # x, r, obj_val, term_status, soln_time, rel_gap, nodes = solve_game(A_vec[i], c_vec[i], B)
                # B_used = r' * c_vec[i]

                # write(f, "$i\t$i\t$B\t$(1 / budget_denominators[k])\t$B_used\t$obj_val\t$term_status\t$soln_time\t$rel_gap\t$nodes\n")
                # flush(f)

                x, r, obj_val, obj_bound, dual_obj, term_status, soln_time, rel_gap, nodes = solve_game(A_vec[i], c_vec[i], B)
                B_used = r' * c_vec[i]

                write(f, "$i\t$i\t$B\t$(1 / budget_denominators[k])\t$B_used\t$obj_val\t$obj_bound\t$dual_obj\t$term_status\t$soln_time\t$rel_gap\t$nodes\n")
                flush(f)
            end

        # end
    end
    close(f)
end


###################################################################################
################### Comparison with Naive and Greedy Approaches ###################

num_rows_vec = [10,100,1000] # [10, 100, 10000]
num_cols_vec = [10,100,1000]
MIPGap = 1E-2
TimeLimit = 300
total_experiments_per_matrix_size = 1
budget_denominators = [1, 4/3, 2, 3, 4, 10]
matrix_entry_range = -10:100
costs_entry_range = 1:10

set_num = 3
subpath = "./Experiments MIP vs Naive vs Greedy/Set $set_num/"
mkpath(subpath)

filenames = String[]
push!(filenames, "Set $set_num summary.txt")
push!(filenames, "Sets summary.txt")
for file in filenames
    if file == "Sets summary.txt"
        filename = "./Experiments MIP vs Naive vs Greedy/" * file
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

test_MIP = false; test_naive =false; test_greedy = true; test_greedy_LP = true
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
        for i in Base.OneTo(total_experiments_per_matrix_size)
        # for (i,A) in enumerate(A_vec)
        #     for (j,c) in enumerate(c_vec)
                
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
            println("Experiment $i of $total_experiments_per_matrix_size")
        # for (i,A) in enumerate(A_vec)
        #     for (j,c) in enumerate(c_vec)
                
                B_vec = [sum(c_vec[i]) / d for d = budget_denominators]
                for (k,B) in enumerate(B_vec)
                    println("Budget fraction: $(1 / budget_denominators[k])")
                    x, r, obj_val, time_achieved, time_elapsed, soln_attempts = solve_game_naive(A_vec[i], c_vec[i], B)
                    B_used = r' * c_vec[i]

                    write(f, "$i\t$i\t$B\t$(1 / budget_denominators[k])\t$B_used\t$obj_val\t$time_achieved\t$time_elapsed\t$soln_attempts\n")
                    flush(f)
                end

            # end
        end
        close(f)
    end

    #### Greedy Algorithm Naive
    if test_greedy_naive
        filename = subpath * "Matrices $num_rows by $num_cols greedy naive.txt"
        f = open(filename, "a")
        write(f, "num_rows = $num_rows\n")
        write(f, "num_cols = $num_cols\n")
        write(f, "TimeLimit = $TimeLimit\n")
        write(f, "MIPGap = $MIPGap\n")
        write(f, "Matrix seed\tCosts seed\tB\tBudget fraction\tBudget spent\tObj val\tTermination status\tSolve time\tNum purchases\tNodes\n")
        flush(f)
        for i in Base.OneTo(total_experiments_per_matrix_size)
            println("Experiment $i of $total_experiments_per_matrix_size")
        # for (i,A) in enumerate(A_vec)
        #     for (j,c) in enumerate(c_vec)
                
                B_vec = [sum(c_vec[i]) / d for d = budget_denominators]
                for (k,B) in enumerate(B_vec)
                    println("Budget fraction: $(1 / budget_denominators[k])")
                    x, r, obj_val, term_status, time_elapsed, num_purchases, _, _, _ = solve_game_greedy_naive(A_vec[i], c_vec[i], B, MIPGap=MIPGap, TimeLimit=TimeLimit)
                    B_used = r' * c_vec[i]

                    write(f, "$i\t$i\t$B\t$(1 / budget_denominators[k])\t$B_used\t$obj_val\t$term_status\t$time_elapsed\t$num_purchases\t$nodes\n")
                    flush(f)
                end

            # end
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
        write(f, "Matrix seed\tCosts seed\tB\tBudget fraction\tBudget spent\tObj val\tTermination status\tSolve time\tNum purchases\tNodes\n")
        flush(f)
        for i in Base.OneTo(total_experiments_per_matrix_size)
            println("Experiment $i of $total_experiments_per_matrix_size")
        # for (i,A) in enumerate(A_vec)
        #     for (j,c) in enumerate(c_vec)
                
                B_vec = [sum(c_vec[i]) / d for d = budget_denominators]
                for (k,B) in enumerate(B_vec)
                    println("Budget fraction: $(1 / budget_denominators[k])")
                    x, r, obj_val, term_status, time_elapsed, num_purchases, _, _, _ = solve_game_greedy_LP(A_vec[i], c_vec[i], B, MIPGap=MIPGap, TimeLimit=TimeLimit)
                    B_used = r' * c_vec[i]

                    write(f, "$i\t$i\t$B\t$(1 / budget_denominators[k])\t$B_used\t$obj_val\t$term_status\t$time_elapsed\t$num_purchases\t$nodes\n")
                    flush(f)
                end

            # end
        end
        close(f)
    end
end


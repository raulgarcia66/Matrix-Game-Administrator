work_dir = pwd()
include(work_dir * "/src/parameters.jl")
include(joinpath(pwd(), "src/solve.jl"))

##### Parameters
# A = create_matrix(-10:10, 6, 7)
A = create_matrix(-2:10, 7, seed=1);
# A = create_matrix_symmetric(-10:10, 90);

entry_range = [i/20 for i = -5:10]
seed = 3
A = create_matrix(entry_range, 50, seed=seed)

num_rows, num_cols = size(A)

c_r = create_cost_vector(2:5, num_rows, seed=seed);
c_s = create_cost_vector(10:12, num_cols, seed=seed);

B = (sum(c_r) + sum(c_s)) * 0.3
# B = sum(c_r) + sum(c_s)


########################################################################
x, r, s, obj_val, obj_bound, dual_obj, term_status, soln_time, rel_gap, nodes = solve_game(A, c_r, c_s, B, TimeLimit=120)
B_used = r' * c_r + s' * c_s
R = filter(i -> r[i] == 1, eachindex(r))
S = filter(j -> s[j] == 1, eachindex(s))

opt_val, x_opt, r_opt, s_opt, R_opt, S_opt, B_used_opt = copy(obj_val), copy(x), copy(r), copy(s), copy(R), copy(S), copy(B_used)

x, r, s, obj_val, dual_var_s, _, _  = solve_game_LP(A, c_r, c_s, B)
x, r, s, obj_val = solve_game_greedy_frequency(A, c_r, c_s, B)

### Containment
B_1 = 20
x, r, s, obj_val, obj_bound, dual_obj, term_status, soln_time, rel_gap, nodes = solve_game(A, c_r, c_s, B_1)
B_used = r' * c_r + s' * c_s
obj_val_1 = copy(obj_val)
R_1 = filter(i -> r[i] == 1, eachindex(r))
S_1 = filter(j -> s[j] == 1, eachindex(s))

B_2 = 25
x, r, s, obj_val, obj_bound, dual_obj, term_status, soln_time, rel_gap, nodes = solve_game(A, c_r, c_s, B_2)
B_used = r' * c_r + s' * c_s
obj_val_2 = copy(obj_val)
R_2 = filter(i -> r[i] == 1, eachindex(r))
S_2 = filter(j -> s[j] == 1, eachindex(s))

########################################################################
r_fix = rand(num_rows) # rand([0,1], num_rows)   # values of 1 and 0 are fixed, all else ignored
s_fix = ones(num_cols) # rand([0,1], num_cols)
x, r, s, obj_val, term_status, soln_time, rel_gap, nodes = solve_game(A, c_r, c_s, sum(c_s) + sum(c_r), r_fix, s_fix)
B_used = r' * c_r + s' * c_s
R = filter(i -> r[i] == 1, eachindex(r))
S = filter(j -> s[j] == 1, eachindex(s))

### Check for tightness of z
tight_cols = []
smaller_cols = []
for j in 1:num_cols
    # println("$(maximum(A) - minimum(A) + A[:,j]' * x)")
    if obj_val == (maximum(A) - minimum(A[:,j]) + A[:,j]' * x)
        push!(tight_cols, j)
    elseif obj_val - 1E-6 > (maximum(A) - minimum(A[:,j]) + A[:,j]' * x)  # small tolerance
        push!(smaller_cols, j)
    end
end
tight_cols
smaller_cols

########################################################################
# Greedy
B = 30
x, r, s, obj_val, term_status, time_elapsed, num_purchases, nodes, R, S, sequence, obj_val_vec, B_used_vec = solve_game_greedy_MILP(A, c_r, c_s, B)
B_used = r' * c_r + s' * c_s

opt_val
obj_val
B_used_opt
B_used
x_opt
x
R_opt
R
S_opt
S
length(R_opt) + length(S_opt)
num_purchases

nodes
sequence
obj_val_vec
B_used_vec

gains = compute_gains(obj_val_vec)


###################################################################################
####################### Comparison with Greedy Approaches #######################
# MILP
# On: Cols 2:5 [10,50,100] × [100];
# Done: Cols 2:5 [10,50,100] × [10,50,100];
# Done: Cols 11:15 [10,50] × [10,50]; [] × []; 
# Done: Cols 21:25 [10,50] × [10,50]; [] × [];
# Greedy
# On: Cols 2:5 [10,50] × [10,50,100];
# Done: Cols 2:5 [] × []; [] × []; 
# Done: Cols 11:15 [] × []; [] × []; 
# Done: Cols 21:25 [] × []; [] × []; 
num_rows_vec = [10,50,100] # [10,50]
num_cols_vec = [10,50,100] # [10,50]
MIPGap = 1E-2
TimeLimit = 1200
total_experiments_per_matrix_size = 5
budget_denominators = [4/3, 2, 4]
matrix_entry_range = [i/20 for i = -5:10]
attack_scale_entry_range = [i * 10 for i = 1:5]
c_r_entry_range = 2:5
c_s_entry_range_vec = [2:5, 11:15, 21:25]
# c_s_entry_range = c_s_entry_ranges[3] # 6:9, 10:15, 15:20, 20:30, 30:35

# Logging
exp_type = "Set MILP"
exp_type = "Set Greedy MILP"
exp_type, ranking, approach = "Set Greedy Freq", "with dual var", "dual"
# exp_type, ranking, approach = "Set Greedy Freq", "with s", "simple"
set_num = 1
subpath = "./Experiments/$exp_type $set_num/"
mkpath(subpath)

filenames = String[]
push!(filenames, "$(exp_type) $set_num summary.txt")
push!(filenames, "$(exp_type) summary.txt")
# push!(filenames, "Set MILP summary.txt")
# push!(filenames, "Set Greedy Freq summary.txt")  # Contains both ranking approaches
# push!(filenames, "Set Greedy MILP summary.txt")
for file in filenames
    if file == "Set MILP summary.txt" || file == "Set Greedy Freq summary.txt" || file == "Set Greedy MILP summary.txt"
        filename = "./Experiments/" * file
        f = open(filename, "a")
    else
        filename = subpath * file
        f = open(filename, "w")  # normally
        # f = open(filename, "a")
    end
    write(f, "$exp_type $set_num\n")
    write(f, "num_rows_vec = $num_rows_vec\n")
    write(f, "num_cols_vec = $num_cols_vec\n")
    write(f, "TimeLimit = $TimeLimit\n")
    write(f, "MIPGap = $MIPGap\n")
    write(f, "total_experiments_per_matrix_size = $total_experiments_per_matrix_size\n")
    write(f, "budget_fraction = $(1 ./ budget_denominators)\n")
    write(f, "matrix_entry_range = $matrix_entry_range\n")
    write(f, "attack_scale_entry_range = $attack_scale_entry_range\n")
    write(f, "c_r_entry_range = $c_r_entry_range\n")
    write(f, "c_s_entry_range_vec = $c_s_entry_range_vec\n")
    write(f, "\n")
    close(f)
end

# Experiments
test_MILP = false; test_greedy_MILP = false; test_greedy_freq = true;
for num_rows in num_rows_vec, num_cols in num_cols_vec, (c_s_ind, c_s_entry_range) in pairs(c_s_entry_range_vec)
    attack_scale_vec = map(seed -> create_attack_scale_vector(attack_scale_entry_range, num_cols, seed=seed), Base.OneTo(total_experiments_per_matrix_size))
    A_vec = map(seed -> create_matrix(matrix_entry_range, attack_scale_vec[seed], num_rows, num_cols, seed=seed), Base.OneTo(total_experiments_per_matrix_size))
    c_r_vec = map(seed -> create_cost_vector(c_r_entry_range, num_rows, seed=seed), Base.OneTo(total_experiments_per_matrix_size))
    c_s_vec = map(seed -> create_cost_vector(c_s_entry_range, num_cols, seed=seed), Base.OneTo(total_experiments_per_matrix_size))

    #### MILP
    if test_MILP
        filename = subpath * "Matrices $num_rows by $num_cols column prices index $(c_s_ind) MILP.txt"
        f = open(filename, "a")
        write(f, "TimeLimit = $TimeLimit\n")
        write(f, "MIPGap = $MIPGap\n")
        write(f, "Matrix seed\tCosts seed\tNum rows\tNum cols\tB\tBudget fraction\tBudget spent\tObj val\tRows purchased\tCols purchased\tObj bound\tDual obj\tTermination status\tSolve time\tRelative gap\tNode count\n")
        flush(f)
        for i in Base.OneTo(total_experiments_per_matrix_size)
                
            B_vec = [(sum(c_r_vec[i]) + sum(c_s_vec[i])) / d for d = budget_denominators]
            for (k,B) in enumerate(B_vec)
                x, r, s, obj_val, obj_bound, dual_obj, term_status, soln_time, rel_gap, nodes = solve_game(A_vec[i], c_r_vec[i], c_s_vec[i], B, MIPGap=MIPGap, TimeLimit=TimeLimit)
                B_used = r' * c_r_vec[i] + s' * c_s_vec[i]
                R = filter(i -> r[i] == 1, eachindex(r))
                S = filter(j -> s[j] == 1, eachindex(s))

                write(f, "$i\t$i\t$num_rows\t$num_cols\t$B\t$(1 / budget_denominators[k])\t$B_used\t$obj_val\t$(length(R))\t$(length(S))\t$obj_bound\t$dual_obj\t$term_status\t$soln_time\t$rel_gap\t$nodes\n")
                flush(f)
            end

            # end
        end
        close(f)
    end

    #### Greedy MIP Algorithm
    if test_greedy_MILP
        filename = subpath * "Matrices $num_rows by $num_cols column prices index $(c_s_ind) greedy MILP.txt"
        f = open(filename, "a")
        write(f, "TimeLimit = $TimeLimit\n")
        write(f, "MIPGap = $MIPGap\n")
        write(f, "Matrix seed\tCosts seed\tNum rows\tNum cols\tB\tBudget fraction\tBudget spent\tObj val\tRows purchased\tCols purchased\tTermination status\tSolve time\tNum purchases\tNodes\n")
        flush(f)
        for i in Base.OneTo(total_experiments_per_matrix_size)
                
            B_vec = [(sum(c_r_vec[i]) + sum(c_s_vec[i])) / d for d = budget_denominators]
            for (k,B) in enumerate(B_vec)
                x, r, s, obj_val, term_status, time_elapsed, num_purchases, nodes, R, S, sequence, obj_val_vec, B_used_vec = solve_game_greedy_MILP(A_vec[i], c_r_vec[i], c_s_vec[i], B, MIPGap=MIPGap, TimeLimit=TimeLimit)
                B_used = r' * c_r_vec[i] + s' * c_s_vec[i]

                write(f, "$i\t$i\t$num_rows\t$num_cols\t$B\t$(1 / budget_denominators[k])\t$B_used\t$obj_val\t$(length(R))\t$(length(S))\t$term_status\t$time_elapsed\t$num_purchases\t$nodes\n")
                flush(f)
            end
        end
        close(f)
    end

    #### Greedy Frequency Algorithm
    if test_greedy_freq
        filename = subpath * "Matrices $num_rows by $num_cols column prices index $(c_s_ind) greedy freq ranking $ranking.txt"
        f = open(filename, "a")
        write(f, "TimeLimit = $TimeLimit\n")
        write(f, "MIPGap = $MIPGap\n")
        # write(f, "Matrix seed\tCosts seed\tNum rows\tNum cols\tB\tBudget fraction\tBudget spent\tRows purchased\tCols purchased\tObj val\tTermination status\tSolve time\n")
        write(f, "Matrix seed\tCosts seed\tNum rows\tNum cols\tB\tBudget fraction\tBudget spent\tObj val\tRows purchased\tCols purchased\n")
        flush(f)
        for i in Base.OneTo(total_experiments_per_matrix_size)
                
            B_vec = [(sum(c_r_vec[i]) + sum(c_s_vec[i])) / d for d = budget_denominators]
            for (k,B) in enumerate(B_vec)
                local r,s
                x, r, s, obj_val = solve_game_greedy_frequency(A_vec[i], c_r_vec[i], c_s_vec[i], B, MIPGap=MIPGap, TimeLimit=TimeLimit, approach=approach)
                B_used = r' * c_r_vec[i] + s' * c_s_vec[i]
                R = filter(i -> r[i] == 1, eachindex(r))
                S = filter(j -> s[j] == 1, eachindex(s))

                # write(f, "$i\t$i\t$num_rows\t$num_cols\t$B\t$(1 / budget_denominators[k])\t$B_used\t$obj_val\t$(length(R))\t$(length(S))\t$term_status\t$time_elapsed\n")
                write(f, "$i\t$i\t$num_rows\t$num_cols\t$B\t$(1 / budget_denominators[k])\t$B_used\t$obj_val\t$(length(R))\t$(length(S))\n")
                flush(f)
            end
        end
        close(f)
    end
end

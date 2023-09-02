work_dir = pwd()
include(work_dir * "/src/parameters.jl")
include(joinpath(pwd(), "src/solve.jl"))

#########################################################################################
################################# Conditional Dominance #################################

matrix_entry_range = [i/20 for i = -5:10]
attack_scale_entry_range = [i * 10 for i = 1:5]

num_rows, num_cols = 50, 50

seed = 2
attack_scales = create_attack_scale_vector(attack_scale_entry_range, num_cols, seed=seed)
A = create_matrix(matrix_entry_range, attack_scales, num_rows, num_cols, seed=seed)
c_r = create_cost_vector(2:5, num_rows, seed=seed);
c_s = create_cost_vector(10:12, num_cols, seed=seed);
# # Make a row conditionally dominated
# cond_dominated_row = rand(1:num_rows)
# cond_dominating_row = cond_dominated_row == num_rows ? cond_dominated_row - 1 : cond_dominated_row + 1
# A[cond_dominated_row,:] = 0.5 * A[cond_dominating_row,:]
# except_col = rand(1:num_cols)
# A[cond_dominated_row, except_col] = A[cond_dominating_row, except_col] + 10

# c_r = create_cost_vector(2:5, num_rows, seed=seed);
# c_s = create_cost_vector(10:12, num_cols, seed=seed);
# c_r[cond_dominated_row] = 1 + c_r[cond_dominating_row]

# Make a number of rows conditionally dominated
num_cond_dom_rows = 4
rows_to_select_from = randperm(num_rows)[1:2*num_cond_dom_rows]
cond_dominated_rows = copy(rows_to_select_from[1:num_cond_dom_rows])
cond_dominating_rows = copy(rows_to_select_from[num_cond_dom_rows+1:end])
except_cols = randperm(num_cols)[1:num_cond_dom_rows]

for (i,k) in zip(cond_dominated_rows, cond_dominating_rows)
    A[i,:] = A[k,:] .- 10
end

for (i,k,j) in zip(cond_dominated_rows, cond_dominating_rows, except_cols)
    A[i,j] = A[k,j] + 10
end

foreach((i,j) -> c_r[i] = 1 + c_r[j], cond_dominated_rows, cond_dominating_rows)

B = (sum(c_r) + sum(c_s)) * 0.5

relax = false

x, r, s, obj_val, obj_bound, dual_obj, term_status, soln_time, _, _ = solve_game(A, c_r, c_s, B, relax=relax, MIPGap=0.01)
R = filter(i -> r[i] == 1, eachindex(r))
S = filter(j -> s[j] == 1, eachindex(s))
b_spent = c_r' * r + c_s' * s

cond_dominated_data = [(i, j) for (i,j) in zip(cond_dominated_rows, except_cols)]
x_cuts, r_cuts, s_cuts, obj_val_cuts, obj_bound_cuts, dual_obj_cuts, term_status_cuts, soln_time_cuts, _, _ = solve_game_with_cuts(A, c_r, c_s, B, cond_dominated_data, relax=relax, MIPGap=0.01)
R_cuts = filter(i -> r_cuts[i] == 1, eachindex(r_cuts))
S_cuts = filter(j -> s_cuts[j] == 1, eachindex(s_cuts))
b_spent_cuts = c_r' * r_cuts + c_s' * s_cuts

obj_val
obj_val_cuts
x
x_cuts
R
R_cuts
S
S_cuts
soln_time
soln_time_cuts

setdiff(R, R_cuts)
setdiff(R_cuts, R)
setdiff(S, S_cuts)
setdiff(S_cuts, S)
intersect(R, cond_dominated_rows)
intersect(R_cuts, cond_dominated_rows)
intersect(R, cond_dominating_rows)
intersect(R_cuts, cond_dominating_rows)
intersect(S, except_cols)
intersect(S_cuts, except_cols)

cond_dominated_rows
cond_dominating_rows
except_cols


#### Experiments
num_rows_vec = [10,50,100]
num_cols_vec = [10,50,100]
MIPGap = 1E-2
TimeLimit = 600
total_experiments_per_matrix_size = 1
budget_denominators = [4/3, 2, 4]
matrix_entry_range = [i/20 for i = -5:10]
attack_scale_entry_range = [i * 10 for i = 1:5]
c_r_entry_range = 2:5
# c_s_entry_range = 2:5
c_s_entry_range_vec = [2:5, 11:15, 21:25]
num_cond_dom_rows_vec = [1,5,10]

# Logging
set_type, relax = "Set MILP Cuts", false
set_type, relax = "Set LP Cuts", true
set_num = 1
subpath = "./Experiments/$set_type $set_num/"
mkpath(subpath)

filenames = String[]
push!(filenames, "$(set_type) $set_num summary.txt")
push!(filenames, "$(set_type)s summary.txt")
# push!(filenames, "$(set_type)s Greedy Freq summary.txt")
for file in filenames
    if file == "$(set_type)s summary.txt"
        filename = "./Experiments/" * file
        f = open(filename, "a")
    else
        filename = subpath * file
        f = open(filename, "w")  # normally
        # f = open(filename, "a")
    end
    write(f, "$set_type $set_num\n")
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
    write(f, "num_cond_dom_rows_vec = $num_cond_dom_rows_vec\n")
    write(f, "relax = $relax\n")
    write(f, "\n")
    close(f)
end

for num_rows in num_rows_vec, num_cols in num_cols_vec, (c_s_ind, c_s_entry_range) in pairs(c_s_entry_range_vec)
    attack_scale_vec = map(seed -> create_attack_scale_vector(attack_scale_entry_range, num_cols, seed=seed), Base.OneTo(total_experiments_per_matrix_size))
    A_vec = map(seed -> create_matrix(matrix_entry_range, attack_scale_vec[seed], num_rows, num_cols, seed=seed), Base.OneTo(total_experiments_per_matrix_size))
    c_r_vec = map(seed -> create_cost_vector(c_r_entry_range, num_rows, seed=seed), Base.OneTo(total_experiments_per_matrix_size))
    c_s_vec = map(seed -> create_cost_vector(c_s_entry_range, num_cols, seed=seed), Base.OneTo(total_experiments_per_matrix_size))

    for exp_type in ["cuts", "no cuts"]
        Random.seed!(1)  # Need randomness to be the same for both methods

        filename = relax ? subpath * "Matrices $num_rows by $num_cols column prices index $(c_s_ind) LP $exp_type.txt" :
                subpath * "Matrices $num_rows by $num_cols column prices index $(c_s_ind) MILP $exp_type.txt"
        # filename = subpath * "Matrices $num_rows by $num_cols column prices index $(c_s_ind) MILP $exp_type.txt"
        f = open(filename, "a")
        write(f, "TimeLimit = $TimeLimit\n")
        write(f, "MIPGap = $MIPGap\n")
        write(f, "Matrix seed\tCosts seed\tNum rows\tNum cols\tc_s range\tNum cond dom rows\tB\tBudget fraction\tBudget spent\tObj val\tRows purchased\tCols purchased\tObj bound\tDual obj\tTermination status\tSolve time\tRelative gap\tNode count\n")
        flush(f)

        for num_cond_dom_rows in num_cond_dom_rows_vec
            if num_cond_dom_rows == 10 && (num_rows == 10 || num_cols == 10)
                continue
            end

            for exp_num in Base.OneTo(total_experiments_per_matrix_size)
                # Make a number of rows conditionally dominated
                rows_to_select_from = randperm(num_rows)[1:2*num_cond_dom_rows]
                cond_dominated_rows = copy(rows_to_select_from[1:num_cond_dom_rows])
                cond_dominating_rows = copy(rows_to_select_from[num_cond_dom_rows+1:end])
                except_cols = randperm(num_cols)[1:num_cond_dom_rows]

                local A, c_r, c_s
                A = deepcopy(A_vec[exp_num])
                c_r = deepcopy(c_r_vec[exp_num])
                c_s = deepcopy(c_s_vec[exp_num])

                for (i,k) in zip(cond_dominated_rows, cond_dominating_rows)
                    A[i,:] = A[k,:] .- 10
                end
                
                for (i,k,j) in zip(cond_dominated_rows, cond_dominating_rows, except_cols)
                    A[i,j] = A[k,j] + 10
                end

                foreach((i,j) -> c_r[i] = 1 + c_r[j], cond_dominated_rows, cond_dominating_rows)

                # foreach((i,j) -> A_vec[exp_num][i,:] = 0.5 * A_vec[exp_num][j,:], cond_dominated_rows, cond_dominating_rows)
                # foreach(i -> foreach(j -> begin 
                #     if A_vec[exp_num][i,j] < 0
                #         A_vec[exp_num][i,j] = -A_vec[exp_num][i,j]
                #     end
                # end, eachindex(A_vec[exp_num][i,:])), cond_dominated_rows)

                # foreach((i,j) -> c_r_vec[exp_num][i] = 1 + c_r_vec[exp_num][j], cond_dominated_rows, cond_dominating_rows)

                B_vec = [(sum(c_r) + sum(c_s)) / d for d = budget_denominators]

                if exp_type == "cuts"
                    
                    cond_dominated_data = [(i, j) for (i,j) in zip(cond_dominated_rows, except_cols)]

                    for (k,B) in enumerate(B_vec)
                        # x, r, s, obj_val, obj_bound, dual_obj, term_status, soln_time, rel_gap, nodes = solve_game_with_cuts(A_vec[exp_num], c_r_vec[exp_num], c_s_vec[exp_num], B, cond_dominated_data, MIPGap=0.01, TimeLimit=TimeLimit, relax=relax)
                        x, r, s, obj_val, obj_bound, dual_obj, term_status, soln_time, rel_gap, nodes = solve_game_with_cuts(A, c_r, c_s, B, cond_dominated_data, MIPGap=0.01, TimeLimit=TimeLimit, relax=relax)
                        R = filter(i -> r[i] == 1, eachindex(r))
                        S = filter(j -> s[j] == 1, eachindex(s))
                        # B_used = r' * c_r_vec[exp_num] + s' * c_s_vec[exp_num]
                        B_used = r' * c_r + s' * c_s

                        write(f, "$exp_num\t$exp_num\t$num_rows\t$num_cols\t$(c_s_entry_range)\t$num_cond_dom_rows\t$B\t$(1 / budget_denominators[k])\t$B_used\t$obj_val\t$(length(R))\t$(length(S))\t$obj_bound\t$dual_obj\t$term_status\t$soln_time\t$rel_gap\t$nodes\n")
                        flush(f)
                    end

                elseif exp_type == "no cuts"
                
                    for (k,B) in enumerate(B_vec)
                        # x, r, s, obj_val, obj_bound, dual_obj, term_status, soln_time, rel_gap, nodes = solve_game(A_vec[exp_num], c_r_vec[exp_num], c_s_vec[exp_num], B, MIPGap=0.01, TimeLimit=TimeLimit, relax=relax)
                        x, r, s, obj_val, obj_bound, dual_obj, term_status, soln_time, rel_gap, nodes = solve_game(A, c_r, c_s, B, MIPGap=0.01, TimeLimit=TimeLimit, relax=relax)
                        R = filter(i -> r[i] == 1, eachindex(r))
                        S = filter(j -> s[j] == 1, eachindex(s))
                        # B_used = r' * c_r_vec[exp_num] + s' * c_s_vec[exp_num]
                        B_used = r' * c_r + s' * c_s

                        write(f, "$exp_num\t$exp_num\t$num_rows\t$num_cols\t$(c_s_entry_range)\t$num_cond_dom_rows\t$B\t$(1 / budget_denominators[k])\t$B_used\t$obj_val\t$(length(R))\t$(length(S))\t$obj_bound\t$dual_obj\t$term_status\t$soln_time\t$rel_gap\t$nodes\n")
                        flush(f)
                    end
                end

            end
        end
        close(f)
    end
end

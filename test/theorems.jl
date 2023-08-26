work_dir = pwd()
include(work_dir * "/src/parameters.jl")
include(joinpath(pwd(), "src/solve.jl"))


##### Conditional Dominance

matrix_entry_range = [i/20 for i = -5:10]
attack_scale_entry_range = [i * 10 for i = 1:5]

num_rows, num_cols = 75, 75

seed = 2
attack_scales = create_attack_scale_vector(attack_scale_entry_range, num_cols, seed=seed)
A = create_matrix(matrix_entry_range, attack_scales, num_rows, num_cols, seed=seed)
# Make a row conditionally dominated
cond_dominated_row = rand(1:num_rows)
cond_dominating_row = cond_dominated_row == num_rows ? cond_dominated_row - 1 : cond_dominated_row + 1
A[cond_dominated_row,:] = 0.5 * A[cond_dominating_row,:]
except_col = rand(1:num_cols)
A[cond_dominated_row, except_col] = A[cond_dominating_row, except_col] + 1

c_r = create_cost_vector(2:5, num_rows, seed=seed);
c_r[cond_dominated_row] = 1 + c_r[cond_dominating_row]
c_s = create_cost_vector(10:12, num_cols, seed=seed);

B = (sum(c_r) + sum(c_s)) * 0.5

x, r, s, obj_val, obj_bound, dual_obj, term_status, soln_time, _, _ = solve_game(A, c_r, c_s, B, relax=false, MIPGap=0.01)
R = filter(i -> r[i] == 1, eachindex(r))
S = filter(j -> s[j] == 1, eachindex(s))

cond_dominated_data = [(cond_dominated_row, except_col)]
x_cuts, r_cuts, s_cuts, obj_val_cuts, obj_bound_cuts, dual_obj_cuts, term_status_cuts, soln_time_cuts, _, _ = solve_game_with_cuts(A, c_r, c_s, B, cond_dominated_data, relax=false, MIPGap=0.01)
R_cuts = filter(i -> r_cuts[i] == 1, eachindex(r_cuts))
S_cuts = filter(j -> s_cuts[j] == 1, eachindex(s_cuts))

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
work_dir = pwd()
include(work_dir * "/src/parameters.jl")
include(joinpath(pwd(), "src/solve.jl"))

##### Parameters
A = create_matrix(-10:10, 6, 7)
A = create_matrix(-2:10, 50);
A = create_matrix_symmetric(-10:10, 90);

num_rows, num_cols = size(A)

c_r = create_cost_vector(2:5, num_rows);
c_s = create_cost_vector(5:10, num_cols);

B = (sum(c_r) + sum(c_s)) * 0.4
# B = sum(c_r) + sum(c_s)

########################
x, r, s, obj_val, obj_bound, dual_obj, term_status, soln_time, rel_gap, nodes = solve_game(A, c_r, c_s, B) #, TimeLimit=30)
B_used = r' * c_r + s' * c_s
R = filter(i -> r[i] == 1, eachindex(r))
S = filter(j -> s[j] == 1, eachindex(s))

opt_val, x_opt, r_opt, s_opt, R_opt, S_opt, B_used_opt = copy(obj_val), copy(x), copy(r), copy(s), copy(R), copy(S), copy(B_used)

########################
r_fix = rand(num_rows) # rand([0,1], num_rows)   # values of 1 and 0 are fixed, all else ignored
s_fix = ones(num_cols) # rand([0,1], num_cols)
x, r, s, obj_val, term_status, soln_time, rel_gap, nodes = solve_game(A, c_r, c_s, sum(c_s) + sum(c_r), r_fix, s_fix)
B_used = r' * c_r + s' * c_s
R = filter(i -> r[i] == 1, eachindex(r))
S = filter(j -> s[j] == 1, eachindex(s))

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

########################
# Greedy
x, r, s, obj_val, term_status, time_elapsed, num_purchases, nodes, R, S, sequence, obj_val_vec, B_spent_vec = solve_game_greedy(A, c_r, c_s, B)
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

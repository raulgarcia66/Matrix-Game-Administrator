using Gurobi
using JuMP
import Random

function compute_big_M_parameters(A::Matrix{T}) where T
    a_max = maximum(A)
    M = map(j -> a_max - minimum(A[:,j]), 1:size(A,2))
    return M
end

"""
Solve matrix game designer (MGD).
"""
function solve_game(A::Matrix{T}, c_r::Vector{U}, c_s::Vector{V}, B::W; relax::Bool=false, TimeLimit::X=300, MIPGap::Y=0.01) where {T,U,V,W,X,Y <: Real}
    
    num_rows, num_cols = size(A)

    #### Initialize model
    # model = Model(Gurobi.Optimizer)
    model = Model(optimizer_with_attributes(Gurobi.Optimizer, "MIPGap"=>MIPGap, "TimeLimit"=>TimeLimit))

    @variable(model, x[1:num_rows] >= 0)   # Don't need upper bound since sum(x) == 1 will enforce it
    @variable(model, r[1:num_rows], Bin)   # 1 means play is available
    @variable(model, s[1:num_cols], Bin)   # 1 means play is available
    @variable(model , z)

    @objective(model, Max, z)

    M = compute_big_M_parameters(A)
    @constraint(model, [j = 1:num_cols], z - A[:,j]' * x <= M[j] * s[j] )
    @constraint(model, sum(c_r' * r) + sum(c_s' * s) <= B)
    @constraint(model, sum(x) == 1)
    @constraint(model, x .<= r)

    #### Solve
    if relax
        relax_integrality(model)
    end
    optimize!(model)

    if termination_status(model) == INFEASIBLE_OR_UNBOUNDED || termination_status(model) == INFEASIBLE
        error("Model is infeasible.")
    elseif termination_status(model) == TIME_LIMIT && !has_values(model)
        error("No primal solution obtained within the time limit.")
    end

    if relax
        rel_gap = 0
        nodes = 0
    else
        rel_gap = relative_gap(model)
        nodes = node_count(model)
    end

    return value.(x), value.(r), value.(s), objective_value(model), objective_bound(model), dual_objective_value(model), termination_status(model), solve_time(model), rel_gap, nodes #, result_count(model)
end

"""
Solve MGD_row with certain rows selections and column removals fixed to purchase and to not purchase (entries of r_fix and s_fix equal to 1 and 0, resp.). 
If a positive k is given, need to purchase k many rows.
"""
function solve_game(A::Matrix{T}, c_r::Vector{U}, c_s::Vector{U}, B::W, r_fix::Vector{V}, s_fix::Vector{V}; 
        relax::Bool=false, TimeLimit::X=300, MIPGap::Y=0.01, k::Int=-1) where {T,U,V,W,X,Y <: Real}
    
    num_rows, num_cols = size(A)

    #### Initialize model
    # model = Model(Gurobi.Optimizer)
    model = Model(optimizer_with_attributes(Gurobi.Optimizer, "MIPGap"=>MIPGap, "TimeLimit"=>TimeLimit))

    @variable(model, x[1:num_rows] >= 0)   # Don't need upper bound since sum(x) == 1 will enforce it
    @variable(model, r[1:num_rows], Bin)   # 1 means play is available
    @variable(model, s[1:num_cols], Bin)   # 1 means play is available
    @variable(model , z)

    @objective(model, Max, z)

    M = compute_big_M_parameters(A)
    @constraint(model, [j = 1:num_cols], z - A[:,j]' * x <= M[j] * s[j] )
    @constraint(model, sum(c_r' * r) + sum(c_s' * s) <= B)
    @constraint(model, sum(x) == 1)
    @constraint(model, x .<= r)
    if k > 0
        @constraint(model, sum(r) + sum(s) == k)
    end

    #### Fix values of r (1 or 0 will be fixed, all other values in r_fix are ignored)
    if any(e -> e == 1 || e == 0, r_fix)
        for i in eachindex(r_fix)
            if r_fix[i] == 1
                @constraint(model, r[i] == 1)
            elseif r_fix[i] == 0
                @constraint(model, r[i] == 0)
            end
        end
    end
    #### Fix values of s (1 or 0 will be fixed, all other values in s_fix are ignored)
    if any(e -> e == 1 || e == 0, s_fix)
        for i in eachindex(s_fix)
            if s_fix[i] == 1
                @constraint(model, s[i] == 1)
            elseif s_fix[i] == 0
                @constraint(model, s[i] == 0)
            end
        end
    end

    #### Solve
    if relax
        relax_integrality(model)
    end
    optimize!(model)

    if termination_status(model) == INFEASIBLE_OR_UNBOUNDED || termination_status(model) == INFEASIBLE
        error("Model is infeasible.")
    elseif termination_status(model) == TIME_LIMIT && !has_values(model)
        error("No primal solution obtained within the time limit.")
    end

    if relax
        rel_gap = 0
        nodes = 0
    else
        rel_gap = relative_gap(model)
        nodes = node_count(model)
    end

    return value.(x), value.(r), value.(s), objective_value(model), termination_status(model), solve_time(model), rel_gap, nodes #, result_count(model)
end

"""
Solve MGD via a Greedy algorithm which purchases the row selection or column removal that provides the greatest increase in value.
The best row selection or column removal is computed by solving an LP for each unpurchased row selection/column removal.
"""
function solve_game_greedy(A::Matrix{T}, c_r::Vector{U}, c_s::Vector{U}, B::V; TimeLimit::W=300, MIPGap::X=0.01) where {T,U,V,W,X <: Real}

    num_rows, num_cols = size(A)

    term_status = "TIME LIMIT"
    obj_val = -Inf
    x = zeros(num_rows)
    r = zeros(num_rows)
    s = zeros(num_cols)

    num_purchases = 0
    nodes = 0

    r_vec = Int[]  # keep track of rows purchased
    s_vec = Int[]  # keep track of columns purchased
    obj_val_vec = Float64[]
    B_spent_vec = Float64[]
    
    start_time = time()

    while time() < start_time + TimeLimit
        r_input = map(e -> abs(e - 0) < 1E-3 ? -1 : e, r)  # r is a solution; need to only fix the 1's (i.e., the purchases)
        s_input = map(e -> abs(e - 0) < 1E-3 ? -1 : e, s)  # similarly for s

        x_sub, r_sub, s_sub, obj_val_sub, term_status_sub, _, nodes_sub = solve_game(A, c_r, c_s, B, r_input, s_input, MIPGap=MIPGap, 
                                                                        TimeLimit=maximum([ceil(Int,TimeLimit-(time()-start_time)), 1]), k=num_purchases+1)
        nodes += nodes_sub
        
        if term_status_sub == INFEASIBLE_OR_UNBOUNDED || term_status_sub == INFEASIBLE || (r_sub == ones(num_rows) && s_sub == ones(num_cols))
            term_status = "FINISHED"
            break
        end

        x = deepcopy(x_sub)
        r = deepcopy(r_sub)
        s = deepcopy(s_sub)
        obj_val = copy(obj_val_sub)
        num_purchases += 1

        row_selected = findfirst(i -> r[i] > 0.5 && !(i in r_vec), eachindex(r))
        column_removed = findfirst(j -> s[j] > 0.5 && !(j in s_vec), eachindex(s))
        if row_selected !== nothing
            push!(r_vec, row_selected)
        elseif column_removed !== nothing
            push!(s_vec, column_removed)
        end
        push!(obj_val_vec, obj_val)
        push!(B_spent_vec, c_r' * r + c_s *s)
    end

    time_elapsed = time() - start_time
    return x, r, s, obj_val, term_status, time_elapsed, num_purchases, nodes, r_vec, s_vec, obj_val_vec, B_spent_vec
end


############################################################################################
############################################################################################
# The following functions pertain to the MGD variant of solely purchasing rows.

"""
Solve matrix game designer (MGD) variant of row purchasing, MGD_row.
"""
function solve_game(A::Matrix{T}, c::Vector{U}, B::V; relax::Bool=false, TimeLimit::W=300, MIPGap::X=0.01) where {T,U,V,W,X <: Real}
    
    num_rows = size(A, 1)

    #### Initialize model
    # model = Model(Gurobi.Optimizer)
    model = Model(optimizer_with_attributes(Gurobi.Optimizer, "MIPGap"=>MIPGap, "TimeLimit"=>TimeLimit))

    @variable(model, x[1:num_rows] >= 0)   # Don't need upper bound since sum(x) == 1 will enforce it
    @variable(model, r[1:num_rows], Bin)   # 1 means play is available
    @variable(model , z)

    @objective(model, Max, z)

    @constraint(model, [j = 1:size(A,2)], z - A[:,j]' * x <= 0 )
    @constraint(model, sum(c' * r) <= B)
    @constraint(model, sum(x) == 1)
    @constraint(model, x .<= r)

    #### Solve
    if relax
        relax_integrality(model)
    end
    optimize!(model)

    if termination_status(model) == INFEASIBLE_OR_UNBOUNDED || termination_status(model) == INFEASIBLE
        error("Model is infeasible.")
    elseif termination_status(model) == TIME_LIMIT && !has_values(model)
        error("No primal solution obtained within the time limit.")
    end

    if relax
        rel_gap = 0
        nodes = 0
    else
        rel_gap = relative_gap(model)
        nodes = node_count(model)
    end

    return value.(x), value.(r), objective_value(model), objective_bound(model), dual_objective_value(model), termination_status(model), solve_time(model), rel_gap, nodes #, result_count(model)
end

"""
Solve MGD_row with certain rows fixed to purchase and to not purchase (entries of r_fix equal to 1 and 0, resp.). 
If a positive k is given, need to purchase k many rows.
"""
function solve_game(A::Matrix{T}, c::Vector{U}, B::V, r_fix::Vector{W}; relax::Bool=false, TimeLimit::X=300, MIPGap::Y=0.01, k::Int=-1) where {T,U,V,W,X,Y <: Real}
    
    num_rows = size(A, 1)

    #### Initialize model
    # model = Model(Gurobi.Optimizer)
    model = Model(optimizer_with_attributes(Gurobi.Optimizer, "MIPGap"=>MIPGap, "TimeLimit"=>TimeLimit, "LogToConsole"=>0))

    @variable(model, x[1:num_rows] >= 0)   # Don't need upper bound since sum(x) == 1 will enforce it
    @variable(model, r[1:num_rows], Bin)   # 1 means play is available
    @variable(model , z)

    @objective(model, Max, z)

    @constraint(model, [j = 1:size(A,2)], z - A[:,j]' * x <= 0 )
    @constraint(model, sum(c' * r) <= B)
    @constraint(model, sum(x) == 1)
    @constraint(model, x .<= r)
    if k > 0
        @constraint(model, sum(r) == k)
    end

    #### Fix values of r (1 or 0 will be fixed, all other values in r_fix are ignored)
    # if r_fix != -ones(Int, length(c))
    if any(e -> e == 1 || e == 0, r_fix)
        for i in eachindex(r_fix)
            if r_fix[i] == 1
                @constraint(model, r[i] == 1)
            elseif r_fix[i] == 0
                @constraint(model, r[i] == 0)
            end
        end
    end

    #### Solve
    if relax
        relax_integrality(model)
    end
    optimize!(model)

    if termination_status(model) == INFEASIBLE_OR_UNBOUNDED || termination_status(model) == INFEASIBLE
        return zeros(num_rows), r_fix, -Inf, termination_status(model), maximum([0.1, solve_time(model)]), 0, 0
        # error("Model is infeasible.")
    elseif termination_status(model) == TIME_LIMIT && !has_values(model)
        return zeros(num_rows), r_fix, -Inf, termination_status(model), TimeLimit, 0, 0
        # error("No primal solution obtained within the time limit.")
    end

    if relax
        rel_gap = 0
        nodes = 0
    else
        rel_gap = relative_gap(model)
        nodes = node_count(model)
    end

    return value.(x), value.(r), objective_value(model), termination_status(model), solve_time(model), rel_gap, nodes
end

"""
Solve MGD_row via a Greedy algorithm which purchases the row that provides the greatest increase in value.
The row is computed by solving an MILP.
"""
function solve_game_greedy(A::Matrix{T}, c::Vector{U}, B::V; TimeLimit::W=300, MIPGap::X=0.01) where {T,U,V,W,X <: Real}

    num_rows = size(A, 1)

    term_status = "TIME LIMIT"
    obj_val = -Inf
    x = zeros(num_rows)
    r = zeros(num_rows)

    num_purchases = 0
    nodes = 0

    r_vec = Int[]  # keep track of rows purchased
    obj_val_vec = Float64[]
    B_spent_vec = Float64[]
    
    start_time = time()

    while time() < start_time + TimeLimit
        r_input = map(e -> abs(e - 0) < 1E-3 ? -1 : e, r)  # r is a solution; need to only fix the 1's (i.e., the purchases)

        x_sub, r_sub, obj_val_sub, term_status_sub, _, _, nodes_sub = solve_game(A, c, B, r_input, MIPGap=MIPGap, TimeLimit=maximum([ceil(Int,TimeLimit-(time()-start_time)), 1]),
                                                                        k=num_purchases+1)
        nodes += nodes_sub
        
        if term_status_sub == INFEASIBLE_OR_UNBOUNDED || term_status_sub == INFEASIBLE || r_sub == ones(num_rows)
            term_status = "FINISHED"
            break
        end

        x = deepcopy(x_sub)
        r = deepcopy(r_sub)
        obj_val = copy(obj_val_sub)
        num_purchases += 1

        row_selected = findfirst(i -> r[i] > 0.5 && !(i in r_vec), eachindex(r))
        if row_selected !== nothing
            push!(r_vec, row_selected)
        end
        push!(obj_val_vec, obj_val)
        push!(B_spent_vec, c' * r)
    end

    time_elapsed = time() - start_time
    return x, r, obj_val, term_status, time_elapsed, num_purchases, nodes, r_vec, obj_val_vec, B_spent_vec
end

"""
Solve MGD_row via a Greedy algorithm which purchases the row that provides the greatest increase in value.
The row is computed by solving an LP for each unpurchased row.
"""
function solve_game_greedy_LP(A::Matrix{T}, c::Vector{U}, B::V; TimeLimit::W=300, MIPGap::X=0.01) where {T,U,V,W,X <: Real}

    num_rows = size(A, 1)

    term_status = "TIME LIMIT"
    obj_val = -Inf
    x = zeros(num_rows)
    r = zeros(num_rows)
    indices_remaining = collect(1:num_rows)
    num_purchases = 0

    r_vec = Int[]  # keep track of rows purchased
    obj_val_vec = Float64[]
    B_spent_vec = Float64[]

    start_time = time()

    while time() < start_time + TimeLimit
        best_obj = -Inf
        best_ind = 0
        x_best = zeros(num_rows)
        r_best = zeros(num_rows)
        local_term_status = "INFEASIBLE"

        for i in eachindex(indices_remaining)
            r_input = copy(r)
            r_input[indices_remaining[i]] = 1

            time_limit = maximum([ceil(Int,TimeLimit-(time()-start_time)), 0])
            if time_limit > 0
                x_sub, r_sub, obj_val_sub, term_status_sub, _, _, _ = solve_game(A, c, B, r_input, MIPGap=MIPGap, relax=true,
                                                                                TimeLimit=time_limit)
            else
                break
            end
            
            if term_status_sub == INFEASIBLE_OR_UNBOUNDED || term_status_sub == INFEASIBLE
                continue
            elseif obj_val_sub > best_obj  # will be triggered at least once if feasible
                local_term_status = "FEASIBLE"
                best_obj = copy(obj_val_sub)
                best_ind = copy(i)
                x_best = deepcopy(x_sub)
                r_best = deepcopy(r_sub)
            end
        end
        
        if local_term_status == "FEASIBLE"
            deleteat!(indices_remaining, best_ind)

            x = deepcopy(x_best)
            r = deepcopy(r_best)
            obj_val = copy(best_obj)
            num_purchases += 1

            row_selected = findfirst(i -> r[i] > 0.5 && !(i in r_vec), eachindex(r))
            if row_selected !== nothing
                push!(r_vec, row_selected)
            end
            push!(obj_val_vec, obj_val)
            push!(B_spent_vec, c' * r)
        else
            term_status = "FINISHED"
            break
        end

        # if local_term_status == "INFEASIBLE"
        #     term_status = "FINISHED"
        #     break
        # end

        # x = deepcopy(x_best)
        # r = deepcopy(r_best)
        # obj_val = copy(best_obj)
        # num_purchases += 1

        # row_selected = findfirst(i -> r[i] > 0.5 && !(i in r_vec), eachindex(r))
        # if row_selected !== nothing
        #     push!(r_vec, row_selected)
        # end
        # push!(obj_val_vec, obj_val)
        # push!(B_spent_vec, c' * r)
    end

    time_elapsed = time() - start_time
    return x, r, obj_val, term_status, time_elapsed, num_purchases, r_vec, obj_val_vec, B_spent_vec
end

"""
Solve MGD_row by choosing random subsets of rows and solving the corresponding LPs for TimeLimit minutes.
"""
function solve_game_naive(A::Matrix{T}, c::Vector{U}, B::V; seed::Int=-1, TimeLimit::W=300) where {T,U,V,W <: Real}

    num_rows = size(A, 1)
    if seed != -1
        Random.seed!(seed)
    end

    time_achieved = 0
    max_obj = -Inf
    soln_attempts = 0
    x = -ones(num_rows)
    r = -ones(num_rows)
    # r_prev = -ones(num_rows)   # would need to keep track of all r's, not just the previous
    start_time = time()
    
    while time() < start_time + TimeLimit
        soln_attempts += 1
        r_fix = compute_r_fix(num_rows)
        if sum(r_fix' * c) > B
            continue
        end

        x_sub, _, obj_val, term_status_sub, _, _, _ = solve_game(A, c, B, r_fix, TimeLimit=maximum([ceil(Int,TimeLimit-(time()-start_time)), 1]), relax=true)

        if term_status_sub != INFEASIBLE_OR_UNBOUNDED && term_status_sub != INFEASIBLE
            println("Current max: $max_obj\nCurrent obj: $obj_val\n")
            if obj_val > max_obj
                # println("Setting max_obj to obj_val")
                max_obj = obj_val
                x = deepcopy(x_sub)
                r = deepcopy(r_fix)
                time_achieved = time() - start_time
            end
        end
    end

    time_elapsed = time() - start_time
    return x, r, max_obj, time_achieved, time_elapsed, soln_attempts
end

"""
Solve MGD_row by choosing random subsets of rows and solve the corresponding LPs for TimeLimit seconds.
If value is within a tolerance of opt_val, the procedure ends.
"""
function solve_game_naive_compare(A::Matrix{T}, c::Vector{U}, B::V, opt_val::W; seed::Int=-1, TimeLimit::X=300) where {T,U,V,W,X <: Real}

    num_rows = size(A, 1)
    if seed != -1
        Random.seed!(seed)
    end

    time_achieved = 0
    max_obj = -Inf
    soln_attempts = 0
    term_status = "FEASIBLE"
    x = -ones(num_rows)
    r = -ones(num_rows)
    gap = Inf
    # r_prev = -ones(num_rows)   # would need to keep track of all r's, not just the previous
    
    start_time = time()

    while time() < start_time + TimeLimit
        soln_attempts += 1
        r_fix = compute_r_fix(num_rows)
        if sum(r_fix' * c) > B
            continue
        end

        x_sub, _, obj_val, term_status_sub, _, _, _ = solve_game(A, c, B, r_fix, TimeLimit=maximum([ceil(Int,TimeLimit-(time()-start_time)), 1]), relax=true)

        if term_status_sub != INFEASIBLE_OR_UNBOUNDED && term_status_sub != INFEASIBLE
            println("Current max: $max_obj\nCurrent obj: $obj_val\n")
            if obj_val > max_obj
                # println("Setting max_obj to obj_val")
                max_obj = obj_val
                x = deepcopy(x_sub)
                r = deepcopy(r_fix)
                time_achieved = time() - start_time
            end

            gap = minimum([gap, abs(obj_val - opt_val)])
            if gap < 1E-2
                term_status = "OPTIMAL"
                break
            end
        end

    end

    time_elapsed = time() - start_time
    return x, r, max_obj, time_achieved, term_status, time_elapsed, gap, soln_attempts
end

"""
Compute a random subset of indices within the range 1:num_rows.
"""
function compute_subset(num_rows)
    # This approach returns a vector with about μ many ones
    μ = rand()
    r = rand(num_rows)
    rows = filter(ind -> r[ind] > μ, eachindex(r))
    if isempty(rows)
        rows = [rand(1:num_rows)]  # vector will have one row
    end
    return rows
end

"""
Compute a binary vector with at least one nonzero entry.
"""
function compute_r_fix(num_rows)
    # # This approach returns a vector with roughly half ones and half zeros
    # return rand([0,1], num_rows)

    # # This approach explicitly computes subsets of rows
    # rows = compute_subset(num_rows)
    # r = [Int(row in rows) for row in 1:num_rows]
    # return r

    # This approach returns a vector with about μ many ones
    μ = rand()
    r_frac = rand(num_rows)
    r = [Int(e > μ) for e in r_frac]
    if sum(r) == 0
        r[rand(1:num_rows)] = 1  # one entry will be 1
    end
    return r
end

# """
# Compute the power set of a given vector x.
# """
# function powerset(x::Vector{T}) where T
#     result = Vector{T}[[]]
#     for elem in x, j in eachindex(result)
#         push!(result, [result[j] ; elem])
#     end
#     return result[2:end]
# end

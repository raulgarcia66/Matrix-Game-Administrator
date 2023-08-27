using Gurobi
using JuMP
using Pipe
import Random

function compute_big_M_parameters(A::Matrix{T}) where T
    a_max = maximum(A)
    M = map(j -> a_max - minimum(A[:,j]), 1:size(A,2))
    return M
end

function compute_gains(obj_val_vec::Vector{T}) where T
    return map(i -> obj_val_vec[i+1] - obj_val_vec[i], collect(eachindex(obj_val_vec))[1:end-1])
end

"""
Solve matrix game designer (MGD).
"""
function solve_game(A::Matrix{T}, c_r::Vector{U}, c_s::Vector{U}, B::W; relax::Bool=false, TimeLimit::X=Inf, MIPGap::Y=0.01) where {T,U,W,X,Y <: Real}
    
    num_rows, num_cols = size(A)

    #### Initialize model
    # model = Model(Gurobi.Optimizer)
    model = Model(optimizer_with_attributes(Gurobi.Optimizer, "MIPGap"=>MIPGap, "TimeLimit"=>TimeLimit))

    @variable(model, x[1:num_rows] >= 0)   # Don't need upper bound since sum(x) == 1 will enforce it
    @variable(model, r[1:num_rows], Bin)   # 1 means row is available
    @variable(model, s[1:num_cols], Bin)   # 1 means column has been removed
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
function solve_game(A::Matrix{T}, c_r::Vector{U}, c_s::Vector{U}, B::W, r_fix::Vector{V}, s_fix::Vector{Z}; 
        relax::Bool=false, TimeLimit::X=Inf, MIPGap::Y=0.01, k::Int=-1) where {T,U,V,W,X,Y,Z <: Real}
    
    num_rows, num_cols = size(A)

    #### Initialize model
    # model = Model(Gurobi.Optimizer)
    model = Model(optimizer_with_attributes(Gurobi.Optimizer, "MIPGap"=>MIPGap, "TimeLimit"=>TimeLimit, "LogToConsole"=>0))

    @variable(model, x[1:num_rows] >= 0)   # Don't need upper bound since sum(x) == 1 will enforce it
    @variable(model, r[1:num_rows], Bin)   # 1 means row is available
    @variable(model, s[1:num_cols], Bin)   # 1 means column has been removed
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

    # if termination_status(model) == INFEASIBLE_OR_UNBOUNDED || termination_status(model) == INFEASIBLE
    #     error("Model is infeasible.")
    # elseif termination_status(model) == TIME_LIMIT && !has_values(model)
    #     error("No primal solution obtained within the time limit.")
    # end

    if termination_status(model) == INFEASIBLE_OR_UNBOUNDED || termination_status(model) == INFEASIBLE
        x, r, s = zeros(num_rows), zeros(num_rows), zeros(num_cols)  # do not return vector of 1's
        obj_val, term_status, soln_time = -Inf, termination_status(model), solve_time(model)
        rel_gap, nodes = 0, 0
    else
        if relax
            rel_gap = 0
            nodes = 0
        else
            rel_gap = relative_gap(model)
            nodes = node_count(model)
        end
        x, r, s = value.(x), value.(r), value.(s)
        obj_val, term_status, soln_time = objective_value(model), termination_status(model), solve_time(model)
    end

    return x, r, s, obj_val, term_status, soln_time, rel_gap, nodes
    # return value.(x), value.(r), value.(s), objective_value(model), termination_status(model), solve_time(model), rel_gap, nodes
end

"""
Solve MGD via a Greedy algorithm which purchases the row selection or column removal that provides the greatest increase in value.
The best row selection or column removal is computed by solving an LP for each unpurchased row selection/column removal.
"""
function solve_game_greedy_MIP(A::Matrix{T}, c_r::Vector{U}, c_s::Vector{U}, B::V; TimeLimit::W=600, MIPGap::X=0.01) where {T,U,V,W,X <: Real}

    num_rows, num_cols = size(A)

    term_status = "TIME LIMIT"
    obj_val = -Inf
    x = zeros(num_rows)
    r = zeros(num_rows)
    s = zeros(num_cols)

    num_purchases = 0
    nodes = 0

    R = Int[]  # keep track of rows purchased
    S = Int[]  # keep track of columns purchased
    sequence = String[]
    obj_val_vec = Float64[]
    B_spent_vec = Float64[]
    
    start_time = time()

    while time() < start_time + TimeLimit
        r_input = map(e -> abs(e - 0) < 1E-3 ? -1 : e, r)  # r is a solution; need to solely fix the 1's (i.e., the purchases)
        s_input = map(e -> abs(e - 0) < 1E-3 ? -1 : e, s)  # similarly for s

        x_sub, r_sub, s_sub, obj_val_sub, term_status_sub, _, nodes_sub = solve_game(A, c_r, c_s, B, r_input, s_input, MIPGap=MIPGap, 
                                                                        TimeLimit=maximum([ceil(Int,TimeLimit-(time()-start_time)), 1]), k=num_purchases+1)
        nodes += nodes_sub
        
        if term_status_sub == INFEASIBLE_OR_UNBOUNDED || term_status_sub == INFEASIBLE || 
                (r_sub == ones(num_rows) && s_sub == ones(num_cols)) || abs(obj_val_sub - obj_val) < 1E-5
            term_status = "FINISHED"
            break
        end

        x = deepcopy(x_sub)
        r = deepcopy(r_sub)
        s = deepcopy(s_sub)
        obj_val = copy(obj_val_sub)
        num_purchases += 1

        row_selected = findfirst(i -> r[i] > 0.5 && !(i in R), eachindex(r))
        column_removed = findfirst(j -> s[j] > 0.5 && !(j in S), eachindex(s))
        if row_selected !== nothing
            push!(R, row_selected)
            push!(sequence, "R")
        elseif column_removed !== nothing
            push!(S, column_removed)
            push!(sequence, "C")
        end
        push!(obj_val_vec, obj_val)
        push!(B_spent_vec, c_r' * r + c_s' * s)
    end

    time_elapsed = time() - start_time
    return x, r, s, obj_val, term_status, time_elapsed, num_purchases, nodes, R, S, sequence, obj_val_vec, B_spent_vec
end

function solve_game_greedy_frequency(A::Matrix{T}, c_r::Vector{U}, c_s::Vector{U}, B::V; TimeLimit::W=Inf, MIPGap::X=0.0001,
    approach::String="simple") where {T,U,V,W,X <: Real}

    if approach == "dual"
        x, r, s, obj_val, dual_var_s, _, _  = solve_game_LP(A, c_r, c_s, B, TimeLimit=TimeLimit, MIPGap=MIPGap)

        # Feed in x instead of r in case x < r
        # purchases = determine_greedy_purchases(x, s, dual_var_s, c_r, c_s, B)  # dual variables of s rankings
    else # "simple"
        x, r, s, obj_val, _, _, term_status, _, _, _ = solve_game(A, c_r, c_s, B, relax=true, TimeLimit=TimeLimit, MIPGap=MIPGap)
        # B_used = r' * c_r + s' * c_s

        # Feed in x instead of r in case x < r
        purchases = determine_greedy_purchases(x, s, c_r, c_s, B) # simple rankings
    end

    # Construct r_fix and s_fix
    R = @pipe filter(tup -> tup[1] == "R", purchases) |> map(tup -> tup[2], _)
    S = @pipe filter(tup -> tup[1] == "S", purchases) |> map(tup -> tup[2], _)
    r_fix = map(i -> i in R ? 1 : 0, 1:length(c_r))
    s_fix = map(j -> j in S ? 1 : 0, 1:length(c_s))

    x, r, s, obj_val, term_status, soln_time, _, _ = solve_game(A, c_r, c_s, B, r_fix, s_fix, TimeLimit=TimeLimit, MIPGap=MIPGap)

    return x, r, s, obj_val
end

function determine_greedy_purchases(r::Vector{T}, s::Vector{T}, c_r::Vector{U}, c_s::Vector{U}, B::V) where {T,U,V <: Real}
    # Ideally r and x from the LP sol'n are the same value. Feed in x for r when calling the function

    # Remove rows with no support and sort remaining
    pos_r = filter(i -> r[i] > 1E-6, eachindex(r))
    pos_s = filter(i -> s[i] > 1E-6, eachindex(s))
    sort!(pos_r, by = t-> r[t], rev=true)
    sort!(pos_s, by = t-> s[t], rev=true)

    purchases = Tuple{String,Int}[]
    spent = 0
    # Need to purchase at least one row. Purchase affordable row with largest value first
    ind = findfirst(i -> r[i] <= B, pos_r)
    if ind !== nothing
        push!(purchases, ("R", pos_r[ind]))
        spent += c_r[pos_r[ind]]
        deleteat!(pos_r, ind)
    else
        # Purchase cheapest row if can't afford any row with support
        val, ind_c_r = findmin(c_r)
        push!(purchases, ("R", ind_c_r))
        spent += c_r[ind_c_r]
    end

    # Sort row and column actions
    candidates = [ [("R", i) for i=pos_r] ; [("S", j) for j=pos_s] ]
    sort!(candidates, by=tup -> tup[2], rev=true)

    for candidate in candidates
        price = 0
        if candidate[1] == "R"
            price = c_r[candidate[2]]
        else
            price = c_s[candidate[2]]
        end 
        if spent + price > B
            break
        end
        push!(purchases, candidate)
        spent += price
    end

    return purchases
end

# function determine_greedy_purchases(r::Vector{T}, s::Vector{T}, dual_var_s::Vector{T}, c_r::Vector{U}, c_s::Vector{U}, B::V) where {T,U,V <: Real}

# end

"""
Purpose of this function is to solve the MGD LP relaxation and also return the dual variables of the z inequalites.
"""
function solve_game_LP(A::Matrix{T}, c_r::Vector{U}, c_s::Vector{U}, B::V; TimeLimit::W=Inf, MIPGap::X=0.0001) where {T,U,V,W,X <: Real}

    num_rows, num_cols = size(A)

    #### Initialize model
    # model = Model(Gurobi.Optimizer)
    model = Model(optimizer_with_attributes(Gurobi.Optimizer, "MIPGap"=>MIPGap, "TimeLimit"=>TimeLimit))

    @variable(model, x[1:num_rows] >= 0)   # Don't need upper bound since sum(x) == 1 will enforce it
    @variable(model, r[1:num_rows], Bin)   # 1 means row is available
    @variable(model, s[1:num_cols], Bin)   # 1 means column has been removed
    @variable(model , z)

    @objective(model, Max, z)

    M = compute_big_M_parameters(A)
    cons_z = @constraint(model, [j = 1:num_cols], z - A[:,j]' * x <= M[j] * s[j] )
    @constraint(model, sum(c_r' * r) + sum(c_s' * s) <= B)
    @constraint(model, sum(x) == 1)
    @constraint(model, x .<= r)

    #### Solve LP
    relax_integrality(model)
    optimize!(model)

    if termination_status(model) == INFEASIBLE_OR_UNBOUNDED || termination_status(model) == INFEASIBLE
        error("Model is infeasible.")
    elseif termination_status(model) == TIME_LIMIT && !has_values(model)
        error("No primal solution obtained within the time limit.")
    elseif !has_duals(model)
        error("No dual solutions available.")
    end

    return value.(x), value.(r), value.(s), objective_value(model), dual.(cons_z), termination_status(model), solve_time(model)
end

"""
Solve MGD with additional inequalities to strengthen relaxation.
"""
function solve_game_with_cuts(A::Matrix{T}, c_r::Vector{U}, c_s::Vector{U}, B::W, cond_dom_data::Vector{NTuple{2,Int}}; relax::Bool=false, TimeLimit::X=Inf, MIPGap::Y=0.0001) where {T,U,W,X,Y <: Real}
    
    num_rows, num_cols = size(A)

    #### Initialize model
    # model = Model(Gurobi.Optimizer)
    model = Model(optimizer_with_attributes(Gurobi.Optimizer, "MIPGap"=>MIPGap, "TimeLimit"=>TimeLimit))

    @variable(model, x[1:num_rows] >= 0)   # Don't need upper bound since sum(x) == 1 will enforce it
    @variable(model, r[1:num_rows], Bin)   # 1 means row is available
    @variable(model, s[1:num_cols], Bin)   # 1 means column has been removed
    @variable(model , z)

    @objective(model, Max, z)

    M = compute_big_M_parameters(A)
    @constraint(model, [j = 1:num_cols], z - A[:,j]' * x <= M[j] * s[j] )
    @constraint(model, sum(c_r' * r) + sum(c_s' * s) <= B)
    @constraint(model, sum(x) == 1)
    @constraint(model, x .<= r)

    # Add conditional dominance cuts 
    # cond_dom_rows, cond_dom_cols = compute_conditionally_dominated_rows(A, c_r, c_s)  # entries coincide
    # @constraint(model, [k in eachindex(cond_dom_rows)], x[cond_dom_rows[k]] + s[cond_dom_cols[k]] <= 1)
    @constraint(model, [k in eachindex(cond_dom_data)], x[cond_dom_data[k][1]] + s[cond_dom_data[k][2]] <= 1)

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

    # println("Cuts added: $(length(cond_dom_rows))")
    return value.(x), value.(r), value.(s), objective_value(model), objective_bound(model), dual_objective_value(model), termination_status(model), solve_time(model), rel_gap, nodes #, result_count(model)
end

# function compute_conditionally_dominated_rows(A::Matrix{T}, c_r::Vector{U}, c_s::Vector{U}) where {T,U <: Real}
#     # This is a big and complex combinatorial problem
    
#     cond_dom_rows = Int[]


#     for i = axes(A,1)
#         for j = axes(A,2)

#         end
#     end
# end

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

    # if termination_status(model) == INFEASIBLE_OR_UNBOUNDED || termination_status(model) == INFEASIBLE
    #     return zeros(num_rows), r_fix, -Inf, termination_status(model), maximum([0.1, solve_time(model)]), 0, 0
    #     # error("Model is infeasible.")
    # elseif termination_status(model) == TIME_LIMIT && !has_values(model)
    #     return zeros(num_rows), r_fix, -Inf, termination_status(model), TimeLimit, 0, 0
    #     # error("No primal solution obtained within the time limit.")
    # end

    # if relax
    #     rel_gap = 0
    #     nodes = 0
    # else
    #     rel_gap = relative_gap(model)
    #     nodes = node_count(model)
    # end

    if termination_status(model) == INFEASIBLE_OR_UNBOUNDED || termination_status(model) == INFEASIBLE
        x = [1.0; zeros(num_rows-1)]
        r = zeros(num_rows)  # do not return vector of 1's
        obj_val, term_status, soln_time = -Inf, termination_status(model), solve_time(model)
        rel_gap, nodes = 0, 0
    else
        if relax
            rel_gap = 0
            nodes = 0
        else
            rel_gap = relative_gap(model)
            nodes = node_count(model)
        end
        x, r = value.(x), value.(r)
        obj_val, term_status, soln_time = objective_value(model), termination_status(model), solve_time(model)
    end

    return x, r, obj_val, term_status, soln_time, rel_gap, nodes
    # return value.(x), value.(r), objective_value(model), termination_status(model), solve_time(model), rel_gap, nodes
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

    R = Int[]  # keep track of rows purchased
    obj_val_vec = Float64[]
    B_spent_vec = Float64[]
    
    start_time = time()

    while time() < start_time + TimeLimit
        r_input = map(e -> abs(e - 0) < 1E-3 ? -1 : e, r)  # r is a solution; need to solely fix the 1's (i.e., the purchases)

        x_sub, r_sub, obj_val_sub, term_status_sub, _, _, nodes_sub = solve_game(A, c, B, r_input, MIPGap=MIPGap, TimeLimit=maximum([ceil(Int,TimeLimit-(time()-start_time)), 1]),
                                                                        k=num_purchases+1)
        nodes += nodes_sub
        
        if term_status_sub == INFEASIBLE_OR_UNBOUNDED || term_status_sub == INFEASIBLE || 
                r_sub == ones(num_rows) || abs(obj_val_sub - obj_val) < 1E-5
            term_status = "FINISHED"
            break
        end

        x = deepcopy(x_sub)
        r = deepcopy(r_sub)
        obj_val = copy(obj_val_sub)
        num_purchases += 1

        row_selected = findfirst(i -> r[i] > 0.5 && !(i in R), eachindex(r))
        if row_selected !== nothing
            push!(R, row_selected)
        end
        push!(obj_val_vec, obj_val)
        push!(B_spent_vec, c' * r)
    end

    time_elapsed = time() - start_time
    return x, r, obj_val, term_status, time_elapsed, num_purchases, nodes, R, obj_val_vec, B_spent_vec
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

    R = Int[]  # keep track of rows purchased
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

            row_selected = findfirst(i -> r[i] > 0.5 && !(i in R), eachindex(r))
            if row_selected !== nothing
                push!(R, row_selected)
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

        # row_selected = findfirst(i -> r[i] > 0.5 && !(i in R), eachindex(r))
        # if row_selected !== nothing
        #     push!(R, row_selected)
        # end
        # push!(obj_val_vec, obj_val)
        # push!(B_spent_vec, c' * r)
    end

    time_elapsed = time() - start_time
    return x, r, obj_val, term_status, time_elapsed, num_purchases, R, obj_val_vec, B_spent_vec
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

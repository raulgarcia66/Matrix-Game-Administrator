using Gurobi
using JuMP
import Random

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


function solve_game(A::Matrix{T}, c::Vector{U}, B::V, r_fix::Vector{W}; relax::Bool=false, TimeLimit::X=300, MIPGap::Y=0.01, greedy_naive::Bool=false, k::Int=-1) where {T,U,V,W,X,Y <: Real}
    
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
    if greedy_naive
        @constraint(model, sum(r) == k)
    end

    #### Fix values of r (1 or 0 will be fixed, all other values in r_fix are ignored)
    if r_fix != -ones(Int, length(c))
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


function solve_game_naive(A::Matrix{T}, c::Vector{U}, B::V; seed::Int=-1, TimeLimit::W=300) where {T,U,V,W <: Real}
    # Choose subsets of rows and solve the corresponding LP for TimeLimit minutes
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


function solve_game_naive_compare(A::Matrix{T}, c::Vector{U}, B::V, opt_val::W; seed::Int=-1, TimeLimit::X=300) where {T,U,V,W,X <: Real}
    #= Choose subsets of rows and solve the corresponding LP for TimeLimit seconds =#

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


function compute_subset(num_rows)
    # This approach returns a vector with about μ many ones
    μ = rand()
    r = rand(num_rows)
    rows = filter(ind -> r[ind] > μ, eachindex(r))
    if isempty(rows)
        rows = [rand(1:num_rows)]
    end
    return rows
end


function compute_r_fix(num_rows)
    # # This will on average return a vector with half ones half zeros
    # return rand([0,1], num_rows)

    # # This approach focuses on computing subsets of rows
    # rows = compute_subset()
    # r = [Int(row in rows) for row in 1:num_rows]
    # return r

    # This approach returns a vector with about μ many ones
    μ = rand()
    r_frac = rand(num_rows)
    r = [Int(e > μ) for e in r_frac]
    if sum(r) == 0
        r[rand(1:num_rows)] = 1
    end
    return r
end


# function powerset(x::Vector{T}) where T
#     result = Vector{T}[[]]
#     for elem in x, j in eachindex(result)
#         push!(result, [result[j] ; elem])
#     end
#     return result[2:end]
# end


function solve_game_greedy_naive(A::Matrix{T}, c::Vector{U}, B::V; TimeLimit::W=300, MIPGap::X=0.01) where {T,U,V,W,X <: Real}

    num_rows = size(A, 1)

    term_status = "TIME LIMIT"
    obj_val = -Inf
    x = zeros(num_rows)
    r = zeros(num_rows)
    # r_prev = -ones(num_rows)   # would need to keep track of all r's, not just the previous
    num_purchases = 0
    nodes = 0

    r_vec = []
    obj_val_vec = []
    B_spent_vec = []
    
    start_time = time()

    while time() < start_time + TimeLimit
        r_input = map(e -> abs(e - 0) < 1E-3 ? -1 : e, r)  # only fix the 1's

        x_sub, r_sub, obj_val_sub, term_status_sub, _, _, nodes_sub = solve_game(A, c, B, r_input, MIPGap=MIPGap, TimeLimit=maximum([ceil(Int,TimeLimit-(time()-start_time)), 1]),
                                                                        greedy_naive=true, k=num_purchases+1)
        nodes += nodes_sub
        
        if term_status_sub == INFEASIBLE_OR_UNBOUNDED || term_status_sub == INFEASIBLE || r_sub == ones(num_rows)
            term_status = "FINISHED"
            break
        end

        x = deepcopy(x_sub)
        r = deepcopy(r_sub)
        obj_val = copy(obj_val_sub)
        num_purchases += 1

        push!(r_vec, r)
        push!(obj_val_vec, obj_val)
        push!(B_spent_vec, c' * r)
    end

    time_elapsed = time() - start_time
    return x, r, obj_val, term_status, time_elapsed, num_purchases, nodes, r_vec, obj_val_vec, B_spent_vec
end


function solve_game_greedy_LP(A::Matrix{T}, c::Vector{U}, B::V; TimeLimit::W=300, MIPGap::X=0.01) where {T,U,V,W,X <: Real}

    num_rows = size(A, 1)

    term_status = "TIME LIMIT"
    obj_val = -Inf
    x = zeros(num_rows)
    r = zeros(num_rows)
    indices_remaining = collect(1:num_rows)
    num_purchases = 0

    r_vec = []
    obj_val_vec = []
    B_spent_vec = []

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
                best_ind = copy(i) # indices_remaining[i]
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

            push!(r_vec, r)
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

        # push!(r_vec, r)
        # push!(obj_val_vec, obj_val)
        # push!(B_spent_vec, c' * r)
    end

    time_elapsed = time() - start_time
    return x, r, obj_val, term_status, time_elapsed, num_purchases, r_vec, obj_val_vec, B_spent_vec
end

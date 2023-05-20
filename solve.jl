using Gurobi
using JuMP
import Random

function solve_game(A::Matrix{T}, c::Vector{U}, B::V; relax::Bool=false, TimeLimit::Int=300, MIPGap::W=0.05) where {T,U,V,W <: Real}
    
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

    # return value.(x), value.(r), objective_value(model), termination_status(model), solve_time(model), relative_gap(model), node_count(model)
    return value.(x), value.(r), objective_value(model), termination_status(model), solve_time(model), rel_gap, nodes
end


function solve_game(A::Matrix{T}, c::Vector{U}, B::V, r_fix::Vector{W}; relax::Bool=false, TimeLimit::Int=300, MIPGap::X=0.05) where {T,U,V,W,X <: Real}
    
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

    # return value.(x), value.(r), objective_value(model), termination_status(model), solve_time(model), relative_gap(model), node_count(model)
    return value.(x), value.(r), objective_value(model), termination_status(model), solve_time(model), rel_gap, nodes
end


function solve_game_naive(A::Matrix{T}, c::Vector{U}, B::V; seed::Int=-1, TimeLimit::Int=300) where {T,U,V <: Real}
    # Choose subsets of rows and solve the corresponding LP for TimeLimit minutes
    num_rows = size(A, 1)
    if seed != -1
        Random.seed!(seed)
    end

    time_elapsed = 0
    max_obj = -Inf
    soln_attempts = 0
    x = -ones(num_rows)
    r = -ones(num_rows)
    # r_prev = -ones(num_rows)   # would need to keep track of all r's, not just the previous
    # TODO: Consider a timer here
    
    while time_elapsed < TimeLimit
        soln_attempts += 1
        r_fix = compute_r_fix(num_rows)
        if sum(r_fix' * c) > B
            continue
        end

        # x_sub, _, obj_val, term_status_sub, soln_time_sub, _, _ = solve_game(A, c, B, r_fix, TimeLimit=minumum([TimeLimit, TimeLimit-time_elapsed]), relax=true)
        x_sub, _, obj_val, term_status_sub, soln_time_sub, _, _ = solve_game(A, c, B, r_fix, TimeLimit=floor(Int,TimeLimit-time_elapsed)+1, relax=true)
        time_elapsed += soln_time_sub

        if term_status_sub != INFEASIBLE_OR_UNBOUNDED && term_status_sub != INFEASIBLE
            println("Current max: $max_obj\nCurrent obj: $obj_val\n")
            if obj_val > max_obj
                # println("Setting max_obj to obj_val")
                max_obj = obj_val
                x = x_sub
                r = r_fix
            end
        end

        # THE FOLLOWING ISN'T VALID BECAUSE THE UPPER BOUND ISN'T COMING FROM THE DUAL
        # EVERY SOLUTION IS FEASIBLE, HENCE A LOWER BOUND
        # if term_status_sub != INFEASIBLE_OR_UNBOUNDED && term_status_sub != INFEASIBLE
        #     println("Current max: $max_obj\nCurrent min: $min_obj\nCurrent obj: $obj_val")
        #     # max_obj = minimum([max_obj, obj_val])
        #     # min_obj = maximum([min_obj, obj_val])
        #     if obj_val < max_obj && obj_val > min_obj
        #         println("Setting max_obj to obj_val")
        #         max_obj = obj_val
        #         x = x_sub
        #         r = r_fix
        #     elseif obj_val > min_obj
        #         println("Setting min_obj to obj_val\n")
        #         min_obj = obj_val
        #         x = x_sub
        #         r = r_fix
        #     end
        # end

        # if max_obj != Inf && min_obj != -Inf
        #     # What if max_obj is 0?
        #     rel_gap = minimum([abs(max_obj - min_obj) / abs(max_obj), rel_gap])
        # end

        # if rel_gap < tol
        #     term_status = "OPTIMAL"
        #     break
        # end
    end

    # return x, r, max_obj, minimum([time_elapsed, TimeLimit]), soln_attempts
    return x, r, max_obj, time_elapsed, soln_attempts
end


function solve_game_naive_test(A::Matrix{T}, c::Vector{U}, B::V, opt_val::W; seed::Int=-1, TimeLimit::Int=300) where {T,U,V,W <: Real}
    # Choose subsets of rows and solve the corresponding LP for TimeLimit minutes
    num_rows = size(A, 1)
    if seed != -1
        Random.seed!(seed)
    end

    time_elapsed = 0
    max_obj = -Inf
    soln_attempts = 0
    term_status = "FEASIBLE"
    x = -ones(num_rows)
    r = -ones(num_rows)
    gap = Inf
    # r_prev = -ones(num_rows)   # would need to keep track of all r's, not just the previous
    # TODO: Consider a timer here
    
    while time_elapsed < TimeLimit
        soln_attempts += 1
        r_fix = compute_r_fix(num_rows)
        if sum(r_fix' * c) > B
            continue
        end

        # x_sub, _, obj_val, term_status_sub, soln_time_sub, _, _ = solve_game(A, c, B, r_fix, TimeLimit=minumum([TimeLimit, TimeLimit-time_elapsed]), relax=true)
        x_sub, _, obj_val, term_status_sub, soln_time_sub, _, _ = solve_game(A, c, B, r_fix, TimeLimit=floor(Int,TimeLimit-time_elapsed)+1, relax=true)
        time_elapsed += soln_time_sub

        if term_status_sub != INFEASIBLE_OR_UNBOUNDED && term_status_sub != INFEASIBLE
            println("Current max: $max_obj\nCurrent obj: $obj_val\n")
            if obj_val > max_obj
                # println("Setting max_obj to obj_val")
                max_obj = obj_val
                x = x_sub
                r = r_fix
            end

            gap = minimum([gap, abs(obj_val - opt_val)])
            if gap < 1E-3
                term_status = "OPTIMAL"
                break
            end
        end
    end

    # return x, r, max_obj, minimum([time_elapsed, TimeLimit]), soln_attempts
    return x, r, max_obj, term_status, time_elapsed, gap, soln_attempts
end


function compute_subset(num_rows)
    μ = rand()
    r = rand([0,1], num_rows)
    rows = filter(ind -> r[ind] > μ, eachindex(r))
    if isempty(rows)
        rows = [rand(1:num_rows)]
    end
    return rows
end


function compute_r_fix(num_rows)
    # # This will on average return a vector with half ones half zeros
    # return rand([0,1], num_rows)

    # rows = compute_subset()
    # r = [Int(row in rows) for row in 1:num_rows]
    # return r

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

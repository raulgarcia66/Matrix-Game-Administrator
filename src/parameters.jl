using LinearAlgebra
using Random

function create_matrix(entry_range, num_rows::Int, num_cols::Int; seed::Int=-1)
    if seed != -1
        Random.seed!(seed)
    end
    A = rand(entry_range, num_rows, num_cols)
    return A
end

function create_matrix(entry_range, R::Vector{T}, num_rows::Int, num_cols::Int; seed::Int=-1) where T
    if length(R) != num_cols
        error("My brother, try again.")
    end
    A = create_matrix(entry_range, num_rows, num_cols, seed=seed)
    foreach(j -> A[:,j] = R[j] * A[:,j], eachindex(R))
    return A
end

function create_matrix(entry_range, num_rows_cols::Int; seed::Int=-1)
    if seed != -1
        Random.seed!(seed)
    end
    A = rand(entry_range, num_rows_cols, num_rows_cols)
    return A
end

function create_matrix(entry_range, R::Vector{T}, num_rows_cols::Int; seed::Int=-1) where T
    if length(R) != num_rows_cols
        error("My brother, try again.")
    end
    A = create_matrix(entry_range, num_rows_cols, seed=seed)
    foreach(j -> A[:,j] = R[j] * A[:,j], eachindex(R))
    return A
end

function create_matrix_symmetric(entry_range, num_rows::Int; seed::Int=-1)
    if seed != -1
        Random.seed!(seed)
    end
    A = rand(entry_range, num_rows, num_rows)

    for i = 1:num_rows
        # Make symmetric
        for j = i+1:num_rows
            A[i,j] = -A[j,i]
        end
        # Diagonal entry of 0 to ensure symmetry
        A[i,i] = 0
    end

    return A
end

function create_cost_vector(entry_range, num_entries::Int; seed::Int=-1)
    if seed != -1
        Random.seed!(seed)
    end
    return rand(entry_range, num_entries)
end

function create_attack_scale_vector(entry_range, num_entries::Int; seed::Int=-1)
    # if seed != -1
    #     Random.seed!(seed)
    # end
    # return rand(entry_range, num_entries)
    return create_cost_vector(entry_range, num_entries, seed=seed)
end

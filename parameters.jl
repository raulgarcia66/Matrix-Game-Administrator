using LinearAlgebra

function create_matrix(entry_range, num_rows::Int, num_cols::Int)
    A = rand(entry_range, num_rows, num_cols)

    ### Anything else? Make value at least a certain amount?

    return A
end

function create_matrix_symmetric(entry_range, num_rows::Int)
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
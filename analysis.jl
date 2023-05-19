using DataFrames
using CSV 
using Pipe
using LinearAlgebra
using Statistics

# set_num = 3
# subpath = "./Experiments/Set $set_num/"
dfs = []

num_rows_vec = [10, 100, 1000]
num_cols_vec = [10, 100, 1000]

for (i,num_rows) in enumerate(num_rows_vec)
    subpath = "./Experiments/Set $i/"

    for num_cols in num_cols_vec
        filename = subpath * "Matrices $num_rows by $num_cols.txt"
        df = CSV.File(filename,
                        delim='\t',
                        ignorerepeated=true,
                        header = 3, # 3 for on line 3 or Vector of the names as strings or symbols
                        skipto = 4,
                        ) |> DataFrame

        df[:,"Num rows"] .= num_rows
        df[:, "Num cols"] .= num_cols

        select!(df, ["Num rows"; "Num cols"; names(df)[3:end-2]; names(df)[1:2]])
        push!(dfs, df)
    end
end

i = 9
fdf = dfs[i][dfs[i][!,"Termination status"] .== "OPTIMAL",:]
using Statistics
combine(fdf, "Solve time" => mean, "Solve time" => std)

cdf_vec = []
for df in dfs
    fdf = df[df[!,"Termination status"] .== "OPTIMAL",:]
    cdf = combine(fdf, "Solve time" => mean, "Solve time" => std, "Solve time" => length)
    cdf[:, "Num rows"] .= df[1, "Num rows"]
    cdf[:, "Num cols"] .= df[1, "Num cols"]
    push!(cdf_vec, cdf)
end

cdfs = vcat(cdf_vec...)

# Choose subsets of rows and solve the corresponding LP for 5 minutes
# Compute the relative gap and compare with that of MIP


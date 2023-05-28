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

##################################################################################
##################################################################################
dfs_naive = []
dfs_MIP = []

num_rows_vec = [10, 100, 1000]
num_cols_vec = [10, 100, 1000]

for (i,num_rows) in enumerate(num_rows_vec)
    subpath = "./Experiments Naive vs MIP/Set $i/"

    for num_cols in num_cols_vec
        # Naive files
        filename_naive = subpath * "Matrices $num_rows by $num_cols naive.txt"        
        df_naive = CSV.File(filename_naive,
                        delim='\t',
                        ignorerepeated=true,
                        header = 3, # 3 for on line 3 or Vector of the names as strings or symbols
                        skipto = 4,
                        ) |> DataFrame

        df_naive[:,"Num rows"] .= num_rows
        df_naive[:, "Num cols"] .= num_cols

        select!(df_naive, ["Num rows"; "Num cols"; names(df_naive)[3:end-2]; names(df_naive)[1:2]])  # Move Num rows and Num cols to front, Matrix seed and Costs seed to end
        push!(dfs_naive, df_naive)

        # MIP files
        filename_MIP = subpath * "Matrices $num_rows by $num_cols MIP.txt"
        df_MIP = CSV.File(filename_MIP,
                        delim='\t',
                        ignorerepeated=true,
                        header = 3, # 3 for on line 3 or Vector of the names as strings or symbols
                        skipto = 4,
                        ) |> DataFrame

        df_MIP[:,"Num rows"] .= num_rows
        df_MIP[:, "Num cols"] .= num_cols

        select!(df_MIP, ["Num rows"; "Num cols"; names(df_MIP)[3:end-2]; names(df_MIP)[1:2]])  # Move Num rows and Num cols to front, Matrix seed and Costs seed to end
        push!(dfs_MIP, df_MIP)
    end
end

# for i in 1:9
#     println("$i")
#     println("$(dfs_naive[i])")
#     println("$(dfs_MIP[i])")
# end

df_naive = deepcopy(dfs_naive[6])
df_MIP = deepcopy(dfs_MIP[6])
df_MIP[df_MIP[:,"Termination status"] .== "OPTIMAL", :]
df_naive[df_MIP[:,"Termination status"] .== "OPTIMAL", :]

df_naive[:,"Gap with MIP dual obj"] = abs.(df_naive[:,"Best obj val"] - df_MIP[:, "Dual obj"]) ./ abs.(df_naive[:, "Best obj val"])
df_naive[:,"Gap with MIP obj"] = abs.(df_naive[:,"Best obj val"] - df_MIP[:, "Obj val"]) ./ abs.(df_MIP[:, "Obj val"])
df_naive[:, "Time excess"] = df_naive[!,"Time achieved"] - df_MIP[!,"Solve time"]
df_naive
describe(df_naive)

df_naive_agg = combine(df_naive, "Gap with MIP dual obj" => mean, "Gap with MIP dual obj" => std, "Gap with MIP dual obj" => minimum, "Gap with MIP dual obj" => maximum,
        "Gap with MIP obj" => mean, "Gap with MIP obj" => std, "Gap with MIP obj" => minimum, "Gap with MIP obj" => maximum,
        "Time excess" => mean, "Time excess" => std, "Time excess" => minimum, "Time excess" => maximum
        )
df_naive_agg[:, "Num rows"] .= 1000
df_naive_agg[:, "Num cols"] .= 1000
select!(df_naive_agg, "Num rows", "Num cols", :)

# df_naive_agg[:,end-3:end]
dfs_naive_agg = []

for (i,num_rows) in enumerate(num_rows_vec)
    for (j,num_cols) in enumerate(num_cols_vec)
        ind = (i-1)*length(num_cols_vec) + j
        df_naive = deepcopy(dfs_naive[ind])
        df_MIP = deepcopy(dfs_MIP[ind])
        df_MIP[df_MIP[:,"Termination status"] .== "OPTIMAL", :]
        df_naive[df_MIP[:,"Termination status"] .== "OPTIMAL", :]

        df_naive[:,"Gap with MIP dual obj"] = abs.(df_naive[:,"Best obj val"] - df_MIP[:, "Dual obj"]) ./ abs.(df_naive[:, "Best obj val"])
        df_naive[:,"Gap with MIP obj"] = abs.(df_naive[:,"Best obj val"] - df_MIP[:, "Obj val"]) ./ abs.(df_MIP[:, "Obj val"])
        df_naive[:, "Time excess"] = df_naive[!,"Time achieved"] - df_MIP[!,"Solve time"]
        df_naive
        describe(df_naive)

        df_naive_agg = combine(df_naive, "Gap with MIP dual obj" => mean, "Gap with MIP dual obj" => std, "Gap with MIP dual obj" => minimum, "Gap with MIP dual obj" => maximum,
                "Gap with MIP obj" => mean, "Gap with MIP obj" => std, "Gap with MIP obj" => minimum, "Gap with MIP obj" => maximum,
                "Time excess" => mean, "Time excess" => std, "Time excess" => minimum, "Time excess" => maximum
                )
        df_naive_agg[:, "Num rows"] .= 1000
        df_naive_agg[:, "Num cols"] .= 1000
        select!(df_naive_agg, "Num rows", "Num cols", :)
        push!(dfs_naive_agg , df_naive_agg)

    end
end

dfs_naive_agg_concat = vcat(dfs_naive_agg...)
dfs_naive_agg_concat[:,1:6]
dfs_naive_agg_concat[:,[collect(1:2);collect(7:10)]]
dfs_naive_agg_concat[:,[collect(1:2);collect(11:14)]]


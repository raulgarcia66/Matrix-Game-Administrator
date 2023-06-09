using DataFrames
using CSV 
using Pipe
using LinearAlgebra
using Statistics

##################################################################################
########################### Load experiments one method ##########################
work_dir = pwd()
set_num = 4
subpath = work_dir * "/Experiments/Set $set_num/"
dfs = []

num_rows_vec = [10, 100, 1000]
num_cols_vec = [10, 100, 1000]

for num_rows in num_rows_vec
    for num_cols in num_cols_vec
        filename = subpath * "Matrices $num_rows by $num_cols MIP.txt"
        df = CSV.File(filename,
                        delim='\t',
                        ignorerepeated=true,
                        header = 5, # 5 for on line 5 or Vector of the names as strings or symbols
                        skipto = 6,
                        ) |> DataFrame

        df[:,"Num rows"] .= num_rows
        df[:, "Num cols"] .= num_cols

        select!(df, ["Num rows"; "Num cols"; names(df)[3:end-2]; names(df)[1:2]])  # Move "Num rows" and "Num cols" to front, "Matrix seed" and "Costs seed" to end
        push!(dfs, df)
    end
end

i = 9
# fdf = dfs[i][dfs[i][!,"Termination status"] .== "OPTIMAL",:]
gdf = groupby(dfs[i], "Termination status")
cdf = combine(gdf, "Solve time" => mean, "Solve time" => std, "Solve time" => length => "Count")

cdf_vec = []
for df in dfs
    # fdf = df[df[!,"Termination status"] .== "OPTIMAL",:]
    gdf = groupby(dfs[i], "Termination status")
    cdf = combine(gdf, "Solve time" => mean, "Solve time" => std, "Solve time" => length => "Count")
    cdf[:, "Num rows"] .= df[1, "Num rows"]
    cdf[:, "Num cols"] .= df[1, "Num cols"]
    push!(cdf_vec, cdf)
end

cdfs = vcat(cdf_vec...)

##################################################################################
################################# Load all files #################################
dfs_MIP = []
# dfs_naive = []
dfs_greedy = []
dfs_greedy_LP = []

num_rows_vec = [10, 100, 1000]
num_cols_vec = [10, 100, 1000]

work_dir = pwd()
set_num = 4
subpath = work_dir * "/Experiments/Set $set_num/"

for num_rows in num_rows_vec
    for num_cols in num_cols_vec

        # MIP files
        filename_MIP = subpath * "Matrices $num_rows by $num_cols MIP.txt"
        df_MIP = CSV.File(filename_MIP,
                        delim='\t',
                        ignorerepeated=true,
                        header = 5,
                        skipto = 6,
                        ) |> DataFrame

        df_MIP[:,"Num rows"] .= num_rows
        df_MIP[:, "Num cols"] .= num_cols

        select!(df_MIP, ["Num rows"; "Num cols"; names(df_MIP)[3:end-2]; names(df_MIP)[1:2]])  # Move "Num rows" and "Num cols" to front, "Matrix seed" and "Costs seed" to end
        push!(dfs_MIP, df_MIP)

        # # Naive enumeration files
        # filename_naive = subpath * "Matrices $num_rows by $num_cols naive.txt"        
        # df_naive = CSV.File(filename_naive,
        #                 delim='\t',
        #                 ignorerepeated=true,
        #                 header = 5,
        #                 skipto = 6,
        #                 ) |> DataFrame

        # df_naive[:,"Num rows"] .= num_rows
        # df_naive[:, "Num cols"] .= num_cols

        # select!(df_naive, ["Num rows"; "Num cols"; names(df_naive)[3:end-2]; names(df_naive)[1:2]])
        # push!(dfs_naive, df_naive)

        # Greedy files
        filename_greedy = subpath * "Matrices $num_rows by $num_cols greedy.txt"
        df_greedy = CSV.File(filename_greedy,
                        delim='\t',
                        ignorerepeated=true,
                        header = 5,
                        skipto = 6,
                        ) |> DataFrame

        df_greedy[:,"Num rows"] .= num_rows
        df_greedy[:, "Num cols"] .= num_cols

        select!(df_greedy, ["Num rows"; "Num cols"; names(df_greedy)[3:end-2]; names(df_greedy)[1:2]])
        push!(dfs_greedy, df_greedy)

        # Greedy LP files
        filename_greedy_LP = subpath * "Matrices $num_rows by $num_cols greedy LP.txt"
        df_greedy_LP = CSV.File(filename_greedy_LP,
                        delim='\t',
                        ignorerepeated=true,
                        header = 5,
                        skipto = 6,
                        ) |> DataFrame

        df_greedy_LP[:,"Num rows"] .= num_rows
        df_greedy_LP[:, "Num cols"] .= num_cols

        select!(df_greedy_LP, ["Num rows"; "Num cols"; names(df_greedy_LP)[3:end-2]; names(df_greedy_LP)[1:2]])
        push!(dfs_greedy_LP, df_greedy_LP)
    end
end

for i in 1:1
    println("\n$i")
    println("$(dfs_MIP[i])")
    # println("$(dfs_naive[i])")
    println("$(dfs_greedy[i])")
    println("$(dfs_greedy_LP[i])")
end

##################################################################################
############################# Analysis individually ##############################
t = 5
df_MIP = deepcopy(dfs_MIP[t])
# df_naive = deepcopy(dfs_naive[t])
df_greedy = deepcopy(dfs_greedy[t])
df_greedy_LP = deepcopy(dfs_greedy_LP[t])

#### Remove row with full budget because its an LP
filter!(row -> row["Budget fraction"] != 1.0, df_MIP)
# filter(row -> row["Budget fraction"] != 1.0, df_naive)
filter!(row -> row["Budget fraction"] != 1.0, df_greedy)
filter!(row -> row["Budget fraction"] != 1.0, df_greedy_LP)

#### If we want to only work with instances with MIP optimality (using groups instead)
# df_MIP = df_MIP[df_MIP[:,"Termination status"] .== "OPTIMAL", :]
# df_naive = df_naive[df_naive[:,"MIP Term status"] .== "OPTIMAL", :]
# df_greedy = df_greedy[df_greedy[:,"MIP Term status"] .== "OPTIMAL", :]
# df_greedy_LP = df_greedy_LP[df_greedy_LP[:,"MIP Term status"] .== "OPTIMAL", :]

#### Add MIP term status for group comparisons
# df_naive[:,"MIP Term status"] = df_MIP[!,"Termination status"]
df_greedy[:,"MIP Term status"] = df_MIP[!,"Termination status"]
df_greedy_LP[:,"MIP Term status"] = df_MIP[!,"Termination status"]

#### Compute comparisons
# df_naive[:,"Gap"] = (df_MIP[:, "Obj val"] - df_naive[:,"Best obj val"])
# df_naive[:,"Rel gap"] = (df_MIP[:, "Obj val"] - df_naive[:,"Best obj val"]) ./ abs.(df_MIP[:, "Obj val"])
# df_naive[:, "Time excess"] = df_naive[!,"Time achieved"] - df_MIP[!,"Solve time"]
# df_naive[:, "Naive max found faster"] = df_naive[!,"Time excess"] .< 0
# # df_naive
# # describe(df_naive)

df_greedy[:,"Gap"] = (df_MIP[:, "Obj val"] - df_greedy[:,"Obj val"])
df_greedy[:,"Rel gap"] = (df_MIP[:, "Obj val"] - df_greedy[:,"Obj val"]) ./ abs.(df_MIP[:, "Obj val"])
df_greedy[:, "Time excess"] = df_greedy[!,"Solve time"] - df_MIP[!,"Solve time"]
df_greedy[:, "Greedy max found faster"] = df_greedy[!,"Time excess"] .< 0
# df_greedy
# describe(df_greedy)

df_greedy_LP[:,"Gap"] = (df_MIP[:, "Obj val"] - df_greedy_LP[:,"Obj val"])
df_greedy_LP[:,"Rel gap"] = (df_MIP[:, "Obj val"] - df_greedy_LP[:,"Obj val"]) ./ abs.(df_MIP[:, "Obj val"])
df_greedy_LP[:, "Time excess"] = df_greedy_LP[!,"Solve time"] - df_MIP[!,"Solve time"]
df_greedy_LP[:, "Greedy max found faster"] = df_greedy_LP[!,"Time excess"] .< 0
# df_greedy
# describe(df_greedy)

#### Group
gdf_MIP = groupby(df_MIP, "Termination status")
# gdf_naive = groupby(df_naive, "MIP Term status")
gdf_greedy = groupby(df_greedy, ["MIP Term status", "Termination status"])
gdf_greedy_LP = groupby(df_greedy_LP, ["MIP Term status", "Termination status"])

#### Compute aggregates
gdf_MIP_agg = combine(gdf_MIP, "Termination status" => length => "Count",
        "Relative gap" => mean, "Relative gap" => minimum, "Relative gap" => maximum, "Relative gap" => std,
        "Solve time" => mean, "Solve time" => minimum, "Solve time" => maximum, "Solve time" => std,
        )
gdf_MIP_agg[:, "Num rows"] .= num_rows_vec[2]
gdf_MIP_agg[:, "Num cols"] .= num_cols_vec[3]
select!(gdf_MIP_agg, "Num rows", "Num cols", :)

# gdf_naive_agg = combine(gdf_naive, "MIP Term status" => length => "MIP status count", "Naive max found faster" => count => "Faster count",
#         "Gap" => mean, "Gap" => minimum, "Gap" => maximum, "Gap" => std,
#         "Rel gap" => mean, "Rel gap" => minimum, "Rel gap" => maximum, "Rel gap" => std,
#         "Time excess" => mean, "Time excess" => minimum, "Time excess" => maximum, "Time excess" => std
#         )
# gdf_naive_agg[:, "Num rows"] .= num_rows_vec[2]
# gdf_naive_agg[:, "Num cols"] .= num_cols_vec[3]
# select!(gdf_naive_agg, "Num rows", "Num cols", :)

gdf_greedy_agg = combine(gdf_greedy, "MIP Term status" => length => "MIP status count", "Termination status" => length => "Greedy status count", "Greedy max found faster" => count => "Faster count",
        "Gap" => mean, "Gap" => minimum, "Gap" => maximum, "Gap" => std,
        "Rel gap" => mean, "Rel gap" => minimum, "Rel gap" => maximum, "Rel gap" => std,
        "Time excess" => mean, "Time excess" => minimum, "Time excess" => maximum, "Time excess" => std
        )
gdf_greedy_agg[:, "Num rows"] .= num_rows_vec[2]
gdf_greedy_agg[:, "Num cols"] .= num_cols_vec[3]
select!(gdf_greedy_agg, "Num rows", "Num cols", :)

gdf_greedy_LP_agg = combine(gdf_greedy_LP, "MIP Term status" => length => "MIP status count", "Termination status" => length => "Greedy LP status count", "Greedy LP max found faster" => count => "Faster count",
        "Gap" => mean, "Gap" => minimum, "Gap" => maximum, "Gap" => std,
        "Rel gap" => mean, "Rel gap" => minimum, "Rel gap" => maximum, "Rel gap" => std,
        "Time excess" => mean, "Time excess" => minimum, "Time excess" => maximum, "Time excess" => std
        )
gdf_greedy_LP_agg[:, "Num rows"] .= num_rows_vec[2]
gdf_greedy_LP_agg[:, "Num cols"] .= num_cols_vec[3]
select!(gdf_greedy_LP_agg, "Num rows", "Num cols", :)

# Observe
gdf_MIP_agg[:,3:7]
gdf_MIP_agg[:,8:end]
gdf_naive_agg[:,3:11]
gdf_naive_agg[:,12:end]
gdf_greedy_agg[:,3:11]
gdf_greedy_agg[:,12:end]
gdf_greedy_LP_agg[:,3:11]
gdf_greedy_LP_agg[:,12:end]

##################################################################################
############################## Analysis all files ################################
gdfs_MIP_agg = []
# gdfs_naive_agg = []
gdfs_greedy_agg = []
gdfs_greedy_LP_agg = []

for (i,num_rows) in enumerate(num_rows_vec)
    for (j,num_cols) in enumerate(num_cols_vec)

        ind = (i-1)*length(num_cols_vec) + j

        df_MIP = deepcopy(dfs_MIP[ind])
        # df_naive = deepcopy(dfs_naive[ind])
        df_greedy = deepcopy(dfs_greedy[ind])
        df_greedy_LP = deepcopy(dfs_greedy_LP[ind])

        #### Remove row with full budget because its an LP
        filter!(row -> row["Budget fraction"] != 1.0, df_MIP)
        # filter(row -> row["Budget fraction"] != 1.0, df_naive)
        filter!(row -> row["Budget fraction"] != 1.0, df_greedy)
        filter!(row -> row["Budget fraction"] != 1.0, df_greedy_LP)

        #### If we want to only work with instances with MIP optimality (using groups instead)
        # df_MIP = df_MIP[df_MIP[:,"Termination status"] .== "OPTIMAL", :]
        # df_naive = df_naive[df_naive[:,"MIP Term status"] .== "OPTIMAL", :]
        # df_greedy = df_greedy[df_greedy[:,"MIP Term status"] .== "OPTIMAL", :]
        # df_greedy_LP = df_greedy_LP[df_greedy_LP[:,"MIP Term status"] .== "OPTIMAL", :]

        #### Add MIP term status for group comparisons
        # df_naive[:,"MIP Term status"] = df_MIP[!,"Termination status"]
        df_greedy[:,"MIP Term status"] = df_MIP[!,"Termination status"]
        df_greedy_LP[:,"MIP Term status"] = df_MIP[!,"Termination status"]

        #### Compute comparisons
        # df_naive[:,"Gap"] = (df_MIP[:, "Obj val"] - df_naive[:,"Best obj val"])
        # df_naive[:,"Rel gap"] = (df_MIP[:, "Obj val"] - df_naive[:,"Best obj val"]) ./ abs.(df_MIP[:, "Obj val"])
        # df_naive[:, "Time excess"] = df_naive[!,"Time achieved"] - df_MIP[!,"Solve time"]
        # df_naive[:, "Naive max found faster"] = df_naive[!,"Time excess"] .< 0

        df_greedy[:,"Gap"] = (df_MIP[:, "Obj val"] - df_greedy[:,"Obj val"])
        df_greedy[:,"Rel gap"] = (df_MIP[:, "Obj val"] - df_greedy[:,"Obj val"]) ./ abs.(df_MIP[:, "Obj val"])
        df_greedy[:, "Time excess"] = df_greedy[!,"Solve time"] - df_MIP[!,"Solve time"]
        df_greedy[:, "Greedy max found faster"] = df_greedy[!,"Time excess"] .< 0

        df_greedy_LP[:,"Gap"] = (df_MIP[:, "Obj val"] - df_greedy_LP[:,"Obj val"])
        df_greedy_LP[:,"Rel gap"] = (df_MIP[:, "Obj val"] - df_greedy_LP[:,"Obj val"]) ./ abs.(df_MIP[:, "Obj val"])
        df_greedy_LP[:, "Time excess"] = df_greedy_LP[!,"Solve time"] - df_MIP[!,"Solve time"]
        df_greedy_LP[:, "Greedy LP max found faster"] = df_greedy_LP[!,"Time excess"] .< 0

        #### Group
        gdf_MIP = groupby(df_MIP, "Termination status")
        # gdf_naive = groupby(df_naive, "MIP Term status")
        gdf_greedy = groupby(df_greedy, ["MIP Term status", "Termination status"])
        gdf_greedy_LP = groupby(df_greedy_LP, ["MIP Term status", "Termination status"])

        #### Compute aggregates
        gdf_MIP_agg = combine(gdf_MIP, "Termination status" => length => "Count",
                "Relative gap" => mean, "Relative gap" => minimum, "Relative gap" => maximum, "Relative gap" => std,
                "Solve time" => mean, "Solve time" => minimum, "Solve time" => maximum, "Solve time" => std,
                )
        gdf_MIP_agg[:, "Num rows"] .= num_rows_vec[2]
        gdf_MIP_agg[:, "Num cols"] .= num_cols_vec[3]
        select!(gdf_MIP_agg, "Num rows", "Num cols", :)

        # gdf_naive_agg = combine(gdf_naive, "MIP Term status" => length => "MIP status count", "Naive max found faster" => count => "Faster count",
        #         "Gap" => mean, "Gap" => minimum, "Gap" => maximum, "Gap" => std,
        #         "Rel gap" => mean, "Rel gap" => minimum, "Rel gap" => maximum, "Rel gap" => std,
        #         "Time excess" => mean, "Time excess" => minimum, "Time excess" => maximum, "Time excess" => std
        #         )
        # gdf_naive_agg[:, "Num rows"] .= num_rows_vec[2]
        # gdf_naive_agg[:, "Num cols"] .= num_cols_vec[3]
        # select!(gdf_naive_agg, "Num rows", "Num cols", :)

        gdf_greedy_agg = combine(gdf_greedy, "MIP Term status" => length => "MIP status count", "Termination status" => length => "Greedy status count", "Greedy max found faster" => count => "Faster count",
                "Gap" => mean, "Gap" => minimum, "Gap" => maximum, "Gap" => std,
                "Rel gap" => mean, "Rel gap" => minimum, "Rel gap" => maximum, "Rel gap" => std,
                "Time excess" => mean, "Time excess" => minimum, "Time excess" => maximum, "Time excess" => std
                )
        gdf_greedy_agg[:, "Num rows"] .= num_rows_vec[2]
        gdf_greedy_agg[:, "Num cols"] .= num_cols_vec[3]
        select!(gdf_greedy_agg, "Num rows", "Num cols", :)

        gdf_greedy_LP_agg = combine(gdf_greedy_LP, "MIP Term status" => length => "MIP status count", "Termination status" => length => "Greedy LP status count", "Greedy LP max found faster" => count => "Faster count",
                "Gap" => mean, "Gap" => minimum, "Gap" => maximum, "Gap" => std,
                "Rel gap" => mean, "Rel gap" => minimum, "Rel gap" => maximum, "Rel gap" => std,
                "Time excess" => mean, "Time excess" => minimum, "Time excess" => maximum, "Time excess" => std
                )
        gdf_greedy_LP_agg[:, "Num rows"] .= num_rows_vec[2]
        gdf_greedy_LP_agg[:, "Num cols"] .= num_cols_vec[3]
        select!(gdf_greedy_LP_agg, "Num rows", "Num cols", :)

        push!(gdfs_MIP_agg , gdf_MIP_agg)
        # push!(gdfs_naive_agg , gdf_naive_agg)
        push!(gdfs_greedy_agg , gdf_greedy_agg)
        push!(gdfs_greedy_LP_agg , gdf_greedy_LP_agg)
    end
end

t = 9
gdfs_MIP_agg[t][:,3:7]
gdfs_MIP_agg[t][:,8:end]
# gdfs_naive_agg[t][:,3:11]
# gdfs_naive_agg[t][:,12:end]
gdfs_greedy_agg[t][:,3:11]
gdfs_greedy_agg[t][:,12:end]
gdfs_greedy_LP_agg[t][:,3:11]
gdfs_greedy_LP_agg[t][:,12:end]

# # Can no longer merge
# dfs_naive_agg_concat = vcat(dfs_naive_agg...)
# dfs_naive_agg_concat[:,1:6]
# dfs_naive_agg_concat[:,[collect(1:2);collect(7:10)]]
# dfs_naive_agg_concat[:,[collect(1:2);collect(11:end)]]


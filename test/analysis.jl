using DataFrames
using CSV 
using Pipe
using LinearAlgebra
using Statistics
using Plots

##################################################################################
############################## Load files one method #############################
work_dir = pwd()
set_num = 5
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
# dfs_greedy_LP = []

num_rows_vec = [10, 100, 1000]
num_cols_vec = [10, 100, 1000]

work_dir = pwd()
set_num = 5
subpath = work_dir * "/Experiments/Set $set_num/"

load_MIP = true; load_naive = false; load_greedy = true; load_greedy_LP = false;
for num_rows in num_rows_vec
    for num_cols in num_cols_vec

        # MIP files
        if load_MIP
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
        end

        # # Naive enumeration files
        # if load_naive
        #     filename_naive = subpath * "Matrices $num_rows by $num_cols naive.txt"        
        #     df_naive = CSV.File(filename_naive,
        #                     delim='\t',
        #                     ignorerepeated=true,
        #                     header = 5,
        #                     skipto = 6,
        #                     ) |> DataFrame

        #     df_naive[:,"Num rows"] .= num_rows
        #     df_naive[:, "Num cols"] .= num_cols

        #     select!(df_naive, ["Num rows"; "Num cols"; names(df_naive)[3:end-2]; names(df_naive)[1:2]])
        #     push!(dfs_naive, df_naive)
        # end

        # Greedy files
        if load_greedy
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
        end

        # Greedy LP files
        if load_greedy_LP
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
end

for i in 1:1
    println("\n$i")
    println("$(dfs_MIP[i])")
    # println("$(dfs_naive[i])")
    println("$(dfs_greedy[i])")
    # println("$(dfs_greedy_LP[i])")
end

##################################################################################
############################# Analysis individually ##############################
t = 5
df_MIP = deepcopy(dfs_MIP[t])
# df_naive = deepcopy(dfs_naive[t])
df_greedy = deepcopy(dfs_greedy[t])
# df_greedy_LP = deepcopy(dfs_greedy_LP[t])

#### Remove row with full budget because its an LP
filter!(row -> row["Budget fraction"] != 1.0, df_MIP)
# filter(row -> row["Budget fraction"] != 1.0, df_naive)
filter!(row -> row["Budget fraction"] != 1.0, df_greedy)
# filter!(row -> row["Budget fraction"] != 1.0, df_greedy_LP)

#### If we want to only work with instances with MIP optimality (using groups instead)
# df_MIP = df_MIP[df_MIP[:,"Termination status"] .== "OPTIMAL", :]
# df_naive = df_naive[df_naive[:,"MIP Term status"] .== "OPTIMAL", :]
# df_greedy = df_greedy[df_greedy[:,"MIP Term status"] .== "OPTIMAL", :]
# df_greedy_LP = df_greedy_LP[df_greedy_LP[:,"MIP Term status"] .== "OPTIMAL", :]

#### Add MIP term status for group comparisons
# df_naive[:,"MIP Term status"] = df_MIP[!,"Termination status"]
df_greedy[:,"MIP Term status"] = df_MIP[!,"Termination status"]
# df_greedy_LP[:,"MIP Term status"] = df_MIP[!,"Termination status"]

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

# df_greedy_LP[:,"Gap"] = (df_MIP[:, "Obj val"] - df_greedy_LP[:,"Obj val"])
# df_greedy_LP[:,"Rel gap"] = (df_MIP[:, "Obj val"] - df_greedy_LP[:,"Obj val"]) ./ abs.(df_MIP[:, "Obj val"])
# df_greedy_LP[:, "Time excess"] = df_greedy_LP[!,"Solve time"] - df_MIP[!,"Solve time"]
# df_greedy_LP[:, "Greedy max found faster"] = df_greedy_LP[!,"Time excess"] .< 0
# # df_greedy
# # describe(df_greedy)

#### Group
gdf_MIP = groupby(df_MIP, ["Num rows", "Num cols", "Budget fraction", "Termination status"])
# gdf_naive = groupby(df_naive, "MIP Term status")
gdf_greedy = groupby(df_greedy, ["Num rows", "Num cols", "Budget fraction", "MIP Term status", "Termination status"])
# gdf_greedy_LP = groupby(df_greedy_LP, ["MIP Term status", "Termination status"])

#### Compute aggregates
gdf_MIP_agg = combine(gdf_MIP, "Termination status" => length => "Status count",
        "Relative gap" => mean, "Relative gap" => minimum, "Relative gap" => maximum, "Relative gap" => std,
        "Solve time" => mean, "Solve time" => minimum, "Solve time" => maximum, "Solve time" => std,
        )
# gdf_MIP_agg[:, "Num rows"] .= num_rows_vec[2]
# gdf_MIP_agg[:, "Num cols"] .= num_cols_vec[3]
# select!(gdf_MIP_agg, "Num rows", "Num cols", :) # reorder

# gdf_naive_agg = combine(gdf_naive, "MIP Term status" => length => "MIP status count", "Naive max found faster" => count => "Faster count",
#         "Gap" => mean, "Gap" => minimum, "Gap" => maximum, "Gap" => std,
#         "Rel gap" => mean, "Rel gap" => minimum, "Rel gap" => maximum, "Rel gap" => std,
#         "Time excess" => mean, "Time excess" => minimum, "Time excess" => maximum, "Time excess" => std
#         )
# # gdf_naive_agg[:, "Num rows"] .= num_rows_vec[2]
# # gdf_naive_agg[:, "Num cols"] .= num_cols_vec[3]
# # select!(gdf_naive_agg, "Num rows", "Num cols", :) # reoder

gdf_greedy_agg = combine(gdf_greedy, "MIP Term status" => length => "MIP status count", "Termination status" => length => "Greedy status count", 
        "Solve time" => mean, "Greedy max found faster" => count => "Faster count",
        "Gap" => mean, "Gap" => minimum, "Gap" => maximum, "Gap" => std,
        "Rel gap" => mean, "Rel gap" => minimum, "Rel gap" => maximum, "Rel gap" => std,
        "Time excess" => mean, "Time excess" => minimum, "Time excess" => maximum, "Time excess" => std
        )
# gdf_greedy_agg[:, "Num rows"] .= num_rows_vec[2]
# gdf_greedy_agg[:, "Num cols"] .= num_cols_vec[3]
# select!(gdf_greedy_agg, "Num rows", "Num cols", :) # reoder

# gdf_greedy_LP_agg = combine(gdf_greedy_LP, "MIP Term status" => length => "MIP status count", "Termination status" => length => "Greedy LP status count", 
#         "Solve time" => mean, "Greedy LP max found faster" => count => "Faster count",
#         "Gap" => mean, "Gap" => minimum, "Gap" => maximum, "Gap" => std,
#         "Rel gap" => mean, "Rel gap" => minimum, "Rel gap" => maximum, "Rel gap" => std,
#         "Time excess" => mean, "Time excess" => minimum, "Time excess" => maximum, "Time excess" => std
#         )
# # gdf_greedy_LP_agg[:, "Num rows"] .= num_rows_vec[2]
# # gdf_greedy_LP_agg[:, "Num cols"] .= num_cols_vec[3]
# # select!(gdf_greedy_LP_agg, "Num rows", "Num cols", :) # reoder

# Observe
gdf_MIP_agg[:,1:6]
gdf_MIP_agg[:,7:end]
# gdf_naive_agg[:,3:11]
# gdf_naive_agg[:,12:end]
gdf_greedy_agg[:,1:8]
gdf_greedy_agg[:,9:end-4]
gdf_greedy_agg[:,end-3:end]
# gdf_greedy_LP_agg[:,3:11]
# gdf_greedy_LP_agg[:,12:end]

##################################################################################
############################## Analysis all files ################################
gdfs_MIP_agg = []
# gdfs_naive_agg = []
gdfs_greedy_agg = []
# gdfs_greedy_LP_agg = []

for (i,num_rows) in enumerate(num_rows_vec)
    for (j,num_cols) in enumerate(num_cols_vec)

        ind = (i-1)*length(num_cols_vec) + j

        df_MIP = deepcopy(dfs_MIP[ind])
        # df_naive = deepcopy(dfs_naive[ind])
        df_greedy = deepcopy(dfs_greedy[ind])
        # df_greedy_LP = deepcopy(dfs_greedy_LP[ind])

        #### Remove row with full budget because its an LP
        filter!(row -> row["Budget fraction"] != 1.0, df_MIP)
        # filter(row -> row["Budget fraction"] != 1.0, df_naive)
        filter!(row -> row["Budget fraction"] != 1.0, df_greedy)
        # filter!(row -> row["Budget fraction"] != 1.0, df_greedy_LP)

        #### If we want to only work with instances with MIP optimality (using groups instead)
        # df_MIP = df_MIP[df_MIP[:,"Termination status"] .== "OPTIMAL", :]
        # df_naive = df_naive[df_naive[:,"MIP Term status"] .== "OPTIMAL", :]
        # df_greedy = df_greedy[df_greedy[:,"MIP Term status"] .== "OPTIMAL", :]
        # df_greedy_LP = df_greedy_LP[df_greedy_LP[:,"MIP Term status"] .== "OPTIMAL", :]

        #### Add MIP term status for group comparisons
        # df_naive[:,"MIP Term status"] = df_MIP[!,"Termination status"]
        df_greedy[:,"MIP Term status"] = df_MIP[!,"Termination status"]
        # df_greedy_LP[:,"MIP Term status"] = df_MIP[!,"Termination status"]

        #### Compute comparisons
        # df_naive[:,"Gap"] = (df_MIP[:, "Obj val"] - df_naive[:,"Best obj val"])
        # df_naive[:,"Rel gap"] = (df_MIP[:, "Obj val"] - df_naive[:,"Best obj val"]) ./ abs.(df_MIP[:, "Obj val"])
        # df_naive[:, "Time excess"] = df_naive[!,"Time achieved"] - df_MIP[!,"Solve time"]
        # df_naive[:, "Naive max found faster"] = df_naive[!,"Time excess"] .< 0

        df_greedy[:,"Gap"] = (df_MIP[:, "Obj val"] - df_greedy[:,"Obj val"])
        df_greedy[:,"Rel gap"] = (df_MIP[:, "Obj val"] - df_greedy[:,"Obj val"]) ./ abs.(df_MIP[:, "Obj val"])
        df_greedy[:, "Time excess"] = df_greedy[!,"Solve time"] - df_MIP[!,"Solve time"]
        df_greedy[:, "Greedy max found faster"] = df_greedy[!,"Time excess"] .< 0

        # df_greedy_LP[:,"Gap"] = (df_MIP[:, "Obj val"] - df_greedy_LP[:,"Obj val"])
        # df_greedy_LP[:,"Rel gap"] = (df_MIP[:, "Obj val"] - df_greedy_LP[:,"Obj val"]) ./ abs.(df_MIP[:, "Obj val"])
        # df_greedy_LP[:, "Time excess"] = df_greedy_LP[!,"Solve time"] - df_MIP[!,"Solve time"]
        # df_greedy_LP[:, "Greedy LP max found faster"] = df_greedy_LP[!,"Time excess"] .< 0

        #### Group
        gdf_MIP = groupby(df_MIP, ["Num rows", "Num cols", "Budget fraction", "Termination status"])
        # gdf_naive = groupby(df_naive, "MIP Term status")
        gdf_greedy = groupby(df_greedy, ["Num rows", "Num cols", "Budget fraction", "MIP Term status", "Termination status"])
        # gdf_greedy_LP = groupby(df_greedy_LP, ["MIP Term status", "Termination status"])

        #### Compute aggregates
        gdf_MIP_agg = combine(gdf_MIP, "Termination status" => length => "Status count",
                "Relative gap" => mean, "Relative gap" => minimum, "Relative gap" => maximum, "Relative gap" => std,
                "Solve time" => mean, "Solve time" => minimum, "Solve time" => maximum, "Solve time" => std,
                )
        # gdf_MIP_agg[:, "Num rows"] .= num_rows_vec[2]
        # gdf_MIP_agg[:, "Num cols"] .= num_cols_vec[3]
        # select!(gdf_MIP_agg, "Num rows", "Num cols", :) # reoder

        # gdf_naive_agg = combine(gdf_naive, "MIP Term status" => length => "MIP status count", "Naive max found faster" => count => "Faster count",
        #         "Gap" => mean, "Gap" => minimum, "Gap" => maximum, "Gap" => std,
        #         "Rel gap" => mean, "Rel gap" => minimum, "Rel gap" => maximum, "Rel gap" => std,
        #         "Time excess" => mean, "Time excess" => minimum, "Time excess" => maximum, "Time excess" => std
        #         )
        # # gdf_naive_agg[:, "Num rows"] .= num_rows_vec[2]
        # # gdf_naive_agg[:, "Num cols"] .= num_cols_vec[3]
        # # select!(gdf_naive_agg, "Num rows", "Num cols", :) # reoder

        gdf_greedy_agg = combine(gdf_greedy, "MIP Term status" => length => "MIP status count",
                "Solve time" => mean, "Termination status" => length => "Greedy status count", "Greedy max found faster" => count => "Faster count",
                "Gap" => mean, "Gap" => minimum, "Gap" => maximum, "Gap" => std,
                "Rel gap" => mean, "Rel gap" => minimum, "Rel gap" => maximum, "Rel gap" => std,
                "Time excess" => mean, "Time excess" => minimum, "Time excess" => maximum, "Time excess" => std
                )
        # gdf_greedy_agg[:, "Num rows"] .= num_rows_vec[2]
        # gdf_greedy_agg[:, "Num cols"] .= num_cols_vec[3]
        # select!(gdf_greedy_agg, "Num rows", "Num cols", :) # reorder

        # gdf_greedy_LP_agg = combine(gdf_greedy_LP, "MIP Term status" => length => "MIP status count",
        #         "Solve time" => mean, "Termination status" => length => "Greedy LP status count", "Greedy LP max found faster" => count => "Faster count",
        #         "Gap" => mean, "Gap" => minimum, "Gap" => maximum, "Gap" => std,
        #         "Rel gap" => mean, "Rel gap" => minimum, "Rel gap" => maximum, "Rel gap" => std,
        #         "Time excess" => mean, "Time excess" => minimum, "Time excess" => maximum, "Time excess" => std
        #         )
        # # gdf_greedy_LP_agg[:, "Num rows"] .= num_rows_vec[2]
        # # gdf_greedy_LP_agg[:, "Num cols"] .= num_cols_vec[3]
        # # select!(gdf_greedy_LP_agg, "Num rows", "Num cols", :) # reoder

        push!(gdfs_MIP_agg , gdf_MIP_agg)
        # push!(gdfs_naive_agg , gdf_naive_agg)
        push!(gdfs_greedy_agg , gdf_greedy_agg)
        # push!(gdfs_greedy_LP_agg , gdf_greedy_LP_agg)
    end
end

t = 4
gdfs_MIP_agg[t][:,1:6]
gdfs_MIP_agg[t][:,7:end]
# gdfs_naive_agg[t][:,3:11]
# gdfs_naive_agg[t][:,12:end]
gdfs_greedy_agg[t][:,1:8]
gdfs_greedy_agg[t][:,9:end-4]
gdfs_greedy_agg[t][:,end-4:end]
# gdfs_greedy_LP_agg[t][:,3:11]
# gdfs_greedy_LP_agg[t][:,12:end]

# # Can no longer merge
# dfs_naive_agg_concat = vcat(dfs_naive_agg...)
# dfs_naive_agg_concat[:,1:6]
# dfs_naive_agg_concat[:,[collect(1:2);collect(7:10)]]
# dfs_naive_agg_concat[:,[collect(1:2);collect(11:end)]]

##################################################################################
############################ Analysis individually 2 #############################
t = 5
df_MIP = deepcopy(dfs_MIP[t])
# df_naive = deepcopy(dfs_naive[t])
df_greedy = deepcopy(dfs_greedy[t])
# df_greedy_LP = deepcopy(dfs_greedy_LP[t])

#### Remove row with full budget because its an LP
filter!(row -> row["Budget fraction"] != 1.0, df_MIP)
# filter(row -> row["Budget fraction"] != 1.0, df_naive)
filter!(row -> row["Budget fraction"] != 1.0, df_greedy)
# filter!(row -> row["Budget fraction"] != 1.0, df_greedy_LP)

#### Merge data DataFrames
df_merged = deepcopy(df_MIP)
select!(df_merged, Not(["Obj bound", "Dual obj", "Node count"]))
rename!(df_merged, "Budget spent" => "Budget spent_MIP",
        "Obj val" => "Obj val_MIP",
        "Termination status" => "Termination status_MIP",
        "Solve time" => "Solve time_MIP",
        "Relative gap" => "Relative gap_MIP"
        )

df_merged[:,"Obj val_G"] = df_greedy[!,"Obj val"]
df_merged[:,"Termination status_G"] = df_greedy[!,"Termination status"]
df_merged[:,"Solve time_G"] = df_greedy[!,"Solve time"]
df_merged[:,"Num purchases_G"] = df_greedy[!,"Num purchases"]

# df_merged[:,"Obj val_GLP"] = df_greedy_LP[!,"Obj val"]
# df_merged[:,"Termination status_GLP"] = df_greedy_LP[!,"Termination status"]
# df_merged[:,"Solve time_GLP"] = df_greedy_LP[!,"Solve time"]
# df_merged[:,"Num purchases_GLP"] = df_greedy_LP[!,"Num purchases"]

#### Compute comparisons
df_merged[:,"Gap_G"] = (df_MIP[:, "Obj val"] - df_greedy[:,"Obj val"])
df_merged[:,"Gap percentage_G"] = (df_MIP[:, "Obj val"] - df_greedy[:,"Obj val"]) ./ abs.(df_MIP[:, "Obj val"]) * 100
df_merged[:, "Time excess_G"] = df_greedy[!,"Solve time"] - df_MIP[!,"Solve time"]
df_merged[:, "Greedy max found faster"] = df_merged[!,"Time excess_G"] .< 0
# df_merged
# describe(df_merged)

# df_merged[:,"Gap_GLP"] = (df_MIP[:, "Obj val"] - df_greedy_LP[:,"Obj val"])
# df_merged[:,"Gap percentage_GLP"] = (df_MIP[:, "Obj val"] - df_greedy_LP[:,"Obj val"]) ./ abs.(df_MIP[:, "Obj val"]) * 100
# df_merged[:, "Time excess_GLP"] = df_greedy_LP[!,"Solve time"] - df_MIP[!,"Solve time"]
# df_merged[:, "Greedy LP max found faster"] = df_merged[!,"Time excess_GLP"] .< 0
# # df_merged
# # describe(df_merged)

#### Group
gdf_merged = groupby(df_merged, ["Num rows", "Num cols", "Budget fraction", "Termination status_MIP", "Termination status_G"])

#### Compute aggregates
gdf_merged_agg = combine(gdf_merged, "Termination status_MIP" => length => "Term status_MIP_count", "Termination status_G" => length => "Term status_G_count",
        "Solve time_MIP" => mean, "Solve time_G" => mean,
        "Relative gap_MIP" => mean, "Gap percentage_G" => mean,
        "Greedy max found faster" => count, "Num purchases_G" => mean, "Time excess_G" => mean,
        "Solve time_MIP" => std, "Solve time_G" => std, "Time excess_G" => std
        )

# Observe
gdf_merged_agg[:,1:7]
gdf_merged_agg[:,8:12]
gdf_merged_agg[:,13:end]

##################################################################################
############################# Analysis all files 2 ###############################
gdfs_merged_agg = []

for (i,num_rows) in enumerate(num_rows_vec)
    for (j,num_cols) in enumerate(num_cols_vec)

        ind = (i-1)*length(num_cols_vec) + j

        df_MIP = deepcopy(dfs_MIP[ind])
        # df_naive = deepcopy(dfs_naive[ind])
        df_greedy = deepcopy(dfs_greedy[ind])
        # df_greedy_LP = deepcopy(dfs_greedy_LP[ind])

        #### Remove row with full budget because its an LP
        filter!(row -> row["Budget fraction"] != 1.0, df_MIP)
        # filter(row -> row["Budget fraction"] != 1.0, df_naive)
        filter!(row -> row["Budget fraction"] != 1.0, df_greedy)
        # filter!(row -> row["Budget fraction"] != 1.0, df_greedy_LP)

        #### Merge data DataFrames
        df_merged = deepcopy(df_MIP)
        select!(df_merged, Not(["Obj bound", "Dual obj", "Node count"]))
        rename!(df_merged, "Budget spent" => "Budget spent_MIP",
                "Obj val" => "Obj val_MIP",
                "Termination status" => "Termination status_MIP",
                "Solve time" => "Solve time_MIP",
                "Relative gap" => "Relative gap_MIP"
                )

        df_merged[:,"Obj val_G"] = df_greedy[!,"Obj val"]
        df_merged[:,"Termination status_G"] = df_greedy[!,"Termination status"]
        df_merged[:,"Solve time_G"] = df_greedy[!,"Solve time"]
        df_merged[:,"Num purchases_G"] = df_greedy[!,"Num purchases"]

        # df_merged[:,"Obj val_GLP"] = df_greedy_LP[!,"Obj val"]
        # df_merged[:,"Termination status_GLP"] = df_greedy_LP[!,"Termination status"]
        # df_merged[:,"Solve time_GLP"] = df_greedy_LP[!,"Solve time"]
        # df_merged[:,"Num purchases_GLP"] = df_greedy_LP[!,"Num purchases"]

        #### Compute comparisons
        df_merged[:,"Gap_G"] = (df_MIP[:, "Obj val"] - df_greedy[:,"Obj val"])
        df_merged[:,"Gap percentage_G"] = (df_MIP[:, "Obj val"] - df_greedy[:,"Obj val"]) ./ abs.(df_MIP[:, "Obj val"]) * 100
        df_merged[:, "Time excess_G"] = df_greedy[!,"Solve time"] - df_MIP[!,"Solve time"]
        df_merged[:, "Greedy max found faster"] = df_merged[!,"Time excess_G"] .< 0
        # df_merged
        # describe(df_merged)

        # df_merged[:,"Gap_GLP"] = (df_MIP[:, "Obj val"] - df_greedy_LP[:,"Obj val"])
        # df_merged[:,"Gap percentage_GLP"] = (df_MIP[:, "Obj val"] - df_greedy_LP[:,"Obj val"]) ./ abs.(df_MIP[:, "Obj val"]) * 100
        # df_merged[:, "Time excess_GLP"] = df_greedy_LP[!,"Solve time"] - df_MIP[!,"Solve time"]
        # df_merged[:, "Greedy LP max found faster"] = df_merged[!,"Time excess_GLP"] .< 0
        # # df_merged
        # # describe(df_merged)

        #### Group
        gdf_merged = groupby(df_merged, ["Num rows", "Num cols", "Budget fraction", "Termination status_MIP", "Termination status_G"])

        #### Compute aggregates
        gdf_merged_agg = combine(gdf_merged, "Termination status_MIP" => length => "Term status_MIP_count", "Termination status_G" => length => "Term status_G_count",
                "Solve time_MIP" => mean, "Solve time_G" => mean,
                "Relative gap_MIP" => mean, "Gap percentage_G" => mean,
                "Greedy max found faster" => count, "Num purchases_G" => mean, "Time excess_G" => mean,
                "Solve time_MIP" => std, "Solve time_G" => std, "Time excess_G" => std
                )

        push!(gdfs_merged_agg, gdf_merged_agg)
    end
end

t = 9
gdfs_merged_agg[t][:,1:7]
gdfs_merged_agg[t][:,8:12]
# gdfs_merged_agg[t][:,13:end]

##################################################################################
##################################### Plots ######################################

foreach(t -> filter!(row -> row["Budget fraction"] != 1.0, dfs_MIP[t]), 1:length(dfs_MIP))

colors = ["maroon", "yellow3","gold", "seagreen", "magenta", "red2", "cyan4", "indigo", "green3"]
plot(xlabel="Budget fraction", ylabel="Solution time (s)", legendtitle="Matrix size", xticks = [0.1, 0.25, 0.33, 0.5, 0.75]);
foreach(t -> plot!(dfs_MIP[t][!,"Budget fraction"], dfs_MIP[t][!,"Solve time"], markershape=:xcross, ms=5, lw=2,
        label="$(dfs_MIP[t][1,"Num rows"])" * "×" * "$(dfs_MIP[t][1,"Num cols"])", color=colors[t]
        ), 1:length(dfs_MIP))
display(plot!())

##################################################################################
################################# Cute MIP Stats #################################
opt_full_values = zeros(length(dfs_MIP))
for t = eachindex(dfs_MIP)
    ind_frac = dfs_MIP[t][:,"Budget fraction"] .== 1
    ind_term = dfs_MIP[t][:,"Termination status"] .== "OPTIMAL"
    row_bools = map((x,y) -> x == true && y == true, ind_frac, ind_term)
    opt_full_values[t] = (dfs_MIP[t][row_bools, "Obj val"])[1]
    # opt_full_values = map(df_MIP -> (df_MIP[row_bools, "Obj val"])[1], dfs_MIP)
end
# Remove full budget case
foreach(t -> filter!(row -> row["Budget fraction"] != 1.0, dfs_MIP[t]), 1:length(dfs_MIP))
# Add columns
foreach(t -> dfs_MIP[t][:, "Matches full v"] = abs.(dfs_MIP[t][!, "Obj val"] .- opt_full_values[t]) .< 1E-2, 1:length(dfs_MIP))

# Concatenate DFs and group
master_df_MIP = vcat(dfs_MIP...)
master_gdf_MIP = groupby(master_df_MIP, "Budget fraction")

combine(master_gdf_MIP, "Matches full v" => count)

master_gdf_MIP_2 = groupby(master_df_MIP, ["Budget fraction", "Termination status"])
combine(master_gdf_MIP_2, "Termination status" => length)

##################################################################################
############################### Greedy vs Greedy LP ##############################
foreach(t -> filter!(row -> row["Budget fraction"] != 1.0, dfs_greedy[t]), 1:length(dfs_greedy))
foreach(t -> filter!(row -> row["Budget fraction"] != 1.0, dfs_greedy_LP[t]), 1:length(dfs_greedy_LP))

for t = eachindex(dfs_greedy)
    dfs_greedy[t][:, "G v > GLP v"] = dfs_greedy[t][!,"Obj val"] .> dfs_greedy_LP[t][!,"Obj val"]
    dfs_greedy[t][:, "G v < GLP v"] = dfs_greedy[t][!,"Obj val"] .< dfs_greedy_LP[t][!,"Obj val"]
    dfs_greedy[t][:, "G faster"] = dfs_greedy[t][!,"Solve time"] .< dfs_greedy_LP[t][!,"Solve time"] .&& dfs_greedy[t][!, "Termination status"] .== "FINISHED"
    dfs_greedy[t][:, "GLP faster"] = dfs_greedy[t][!,"Solve time"] .> dfs_greedy_LP[t][!,"Solve time"] .&& dfs_greedy_LP[t][!, "Termination status"] .== "FINISHED"
end

master_df_greedy = vcat(dfs_greedy...)
master_gdf_greedy = groupby(master_df_greedy, "Budget fraction")

combine(master_gdf_greedy, "G v > GLP v" => count, "G v < GLP v" => count, "G faster" => count, "GLP faster" => count)

##
master_df_greedy = vcat(dfs_greedy...)
master_gdf_greedy = groupby(master_df_greedy, ["Budget fraction", "Termination status"])

master_df_greedy_LP = vcat(dfs_greedy_LP...)
master_gdf_greedy_LP = groupby(master_df_greedy_LP, ["Budget fraction", "Termination status"])

combine(master_gdf_greedy, "Solve time" => mean)
combine(master_gdf_greedy_LP, "Solve time" => mean)

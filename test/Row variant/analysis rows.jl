using DataFrames
using CSV 
using Pipe
using LinearAlgebra
using Statistics
using Plots
using StatsPlots

# TODO: Need to correct directories

##################################################################################
############################## Load files one method #############################
work_dir = pwd()
set_num = 5
subpath = work_dir * "/Experiments/Set $set_num/"
dfs = DataFrame[]

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

cdf_vec = DataFrame[]
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
dfs_MIP = DataFrame[]
# dfs_naive = DataFrame[]
dfs_greedy = DataFrame[]
# dfs_greedy_LP = DataFrame[]

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
gdfs_MIP_agg = DataFrame[]
# gdfs_naive_agg = DataFrame[]
gdfs_greedy_agg = DataFrame[]
# gdfs_greedy_LP_agg = DataFrame[]

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
gdfs_merged_agg = DataFrame[]

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

# #### Save DFs
# work_dir = pwd()
# subpath = work_dir * "/test/"
# for i = 1:length(num_rows_vec), j = 1:length(num_cols_vec)
#     CSV.write(subpath * "Grouped data frame $(num_rows_vec[i]) by $(num_cols_vec[j]).csv", gdfs_merged_agg[(i-1)*length(num_cols_vec) + j])
# end


#####################################################################################
#####################################################################################
#####################################################################################
############################### Load saved DataFrames ###############################
gdfs_merged_agg = DataFrame[]
num_rows_vec = [10, 100, 1000]
num_cols_vec = [10, 100, 1000]

work_dir = pwd()
subpath = work_dir * "/test/"
for i = 1:length(num_rows_vec), j = 1:length(num_cols_vec)
    push!(gdfs_merged_agg, CSV.File(subpath * "Grouped data frame $(num_rows_vec[i]) by $(num_cols_vec[j]).csv") |> DataFrame)
end

names(gdfs_merged_agg[1])
gdfs_merged_agg[5][:,"Gap percentage_G_mean"]

# GroupedDataFrames are the same when loaded but comparisons say otherwise when DataFrames contain NaN's
# foreach(i -> println("$(gdfs_merged_agg_2[i] == gdfs_merged_agg[i])"), 1:9)
# foreach((col, col_2) -> println("$(gdfs_merged_agg_2[8][!,col] == gdfs_merged_agg[8][!,col_2])"), names(gdfs_merged_agg[2]), names(gdfs_merged_agg_2[2]))

##################################################################################
##################################### Plots ######################################
budgets = [1/10, 1/4, 1/3, 1/2, 3/4]
t = 4
gdfs_merged_agg[t]

# foreach(t -> filter!(row -> row["Budget fraction"] != 1.0, dfs_MIP[t]), 1:length(dfs_MIP))
plot_gdfs = map(t -> filter(row -> row["Termination status_MIP"] == "OPTIMAL" && row["Termination status_G"] == "FINISHED", gdfs_merged_agg[t]), 1:length(gdfs_merged_agg))
# foreach(t -> println("$(size(plot_gdfs[t]))"), 1:length(plot_gdfs))
plot_gdfs[t]
master_plot_gdfs = vcat(plot_gdfs[4:end]...)   # ignore DFs with 10 row matrix games

#### Bar plot: matrix size vs time
master_plot_grouped_gdfs = groupby(master_plot_gdfs, "Budget fraction")
master_plot_grouped_gdfs
master_plot_grouped_gdfs[(0.1,)]  # takes rows where the budget fraction is b in (b,)
# ticks = map(row -> string(row["Num rows"]) * " × " * string(row["Num cols"]), eachrow(DataFrame(master_plot_grouped_gdfs[(0.1,)])))
# colors = ["", ""]
for b in budgets
    gdf = master_plot_grouped_gdfs[(b,)]
    ticklabel = @pipe map(row -> string(row["Num rows"]) * " × " * string(row["Num cols"]), eachrow(gdf)) |> reshape(_, 1, length(_))
    display(
    groupedbar([gdf[:,"Solve time_MIP_mean"] gdf[:,"Solve time_G_mean"]],
            xlabel="Matrix size", ylabel="Solution time (s)", labels = ["MIP" "Greedy"], title = "Budget proportion: $(round(b, digits=2))",
            bar_position=:dodge, yaxis=:log, legend=:topleft, # bar_width=0.5, # color=colors,
            xticks=(1:length(ticklabel), ticklabel)
            )
    )
    png(work_dir * "/test/Bar plot matrix size vs time log scale budget $(round(b,digits=2))")
end

#### Bar plot: matrix size vs relative G gap
master_plot_grouped_gdfs = groupby(master_plot_gdfs, "Budget fraction")
master_plot_grouped_gdfs[(0.75,)][:, "Gap percentage_G_mean"]
master_plot_grouped_gdfs[(0.75,)]  # takes rows where the budget fraction is b in (b,)
# ticks = map(row -> string(row["Num rows"]) * " × " * string(row["Num cols"]), eachrow(DataFrame(master_plot_grouped_gdfs[(0.1,)])))
# colors = ["", ""]
for b in budgets
    gdf = master_plot_grouped_gdfs[(b,)]
    ticklabel = @pipe map(row -> string(row["Num rows"]) * " × " * string(row["Num cols"]), eachrow(gdf)) |> reshape(_, 1, length(_))
    display(
    groupedbar(reshape((100 .- gdf[:,"Gap percentage_G_mean"]) / 100, length(gdf[:,"Gap percentage_G_mean"]),1),
            xlabel="Matrix size", ylabel="Relative gap", title = "Budget proportion: $(round(b, digits=2))",
            bar_position=:dodge, legend=:none, # labels = ["MIP" "Greedy"], bar_width=0.5, # color=colors,
            xticks=(1:length(ticklabel), ticklabel),
            ylim = (0,1.05)
            )
    )
    png(work_dir * "/test/Bar plot matrix size vs relative gap budget $(round(b,digits=2))")
end

# Bar plot series
# p_vec = [];
# for b in budgets
#     gdf = master_plot_grouped_gdfs[(b,)]
#     ticklabel = @pipe map(row -> string(row["Num rows"]) * " × " * string(row["Num cols"]), eachrow(gdf)) |> reshape(_, 1, length(_))
#     p = groupedbar([gdf[:,"Solve time_MIP_mean"] gdf[:,"Solve time_G_mean"]],
#             xlabel="Matrix size", ylabel="Sol'n time", labels = ["MIP" "Greedy"], title = "Budget proportion: $(round(b, digits=2))",
#             bar_position=:dodge, yaxis=:log, # bar_width=0.5, # color=colors,
#             xticks=(1:length(ticklabel), ticklabel)
#             )
#     push!(p_vec, p)
#     # png(work_dir * "/test/Bar plot matrix size vs time log scale budget $(round(b,digits=2))")
# end
# plt = plot(p_vec..., layout=@layout([° ° °; ° ° _]))
# display(plt)

#### Plot: budget vs time, all methods
master_plot_grouped_gdfs = groupby(master_plot_gdfs, ["Num rows", "Num cols"]) #, "Budget fraction"])
master_plot_grouped_gdfs[3]
@pipe master_plot_grouped_gdfs[(100,100)][:,"Budget fraction"] |> round.(_, digits=2) |> string.(_) |> reshape(_, 1, length(_))

matrix_dims = [(100,10), (100,100), (100,1000), (1000,10), (1000,100)]
# colors = ["red2", "green3", "magenta", "indigo", "turquoise3"]
colors = ["green3", "red2", "magenta", "gold2", "turquoise3"] # "gold2"

# Plots individually
plt = plot(legend=:outertopright);
for (counter, (num_rows, num_cols)) in enumerate(matrix_dims)
# for num_rows in num_rows_vec[2:end], num_cols in num_cols_vec[2:end]
    gdf = master_plot_grouped_gdfs[(num_rows, num_cols)]
    ticklabel = @pipe gdf[:,"Budget fraction"] |> round.(_, digits=2) |> string.(_) |> reshape(_, 1, length(_)) # @pipe map(row -> string(round(row["Budget fraction"],digits=2)), eachrow(gdf)) |> reshape(_, 1, length(_))
    p = plot(gdf[:,"Budget fraction"], gdf[:,"Solve time_MIP_mean"],
            # [gdf[:,"Solve time_MIP_mean"] gdf[:,"Solve time_G_mean"]], 
            xlabel="Budget proportion", ylabel="Solution time (s)", labels = "MIP $(num_rows) × $(num_cols)",
            color=colors[counter], lw=2, yaxis=:log,
            # xticks=(1:length(ticklabel), ticklabel)
            xticks = (gdf[:,"Budget fraction"], ticklabel), markershape=:xcross, ms=5
            )
    plot!(gdf[:,"Budget fraction"], gdf[:,"Solve time_G_mean"], label="G $(num_rows) × $(num_cols)", 
            color=colors[counter],lw=2,markershape=:circle,ms=5,ls=:dash)
    push!(plots, p)
end
display(plt)
# png(work_dir * "/test/Budget vs time log scale all methods ms 5")

# Subplots
# plots = []
for (counter, (num_rows, num_cols)) in enumerate(matrix_dims)
    if (num_rows, num_cols) == (100, 1000)
        continue
    end

    gdf = master_plot_grouped_gdfs[(num_rows, num_cols)]
    ticklabel = @pipe gdf[:,"Budget fraction"] |> round.(_, digits=2) |> string.(_) |> reshape(_, 1, length(_)) # @pipe map(row -> string(round(row["Budget fraction"],digits=2)), eachrow(gdf)) |> reshape(_, 1, length(_))
    p = plot(gdf[:,"Budget fraction"], gdf[:,"Solve time_MIP_mean"], label=:none,
            xlabel="Budget proportion", ylabel="Solution time (s)", #labels = "MIP $(num_rows) × $(num_cols)",
            color=colors[counter], lw=2, yaxis=:log, title = "Matrix size: $(num_rows) × $(num_cols)",
            xticks = (gdf[:,"Budget fraction"], ticklabel), markershape=:xcross, #ms=5
            )
    plot!(gdf[:,"Budget fraction"], gdf[:,"Solve time_G_mean"], label=:none, #label="G $(num_rows) × $(num_cols)", 
            color=colors[counter],lw=2,ls=:dash,markershape=:circle, #ms=5
            )
    push!(plots, p)
    display(p)
    png(work_dir * "/test/Budget vs time log scale all methods matrix size $(num_rows) by $(num_cols)")
end
# plt = plot(plots..., layout = (2,2), legend=:outertopright)
# png(work_dir * "/test/Budget vs time log scale all methods ms 5")


#### Heatmap: Sol'n time as function of budget and matrix size (fixed row size)
# # TODO: Need x-,y-,z- vectors to correspond to (x,y,z) values. Further issue is some (x,y) don't have z values. Use time limits?
# master_plot_grouped_gdfs = groupby(master_plot_gdfs, ["Num rows"])
# master_plot = master_plot_grouped_gdfs[(100,)]
# unique(master_plot[:,"Budget fraction"])
# unique(master_plot[:,"Num cols"])

# heatmap(unique(master_plot[:,"Num cols"]), unique(master_plot[:,"Budget fraction"]), 
#     master_plot[:,"Solve time_MIP_mean"],
#     # c=cgrad([:blue, :white,:red, :yellow]),
#     xlabel="x values", ylabel="y values", # zlabel="z range",
#     title="My title")


#########
# colors = ["maroon", "yellow3","gold", "seagreen", "magenta", "red2", "cyan4", "indigo", "green3"]
# plot(xlabel="Budget fraction", ylabel="Solution time (s)", legendtitle="Matrix size", xticks = [0.1, 0.25, 0.33, 0.5, 0.75]);
# foreach(t -> plot!(dfs_MIP[t][!,"Budget fraction"], dfs_MIP[t][!,"Solve time"], markershape=:xcross, ms=5, lw=2,
#         label="$(dfs_MIP[t][1,"Num rows"])" * "×" * "$(dfs_MIP[t][1,"Num cols"])", color=colors[t]
#         ), 1:length(dfs_MIP))
# display(plot!())

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


using DataFrames
using CSV 
using Pipe
using LinearAlgebra
using Statistics
using Plots
using StatsPlots

##################################################################################
############################## Load files one method #############################
work_dir = pwd()
# Analysis does one of either MILP or LP at a time. We compare cuts and no cuts, for a fixed set type
# set_type, relax = "Set MILP Cuts", false
set_type, relax = "Set LP Cuts", true
set_num = 1
subpath = work_dir * "/Experiments/$set_type $set_num/"

num_rows_vec = [10,25,50,100]
num_cols_vec = [10,25,50,100]
c_r_entry_range = 2:5
c_s_entry_range_vec = [2:5, 11:15, 21:25]
num_cond_dom_rows_vec = [1,5,10]
exp_type_vec = ["cuts", "no cuts"]

# Individual
dfs = DataFrame[]
for num_rows in num_rows_vec, num_cols in num_cols_vec, (c_s_ind, c_s_entry_range) in pairs(c_s_entry_range_vec), exp_type in exp_type_vec

    filename = relax ? subpath * "Matrices $num_rows by $num_cols column prices index $(c_s_ind) LP $exp_type.txt" :
            subpath * "Matrices $num_rows by $num_cols column prices index $(c_s_ind) MILP $exp_type.txt"

    df = CSV.File(filename,
                    delim='\t',
                    ignorerepeated=true,
                    header = 3, # 3 for on line 3 or Vector of the names as strings or symbols
                    skipto = 4,
                    ) |> DataFrame

    select!(df, Not("c_s range"))
    df[:, "c_r range"] .= "$c_r_entry_range"
    df[:, "c_s range"] .= "$c_s_entry_range"

    select!(df, [names(df)[3:end]; names(df)[1:2]])  # Move "Matrix seed" and "Costs seed" to end
    push!(dfs, df)
end

dfs[1][:,end-7:end]

# Load cuts and no cuts into Dictionary
master_dict = Dict{String,DataFrame}()
# DFs defer only by cuts and no cuts

for exp_type in exp_type_vec
    dfs = DataFrame[]
    for num_rows in num_rows_vec, num_cols in num_cols_vec, (c_s_ind, c_s_entry_range) in pairs(c_s_entry_range_vec)

        filename = relax ? subpath * "Matrices $num_rows by $num_cols column prices index $(c_s_ind) LP $exp_type.txt" :
                subpath * "Matrices $num_rows by $num_cols column prices index $(c_s_ind) MILP $exp_type.txt"

        df = CSV.File(filename,
                        delim='\t',
                        ignorerepeated=true,
                        header = 3, # 5 for on line 5 or Vector of the names as strings or symbols
                        skipto = 4,
                        ) |> DataFrame

        select!(df, Not("c_s range"))
        df[:, "c_r range"] .= "$c_r_entry_range"
        df[:, "c_s range"] .= "$c_s_entry_range"

        select!(df, [names(df)[3:end]; names(df)[1:2]])  # Move "Matrix seed" and "Costs seed" to end
        push!(dfs, df)
    end
    # Merge DFs of different matrix sizes, c_s range so that they defer only by cuts and no cuts
    df_stacked = vcat(dfs...)
    master_dict["$exp_type"] = df_stacked
end

master_dict["cuts"]
master_dict["no cuts"]

# Quick sanity check
obj_val_diff = master_dict["cuts"][!,"Obj val"] - master_dict["no cuts"][!,"Obj val"]
max_val, max_ind = findmax(obj_val_diff)
master_dict["cuts"][max_ind,:]
master_dict["no cuts"][max_ind,:]
min_val, min_ind = findmin(obj_val_diff)
master_dict["cuts"][min_ind,end-8:end]
master_dict["no cuts"][min_ind,end-8:end]


##################################################################################
##### Analysis (two separate DFs)
gdf_cuts = groupby(master_dict["cuts"], ["c_s range", "Budget fraction", "Num cond dom rows", "Termination status"])
gdf_no_cuts = groupby(master_dict["no cuts"], ["c_s range", "Budget fraction", "Num cond dom rows", "Termination status"])
gdf_cuts[("2:5",0.75,1,"OPTIMAL")]

# gdf_agg_cuts = combine(gdf_cuts, "Termination status" => length => "Term status_count",
#         "Solve time" => mean, "Solve time" => std,
#         "Rows purchased" => mean, "Cols purchased" => mean,
#         "Relative gap" => mean,
#         # "Methods gap percentage" => mean, "Cuts faster" => count,
#         # "Time excess no cuts over cuts" => mean, "Time excess no cuts over cuts" => std,
#         )
# gdf_agg_cuts

# gdf_agg_no_cuts = combine(gdf_no_cuts, "Termination status" => length => "Term status_count",
#         "Solve time" => mean, "Solve time" => std,
#         "Rows purchased" => mean, "Cols purchased" => mean,
#         "Relative gap" => mean,
#         # "Methods gap percentage" => mean, "Cuts faster" => count,
#         # "Time excess no cuts over cuts" => mean, "Time excess no cuts over cuts" => std,
#         )
# gdf_agg_no_cuts


##################################################################################
##### Analysis (one stacked DF, cuts specifier added as a column)
df_temp_cuts = deepcopy(master_dict["cuts"])
df_temp_no_cuts = deepcopy(master_dict["no cuts"])
df_temp_cuts[:,"Cuts"] .= "cuts"
df_temp_no_cuts[:,"Cuts"] .= "no cuts"
master_df = vcat(df_temp_cuts, df_temp_no_cuts)

gdf = groupby(master_df, ["Cuts","c_s range", "Budget fraction", "Num cond dom rows", "Termination status"])
gdf[(true,"2:5",0.75,1,"OPTIMAL")]

gdf_agg = combine(gdf, "Termination status" => length => "Term status_count",
        "Solve time" => mean, "Solve time" => std,
        "Rows purchased" => mean, "Cols purchased" => mean,
        "Relative gap" => mean,
        # "Methods gap percentage" => mean, "Cuts faster" => count,
        # "Time excess no cuts over cuts" => mean, "Time excess no cuts over cuts" => std,
        )
gdf_agg


##################################################################################
##### Analysis merged DF
df_temp_cuts = deepcopy(master_dict["cuts"])
df_temp_no_cuts = deepcopy(master_dict["no cuts"])
select!(df_temp_cuts, Not(["Matrix seed", "Costs seed"]))
select!(df_temp_no_cuts, Not(["Matrix seed", "Costs seed"]))

df_merged = deepcopy(df_temp_cuts)
rename!(df_merged, "Budget spent" => "Budget spent_cuts",
        "Obj val" => "Obj val_cuts", "Obj bound" => "Obj bound_cuts", "Dual obj" => "Dual obj_cuts",
        "Rows purchased" => "Rows purchased_cuts", "Cols purchased" => "Cols purchased_cuts",
        "Termination status" => "Term status_cuts", "Solve time" => "Solve time_cuts",
        "Relative gap" => "Rel gap_cuts", "Node count" => "Node count_cuts"
        )

# Add cols from DF with no cuts
df_merged[:,"Budget spent_no_cuts"] = df_temp_no_cuts[!,"Budget spent"]
df_merged[:,"Obj val_no_cuts"] = df_temp_no_cuts[!,"Obj val"]
df_merged[:,"Obj bound_no_cuts"] = df_temp_no_cuts[!,"Obj bound"]
df_merged[:,"Dual obj_no_cuts"] = df_temp_no_cuts[!,"Dual obj"]
df_merged[:,"Rows purchased_no_cuts"] = df_temp_no_cuts[!,"Rows purchased"]
df_merged[:,"Cols purchased_no_cuts"] = df_temp_no_cuts[!,"Cols purchased"]
df_merged[:,"Term status_no_cuts"] = df_temp_no_cuts[!,"Termination status"]
df_merged[:,"Solve time_no_cuts"] = df_temp_no_cuts[!,"Solve time"]
df_merged[:,"Rel gap_no_cuts"] = df_temp_no_cuts[!,"Relative gap"]
df_merged[:,"Node count_no_cuts"] = df_temp_no_cuts[!,"Node count"]

##### Compute comparison columns
# Counting
df_merged[:,"Faster_cuts"] = df_merged[:,"Solve time_cuts"] .< df_merged[:,"Solve time_no_cuts"]
df_merged[:,"Faster_no_cuts"] = df_merged[:,"Solve time_cuts"] .> df_merged[:,"Solve time_no_cuts"]
# df_merged[:,"Rel gap better_cuts"] = df_merged[:,"Rel gap_cuts"] .< df_merged[:,"Rel gap_no_cuts"]
# df_merged[:,"Node count better_cuts"] = df_merged[:,"Node count_cuts"] .< df_merged[:,"Node count_no_cuts"]
# Quantifying
df_merged[:,"Time gap_cuts"] = df_merged[:,"Solve time_cuts"] - df_merged[:,"Solve time_no_cuts"]
df_merged[:,"Time rel gap_cuts"] = (df_merged[:,"Solve time_cuts"] - df_merged[:,"Solve time_no_cuts"]) ./ abs.(df_merged[:,"Solve time_no_cuts"])
# df_merged[:,"Node count excess_cuts"] = df_merged[:,"Node count_cuts"] - df_merged[:,"Node count_no_cuts"]
# Understand suboptimal solutions and also verify problems match
df_merged[:,"Obj val gap LP"] = df_merged[:,"Obj val_no_cuts"] - df_merged[:,"Obj val_cuts"]
df_merged[:,"Obj val rel gap LP"] = (df_merged[:,"Obj val_no_cuts"] - df_merged[:,"Obj val_cuts"] ) ./ df_merged[:,"Obj val_no_cuts"]
df_merged[:,"Obj val larger_cuts"] = df_merged[:,"Obj val_cuts"] .> df_merged[:,"Obj val_no_cuts"]
df_merged[:,"Budget spent gap_cuts"] = df_merged[:,"Budget spent_cuts"] - df_merged[:,"Budget spent_no_cuts"]
df_merged[:,"Budget spent rel gap_cuts"] = (df_merged[:,"Budget spent_cuts"] - df_merged[:,"Budget spent_no_cuts"]) ./ abs.(df_merged[:,"Budget spent_no_cuts"])
df_merged[:,"Budget spent larger_cuts"] = df_merged[:,"Budget spent_cuts"] .> df_merged[:,"Budget spent_no_cuts"]
# df_merged[:,"Rows purchased_cuts - Rows purchased_no_cuts"] = df_merged[:,"Rows purchased_cuts"] - df_merged[:,"Rows purchased_no_cuts"]
# df_merged[:,"Cols purchased_cuts - Cols purchased_no_cuts"] = df_merged[:,"Cols purchased_cuts"] - df_merged[:,"Cols purchased_no_cuts"]


# TODO: Careful with analysis among different matrix sizes (e.g., solve time average)
gdf = groupby(df_merged, ["c_s range", "Num cond dom rows", "Budget fraction"])
# gdf = groupby(df_merged, ["c_s range", "Num cond dom rows", "Budget fraction", "Term status_cuts", "Term status_no_cuts"])

gdf_agg = combine(gdf, nrow => "Group size",
        "Faster_cuts" => count, "Faster_no_cuts" => count,
        "Solve time_cuts" => mean, "Solve time_no_cuts" => mean, 
        "Time gap_cuts" => mean, "Time rel gap_cuts" => mean,
        "Obj val larger_cuts" => count, "Obj val gap LP" => mean, "Obj val rel gap LP" => mean,
        # "Rel gap better_cuts" => count,
        # "Node count better_cuts" => count, "Node count excess_cuts" => mean,
        # "Node count_cuts" => mean, "Node count_no_cuts" => mean,
        # "Rows purchased_cuts - Rows purchased_no_cuts" => mean, "Cols purchased_cuts - Cols purchased_no_cuts" => mean,
        "Budget spent larger_cuts" => count, "Budget spent gap_cuts" => mean, "Budget spent rel gap_cuts" => mean,
        )

gdf_agg[:,1:10]
gdf_agg[:,11:end]
# gdf_agg[:,13:end]

# TODO: Characterize the LP relaxation. If r_i or s_j is 0 in an LP solution, will it be 0 in the MILP solution? i.e, Can we solve the LP and fix those zeros in the MILP?



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
set_type, relax = "Set MILP Cuts", false
# set_type, relax = "Set LP Cuts", true
set_num = 1
subpath = work_dir * "/Experiments/$set_type $set_num/"

# num_rows_vec = [10,25,50,100]  # LP
# num_cols_vec = [10,25,50,100]  # LP
num_rows_vec = [10,50,100]  # MILP
num_cols_vec = [10,50,100]  # MILP
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
# Differences in MILP values are due to having hit the TIME LIMIT
# For matrix size {10,50,100}^2, differences in LP objective occur more often at matrix sizes with one dimension being 10 
# i.e., for an LP analysis, analyze all sizes {10,25,50,100}^2
obj_val_diff = master_dict["cuts"][!,"Obj val"] - master_dict["no cuts"][!,"Obj val"]
max_val, max_ind = findmax(obj_val_diff)
master_dict["cuts"][max_ind,:]
master_dict["no cuts"][max_ind,:]
master_dict["cuts"][max_ind,"Obj val"] - master_dict["no cuts"][max_ind,"Obj val"] / master_dict["no cuts"][max_ind,"Obj val"]
min_val, min_ind = findmin(obj_val_diff)
master_dict["cuts"][min_ind,:]
master_dict["no cuts"][min_ind,:]


# ##################################################################################
# ##### Analysis (one stacked DF, cuts specifier added as a column)
# # Can only do this if columns are the same and it makes it hard to compare methods side by side
# df_temp_cuts = deepcopy(master_dict["cuts"])
# df_temp_no_cuts = deepcopy(master_dict["no cuts"])
# df_temp_cuts[:,"Cuts"] .= "cuts"
# df_temp_no_cuts[:,"Cuts"] .= "no cuts"
# master_df = vcat(df_temp_cuts, df_temp_no_cuts)

# gdf = groupby(master_df, ["Cuts","c_s range", "Budget fraction", "Num cond dom rows", "Termination status"])
# gdf[(true,"2:5",0.75,1,"OPTIMAL")]

# gdf_agg = combine(gdf, "Termination status" => length => "Term status_count",
#         "Solve time" => mean, "Solve time" => std,
#         "Rows purchased" => mean, "Cols purchased" => mean,
#         "Relative gap" => mean,
#         # "Methods gap percentage" => mean, "Cuts faster" => count,
#         # "Time excess no cuts over cuts" => mean, "Time excess no cuts over cuts" => std,
#         )
# gdf_agg


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
        "Relative gap" => "Rel gap_cuts",  # irrelevant for LP
        "Node count" => "Node count_cuts",  # irrelevant for LP
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
df_merged[:,"Time gap_cuts"] = df_merged[:,"Solve time_cuts"] - df_merged[:,"Solve time_no_cuts"]
df_merged[:,"Time rel gap_cuts"] = (df_merged[:,"Solve time_cuts"] - df_merged[:,"Solve time_no_cuts"]) ./ df_merged[:,"Solve time_no_cuts"]
df_merged[:,"Rel gap better_cuts"] = df_merged[:,"Rel gap_cuts"] .< df_merged[:,"Rel gap_no_cuts"]
df_merged[:,"Rel gap gap_cuts"] = df_merged[:,"Rel gap_cuts"] - df_merged[:,"Rel gap_no_cuts"]
df_merged[:,"Rel gap rel gap_cuts"] = df_merged[:,"Rel gap_cuts"] - df_merged[:,"Rel gap_no_cuts"] ./ abs.(df_merged[:,"Rel gap_no_cuts"])
df_merged[:,"Node count better_cuts"] = df_merged[:,"Node count_cuts"] .< df_merged[:,"Node count_no_cuts"]
df_merged[:,"Node count gap_cuts"] = df_merged[:,"Node count_cuts"] - df_merged[:,"Node count_no_cuts"]
df_merged[:,"Node count rel gap_cuts"] = df_merged[:,"Node count_cuts"] - df_merged[:,"Node count_no_cuts"] ./ df_merged[:,"Node count_no_cuts"]

df_merged[:,"Obj val larger_cuts"] = df_merged[:,"Obj val_cuts"] .> df_merged[:,"Obj val_no_cuts"]
df_merged[:,"Obj val larger_no_cuts"] = df_merged[:,"Obj val_cuts"] .< df_merged[:,"Obj val_no_cuts"]
df_merged[:,"Obj val gap_cuts"] = df_merged[:,"Obj val_no_cuts"] - df_merged[:,"Obj val_cuts"]
df_merged[:,"Obj val rel gap_cuts"] = (df_merged[:,"Obj val_no_cuts"] - df_merged[:,"Obj val_cuts"] ) ./ abs.(df_merged[:,"Obj val_no_cuts"])  # What does this mean for negative numbers?
df_merged[:,"Budget spent larger_cuts"] = df_merged[:,"Budget spent_cuts"] .> df_merged[:,"Budget spent_no_cuts"]
df_merged[:,"Budget spent gap_cuts"] = df_merged[:,"Budget spent_cuts"] - df_merged[:,"Budget spent_no_cuts"]
df_merged[:,"Budget spent rel gap_cuts"] = (df_merged[:,"Budget spent_cuts"] - df_merged[:,"Budget spent_no_cuts"]) ./ df_merged[:,"Budget spent_no_cuts"]
df_merged[:,"R_cuts > R_no_cuts"] = df_merged[:,"Rows purchased_cuts"] .> df_merged[:,"Rows purchased_no_cuts"]
df_merged[:,"S_cuts > S_no_cuts"] = df_merged[:,"Cols purchased_cuts"] .> df_merged[:,"Cols purchased_no_cuts"]


# TODO: Careful with analysis among different matrix sizes (small problems solved without branch and bound)
# All problem sizes
# gdf = groupby(df_merged, ["c_s range", "Num cond dom rows", "Budget fraction"])
# gdf = groupby(df_merged, ["c_s range", "Num cond dom rows", "Budget fraction", "Term status_cuts", "Term status_no_cuts"])
gdf = groupby(df_merged, ["Term status_cuts", "Term status_no_cuts", "Num rows", "Num cols", "c_s range", "Budget fraction"])
gdf = groupby(df_merged, ["Term status_cuts", "Term status_no_cuts", "c_s range", "Budget fraction"])

# Remove small problems
df_merged_temp = deepcopy(df_merged)
# df_merged_temp = df_merged_temp[df_merged_temp[:,"Num rows"] .!= 10 .|| df_merged_temp[:,"Num cols"] .!= 10,:]
df_merged_temp = df_merged_temp[df_merged_temp[:,"Num rows"] .!= 10 .&& df_merged_temp[:,"Num cols"] .!= 10,:]
gdf = groupby(df_merged_temp, ["Term status_cuts", "Term status_no_cuts", "c_s range", "Budget fraction"])

gdf_agg = combine(gdf, nrow => "Group size",
        # "Faster_cuts" => count, "Faster_no_cuts" => count,
        "Solve time_cuts" => mean, "Solve time_no_cuts" => mean, 
        # "Time gap_cuts" => mean, "Time rel gap_cuts" => mean,
        # "Obj val larger_cuts" => count, "Obj val larger_no_cuts" => count, "Obj val gap_cuts" => mean, "Obj val rel gap_cuts" => mean,
        # "Budget spent larger_cuts" => count, "Budget spent gap_cuts" => mean, "Budget spent rel gap_cuts" => mean,
        # "Rel gap better_cuts" => count, "Rel gap gap_cuts" => mean, "Rel gap rel gap_cuts" => mean,
        # "Node count better_cuts" => count, "Node count gap_cuts" => maximum, "Node count gap_cuts" => minimum, "Node count gap_cuts" => mean, "Node count rel gap_cuts" => mean,
        # "Node count_cuts" => mean, "Node count_no_cuts" => mean,
        # "R_cuts > R_no_cuts" => count, "S_cuts > S_no_cuts" => count
        )

gdf_agg[:,3:end]

gdf_agg2 = gdf_agg[:, 
    ["Term status_cuts", "Term status_no_cuts", "c_s range", "Budget fraction", "Group size",
    "Faster_cuts_count", "Faster_no_cuts_count",
    # "Solve time_cuts_mean", "Solve time_no_cuts_mean",  # Are similar
    # "Time rel gap_cuts_mean",  # Solve time average for cuts slightly larger, even if it is faster on more instances (i.e., competitive)
    "Obj val larger_cuts_count",  # For time outs
    # "Budget spent larger_cuts_count",
    # "Budget spent gap_cuts_mean",
    # "Rel gap better_cuts_count",  # Since same obj val, this is redundant with "Obj val larger_cuts_count"
    # "Rel gap rel gap_cuts_mean",  # Can be uninformative if values are small
    # "Node count better_cuts_count",  # Sometime competitive
    # "Node count rel gap_cuts_mean",  # Average can be uninformative if values are much larger for some problems
    "Node count gap_cuts_minimum", "Node count gap_cuts_maximum",
    "Node count_cuts_mean", "Node count_no_cuts_mean"
    # "R_cuts > R_no_cuts_count", "S_cuts > S_no_cuts_count"  # Are similar
    ]]

gdf_agg2[:,1:7]
gdf_agg2[:,8:11]
gdf_agg2[:,12:16]
gdf_agg2[:, 17:19]

gdf_agg2
gdf_agg2[:,1:8]
gdf_agg2[:,9:end]

# Do analysis for LP? Would want to do for {10,25,50,100} to get more examples of cutting solution off

# TODO: Characterize the LP relaxation. If r_i or s_j is 0 in an LP solution, will it be 0 in the MILP solution? i.e, Can we solve the LP and fix those zeros in the MILP?



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
set_type = "Set MILP"
set_type = "Set Greedy MILP"
set_num = 1
subpath = work_dir * "/Experiments/$set_type $set_num/"

num_rows_vec = [10,50,100]
num_cols_vec = [10,50,100]
c_r_entry_range = 2:5
c_s_entry_range_vec = [2:5, 11:15, 21:25]
set_type_vec = ["Set MILP", "Set Greedy MILP"]

# Individual
dfs = DataFrame[]
for num_rows in num_rows_vec, num_cols in num_cols_vec, (c_s_ind, c_s_entry_range) in pairs(c_s_entry_range_vec)

    filename = subpath * "Matrices $num_rows by $num_cols column prices index $(c_s_ind) MILP.txt"

    df = CSV.File(filename,
                    delim='\t',
                    ignorerepeated=true,
                    header = 3, # 3 for on line 3 or Vector of the names as strings or symbols
                    skipto = 4,
                    ) |> DataFrame

    # select!(df, Not("c_s range"))  # c_s range not logged in these files
    df[:, "c_r range"] .= "$c_r_entry_range"
    df[:, "c_s range"] .= "$c_s_entry_range"

    select!(df, [names(df)[3:end]; names(df)[1:2]])  # Move "Matrix seed" and "Costs seed" to end
    push!(dfs, df)
end

dfs[1]
df_stacked = vcat(dfs...)

# Load MILP and Greedy MILP into a Dictionary
master_dict = Dict{String,DataFrame}()

for set_type in set_type_vec
    dfs = DataFrame[]
    for num_rows in num_rows_vec, num_cols in num_cols_vec, (c_s_ind, c_s_entry_range) in pairs(c_s_entry_range_vec)

        subpath = work_dir * "/Experiments/$set_type $set_num/"
        filename = set_type == "Set MILP" ? subpath * "Matrices $num_rows by $num_cols column prices index $(c_s_ind) MILP.txt" :
                subpath * "Matrices $num_rows by $num_cols column prices index $(c_s_ind) greedy MILP.txt"

        df = CSV.File(filename,
                        delim='\t',
                        ignorerepeated=true,
                        header = 3, # 3 for on line 3 or Vector of the names as strings or symbols
                        skipto = 4,
                        ) |> DataFrame

        # select!(df, Not("c_s range"))  # c_s range not logged in these files
        df[:, "c_r range"] .= "$c_r_entry_range"
        df[:, "c_s range"] .= "$c_s_entry_range"

        select!(df, [names(df)[3:end]; names(df)[1:2]])  # Move "Matrix seed" and "Costs seed" to end
        push!(dfs, df)
    end
    df_stacked = vcat(dfs...)
    if set_type == "Set MILP"
        master_dict["MILP"] = df_stacked
    elseif set_type == "Set Greedy MILP"
        master_dict["greedy"] = df_stacked
    end
end

master_dict["MILP"]
master_dict["greedy"]

# Quick sanity check
obj_val_diff = master_dict["MILP"][!,"Obj val"] - master_dict["greedy"][!,"Obj val"]
max_val, max_ind = findmax(obj_val_diff)
master_dict["greedy"][max_ind,:]
master_dict["MILP"][max_ind,:]
min_val, min_ind = findmin(obj_val_diff)
master_dict["greedy"][min_ind,:]
master_dict["MILP"][min_ind,:]
# The optimal objectives are within 1% optimality gap

##################################################################################
##### Analysis merged DF
df_temp_MILP = deepcopy(master_dict["MILP"])
df_temp_greedy = deepcopy(master_dict["greedy"])
select!(df_temp_MILP, Not(["Matrix seed", "Costs seed"]))
select!(df_temp_greedy, Not(["Matrix seed", "Costs seed"]))

df_merged = deepcopy(df_temp_MILP)
rename!(df_merged, "Budget spent" => "Budget spent_MILP",
        "Obj val" => "Obj val_MILP", "Obj bound" => "Obj bound_MILP", "Dual obj" => "Dual obj_MILP",
        "Rows purchased" => "Rows purchased_MILP", "Cols purchased" => "Cols purchased_MILP",
        "Termination status" => "Term status_MILP", "Solve time" => "Solve time_MILP",
        "Relative gap" => "Rel gap_MILP", "Node count" => "Node count_MILP"
        )

# Add cols from Greedy MILP
df_merged[:,"Budget spent_greedy"] = df_temp_greedy[!,"Budget spent"]
df_merged[:,"Obj val_greedy"] = df_temp_greedy[!,"Obj val"]
df_merged[:,"Rows purchased_greedy"] = df_temp_greedy[!,"Rows purchased"]
df_merged[:,"Cols purchased_greedy"] = df_temp_greedy[!,"Cols purchased"]
df_merged[:,"Term status_greedy"] = df_temp_greedy[!,"Termination status"]
df_merged[:,"Solve time_greedy"] = df_temp_greedy[!,"Solve time"]
# df_merged[:,"Num purchases_greedy"] = df_temp_greedy[!,"Num purchases"]
# df_merged[:,"Node count_greedy"] = df_temp_greedy[!,"Node count"]  # Not actual node countno_MILP

##### Compute comparison columns
df_merged[:,"Faster_MILP"] = df_merged[:,"Solve time_MILP"] .< df_merged[:,"Solve time_greedy"]
df_merged[:,"Faster_greedy"] = df_merged[:,"Solve time_MILP"] .> df_merged[:,"Solve time_greedy"]
df_merged[:,"Time gap_MILP"] = df_merged[:,"Solve time_MILP"] - df_merged[:,"Solve time_greedy"]
df_merged[:,"Time rel gap_MILP"] = (df_merged[:,"Solve time_MILP"] - df_merged[:,"Solve time_greedy"]) ./ df_merged[:,"Solve time_greedy"]
# df_merged[:,"Node count better_MILP"] = df_merged[:,"Node count_MILP"] .< df_merged[:,"Node count_greedy"]
# df_merged[:,"Node count gap_MILP"] = df_merged[:,"Node count_MILP"] - df_merged[:,"Node count_greedy"]
df_merged[:,"Obj val larger_greedy"] = df_merged[:,"Obj val_greedy"] .> df_merged[:,"Obj val_MILP"]
df_merged[:,"Obj val gap_greedy"] = df_merged[:,"Obj val_MILP"] - df_merged[:,"Obj val_greedy"]
df_merged[:,"Obj val rel gap_greedy"] = (df_merged[:,"Obj val_MILP"] - df_merged[:,"Obj val_greedy"] ) ./ abs.(df_merged[:,"Obj val_MILP"])  # What does this mean for negative numbers?
df_merged[:,"Budget spent larger_greedy"] = df_merged[:,"Budget spent_greedy"] .> df_merged[:,"Budget spent_MILP"]
df_merged[:,"Budget spent larger_MILP"] = df_merged[:,"Budget spent_greedy"] .< df_merged[:,"Budget spent_MILP"]
df_merged[:,"Budget spent gap_greedy"] = df_merged[:,"Budget spent_greedy"] - df_merged[:,"Budget spent_MILP"]
df_merged[:,"Budget spent rel gap_greedy"] = (df_merged[:,"Budget spent_greedy"] - df_merged[:,"Budget spent_MILP"]) ./ df_merged[:,"Budget spent_MILP"]
# TODO: Compare num row purchases and num col purchases (on 2nd thought, the size of the matrix and the budget means the averages don't summarize well)
# df_merged[:,"Rows purchased_MILP - Rows purchased_greedy"] = df_merged[:,"Rows purchased_MILP"] - df_merged[:,"Rows purchased_greedy"]
# df_merged[:,"Cols purchased_MILP - Cols purchased_greedy"] = df_merged[:,"Cols purchased_MILP"] - df_merged[:,"Cols purchased_greedy"]

# Group
gdf = groupby(df_merged, ["Term status_MILP", "c_s range", "Budget fraction"])
gdf = groupby(df_merged, ["c_s range", "Term status_MILP"])
gdf = groupby(df_merged, ["Num rows", "Num cols", "Term status_MILP"])
gdf = groupby(df_merged, ["Budget fraction", "Term status_MILP"])
gdf = groupby(df_merged, ["Term status_MILP", "Budget fraction"])

gdf_agg = combine(gdf, nrow => "Group size",
        "Faster_MILP" => count, "Faster_greedy" => count,
        "Solve time_MILP" => mean, "Solve time_greedy" => mean, 
        "Time gap_MILP" => mean, "Time rel gap_MILP" => mean,
        "Obj val larger_greedy" => count, "Obj val gap_greedy" => mean, "Obj val rel gap_greedy" => mean,
        "Rel gap_MILP" => mean,
        "Budget spent larger_greedy" => count, "Budget spent larger_MILP" => count, "Budget spent gap_greedy" => mean, "Budget spent rel gap_greedy" => mean,
        # "Node count better_MILP" => count, "Node count gap_MILP" => mean,
        # "Node count_MILP" => mean, "Node count_greedy" => mean,
        # "Rows purchased_cuts - Rows purchased_greedy" => mean, "Cols purchased_cuts - Cols purchased_greedy" => mean,
        )

gdf_agg[:,1:end-7]
gdf_agg[:,end-6:end]

gdf_agg[:, ["Term status_MILP", "c_s range", "Budget fraction", "Group size", "Faster_greedy_count", "Obj val rel gap_greedy_mean", "Rel gap_MILP_mean"]]

# TODO: Come up with MILP performance analysis, e.g., time vs budget, time vs column Costs
# TODO: Come up with analysis to showcase preference for removing columns
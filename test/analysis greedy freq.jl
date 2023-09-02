
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
set_type, ranking = "Set Greedy Freq", "with dual var"
set_type, ranking = "Set Greedy Freq", "with s"
set_num = 1
subpath = work_dir * "/Experiments/$set_type $set_num/"

num_rows_vec = [10,50,100]
num_cols_vec = [10,50,100]
c_r_entry_range = 2:5
c_s_entry_range_vec = [2:5, 11:15, 21:25]
set_type_vec = ["Set MILP", "Set Greedy Freq"]
ranking_vec = ["with dual var", "with s"]

# Individual
dfs = DataFrame[]
for num_rows in num_rows_vec, num_cols in num_cols_vec, (c_s_ind, c_s_entry_range) in pairs(c_s_entry_range_vec)

    filename = set_type == "Set MILP" ? subpath * "Matrices $num_rows by $num_cols column prices index $(c_s_ind) MILP.txt" : 
            subpath * "Matrices $num_rows by $num_cols column prices index $(c_s_ind) greedy freq ranking $ranking.txt"

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

# Load MILP and Greedy Freqs into a Dictionary
master_dict = Dict{String,DataFrame}()

for set_type in set_type_vec, ranking in ranking_vec
    if set_type == "Set Greedy Freq"
        for ranking in ranking_vec
            dfs = DataFrame[]
            for num_rows in num_rows_vec, num_cols in num_cols_vec, (c_s_ind, c_s_entry_range) in pairs(c_s_entry_range_vec)

                subpath = work_dir * "/Experiments/$set_type $set_num/"
                filename = subpath * "Matrices $num_rows by $num_cols column prices index $(c_s_ind) greedy freq ranking $ranking.txt"

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
            master_dict["Greedy $ranking"] = df_stacked
        end
    else
        dfs = DataFrame[]
        for num_rows in num_rows_vec, num_cols in num_cols_vec, (c_s_ind, c_s_entry_range) in pairs(c_s_entry_range_vec)

            subpath = work_dir * "/Experiments/$set_type $set_num/"
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
        df_stacked = vcat(dfs...)
        master_dict["MILP"] = df_stacked
    end
end

master_dict["MILP"]
master_dict["Greedy with s"]
master_dict["Greedy with dual var"]


##################################################################################
##### Analysis merged DF (MILP vs one Greedy freq ranking)
df_temp_MILP = deepcopy(master_dict["MILP"])
# df_temp_greedy, g_suffix = deepcopy(master_dict["Greedy with s"]), "s"
df_temp_greedy, g_suffix = deepcopy(master_dict["Greedy with dual var"]), "dual"
select!(df_temp_MILP, Not(["Matrix seed", "Costs seed", "Node count"]))
select!(df_temp_greedy, Not(["Matrix seed", "Costs seed"]))

df_merged = deepcopy(df_temp_MILP)
rename!(df_merged, "Budget spent" => "Budget spent_MILP",
        "Obj val" => "Obj val_MILP", "Obj bound" => "Obj bound_MILP", "Dual obj" => "Dual obj_MILP",
        "Rows purchased" => "Rows purchased_MILP", "Cols purchased" => "Cols purchased_MILP",
        "Termination status" => "Term status_MILP", "Solve time" => "Solve time_MILP",
        "Relative gap" => "Rel gap_MILP"
        )

# Add cols from greedy freq
df_merged[:, "Budget spent_$g_suffix"] = df_temp_greedy[!, "Budget spent"]
df_merged[:, "Obj val_$g_suffix"] = df_temp_greedy[!, "Obj val"]
df_merged[:, "Rows purchased_$g_suffix"] = df_temp_greedy[!, "Rows purchased"]
df_merged[:, "Cols purchased_$g_suffix"] = df_temp_greedy[!, "Cols purchased"]

#### Compute comparison columns
# Counting
df_merged[:,"Obj val larger_MILP"] = df_merged[:,"Obj val_MILP"] .> df_merged[:,"Obj val_$g_suffix"]  # should always be true
df_merged[:,"Obj val larger or equal_$g_suffix"] = df_merged[:,"Obj val_MILP"] .<= df_merged[:,"Obj val_$g_suffix"]  # could be due to rel_gap solver setting
df_merged[:,"Obj val_MILP - Obj val_$g_suffix"] = df_merged[:,"Obj val_MILP"] - df_merged[:,"Obj val_$g_suffix"]
# TODO: Normalize difference in obj values
df_merged[:,"Budget spent smaller_MILP"] = df_merged[:,"Budget spent_MILP"] .< df_merged[:,"Budget spent_$g_suffix"]
# TODO: Add columns that count how many optimal rows and cols the Greedy approach purchased

minimum(df_merged[:,"Obj val_MILP - Obj val_$g_suffix"])

gdf = groupby(df_merged, ["c_s range", "Budget fraction", "Term status_MILP"])

gdf_agg = combine(gdf, nrow => "Group size", "Obj val larger_MILP" => count, "Budget spent smaller_MILP" => count,
        "Solve time_MILP" => mean, "Rel gap_MILP"
        )

gdf[("11:15", 0.5, "TIME_LIMIT")][:,"Obj val_MILP - Obj val_$g_suffix"]


##################################################################################
##### Analysis merged DF (Greedy freq s vs Greedy freq dual)
# df_temp_MILP = deepcopy(master_dict["MILP"])
df_temp_greedy_s = deepcopy(master_dict["Greedy with s"])
df_temp_greedy_dual = deepcopy(master_dict["Greedy with dual var"])
select!(df_temp_greedy_s, Not(["Matrix seed", "Costs seed"]))
select!(df_temp_greedy_dual, Not(["Matrix seed", "Costs seed"]))

df_merged = deepcopy(df_temp_greedy_s)
rename!(df_merged, "Obj val" => "Obj val_s", "Budget spent" => "Budget spent_s",
        "Rows purchased" => "Rows purchased_s", "Cols purchased" => "Cols purchased_s",
        )

# Add cols from greedy freq
df_merged[:,"Obj val_dual"] = df_temp_greedy_dual[!,"Obj val"]
df_merged[:,"Budget spent_dual"] = df_temp_greedy_dual[!,"Budget spent"]
df_merged[:,"Rows purchased_dual"] = df_temp_greedy_dual[!,"Rows purchased"]
df_merged[:,"Cols purchased_dual"] = df_temp_greedy_dual[!,"Cols purchased"]
#### Compute comparisons
df_merged[:,"Obj val larger_dual"] = df_merged[:,"Obj val_dual"] .> df_merged[:,"Obj val_s"]
df_merged[:,"Obj val larger_s"] = df_merged[:,"Obj val_dual"] .< df_merged[:,"Obj val_s"]
df_merged[:,"Obj vals equal"] = df_merged[:,"Obj val_dual"] .== df_merged[:,"Obj val_s"]
df_merged[:,"Obj val_dual - Obj val_s"] = df_merged[:,"Obj val_dual"] - df_merged[:,"Obj val_s"]
# TODO: Normalize difference in obj values
# TODO: Add columns that count how rows and cols each Greedy freq approach purchased that the other didn'table

# TODO: Careful with analysis among different matrix sizes (e.g., solve time average)
gdf = groupby(df_merged, ["c_s range", "Budget fraction"])

gdf_agg = combine(gdf, nrow => "Group size", "Obj val larger_dual" => count => "Dual obj larger", "Obj val larger_s" => count => "s obj larger",
        "Obj vals equal" => count => "Obj vals equal", "Obj val_dual - Obj val_s" => mean, "Budget spent_dual" => mean, "Budget spent_s" => mean,
        "Rows purchased_dual" => mean, "Rows purchased_s" => mean, "Cols purchased_dual" => mean, "Cols purchased_s" => mean,
        )

gdf_agg[:,1:7]
gdf_agg[:,8:end]


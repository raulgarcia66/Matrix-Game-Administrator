
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
        master_dict["Greedy MILP"] = df_stacked
    end
end

master_dict["MILP"]
master_dict["Greedy"]


##################################################################################
##### Analysis merged DF
df_temp_MILP = deepcopy(master_dict["MILP"])
df_temp_greedy = deepcopy(master_dict["Greedy"])
select!(df_temp_MILP, Not(["Matrix seed", "Costs seed"]))
select!(df_temp_greedy, Not(["Matrix seed", "Costs seed"]))

# df_merged = deepcopy(df_temp_cuts)
# rename!(df_merged, "Budget spent" => "Budget spent_cuts",
#         "Obj val" => "Obj val_cuts", "Obj bound" => "Obj bound_cuts", "Dual obj" => "Dual obj_cuts",
#         "Rows purchased" => "Rows purchased_cuts", "Cols purchased" => "Cols purchased_cuts",
#         "Termination status" => "Term status_cuts", "Solve time" => "Solve time_cuts",
#         "Relative gap" => "Rel gap_cuts", "Node count" => "Node count_cuts"
#         )

# Add cols from DF with no cuts
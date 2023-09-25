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
set_num = 1
subpath = work_dir * "/Experiments/$set_type $set_num/"

num_rows_vec = [10,50,100]
num_cols_vec = [10,50,100]
c_r_entry_range = 2:5
c_s_entry_range_vec = [2:5, 11:15, 21:25]

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

df = vcat(dfs...)
select!(df, Not(["Matrix seed", "Costs seed", "Obj bound", "Dual obj"]))
# gdf = groupby(df, ["Termination status", "c_s range", "Budget fraction"])  # 14 groups
# gdf = groupby(df, ["Termination status", "Num rows", "Num cols", "c_s range", "Budget fraction"])  # 83 groups
gdf = groupby(df, ["Termination status", "Num rows", "Num cols"])  # 12 groups
cdf = combine(gdf, nrow => "Group size", "Solve time" => mean, "Solve time" => std)



##################################################################################
##################################### Plots ######################################

#### Solution time vs budget fraction. Include a curve for each matrix size, column prices will be accumulated
gdf = groupby(df, ["Num rows", "Num cols", "Budget fraction"])  # "Termination status",
cdf = combine(gdf, nrow => "Group size", "Solve time" => mean, "Solve time" => std)
ggdf = groupby(cdf, ["Num rows", "Num cols"])  # "Termination status",
# ggdf[(100,100)]  # "OPTIMAL"

# For labels
budgets = [0.25, 0.50, 0.75]
ticklabel = ["0.25", "0.50", "0.75"]
# colors = ["green3", "red2", "chocolate", "gold2", "turquoise3", "magenta", "darkorchid4", "red4", "blue"]
colors = ["green3", "red2", "chocolate", "gold2", "turquoise3", "magenta", "darkorchid3", "red4", "blue"]

p = plot(xlabel="Budget proportion", ylabel="Solution time (s)",legend=:outertopright); # title = "MGD MILP Formulation"
counter = 1
for num_rows in num_rows_vec, num_cols in num_cols_vec

    pdf = ggdf[(num_rows,num_cols)]  # "OPTIMAL"
    plot!(pdf[:, "Budget fraction"], pdf[:,"Solve time_mean"], xticks = budgets, # (budgets, ticklabel),
        xlabel="Budget proportion", ylabel="Solution time (s)", labels = "$(num_rows) × $(num_cols)",
        lw=2, yaxis=:log, markershape=:circle, ms=5, color=colors[counter]
    )
    counter += 1
end
display(p)
png(joinpath(work_dir, "Experiments", "Plots", "MGD MILP Formulation Budget Fraction vs Solve Time per Matrix Size"))
# savefig(p, joinpath(work_dir, "Experiments", "Plots", "MGD MILP Formulation Budget Fraction vs Solve Time per Matrix Size.png"))

#### Another plot showing effect of column prices on solution time, per budget fraction
gdf = groupby(df, ["Num rows", "Num cols", "Budget fraction", "c_s range"])  # "Termination status"
cdf = combine(gdf, nrow => "Group size", "Solve time" => mean, "Solve time" => std)
ggdf = groupby(cdf, ["Num rows", "Num cols", "Budget fraction"])  # "Termination status"
# ggdf[(100,10,0.75)]  # "OPTIMAL"

# # Extract data matrix
# data_matrix = zeros(length(num_rows_vec) * length(num_cols_vec), 3)
# # for b in budgets
# dim_counter = 0
# for num_rows in num_rows_vec, num_cols in num_cols_vec
#     dim_counter += 1

#     b = 0.50
#     pdf = ggdf[(num_rows, num_cols, b)]

#     for (i, c_s) in enumerate(pdf[:,"c_s range"])
#         if c_s == "2:5"
#             data_matrix[dim_counter,1] = pdf[i, "Solve time_mean"]
#         elseif c_s == "11:15"
#             data_matrix[dim_counter,2] = pdf[i, "Solve time_mean"]
#         elseif c_s == "21:25"
#             data_matrix[dim_counter,3] = pdf[i, "Solve time_mean"]
#         end
#     end
# end
# # end
# data_matrix


for b in budgets

    data_matrix = zeros(length(num_rows_vec) * length(num_cols_vec), 3)  # 3 is for number of c_s prices

    ticklabel = []
    dim_counter = 0
    for num_rows in num_rows_vec, num_cols in num_cols_vec
        dim_counter += 1
        push!(ticklabel, "$num_rows×$num_cols")

        pdf = try
            ggdf[(num_rows, num_cols, b)]  # "OPTIMAL"
        catch e
            println("$e")
            continue
        end

        for (i, c_s) in enumerate(pdf[:,"c_s range"])
            if c_s == "2:5"
                data_matrix[dim_counter,1] = pdf[i, "Solve time_mean"]
            elseif c_s == "11:15"
                data_matrix[dim_counter,2] = pdf[i, "Solve time_mean"]
            elseif c_s == "21:25"
                data_matrix[dim_counter,3] = pdf[i, "Solve time_mean"]
            end
        end
    end

    ctg = repeat(["1-LOW", "2-MEDIUM", "3-HIGH"], inner = 9)  # Bars are arranged in alphabetical order

    p = groupedbar(data_matrix, bar_position = :dodge, bar_width=0.8,
        title = "Budget Proportion: $b", xlabel="Matrix sizes", ylabel="Solution time (s)", legend=:topleft,
        group = ctg, # labels =["LOW", "MEDIUM", "HIGH"], 
        yaxis=:log, xticks=(1:length(ticklabel), ticklabel), lw = 0
    )
    display(p)
    png(joinpath(work_dir, "Experiments", "Plots", "MGD MILP Formulation Matrix Size vs Solution Time per Column Price with Budget $(round(b, digits=2))"))
end

##################################################################################
############################### Greedy vs Rank ###################################

work_dir = pwd()
# set_type = "Set Greedy MILP"
set_num = 1
# subpath = work_dir * "/Experiments/$set_type $set_num/"

num_rows_vec = [10,50,100]
num_cols_vec = [10,50,100]
c_r_entry_range = 2:5
c_s_entry_range_vec = [2:5, 11:15, 21:25]
set_type_vec = ["Set MILP", "Set Greedy MILP", "Set Greedy Freq"]

# Load Greedy MILP and Greedy Freq Dual into a Dictionary
master_dict = Dict{String,DataFrame}()

for set_type in set_type_vec
    dfs = DataFrame[]
    for num_rows in num_rows_vec, num_cols in num_cols_vec, (c_s_ind, c_s_entry_range) in pairs(c_s_entry_range_vec)

        subpath = work_dir * "/Experiments/$set_type $set_num/"
        filename = ""
        if set_type == "Set MILP"
            filename = subpath * "Matrices $num_rows by $num_cols column prices index $(c_s_ind) MILP.txt"
        elseif set_type == "Set Greedy MILP"
            filename = subpath * "Matrices $num_rows by $num_cols column prices index $(c_s_ind) greedy MILP.txt"
        else
            filename = subpath * "Matrices $num_rows by $num_cols column prices index $(c_s_ind) greedy freq ranking with dual var.txt"
        end

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
        master_dict["M"] = df_stacked
    elseif set_type == "Set Greedy MILP"
        master_dict["G"] = df_stacked
    elseif set_type == "Set Greedy Freq"
        master_dict["R"] = df_stacked
    end
end

master_dict["M"]
master_dict["G"]
master_dict["R"]

#### Analysis merged DF
df_temp_m = deepcopy(master_dict["M"])
df_temp_g = deepcopy(master_dict["G"])
df_temp_r = deepcopy(master_dict["R"])
select!(df_temp_m, Not(["Matrix seed", "Costs seed", "Obj bound", "Dual obj"]))

df_merged = deepcopy(df_temp_m)
rename!(df_merged, "Budget spent" => "Budget spent_m",
        "Obj val" => "Obj val_m", "Termination status" => "Term status_m", "Solve time" => "Solve time_m",
        "Rows purchased" => "Rows purchased_m", "Cols purchased" => "Cols purchased_m",
        "Relative gap" => "Rel gap_m", "Node count" => "Node count_m"
        )

# Add cols from Greedy MILP
df_merged[:,"Budget spent_g"] = df_temp_g[!,"Budget spent"]
df_merged[:,"Obj val_g"] = df_temp_g[!,"Obj val"]
df_merged[:,"Rows purchased_g"] = df_temp_g[!,"Rows purchased"]
df_merged[:,"Cols purchased_g"] = df_temp_g[!,"Cols purchased"]
df_merged[:,"Term status_g"] = df_temp_g[!,"Termination status"]
df_merged[:,"Solve time_g"] = df_temp_g[!,"Solve time"]
# Add cols from Greedy Freq Dual
df_merged[:, "Budget spent_r"] = df_temp_r[!, "Budget spent"]
df_merged[:, "Obj val_r"] = df_temp_r[!, "Obj val"]
df_merged[:, "Rows purchased_r"] = df_temp_r[!, "Rows purchased"]
df_merged[:, "Cols purchased_r"] = df_temp_r[!, "Cols purchased"]

#### Compute comparison columns
# Greedy MILP
# df_merged[:,"Obj val rel gap_g"] = 1 .- ((df_merged[:,"Obj val_m"] - df_merged[:,"Obj val_g"] ) ./ abs.(df_merged[:,"Obj val_m"]))
df_merged[:,"Obj val prop_g"] = df_merged[:,"Obj val_g"] ./ df_merged[:,"Obj val_m"]
# Greedy Freq Dual
# df_merged[:,"Obj val rel gap_r"] = 1 .- ((df_merged[:,"Obj val_m"] - df_merged[:,"Obj val_r"] ) ./ abs.(df_merged[:,"Obj val_m"]))
df_merged[:,"Obj val prop_r"] = df_merged[:,"Obj val_r"] ./ df_merged[:,"Obj val_m"]

val_max_g, i_max_g = findmax(df_merged[:,"Obj val prop_g"])
val_max_r, i_max_r = findmax(df_merged[:,"Obj val prop_r"])

val_min_g, i_min_g = findmin(df_merged[:,"Obj val prop_g"])
val_min_r, i_min_r = findmin(df_merged[:,"Obj val prop_r"])

df_merged[i_min_g,"Obj val_m"]
df_merged[i_min_g,"Obj val_g"]
df_merged[i_max_g,"Obj val_m"]
df_merged[i_max_g,"Obj val_g"]

df_merged[i_min_r,"Obj val_m"]
df_merged[i_min_r,"Obj val_r"]
df_merged[i_max_r,"Obj val_m"]
df_merged[i_max_r,"Obj val_r"]

# count(e -> e < 0, df_merged[!,"Obj val prop_r"])  # Just one entry is negative (at 51)
# count(e -> e > 1, df_merged[!,"Obj val prop_g"])  # Two entries 

# Fix negative entries
for (i, e) in enumerate(df_merged[:,"Obj val prop_r"])
    if e < 0
        # If e is < 0, the objs must be opposite signs and we know MILP is larger
        df_merged[i,"Obj val prop_r"] = abs(df_merged[i, "Obj val_r"]) / (df_merged[i,"Obj val_m"] + abs(df_merged[i, "Obj val_r"]))
    end
end

# Fix values > 1
for (i, e) in enumerate(df_merged[:,"Obj val prop_g"])
    if e > 1
        # Fix MIPGap + numerical error
        df_merged[i,"Obj val prop_g"] = 1.0
    end
end

#### Rel gap of Greedy MILP and Greedy Freq Dual, per price level and budget fraction
gdf = groupby(df_merged, ["c_s range", "Budget fraction"])
cdf = combine(gdf, nrow => "Group size", "Obj val prop_g" => mean, "Obj val prop_r" => mean)
ggdf = groupby(cdf, ["c_s range", "Budget fraction"])  # "Termination status"
# ggdf[(100,10,0.75)]  # "OPTIMAL"

# For labels
budgets = [0.25, 0.50, 0.75]

data_matrix = zeros(length(num_rows_vec) * length(num_cols_vec), 2)  # 2 is for the two methods, G and R

ticklabel = Vector{String}(undef,9)
dim_counter = 10
for c_s in ["2:5", "11:15", "21:25"], b in budgets
    dim_counter -= 1

    data_matrix[dim_counter,1] = ggdf[(c_s, b)][1,"Obj val prop_g_mean"]
    data_matrix[dim_counter,2] = ggdf[(c_s, b)][1,"Obj val prop_r_mean"]

    if c_s == "2:5"
        # push!(ticklabel, "L×$(round(b, digits=2))")
        ticklabel[dim_counter] = "L×$(round(b, digits=2))"
    elseif c_s == "11:15"
        # push!(ticklabel, "M×$(round(b, digits=2))")
        ticklabel[dim_counter] = "M×$(round(b, digits=2))"
    elseif c_s == "21:25"
        # push!(ticklabel, "H×$(round(b, digits=2))")
        ticklabel[dim_counter] = "H×$(round(b, digits=2))"
    end
end

ctg = repeat(["H_G", "H_R"], inner = 9)  # Bars are arranged in alphabetical order

p = groupedbar(data_matrix, bar_position = :dodge, bar_width=0.6,
    guidefont = font(11,"Dejavu Sans"), # title = "Ratio with MILP objective", 
    xlabel="Column Price Level × Budget Proportion", ylabel="Ratio with MILP Obj",
    group = ctg, ylims = (0,1), legend=:topleft,
    xticks=(1:length(ticklabel), ticklabel), lw = 1.4,
    color = ["grey" "white"]
)
display(p)
png(joinpath(work_dir, "test", "Greedy vs Ranking Ratios With MILP Objective H_G H_R black and white"))
# png(joinpath(work_dir, "Experiments", "Plots", "Greedy vs Ranking Ratios With MILP Objective H_G H_R"))
# png(joinpath(work_dir, "Experiments", "Plots", "Greedy vs Ranking Ratios With MILP Objective H_G H_R with border"))

# data = [[0.745863, 0.222339]
#     [0.808621, 0.343113]
#     [0.785822, 0.488157]
#     [0.721516, 0.201917]
#     [0.735018, 0.301307]
#     [0.772852, 0.47571]
#     [0.674362, 0.179662]
#     [0.626785, 0.213839]
#     [0.623989, 0.330781]]
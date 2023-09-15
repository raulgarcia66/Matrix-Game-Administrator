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
ggdf[(100,100)]  # "OPTIMAL"

# For labels
budgets = [0.25, 0.50, 0.75]
ticklabel = ["0.25", "0.50", "0.75"]
# colors = ["green3", "red2", "magenta", "gold2", "turquoise3", "red4", "darkorchid4", "chocolate", "blue"]
colors = ["green3", "red2", "chocolate", "gold2", "turquoise3", "magenta", "darkorchid4", "red4", "blue"]

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
png("MGD MILP Formulation Budget Fraction vs Solve Time per Matrix Size")

#### Another plot showing effect of column prices on solution time, per budget fraction
gdf = groupby(df, ["Num rows", "Num cols", "Budget fraction", "c_s range"])  # "Termination status"
cdf = combine(gdf, nrow => "Group size", "Solve time" => mean, "Solve time" => std)
ggdf = groupby(cdf, ["Num rows", "Num cols", "Budget fraction"])  # "Termination status"
ggdf[(100,10,0.75)]  # "OPTIMAL"

for b in budgets

    data_matrix = zeros(length(num_rows_vec) * length(num_cols_vec), 3)  # 3 is for number of c_s prices

    ticklabel = []
    dim_counter = 0
    for num_rows in num_rows_vec, num_cols in num_cols_vec
        dim_counter += 1
        push!(ticklabel, "$num_rows × $num_cols")

        try
            pdf = ggdf[(num_rows, num_cols, b)]  # "OPTIMAL"
        catch
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
    png("MGD MILP Formulation Matrix Size vs Solution Time per Column Price with Budget $(round(b, digits=2))")
end

# TODO: Updates colors

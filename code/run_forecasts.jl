# ==============================================================================
# Forecasting Script: run_forecasts.jl
# Implements and compares methods 1a, 1b, 1c, 2, and 3 for parameter forecasting.
# Historical region: piecewise constant parameters
# Forecast region: extrapolation strategy
# ==============================================================================

include("opt_vars.jl")
include("get_data_new.jl")

using DelimitedFiles
using LinearAlgebra
using Statistics
using PositiveIntegrators
using PyPlot
using Trapz
using Printf

println("--- Starting Forecast Pipeline ---")

_local_suffix = isdefined(Main, :forecast_suffix) ? Main.forecast_suffix : ""

# ==============================================================================
# Setup forecasting horizon
# ==============================================================================
base_p = all_models[Model_no]
if !@isdefined(T_forecast)
    T_forecast = 7.0
end
if !@isdefined(data_bound)
    data_bound = 74.0
end

plot_top_M = 3 # Number of best methods to plot (<= 0 for all)
# ==============================================================================
# 1. Load optimized parameters from CSV backups
# ==============================================================================

println("Loading optimization data from backups...")

if !@isdefined(BACKUP_ID)
    BACKUP_ID = "backupI"
end

midpoints_data = vec(DelimitedFiles.readdlm("midpoints_$(BACKUP_ID).csv", ',', Float64))
param_matrix_data = DelimitedFiles.readdlm("param_matrix_$(BACKUP_ID).csv", ',', Float64)

# Store full (unfiltered) data for the piecewise constant reference solution
midpoints_full = copy(midpoints_data)
param_matrix_full = copy(param_matrix_data)

N_cells_full = length(midpoints_full)
X_interfaces_full = zeros(N_cells_full + 1)
X_interfaces_full[1] = 0.0
for i in 1:N_cells_full
    X_interfaces_full[i+1] = 2.0 * midpoints_full[i] - X_interfaces_full[i]
end
h_widths_full = diff(X_interfaces_full)

# Find start cell using data_bound as the start time
start_cell = findfirst(x -> x >= data_bound, X_interfaces_full)
if start_cell === nothing || start_cell == 1
    start_cell = N_cells_full + 1
end

t_end = X_interfaces_full[start_cell]
t_max_sim = t_end + T_forecast

num_opt_vars = size(param_matrix_full, 1)

println("Optimization period ends at t = $t_end")
println("Forecasting up to t = $t_max_sim")

# ==============================================================================
# 3. Import WENO and Predictor
# ==============================================================================

include("weno.jl")
using .WenoInterpolation: weno_evaluate_non_uniform

include("ParameterPredictor.jl")
using .ParameterPredictor: predict_single_step

# ==============================================================================
# 4. Reused functions from run_rolling_forecast.jl
# ==============================================================================

function p_eval_piecewise_constant(
    t,
    param_matrix,
    X_interfaces
)
    N_cells = size(param_matrix, 2)
    idx = searchsortedlast(X_interfaces, t)
    idx = clamp(idx, 1, N_cells)
    return param_matrix[:, idx]
end

function build_forecast_methods(
    start_cell,
    t_start,
    t_forecast_end,
    param_matrix,
    X_interfaces,
    h_widths,
    num_opt_vars
)
    history_idx = 1:(start_cell-1)
    interface_idx = 1:start_cell

    X_interfaces_hist = X_interfaces[interface_idx]
    h_widths_hist = h_widths[history_idx]

    local_itp_weno_k2_const = []
    local_itp_weno_k2_lin = []
    local_itp_weno_k3_const = []
    local_itp_weno_k3_lin = []

    for j in 1:num_opt_vars
        y_values_hist = param_matrix[j, history_idx]
        ub_val = isdefined(Main, :bounds) ? Main.bounds[2, j] : Inf

        push!(local_itp_weno_k2_const, t -> weno_evaluate_non_uniform(t, y_values_hist, h_widths_hist, X_interfaces_hist, 2; apply_zhang_shu=true, lower_bound=0.0, upper_bound=ub_val, extrapolation=:constant))
        push!(local_itp_weno_k2_lin,   t -> weno_evaluate_non_uniform(t, y_values_hist, h_widths_hist, X_interfaces_hist, 2; apply_zhang_shu=true, lower_bound=0.0, upper_bound=ub_val, extrapolation=:linear))
        push!(local_itp_weno_k3_const, t -> weno_evaluate_non_uniform(t, y_values_hist, h_widths_hist, X_interfaces_hist, 3; apply_zhang_shu=true, lower_bound=0.0, upper_bound=ub_val, extrapolation=:constant))
        push!(local_itp_weno_k3_lin,   t -> weno_evaluate_non_uniform(t, y_values_hist, h_widths_hist, X_interfaces_hist, 3; apply_zhang_shu=true, lower_bound=0.0, upper_bound=ub_val, extrapolation=:linear))
    end

    params_1b = Dict{Int, Dict{String, Vector{Float64}}}()
    for d in 1:Int(T_forecast)
        t_upper = min(t_start + Float64(d), t_forecast_end)
        
        p_k2_c = zeros(num_opt_vars)
        p_k2_l = zeros(num_opt_vars)
        p_k3_c = zeros(num_opt_vars)
        p_k3_l = zeros(num_opt_vars)

        for j in 1:num_opt_vars
            p_k2_c[j] = mean([local_itp_weno_k2_const[j](τ) for τ in range(t_start, t_upper, length=50)])
            p_k2_l[j] = mean([local_itp_weno_k2_lin[j](τ)   for τ in range(t_start, t_upper, length=50)])
            p_k3_c[j] = mean([local_itp_weno_k3_const[j](τ) for τ in range(t_start, t_upper, length=50)])
            p_k3_l[j] = mean([local_itp_weno_k3_lin[j](τ)   for τ in range(t_start, t_upper, length=50)])
        end
        
        params_1b[d] = Dict("k2_const" => p_k2_c, "k2_lin" => p_k2_l, "k3_const" => p_k3_c, "k3_lin" => p_k3_l)
    end

    params_2 = p_eval_piecewise_constant(t_start, param_matrix, X_interfaces)

    params_3 = try
        p_3, _ = predict_single_step(param_matrix, h_widths, start_cell - 1, 1, 1.0, 1.0)
        p_3
    catch
        params_2
    end

    methods_dict = Dict{String, Any}(
        "1a) k=2 const" => t -> t <= t_start ? p_eval_piecewise_constant(t, param_matrix, X_interfaces) : [local_itp_weno_k2_const[j](t) for j in 1:num_opt_vars],
        "1a) k=2 lin"   => t -> t <= t_start ? p_eval_piecewise_constant(t, param_matrix, X_interfaces) : [local_itp_weno_k2_lin[j](t) for j in 1:num_opt_vars],
        "1a) k=3 const" => t -> t <= t_start ? p_eval_piecewise_constant(t, param_matrix, X_interfaces) : [local_itp_weno_k3_const[j](t) for j in 1:num_opt_vars],
        "1a) k=3 lin"   => t -> t <= t_start ? p_eval_piecewise_constant(t, param_matrix, X_interfaces) : [local_itp_weno_k3_lin[j](t) for j in 1:num_opt_vars],
        "2) Current Cell Avg" => t -> t <= t_start ? p_eval_piecewise_constant(t, param_matrix, X_interfaces) : params_2,
        "3) Hist. Predictor" => t -> t <= t_start ? p_eval_piecewise_constant(t, param_matrix, X_interfaces) : params_3
    )

    for d in 1:Int(T_forecast)
        methods_dict["1b) $(d)-Day Avg k=2 const"] = t -> t <= t_start ? p_eval_piecewise_constant(t, param_matrix, X_interfaces) : params_1b[d]["k2_const"]
        methods_dict["1b) $(d)-Day Avg k=2 lin"]   = t -> t <= t_start ? p_eval_piecewise_constant(t, param_matrix, X_interfaces) : params_1b[d]["k2_lin"]
        methods_dict["1b) $(d)-Day Avg k=3 const"] = t -> t <= t_start ? p_eval_piecewise_constant(t, param_matrix, X_interfaces) : params_1b[d]["k3_const"]
        methods_dict["1b) $(d)-Day Avg k=3 lin"]   = t -> t <= t_start ? p_eval_piecewise_constant(t, param_matrix, X_interfaces) : params_1b[d]["k3_lin"]
    end

    return methods_dict
end

function P_non_autonomous(u, p_static, t, p_eval)
    p_dynamic = copy(p_static)
    p_vals = p_eval(t)
    for j in 1:num_opt_vars
        p_dynamic[opt_positions_filtered[j]] = p_vals[j]
    end
    return P(u, Tuple(p_dynamic), t)
end

function d_non_autonomous(u, p_static, t, p_eval)
    p_dynamic = copy(p_static)
    p_vals = p_eval(t)
    for j in 1:num_opt_vars
        p_dynamic[opt_positions_filtered[j]] = p_vals[j]
    end
    return d(u, Tuple(p_dynamic), t)
end

# ==============================================================================
# 5. Solve with full piecewise constant parameters (numerical reference)
# ==============================================================================

println("Solving ODE with full piecewise constant parameters (numerical reference)...")

p_eval_full_piecewise(t) = p_eval_piecewise_constant(t, param_matrix_full, X_interfaces_full)

P_const_pw(u, p_static, t) = P_non_autonomous(u, p_static, t, p_eval_full_piecewise)
d_const_pw(u, p_static, t) = d_non_autonomous(u, p_static, t, p_eval_full_piecewise)

prob_const_pw = PDSProblem(P_const_pw, d_const_pw, original_u0, (0.0, t_max_sim), base_p)
sol_const_pw = solve(prob_const_pw, method)

println("Full piecewise constant reference solution complete.")

# ==============================================================================
# 6. Build methods and solve forecast PDS problems from t_end
# ==============================================================================

methods_dict = build_forecast_methods(
    start_cell,
    t_end,
    t_max_sim,
    param_matrix_full,
    X_interfaces_full,
    h_widths_full,
    num_opt_vars
)


solutions = Dict()

# Extract initial condition for the forecast period from the reference solution
u0_forecast = sol_const_pw(t_end)

for (label, p_eval) in methods_dict
    println("Solving ODE using $label...")
    P_temp(u, p_static, t) = P_non_autonomous(u, p_static, t, p_eval)
    d_temp(u, p_static, t) = d_non_autonomous(u, p_static, t, p_eval)

    prob = PDSProblem(P_temp, d_temp, u0_forecast, (t_end, t_max_sim), base_p)
    sol = solve(prob, method)
    solutions[label] = sol
end

println("ODE solving complete!")

# 7. Visualization: Compare compartment dynamics (Infected I and Deaths D_covid)
println("Plotting compartment comparisons...")

comp_labels = ["S", "V", "E", "L", "I", "H", "R", "Q", "D", "Ve", "D_covid"]
idx_I = 5
idx_D = 11

# Helper function to get compartment trajectory from solution
function get_compartment_traj(sol, p_eval, idx)
    t_vals = sol.t
    y_vals = Float64[]
    for k in eachindex(t_vals)
        tk = t_vals[k]
        uk = sol.u[k]
        
        p_dynamic = copy(base_p)
        p_vals = p_eval(tk)
        for j in 1:num_opt_vars
            p_dynamic[opt_positions_filtered[j]] = p_vals[j]
        end
        
        full_res = Main.getplotvar([tk], [uk], p_dynamic, idx)
        push!(y_vals, full_res[1])
    end
    return y_vals
end

# Reference data from Main session
tref_I = isdefined(Main, :tref_I) ? Main.tref_I : []
refI   = isdefined(Main, :refI) ? Main.refI : []
tref_D = isdefined(Main, :tref_D) ? Main.tref_D : []
refD   = isdefined(Main, :refD) ? Main.refD : []

function compute_forecast_cost(sol, p_eval, t_start_eval::Float64, t_end_eval::Float64, base_p, opt_positions_filtered, num_opt_vars)
    tref_I = isdefined(Main, :tref_I) ? Main.tref_I : []
    refI   = isdefined(Main, :refI) ? Main.refI : []
    tref_D = isdefined(Main, :tref_D) ? Main.tref_D : []
    refD   = isdefined(Main, :refD) ? Main.refD : []
    
    flag_I = isdefined(Main, :flag_I) ? Main.flag_I : true
    flag_D = isdefined(Main, :flag_D) ? Main.flag_D : true
    
    weights = ones(11)
    err = zeros(11)
    
    for i in 1:11
        tref = sol.t
        solref = getindex.(sol.u, i)
        
        if i == 5 && flag_I
            tref = tref_I; solref = refI
        elseif i == 11 && flag_D
            tref = tref_D; solref = refD
        else
            continue
        end
        
        valid_indices = findall(x -> x >= t_start_eval && x <= t_end_eval, tref)
        if isempty(valid_indices)
            continue
        end
        tref_filtered = tref[valid_indices]
        solref_filtered = solref[valid_indices]
        
        num_sol = Float64[]
        for tk in tref_filtered
            uk = sol(tk)
            p_dynamic = copy(base_p)
            p_vals = p_eval(tk)
            for j in 1:num_opt_vars
                p_dynamic[opt_positions_filtered[j]] = p_vals[j]
            end
            val = Main.getplotvar([tk], [uk], p_dynamic, i)[1]
            push!(num_sol, val)
        end 
        
        aux = num_sol .- solref_filtered
        denominator = max.(solref_filtered, 1.0)
        err[i] = weights[i] * (norm(aux) / norm(denominator))
        if isnan(err[i]) || isinf(err[i])
            err[i] = 1e16
        end
    end
    return sum(err)
end

if (plot_top_M > 0)
    println("Determining top $plot_top_M methods based on average cost...")
    method_avg_costs = Dict{String, Float64}()
    num_days = Int(T_forecast)
    for (label, p_eval) in methods_dict
        sol = solutions[label]
        day_costs = Float64[]
        for d in 1:num_days
            push!(day_costs, compute_forecast_cost(sol, p_eval, t_end, t_end + Float64(d), base_p, opt_positions_filtered, num_opt_vars))
        end
        method_avg_costs[label] = mean(day_costs)
    end
    sorted_methods = sort(collect(method_avg_costs), by = x -> x[2])
    methods_to_plot = [x[1] for x in sorted_methods[1:min(plot_top_M, length(sorted_methods))]]
    println("Top $plot_top_M methods to plot: ", methods_to_plot)
else
    methods_to_plot = collect(keys(methods_dict))
end

# Create Figure 1
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))

# ==============================================================================
# Rank based coloring and styling for forecast methods
# ==============================================================================

using PyPlot

styles = Dict{String, Tuple{Any,String,String}}()

rank_cmap = plt.get_cmap("tab10")

rank_linestyles = ["-", "--", "-.", ":"]

rank_markers = [
    "o", "s", "^", "D",
    "v", "P", "X", "*",
    "h", "p"
]

for (rank, label) in enumerate(methods_to_plot)

    c =
        rank_cmap(
            (rank - 1) /
            max(length(methods_to_plot) - 1, 1)
        )

    ls =
        rank_linestyles[
            mod1(rank, length(rank_linestyles))
        ]

    mk =
        rank_markers[
            mod1(rank, length(rank_markers))
        ]

    styles[label] = (c, ls, mk)

end

# Plot full piecewise constant numerical reference solution (plotted first/underneath with zorder=2)
y_I_const = get_compartment_traj(sol_const_pw, p_eval_full_piecewise, idx_I)
y_D_const = get_compartment_traj(sol_const_pw, p_eval_full_piecewise, idx_D)
ax1.plot(sol_const_pw.t, y_I_const, color="black", linestyle="-", linewidth=2.5, label="Numerical (pw. const.)", zorder=2)
ax2.plot(sol_const_pw.t, y_D_const, color="black", linestyle="-", linewidth=2.5, label="Numerical (pw. const.)", zorder=2)

for label in methods_to_plot
    sol = solutions[label]
    t = sol.t
    p_eval = methods_dict[label]
    y_I = get_compartment_traj(sol, p_eval, idx_I)
    y_D = get_compartment_traj(sol, p_eval, idx_D)
    
    col, ls, mk = styles[label]
    me = max(1, div(length(t), 20))  # show ~20 markers along the line
    
    ax1.plot(t, y_I, color=col, linestyle=ls, linewidth=2.0, label=label,
             marker=mk, markevery=me, markersize=5, markeredgewidth=0.5, markeredgecolor="white", zorder=3)
    ax2.plot(t, y_D, color=col, linestyle=ls, linewidth=2.0, label=label,
             marker=mk, markevery=me, markersize=5, markeredgewidth=0.5, markeredgecolor="white", zorder=3)
end


# Plot Reference Data if available (split into optimization period and forecast period)
if !isempty(tref_I) && !isempty(refI)
    valid_I_hist = findall(x -> x <= t_end, tref_I)
    valid_I_fore = findall(x -> t_end < x <= t_max_sim, tref_I)
    
    if !isempty(valid_I_hist)
        ax1.scatter(tref_I[valid_I_hist], refI[valid_I_hist], color="#546e7a", s=20, alpha=0.5, label="Data (Historical)", zorder=5)
    end
    if !isempty(valid_I_fore)
        ax1.scatter(tref_I[valid_I_fore], refI[valid_I_fore], color="#d32f2f", marker="o", s=35, label="Data (Forecast)", zorder=5)
    end
end

if !isempty(tref_D) && !isempty(refD)
    valid_D_hist = findall(x -> x <= t_end, tref_D)
    valid_D_fore = findall(x -> t_end < x <= t_max_sim, tref_D)
    
    if !isempty(valid_D_hist)
        ax2.scatter(tref_D[valid_D_hist], refD[valid_D_hist], color="#546e7a", s=20, alpha=0.5, label="Data (Historical)", zorder=5)
    end
    if !isempty(valid_D_fore)
        ax2.scatter(tref_D[valid_D_fore], refD[valid_D_fore], color="#d32f2f", marker="o", s=35, label="Data (Forecast)", zorder=5)
    end
end

# Draw vertical forecast threshold line
ax1.axvline(x=t_end, color="black", linestyle="--", linewidth=1.5, alpha=0.8)
ax1.text(t_end + 1.0, ax1.get_ylim()[2] * 0.9, " ", fontsize=10, fontweight="bold", zorder=10)

ax2.axvline(x=t_end, color="black", linestyle="--", linewidth=1.5, alpha=0.8)
ax2.text(t_end + 1.0, ax2.get_ylim()[2] * 0.9, " ", fontsize=10, fontweight="bold", zorder=10)

# Styling ax1
ax1.set_title("Infected Population (I) Forecast Comparison", fontsize=14, fontweight="bold")
ax1.set_xlabel("Days t", fontsize=12)
ax1.set_ylabel("Infected", fontsize=12)
ax1.grid(true, linestyle=":", alpha=0.5)
ax1.legend(loc="upper left")

# Styling ax2
ax2.set_title("COVID Deaths (D_covid) Forecast Comparison", fontsize=14, fontweight="bold")
ax2.set_xlabel("Days t", fontsize=12)
ax2.set_ylabel("Deaths", fontsize=12)
ax2.grid(true, linestyle=":", alpha=0.5)
ax2.legend(loc="upper left")

plt.tight_layout()
_fc_base1 = "forecast_comparison_$(BACKUP_ID)_db$(data_bound)_T$(T_forecast)$(_local_suffix)"
plt.savefig("$(_fc_base1).png", dpi=300)
plt.savefig("$(_fc_base1).eps", format="eps")
plt.savefig("$(_fc_base1).pdf", format="pdf")
if isdefined(Main, :AntigravityHeadless)
    println("Headless: Saved forecast_comparison_$(BACKUP_ID)_db$(data_bound)_T$(T_forecast).png")
else
    plt.show()
end


    
println("--- Forecast Pipeline Complete ---")

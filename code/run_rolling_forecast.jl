# ==============================================================================
# Rolling Forecast Evaluation
# ==============================================================================

include("opt_vars.jl")
include("get_data_new.jl")

using DelimitedFiles
using LinearAlgebra
using Statistics
using PositiveIntegrators
using PyPlot

println("------------------------------------------------------------")
println("Starting Rolling Forecast Evaluation")
println("------------------------------------------------------------")

# ==============================================================================
# SETTINGS
# ==============================================================================
T_forecasts = [2.0, 3.0, 5.0, 7.0, 10.0] # forecast horizons in days
t_final = 180.0 # total time span for forecasts (should be covered by data)
plot_top_M = 3 # Number of best methods to plot (set to <= 0 to plot all)

if !@isdefined(BACKUP_ID)
    BACKUP_ID = "backupI"
end

# Collect all (data_bound, T_forecast, cost, method) tuples across every T and method
all_real_data_points = Vector{NamedTuple{(:data_bound, :T_forecast, :cost, :method), Tuple{Float64, Float64, Float64, String}}}()

for T_forecast in T_forecasts
    println("Running forecasts with horizon = ", T_forecast, " days")


    # earliest forecast start:
    min_forecast_cell = 5

    # ==============================================================================
    # LOAD DATA
    # ==============================================================================

    midpoints = vec(
        DelimitedFiles.readdlm(
            "midpoints_$(BACKUP_ID).csv",
            ',',
            Float64
        )
    )

    param_matrix = DelimitedFiles.readdlm(
        "param_matrix_$(BACKUP_ID).csv",
        ',',
        Float64
    )

    N_cells = length(midpoints)
    num_opt_vars = size(param_matrix, 1)

    # ==============================================================================
    # RECONSTRUCT GRID
    # ==============================================================================

    if N_cells == 1
        error(
            "N_cells == 1: rolling forecast requires at least 2 optimisation " *
            "intervals to have a meaningful history. Aborting."
        )
    end

    X_interfaces = zeros(N_cells + 1)

    X_interfaces[1] = 0.0

    for i in 1:N_cells
        X_interfaces[i+1] =
            2.0 * midpoints[i] - X_interfaces[i]
    end


    h_widths = diff(X_interfaces)

  

   

    println("Total optimization horizon = ", t_final)

    # ==============================================================================
    # IMPORT WENO + PREDICTOR
    # ==============================================================================

    include("weno.jl")
    using .WenoInterpolation: weno_evaluate_non_uniform

    include("ParameterPredictor.jl")
    using .ParameterPredictor: predict_single_step

    # ==============================================================================
    # PIECEWISE CONSTANT PARAMETER EVALUATION
    # ==============================================================================

    function p_eval_piecewise_constant(
        t,
        p_matrix,
        x_interfaces
    )

        n_cells = size(p_matrix, 2)

        idx = searchsortedlast(x_interfaces, t)

        idx = clamp(idx, 1, n_cells)

        return p_matrix[:, idx]

    end

    # ==============================================================================
    # PARAMETER FORECAST METHODS (DYNAMIC HISTORICAL WINDOW)
    # ==============================================================================

    function build_forecast_methods(
        start_cell,
        t_start,
        t_forecast_end,
        param_matrix,
        X_interfaces,
        h_widths,
        num_opt_vars
    )
        # Restrict grid and boundaries to the current rolling forecast start point
        # so that evaluation times t > t_start correctly trigger extrapolation.
        interface_idx = 1:start_cell
        history_idx = 1:(start_cell-1)
        

        X_interfaces_hist = X_interfaces[interface_idx]
        h_widths_hist = h_widths[history_idx]

        local_itp_weno_k2_const = []
        local_itp_weno_k2_lin = []
        local_itp_weno_k3_const = []
        local_itp_weno_k3_lin = []

        for j in 1:num_opt_vars

            y_values_hist = param_matrix[j, history_idx]

            ub_val =
                isdefined(Main, :bounds) ?
                Main.bounds[2, j] :
                Inf

            push!(
                local_itp_weno_k2_const,
                t -> weno_evaluate_non_uniform(
                    t,
                    y_values_hist,
                    h_widths_hist,
                    X_interfaces_hist,
                    2;
                    apply_zhang_shu=true,
                    lower_bound=0.0,
                    upper_bound=ub_val,
                    extrapolation=:constant
                )
            )

            push!(
                local_itp_weno_k2_lin,
                t -> weno_evaluate_non_uniform(
                    t,
                    y_values_hist,
                    h_widths_hist,
                    X_interfaces_hist,
                    2;
                    apply_zhang_shu=true,
                    lower_bound=0.0,
                    upper_bound=ub_val,
                    extrapolation=:linear
                )
            )

            push!(
                local_itp_weno_k3_const,
                t -> weno_evaluate_non_uniform(
                    t,
                    y_values_hist,
                    h_widths_hist,
                    X_interfaces_hist,
                    3;
                    apply_zhang_shu=true,
                    lower_bound=0.0,
                    upper_bound=ub_val,
                    extrapolation=:constant
                )
            )

            push!(
                local_itp_weno_k3_lin,
                t -> weno_evaluate_non_uniform(
                    t,
                    y_values_hist,
                    h_widths_hist,
                    X_interfaces_hist,
                    3;
                    apply_zhang_shu=true,
                    lower_bound=0.0,
                    upper_bound=ub_val,
                    extrapolation=:linear
                )
            )

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
            
            params_1b[d] = Dict(
                "k2_const" => p_k2_c,
                "k2_lin" => p_k2_l,
                "k3_const" => p_k3_c,
                "k3_lin" => p_k3_l
            )
        end

        params_2 =
            p_eval_piecewise_constant(
                t_start,
                param_matrix,
                X_interfaces
            )

        params_3 = try

            p_3, _ = predict_single_step(
                param_matrix,
                h_widths,
                start_cell - 1, # Use available historic cell count
                1,
                1.0,
                1.0
            )

            p_3

        catch

            params_2

        end

        methods_dict = Dict{String, Any}(
            "1a) k=2 const" => t -> t <= t_start ? p_eval_piecewise_constant(t, param_matrix, X_interfaces) : [local_itp_weno_k2_const[j](t) for j in 1:num_opt_vars],
            "1a) k=2 lin"   => t -> t <= t_start ? p_eval_piecewise_constant(t, param_matrix, X_interfaces) : [local_itp_weno_k2_lin[j](t) for j in 1:num_opt_vars],
            "1a) k=3 const" => t -> t <= t_start ? p_eval_piecewise_constant(t, param_matrix, X_interfaces) : [local_itp_weno_k3_const[j](t) for j in 1:num_opt_vars],
            "1a) k=3 lin"   => t -> t <= t_start ? p_eval_piecewise_constant(t, param_matrix, X_interfaces) : [local_itp_weno_k3_lin[j](t) for j in 1:num_opt_vars],
            "2) current" => t -> t <= t_start ? p_eval_piecewise_constant(t, param_matrix, X_interfaces) : params_2,
            "3) predictor" => t -> t <= t_start ? p_eval_piecewise_constant(t, param_matrix, X_interfaces) : params_3
        )

        for d in 1:Int(T_forecast)
            methods_dict["1b) $(d)-Day Avg k=2 const"] = t -> t <= t_start ? p_eval_piecewise_constant(t, param_matrix, X_interfaces) : params_1b[d]["k2_const"]
            methods_dict["1b) $(d)-Day Avg k=2 lin"]   = t -> t <= t_start ? p_eval_piecewise_constant(t, param_matrix, X_interfaces) : params_1b[d]["k2_lin"]
            methods_dict["1b) $(d)-Day Avg k=3 const"] = t -> t <= t_start ? p_eval_piecewise_constant(t, param_matrix, X_interfaces) : params_1b[d]["k3_const"]
            methods_dict["1b) $(d)-Day Avg k=3 lin"]   = t -> t <= t_start ? p_eval_piecewise_constant(t, param_matrix, X_interfaces) : params_1b[d]["k3_lin"]
        end

        return methods_dict

    end

    # ==============================================================================
    # NONAUTONOMOUS WRAPPERS
    # ==============================================================================

    function P_non_autonomous(
        u,
        p_static,
        t,
        p_eval
    )

        p_dynamic = copy(p_static)

        p_vals = p_eval(t)

        for j in 1:num_opt_vars

            p_dynamic[
                opt_positions_filtered[j]
            ] = p_vals[j]

        end

        return P(
            u,
            Tuple(p_dynamic),
            t
        )

    end

    function d_non_autonomous(
        u,
        p_static,
        t,
        p_eval
    )

        p_dynamic = copy(p_static)

        p_vals = p_eval(t)

        for j in 1:num_opt_vars

            p_dynamic[
                opt_positions_filtered[j]
            ] = p_vals[j]

        end

        return d(
            u,
            Tuple(p_dynamic),
            t
        )

    end

    # ==============================================================================
    # BASE REFERENCE SOLUTION
    # ==============================================================================

    println("Computing reference solution...")

    base_p = all_models[Model_no]

    function p_eval_reference(t)

        return p_eval_piecewise_constant(
            t,
            param_matrix,
            X_interfaces
        )

    end

    P_ref(u, p, t) =
        P_non_autonomous(
            u,
            p,
            t,
            p_eval_reference
        )

    d_ref(u, p, t) =
        d_non_autonomous(
            u,
            p,
            t,
            p_eval_reference
        )

    prob_ref = PDSProblem(
        P_ref,
        d_ref,
        original_u0,
        (0.0, t_final),
        base_p
    )

    sol_ref = solve(prob_ref, method)

    println("Reference solution complete.")

    # ==============================================================================
    # COST FUNCTION
    # ==============================================================================

    function compute_simple_cost(
        sol_forecast,
        sol_reference,
        t0,
        t1
    )

        ts = range(t0, t1, length=100)

        err = 0.0
        ref_norm = 0.0

        for t in ts

            u_ref = sol_reference(t)
            u_fc = sol_forecast(t)
            err += norm(u_ref - u_fc)
            ref_norm += norm(max.(u_ref, 1.0))

        end

        return err / ref_norm

    end

    function compute_numerical_ID_cost(
        sol_forecast,
        sol_reference,
        t0,
        t1
    )
        ts = range(t0, t1, length=100)
        
        ref_I = [sol_reference(t)[5] for t in ts]
        fc_I = [sol_forecast(t)[5] for t in ts]
        denom_I = max.(ref_I, 1.0)
        cost_I = norm(fc_I .- ref_I) / norm(denom_I)
        
        ref_D = [sol_reference(t)[11] for t in ts]
        fc_D = [sol_forecast(t)[11] for t in ts]
        denom_D = max.(ref_D, 1.0)
        cost_D = norm(fc_D .- ref_D) / norm(denom_D)
        
        return cost_I + cost_D
    end

    # ==============================================================================
    # REAL-DATA COST (time-dependent parameters, direct comparison to observations)
    # ==============================================================================
    # Evaluates how well sol_forecast (already solved with time-dependent p_eval(t))
    # fits the available real data on the window [t0, t1].
    # Uses the same weighted relative-error logic as cost() in main.jl.

    function compute_real_data_cost(
        sol_forecast,
        t0,
        t1
    )

        weights = ones(11)

        err = zeros(11)

        # Compartment 5 (I – currently infected)
        if isdefined(Main, :tref_I) && isdefined(Main, :refI)
            valid = findall(x -> x >= t0 && x <= t1, Main.tref_I)
            if !isempty(valid)
                ts = Main.tref_I[valid]
                data = Float64.(Main.refI[valid])
                sim = [sol_forecast(t)[5] for t in ts]
                denom = max.(data, 1.0)
                err[5] = weights[5] * norm(sim .- data) / norm(denom)
            end
        end

        # Compartment 11 (D_covid – cumulative COVID deaths)
        if isdefined(Main, :tref_D) && isdefined(Main, :refD)
            valid = findall(x -> x >= t0 && x <= t1, Main.tref_D)
            if !isempty(valid)
                ts = Main.tref_D[valid]
                data = Float64.(Main.refD[valid])
                sim = [sol_forecast(t)[11] for t in ts]
                denom = max.(data, 1.0)
                err[11] = weights[11] * norm(sim .- data) / norm(denom)
            end
        end

        total = sum(err)
        return isnan(total) || isinf(total) ? 1e16 : total

    end

    # ==============================================================================
    # BUILD SHIFTED INITIAL CONDITION (mass-conserving, R absorbs difference)
    # ==============================================================================
    # Constructs a modified u0 where compartments with real data at t_start are
    # snapped to the observed values, and compartment S (index 1) compensates
    # for the change in total mass so that sum(u0_shifted) ≈ sum(u0_forecast).

    function build_shifted_u0(
        u0_forecast,
        t_start
    )

        u0_shifted = copy(u0_forecast)

        delta_mass = 0.0   # total mass injected into / removed from other compartments

        # Compartment 5: currently infected
        if isdefined(Main, :tref_I) && isdefined(Main, :refI)
            idx_near = argmin(abs.(Main.tref_I .- t_start))
            if abs(Main.tref_I[idx_near] - t_start) <= 2.0  # within 2 days
                obs_I = Float64(Main.refI[idx_near])
                delta_mass += obs_I - u0_shifted[5]
                u0_shifted[5] = obs_I
            end
        end

        # Compartment 11: cumulative COVID deaths
        if isdefined(Main, :tref_D) && isdefined(Main, :refD)
            idx_near = argmin(abs.(Main.tref_D .- t_start))
            if abs(Main.tref_D[idx_near] - t_start) <= 2.0
                obs_D = Float64(Main.refD[idx_near])
                delta_mass += obs_D - u0_shifted[11]
                u0_shifted[11] = obs_D
            end
        end

        # Compensate in compartment S (index 1) to conserve total mass
        u0_shifted[1] = max(0.0, u0_shifted[1] - delta_mass)

        return u0_shifted

    end

    # ==============================================================================
    # ROLLING FORECAST LOOP
    # ==============================================================================

    forecast_times = Float64[]

    method_costs = Dict{String,Vector{Float64}}()
    method_real_data_costs = Dict{String,Vector{Float64}}()
    method_shifted_u0_costs = Dict{String,Vector{Float64}}()
    method_numerical_ID_costs = Dict{String,Vector{Float64}}()

    method_names = [
        "1a) k=2 const",
        "1a) k=2 lin",
        "1a) k=3 const",
        "1a) k=3 lin",
        "2) current",
        "3) predictor"
    ]
    for d in 1:Int(T_forecast)
        push!(method_names, "1b) $(d)-Day Avg k=2 const")
        push!(method_names, "1b) $(d)-Day Avg k=2 lin")
        push!(method_names, "1b) $(d)-Day Avg k=3 const")
        push!(method_names, "1b) $(d)-Day Avg k=3 lin")
    end

    for m in method_names
        method_costs[m] = Float64[]
        method_real_data_costs[m] = Float64[]
        method_shifted_u0_costs[m] = Float64[]
        method_numerical_ID_costs[m] = Float64[]
    end

    println("Starting rolling forecasts...")

    # Last interface that still has a full T_forecast window within the available data
    t_data_end = min(X_interfaces[end], t_final)

    max_start_cell = something(
        findlast(i -> X_interfaces[i] + T_forecast <= t_data_end, 1:length(X_interfaces)),
        min_forecast_cell
    )

    println("Data range: [0.0, ", t_data_end, "]  |  Forecast window: ", T_forecast, " days")
    println("start_cell range: ", min_forecast_cell, " to ", max_start_cell)

    for start_cell in min_forecast_cell:max_start_cell

        t_start = X_interfaces[start_cell]
        t_end_forecast = t_start + T_forecast   # always exactly T_forecast days; within data range by construction

        println("--------------------------------")
        println("Forecast start = ", t_start)
        println("--------------------------------")

        push!(forecast_times, t_start)

        u0_forecast = sol_ref(t_start)

        # Dynamic slice built contextually per loop step
        methods_dict =
            build_forecast_methods(
                start_cell,
                t_start,
                t_end_forecast,
                param_matrix,
                X_interfaces,
                h_widths,
                num_opt_vars
            )

        for method_name in method_names

            println("Running ", method_name)

            p_eval = methods_dict[method_name]

            P_temp(u, p, t) =
                P_non_autonomous(
                    u,
                    p,
                    t,
                    p_eval
                )

            d_temp(u, p, t) =
                d_non_autonomous(
                    u,
                    p,
                    t,
                    p_eval
                )

            prob = PDSProblem(
                P_temp,
                d_temp,
                u0_forecast,
                (t_start, t_end_forecast),
                base_p
            )

            sol_forecast =
                solve(prob, method)

            # 1. Numerische Abweichung zur Referenzlösung berechnen
            ref_cost =
                compute_simple_cost(
                    sol_forecast,
                    sol_ref,
                    t_start,
                    t_end_forecast
                )

            push!(
                method_costs[method_name],
                ref_cost
            )

            id_cost = compute_numerical_ID_cost(
                sol_forecast,
                sol_ref,
                t_start,
                t_end_forecast
            )
            push!(
                method_numerical_ID_costs[method_name],
                id_cost
            )

            # 2. Abweichung zu den realen Daten – mit zeitabhängigen Parametern p_eval(t)
            # Wir lösen das nicht-autonome Problem erneut (wie oben) und vergleichen
            # die Lösung direkt mit den realen Daten (analog zu cost() in main.jl).
            real_data_cost = try
                compute_real_data_cost(
                    sol_forecast,
                    t_start,
                    t_end_forecast
                )
            catch e
                println("Warning in real data cost calculation for ", method_name, ": ", e)
                NaN
            end

            push!(
                method_real_data_costs[method_name],
                real_data_cost
            )

            # 3. Abweichung zu den realen Daten mit verschobenem Anfangswert
            # u0_shifted passt den Startzustand an die realen Daten bei t_start an;
            # Kompartiment R (Index 7) kompensiert die Differenz für Massenerhaltung.
            u0_shifted = build_shifted_u0(
                u0_forecast,
                t_start
            )

            prob_shifted = PDSProblem(
                P_temp,
                d_temp,
                u0_shifted,
                (t_start, t_end_forecast),
                base_p
            )

            sol_shifted = solve(prob_shifted, method)

            shifted_cost = try
                compute_real_data_cost(
                    sol_shifted,
                    t_start,
                    t_end_forecast
                )
            catch e
                println("Warning in shifted-u0 cost calculation for ", method_name, ": ", e)
                NaN
            end

            push!(
                method_shifted_u0_costs[method_name],
                shifted_cost
            )

        end

    end

    # ==============================================================================
    # DETERMINE TOP 3 METHODS
    # ==============================================================================

    function get_top_M_func(cost_dict, M)
        means = [(m, mean(cost_dict[m])) for m in method_names if !isempty(cost_dict[m])]
        sort!(means, by = x -> x[2])
        return [x[1] for x in means[1:min(M, length(means))]]
    end

    if plot_top_M > 0
        println("Filtering to top $plot_top_M methods per graph (independently)...")
        methods_to_plot_costs        = get_top_M_func(method_costs, plot_top_M)
        methods_to_plot_real_data    = get_top_M_func(method_real_data_costs, plot_top_M)
        methods_to_plot_shifted_u0   = get_top_M_func(method_shifted_u0_costs, plot_top_M)
        methods_to_plot_numerical_ID = get_top_M_func(method_numerical_ID_costs, plot_top_M)
        
        println("Top $plot_top_M for reference cost:      ", methods_to_plot_costs)
        println("Top $plot_top_M for real data cost:       ", methods_to_plot_real_data)
        println("Top $plot_top_M for shifted u0 cost:      ", methods_to_plot_shifted_u0)
        println("Top $plot_top_M for numerical I&D cost:   ", methods_to_plot_numerical_ID)
    else
        methods_to_plot_costs        = method_names
        methods_to_plot_real_data    = method_names
        methods_to_plot_shifted_u0   = method_names
        methods_to_plot_numerical_ID = method_names
    end

    # ==============================================================================
    # PLOT COST EVOLUTION (NUMERICAL REFERENCE)
    # ==============================================================================

    println("Plotting rolling forecast costs ($T_forecast Days vs Reference)...")

    fig, ax = plt.subplots(
        figsize=(12, 7)
    )

    # ── Distinguishable rank-based styling helper function ──────────────────────
    function generate_rank_styles(methods_list)
        styles_dict = Dict{String, Tuple{Any, String, String}}()
        cmap = plt.get_cmap("tab10")
        linestyles = ["-", "--", "-.", ":"]
        markers = [
            "o", "s", "^", "D",
            "v", "P", "X", "*",
            "h", "p"
        ]
        n_methods = length(methods_list)
        for (rank, label) in enumerate(methods_list)
            c = cmap((rank - 1) / max(n_methods - 1, 1))
            ls = linestyles[mod1(rank, length(linestyles))]
            mk = markers[mod1(rank, length(markers))]
            styles_dict[label] = (c, ls, mk)
        end
        return styles_dict
    end

    _rf_n_pts = length(forecast_times)
    _rf_me    = max(1, div(_rf_n_pts, 15))  # ~15 markers per curve

    styles_costs = generate_rank_styles(methods_to_plot_costs)
    for method_name in methods_to_plot_costs

        col, ls, mk = styles_costs[method_name]

        ax.semilogy(
            forecast_times,
            method_costs[method_name],
            color=col,
            linestyle=ls,
            linewidth=2,
            label=method_name,
            marker=mk,
            markevery=_rf_me,
            markersize=5,
            markeredgewidth=0.5,
            markeredgecolor="white"
        )

    end

    ax.set_xlabel(
        "Forecast start time",
        fontsize=12,
        fontweight="bold"
    )

    ax.set_ylabel(
        "Relative forecast cost",
        fontsize=12,
        fontweight="bold"
    )

    ax.set_title(
        "Rolling Forecast Performance ($T_forecast Days vs Reference Solution)",
        fontsize=14,
        fontweight="bold"
    )

    ax.grid(true, linestyle=":")
    ax.legend()

    plt.tight_layout()

    _rf_base1 = "rolling_forecast_costs_$(BACKUP_ID)_tf$(t_final)_T$(T_forecast)"
    plt.savefig("$(_rf_base1).png", dpi=300)
    plt.savefig("$(_rf_base1).eps", format="eps")
    plt.savefig("$(_rf_base1).pdf", format="pdf")

    if !isdefined(Main, :AntigravityHeadless)
        plt.show()
    end

    # ==============================================================================
    # PLOT COST EVOLUTION (REAL DATA)
    # ==============================================================================

    println("Plotting rolling forecast costs ($T_forecast Days vs Real Data)...")

    fig2, ax2 = plt.subplots(
        figsize=(12, 7)
    )

    styles_real_data = generate_rank_styles(methods_to_plot_real_data)
    for method_name in methods_to_plot_real_data

        col, ls, mk = styles_real_data[method_name]

        ax2.semilogy(
            forecast_times,
            method_real_data_costs[method_name],
            color=col,
            linestyle=ls,
            linewidth=2,
            label=method_name,
            marker=mk,
            markevery=_rf_me,
            markersize=5,
            markeredgewidth=0.5,
            markeredgecolor="white"
        )

    end

    ax2.set_xlabel(
        "Forecast start time",
        fontsize=12,
        fontweight="bold"
    )

    ax2.set_ylabel(
        "Cost Value ($T_forecast Days Error vs Real Data)",
        fontsize=12,
        fontweight="bold"
    )

    ax2.set_title(
        "Rolling Forecast Performance ($T_forecast Days vs Real Data)",
        fontsize=14,
        fontweight="bold"
    )

    ax2.grid(true, linestyle=":")
    ax2.legend()

    plt.tight_layout()

    _rf_base2 = "rolling_forecast_real_data_costs_$(BACKUP_ID)_tf$(t_final)_T$(T_forecast)"
    plt.savefig("$(_rf_base2).png", dpi=300)
    plt.savefig("$(_rf_base2).eps", format="eps")
    plt.savefig("$(_rf_base2).pdf", format="pdf")

    if !isdefined(Main, :AntigravityHeadless)
        plt.show()
    end

    # ==============================================================================
    # PLOT COST EVOLUTION (REAL DATA, SHIFTED U0)
    # ==============================================================================

    println("Plotting rolling forecast costs ($T_forecast Days vs Real Data, Shifted u0)...")

    fig3, ax3 = plt.subplots(
        figsize=(12, 7)
    )

    styles_shifted_u0 = generate_rank_styles(methods_to_plot_shifted_u0)
    for method_name in methods_to_plot_shifted_u0

        col, ls, mk = styles_shifted_u0[method_name]

        ax3.semilogy(
            forecast_times,
            method_shifted_u0_costs[method_name],
            color=col,
            linestyle=ls,
            linewidth=2,
            label=method_name,
            marker=mk,
            markevery=_rf_me,
            markersize=5,
            markeredgewidth=0.5,
            markeredgecolor="white"
        )

    end

    ax3.set_xlabel(
        "Forecast start time",
        fontsize=12,
        fontweight="bold"
    )

    ax3.set_ylabel(
        "Cost Value ($T_forecast Days Error vs Real Data)",
        fontsize=12,
        fontweight="bold"
    )

    ax3.set_title(
        "Rolling Forecast Performance ($T_forecast Days, u\u2080 shifted to match real data, S compensates)",
        fontsize=14,
        fontweight="bold"
    )

    ax3.grid(true, linestyle=":")
    ax3.legend()

    plt.tight_layout()

    _rf_base3 = "rolling_forecast_shifted_u0_costs_$(BACKUP_ID)_tf$(t_final)_T$(T_forecast)"
    plt.savefig("$(_rf_base3).png", dpi=300)
    plt.savefig("$(_rf_base3).eps", format="eps")
    plt.savefig("$(_rf_base3).pdf", format="pdf")

    if !isdefined(Main, :AntigravityHeadless)
        plt.show()
    end

    # ==============================================================================
    # PLOT COST EVOLUTION (NUMERICAL REFERENCE, ONLY I AND D_COVID)
    # ==============================================================================

    println("Plotting rolling forecast costs ($T_forecast Days vs Numerical Reference, I and D_covid only)...")

    fig4, ax4 = plt.subplots(
        figsize=(12, 7)
    )

    styles_numerical_ID = generate_rank_styles(methods_to_plot_numerical_ID)
    for method_name in methods_to_plot_numerical_ID

        col, ls, mk = styles_numerical_ID[method_name]

        ax4.semilogy(
            forecast_times,
            method_numerical_ID_costs[method_name],
            color=col,
            linestyle=ls,
            linewidth=2,
            label=method_name,
            marker=mk,
            markevery=_rf_me,
            markersize=5,
            markeredgewidth=0.5,
            markeredgecolor="white"
        )

    end

    ax4.set_xlabel(
        "Forecast start time",
        fontsize=12,
        fontweight="bold"
    )

    ax4.set_ylabel(
        "Relative Error (I & D_covid)",
        fontsize=12,
        fontweight="bold"
    )

    ax4.set_title(
        "Rolling Forecast Performance ($T_forecast Days, I and D_covid vs Numerical Reference)",
        fontsize=14,
        fontweight="bold"
    )

    ax4.grid(true, linestyle=":")
    ax4.legend()

    plt.tight_layout()

    _rf_base4 = "rolling_forecast_numerical_ID_costs_$(BACKUP_ID)_tf$(t_final)_T$(T_forecast)"
    plt.savefig("$(_rf_base4).png", dpi=300)
    plt.savefig("$(_rf_base4).eps", format="eps")
    plt.savefig("$(_rf_base4).pdf", format="pdf")

    if !isdefined(Main, :AntigravityHeadless)
        plt.show()
    end

    for method_name in methods_to_plot_real_data
        if haskey(method_real_data_costs, method_name) && !isempty(method_real_data_costs[method_name])
            costs_vec = method_real_data_costs[method_name]
            for (i, c) in enumerate(costs_vec)
                push!(all_real_data_points, (
                    data_bound = forecast_times[i],
                    T_forecast = T_forecast,
                    cost = c,
                    method = method_name
                ))
            end
        end
    end

    println("------------------------------------------------------------")
    println("Rolling Forecast Evaluation Complete for T = $T_forecast")
    println("------------------------------------------------------------")
end

# ==============================================================================
# RUN BEST FORECASTS
# ==============================================================================
if !isempty(all_real_data_points)
    for T_cur in T_forecasts
        points_for_T = filter(p -> p.T_forecast == T_cur, all_real_data_points)
        if isempty(points_for_T)
            continue
        end
        
        println("============================================================")
        println("Finding BEST and WORST for T_forecast = ", T_cur)
        println("============================================================")

        # 1. Find the "best curve" (the method that contains the global minimum cost for this T)
        best_point_for_T = argmin(p -> p.cost, points_for_T)
        best_method_name = best_point_for_T.method
        
        # 2. Filter all points belonging to this specific method for this T
        method_points = filter(p -> p.method == best_method_name, points_for_T)
        
        # 3. Find the best (min cost) and worst (max cost) points for this method
        pair1 = argmin(p -> p.cost, method_points)
        pair2 = argmax(p -> p.cost, method_points)
        
        println("Best method (curve) identified: ", best_method_name)
        println("Pair 1 (min cost on this curve): data_bound = ", pair1.data_bound,
                ", T_forecast = ", pair1.T_forecast,
                ", cost = ", pair1.cost)
        println("Pair 2 (max cost on this curve): data_bound = ", pair2.data_bound,
                ", T_forecast = ", pair2.T_forecast,
                ", cost = ", pair2.cost)
        
        println("Running forecast for Pair 1...")
        global data_bound = pair1.data_bound
        global T_forecast = pair1.T_forecast
        global forecast_suffix = "_best"
        include("run_forecasts.jl")
        
        println("Running forecast for Pair 2...")
        global data_bound = pair2.data_bound
        global T_forecast = pair2.T_forecast
        global forecast_suffix = "_worst"
        include("run_forecasts.jl")
    end
end

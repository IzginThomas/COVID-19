module ParameterPredictor

using LinearAlgebra
using Statistics

export predict_until_time

function predict_single_step(pred_matrix::Matrix{Float64}, pred_h::Vector{Float64}, max_search_idx::Int, L::Int, min_duration::Float64, gamma::Float64)
    num_vars, N_total = size(pred_matrix)
    
    current_state = pred_matrix[:, end]
    current_trend = pred_matrix[:, end] .- pred_matrix[:, end-1]
    
    distances = Float64[]
    valid_indices = Int[]
    
    for t_idx in 2:max_search_idx
        historical_state = pred_matrix[:, t_idx]
        historical_trend = pred_matrix[:, t_idx] .- pred_matrix[:, t_idx-1]
        
        # 1. EUCLIDEAN PARAMETER DISTANCE (Sum of squares)
        # Every parameter contributes to the distance metric.
        sum_sq_dist = sum((current_state .- historical_state).^2)
        dist_abs = sqrt(sum_sq_dist)
        
        # 2. Euclidean Trend Distance
        dist_trend = sqrt(sum((current_trend .- historical_trend).^2))
        
        # 3. Dynamic Direction Protection
        direction_penalty = 0.0
        for v in 1:num_vars
            if current_trend[v] * historical_trend[v] < 0.0
                if abs(current_trend[v]) > 1e-4
                    importance = max(0.0, current_state[v])
                    direction_penalty += 5.0 * importance
                end
            end
        end
        
        # Total distance calculation balanced by gamma weighting
        total_dist = dist_abs + gamma * dist_trend + direction_penalty
        
        push!(distances, total_dist)
        push!(valid_indices, t_idx)
    end
    
    effective_L = min(L, length(valid_indices))
    if effective_L == 0
        return pred_matrix[:, end], mean(pred_h[1:max_search_idx])
    end
    
    sorted_indices = sortperm(distances)
    L_nearest_indices = valid_indices[sorted_indices[1:effective_L]]
    L_nearest_distances = distances[sorted_indices[1:effective_L]]
    
    weights = zeros(effective_L)
    twin_count = sum(L_nearest_distances .< 1e-12)
    
    if twin_count > 0
        for k in 1:effective_L
            if L_nearest_distances[k] < 1e-12
                weights[k] = 1.0 / twin_count
            else
                weights[k] = 0.0
            end
        end
    else
        shift = mean(L_nearest_distances) * 0.1
        inverse_dist = 1.0 ./ (L_nearest_distances .+ shift)
        weights = inverse_dist / sum(inverse_dist)
    end
    
    next_state = zeros(num_vars)
    next_duration = 0.0
    
    for k in 1:effective_L
        hist_idx = L_nearest_indices[k]
        successor_idx = hist_idx + 1
        w = weights[k]
        
        next_state .+= w .* pred_matrix[:, successor_idx]
        
        # Enforce non-negativity constraint
        for v in 1:num_vars
            if next_state[v] < 0.0
                next_state[v] = 0.0
            end
        end
        
        next_duration += w * pred_h[successor_idx]
    end
    
    if next_duration < min_duration
        next_duration = min_duration
    end
    
    return next_state, next_duration
end

function predict_until_time(param_matrix_base::Matrix{Float64}, h_widths_base::Vector{Float64}, X_interfaces_base::Vector{Float64}, target_end_time::Float64; L::Int=3, min_duration::Float64=7.0, gamma::Float64=1.0)
    N_real = size(param_matrix_base, 2)
    pred_matrix = copy(param_matrix_base)
    pred_h = copy(h_widths_base)
    pred_X = copy(X_interfaces_base)
    current_time = pred_X[end]
    
    if current_time >= target_end_time
        return pred_matrix, pred_h, pred_X
    end
    
    println("\n=== Starting Iterative Forecasting (Euclidian Balanced) ===")
    step_counter = 0
    
    while current_time < target_end_time
        step_counter += 1
        max_search_idx = (N_real - 1) - (step_counter - 1)
        
        if max_search_idx < 3
            remaining_duration = target_end_time - current_time
            current_time = target_end_time
            pred_h[end] += remaining_duration
            pred_X[end] = target_end_time
            break
        end
        
        next_params, next_duration = predict_single_step(pred_matrix, pred_h, max_search_idx, L, min_duration, gamma)
        
        if current_time + next_duration > target_end_time
            remaining = target_end_time - current_time
            if remaining < min_duration && step_counter > 1
                pred_h[end] += remaining
                pred_X[end] = target_end_time
                current_time = target_end_time
                break
            else
                next_duration = remaining
            end
        end
        
        current_time += next_duration
        push!(pred_h, next_duration)
        push!(pred_X, current_time)
        pred_matrix = hcat(pred_matrix, next_params)
    end
    
    return pred_matrix, pred_h, pred_X
end

end # module ParameterPredictor


# =====================================================================
# MAIN EXECUTION SCRIPT
# =====================================================================

using DelimitedFiles
using LinearAlgebra # Important for pinv() in the new WENO function

if !@isdefined(BACKUP_ID)
    BACKUP_ID = "backupI"
end

base_p = all_models[Model_no];
println("--- Loading optimization data from CSV backups (BACKUP_ID: $BACKUP_ID)... ---")
midpoints = vec(DelimitedFiles.readdlm("midpoints_$(BACKUP_ID).csv", ',', Float64))
param_matrix = DelimitedFiles.readdlm("param_matrix_$(BACKUP_ID).csv", ',', Float64)
println("Data successfully loaded!")

println("--- Reconstructing time interpolation functions ---")
interpolated_funcs = []
num_opt_vars = length(opt_positions_filtered)

# --------------------------------===================================
# NEW: Derive grid metrics from cell midpoints for the WENO module
# --------------------------------===================================
# Since we only saved the cell midpoints, we reconstruct the cell boundaries (interfaces).
# For a uniform or near-uniform optimization time-grid:
N_cells = length(midpoints)

if N_cells > 1
    # Calculate the step sizes between midpoints
    dt_mid = diff(midpoints)
    # Use the first difference as an approximation for boundary cell widths
    dt_start = dt_mid[1]
    dt_end = dt_mid[end]
    
    # Reconstruct interfaces (cell boundaries)
    X_interfaces = zeros(N_cells + 1)
    X_interfaces[1] = 0.0
    for i in 1:N_cells
        X_interfaces[i+1] = 2.0 * midpoints[i] - X_interfaces[i]
    end
    
    # Compute exact cell widths (delta_x)
    h_widths = diff(X_interfaces)
else
    X_interfaces = [midpoints[1] - 0.5, midpoints[1] + 0.5]
    h_widths = [1.0]
end

if !@isdefined(end_time)
    end_time = length(X_interfaces) > 1 ? X_interfaces[end] : 200.0
end





using Interpolations # If desired for spline-based approaches

function evaluate_piecewise(t::Float64, midpoints::Vector{Float64}, h_widths::Vector{Float64}, 
                            X_interfaces::Vector{Float64}, y_values::Vector{Float64}, k::Int; 
                            method_type::Symbol=:weno, apply_zhang_shu::Bool=true, 
                            lower_bound::Float64=0.0)
    
    N = length(y_values)
    
    # 1. CENTRAL INDEX DETERMINATION
    if t <= X_interfaces[1]
        i = 1
    elseif t >= X_interfaces[end]
        i = N
    else
        i = searchsortedlast(X_interfaces, t)
        i = max(1, min(i, N))
    end
    
    v_raw = 0.0
    
    if method_type == :weno
        v_raw = weno_evaluate_non_uniform(
            t, y_values, h_widths, X_interfaces, k; 
            apply_zhang_shu=false
        )
        
    elseif method_type == :linear
        if N == 1
            v_raw = y_values[1]
        elseif t <= X_interfaces[1]
            v_raw = y_values[1]
        elseif t >= X_interfaces[end]
            v_raw = y_values[end]
        else
            # --- MATHEMATICALLY CORRECT CONSERVATIVE RECONSTRUCTION ---
            # We calculate interface values V_face such that continuity 
            # AND mean-value preservation are guaranteed across the equation system.
            V_face = zeros(N + 1)
            
            # Since we have a non-uniform grid, we use the standard FVD
            # (Finite Volume Reconstruction) derivation for continuous linear profiles:
            for j in 2:N
                # Distance from the center of the left cell to the center of the right cell
                dx_m2m = midpoints[j] - midpoints[j-1]
                # Weighting based on the distance from midpoints to interface surfaces
                V_face[j] = ((midpoints[j] - X_interfaces[j]) * y_values[j-1] + 
                             (X_interfaces[j] - midpoints[j-1]) * y_values[j]) / dx_m2m
            end
            
            # Enforce exact mean-value preservation for boundary cells
            V_face[1]   = 2.0 * y_values[1] - V_face[2]
            V_face[end] = 2.0 * y_values[end] - V_face[end-1]
            
            # Correction step for the interior cells:
            # To ensure the cell mean value is EXACTLY accurate in cell i, the profile slope must be corrected.
            # A purely interpolated V_face ensures continuity but slightly biases the mean value.
            # We shift the line parallelly to guarantee the exact mean value:
            v_L_temp = V_face[i]
            v_R_temp = V_face[i+1]
            temp_avg = 0.5 * (v_L_temp + v_R_temp)
            
            # Error discrepancy relative to the target cell mean
            delta_avg = y_values[i] - temp_avg
            
            # Corrected, exact mean-preserving endpoints for this specific cell
            v_L = v_L_temp + delta_avg
            v_R = v_R_temp + delta_avg
            
            # Linear function implementation inside the interval
            t_L = X_interfaces[i]
            v_raw = v_L + (v_R - v_L) * (t - t_L) / h_widths[i]
        end

    elseif method_type == :constant
        if t <= X_interfaces[1]
            v_raw = y_values[1]
        elseif t >= X_interfaces[end]
            # ----------------------------------------------------------------
            # FUTURE EXTRAPOLATION (t > end_time)
            # ----------------------------------------------------------------
            # Option A: Maintain status quo (value of the last optimized cell)
            v_raw = y_values[end]
            
            # Option B (Alternative): Return to static base value from base_p
            # For this, you would pass the base value as 'default_value' to the function:
            # v_raw = default_value 
        else
            # Normal range inside the data grid bounds
            v_raw = y_values[i]
        end
    else
        error("Unknown interpolation type: :$(method_type)")
    end

    # 3. CENTRAL ZHANG-SHU LIMITER
    if apply_zhang_shu
        v_avg = y_values[i]
        if v_raw < lower_bound
            denom = v_avg - v_raw
            theta = denom > 1e-14 ? min(1.0, (v_avg - lower_bound) / denom) : 1.0
            return theta * (v_raw - v_avg) + v_avg
        end
    end

    return v_raw
end




println("--- Reconstructing time interpolation functions ---")
interpolated_funcs = []
num_opt_vars = length(opt_positions_filtered)

# --- CONFIGURATION REGION FOR RECONSTRUCTION ---
chosen_method    = :weno   # Options: :weno, :linear, :constant
weno_k           = 2       # k=2 for WENO-3, k=3 for WENO-5
limit_positivity = true    # Activate Zhang-Shu Limiter

for j in 1:num_opt_vars
    y_values = param_matrix[j, :]
    
    if length(midpoints) > 1
        # The anonymous function now dynamically accesses our all-round wrapper
        itp = t -> evaluate_piecewise(
            t, 
            midpoints, 
            h_widths, 
            X_interfaces, 
            y_values, 
            weno_k; 
            method_type = chosen_method,
            apply_zhang_shu = limit_positivity, 
            lower_bound = 0.0
        )
    else
        itp = t -> y_values[1]
    end
    push!(interpolated_funcs, itp)
end
println("Interpolation functions successfully generated using method [:", chosen_method, "]!")

# Local wrappers for the non-autonomous dynamics
function P_non_autonomous_load(u, p_static, t)
    p_dynamic = copy(p_static)
    for j in 1:num_opt_vars
        p_dynamic[opt_positions_filtered[j]] = interpolated_funcs[j](t)
    end
    return P(u, Tuple(p_dynamic), t)
end

function d_non_autonomous_load(u, p_static, t)
    p_dynamic = copy(p_static)
    for j in 1:num_opt_vars
        p_dynamic[opt_positions_filtered[j]] = interpolated_funcs[j](t)
    end
    return d(u, Tuple(p_dynamic), t)
end

println("--- Solving the final non-autonomous differential equation from backup ---")
if !@isdefined(original_u0)
    original_u0 = [30416000.0, 0.0, 5.0, 5.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
end
if !@isdefined(method)
    method = MPRK22(1.0)
end

final_prob_load = PDSProblem(P_non_autonomous_load, d_non_autonomous_load, original_u0, (0.0, end_time), base_p)
final_sol_load = solve(final_prob_load, method)

# Call plot execution (uses plot_res defined in test_modul_bayesian_tools.jl)
local_sols = @isdefined(local_solutions_history) ? local_solutions_history : nothing
plot_res(final_sol_load, base_p, interpolated_funcs, local_sols=local_sols)

# Save the final non-autonomous plot
fig_name = "final_non_autonomous_$(BACKUP_ID)_$(end_time).png"
plt.savefig(fig_name, dpi=300)
println("Saved final non-autonomous plot to $fig_name")
if !isdefined(Main, :AntigravityHeadless)
    plt.show()
end


# # 1. Helper function to generate valid LaTeX strings for Matplotlib
function format_param_latex(name::String)
    chars = collect(name)
    len = length(chars)

    # NEW RULE: Handle 2-character parameters with subscripts (e.g., "k1" -> k_1, "αI" -> α_I)
    if len == 2
        # Check if the second character is a digit or an uppercase letter
        if isdigit(chars[2]) || isuppercase(chars[2])
            return "\$$(chars[1])_$(chars[2])\$"
        end
    end

    # If parameter name is short (e.g., "μ", "γ", "x")
    if len < 3
        return "\$$(name)\$"
    end

    name_str = String(chars)

    # Rule for variables ending with "Ve" (e.g., "rIVe", "rLVe", "μVe")
    if endswith(name_str, "Ve")
        if len >= 4
            base = string(chars[1])            # e.g., 'r'
            sup = String(chars[2:(end-2)])     # e.g., 'I' or 'L'
            return "\$$(base)^$(sup)_{V_E}\$"
        else
            # For 3-character parameters like "μVe" -> \mu_{V_E}
            return "\$$(String(chars[1:(end-2)]))_{V_E}\$"
        end
    end

    # Rule for standard 3+ character variables (e.g., "aEL" -> a^E_L)
    m = match(r"^([a-z])([A-Z])([A-Z])$", name_str)
    if m !== nothing
        return "\$$(m.captures[1])^$(m.captures[2])_$(m.captures[3])\$"
    end

    return "\$$(name_str)\$"
end

# 2. Define the original list of all parameter names
parameters = [
    "k1", "k2", "cS", "cV", "cE", "cL", "cI", "cH", "cR", "cQ",
    "cVe", "δV", "δE", "δL", "δH", "δQ", "pv", "Λ", "μ", "φ",
    "ψ", "γ", "μVe", "αI", "αH", "aSV", "aVS", "aVE", "aEI", "aES",
    "aEQ", "aEL", "aLI", "aLR", "aLQ", "rLVe", "aIR", "aIQ", "aIH", "rIVe",
    "aHR", "aRS", "aQR", "b01","b02","b03","b04","b11","b12","b13", "b14","ω1", "ω2"
]

# 3. Extract the original numerical values from x_model for optimized indices
original_vals = x_model

# 4. Fetch text labels from names array and transform to LaTeX formats
raw_opt_names = parameters[opt_positions_filtered]
formatted_opt_names = format_param_latex.(raw_opt_names)

# # 5. Include module and trigger plotter execution
 include("WenoPlotter.jl")
 using .WenoPlotter: plot_weno_parameters

WenoPlotter.plot_weno_parameters(
    midpoints, 
    param_matrix, 
    interpolated_funcs,
    weno_k; 
    opt_positions_filtered = formatted_opt_names, # Passing pure LaTeX strings for titles
    original_values = original_vals,               # Passing constants for horizontal reference lines (axhline)
    backup_id = BACKUP_ID,
    end_time = end_time
)

##

# ... (Loading backups remains as before) ...

# Include the new prediction module defined above
# include("ParameterPredictor.jl") # (Already present within execution memory block)
using .ParameterPredictor: predict_until_time

# --- YOUR TARGET END TIME FOR FORECASTING ---
target_forecast_end = 223.0  # How long should the total simulation run? (e.g., 120 days)

# Execute iterative prediction steps
L_neighbors          = 1      # How many historical patterns to blend together
min_time_window      = 1.0    # Every future step interval spans at least 1.0 day
trend_weight         = 1.5

expanded_param_matrix, expanded_h_widths, expanded_X_interfaces = predict_until_time(
    param_matrix, 
    h_widths, 
    X_interfaces, 
    target_forecast_end; # <-- Watch the semicolon placement!
    L = L_neighbors, 
    min_duration = min_time_window, 
    gamma = trend_weight
)

# Since we now have an extended true grid layout, compute the expanded midpoints
# for internal profile processing evaluation (if needed, e.g., for plotting tasks)
expanded_midpoints = [0.5 * (expanded_X_interfaces[i] + expanded_X_interfaces[i+1]) for i in 1:length(expanded_h_widths)]


# --- GENERATION OF WRAPPER FUNCTIONS (Now parsing the expanded structural matrices!) ---
println("--- Creating time-continuous wrapper pipelines for the ODE solver ---")
interpolated_funcs = []
num_opt_vars = length(opt_positions_filtered)

for j in 1:num_opt_vars
    # We pass the EXPANDED structural parameter arrays from the Predictor engine into evaluate_piecewise
    y_values_expanded = expanded_param_matrix[j, :]
    
    itp = t -> evaluate_piecewise(
        t, 
        expanded_midpoints, 
        expanded_h_widths, 
        expanded_X_interfaces, 
        y_values_expanded, 
        weno_k; 
        method_type = :constant, # Piecewise constant behavior performs best here
        apply_zhang_shu = limit_positivity, 
        lower_bound = 0.0
    )
    push!(interpolated_funcs, itp)
end


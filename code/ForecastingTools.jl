module ForecastingTools

export test_prediction_horizon, test_prediction_horizon_weno_extrapolation

"""
    test_prediction_horizon(cost_func, original_u0, all_models, Model_no, opt_positions_filtered, method, P, d, t_prediction_start::Float64, end_time::Float64, max_cost_threshold::Float64=0.25)

FIRST TEST: Evaluates how many days starting from `t_prediction_start` can be simulated purely 
predictively using frozen parameters before the cost exceeds `max_cost_threshold`.
Based on the sequential interval reconstruction logic (Version 1).
"""
function test_prediction_horizon(
    cost_func, 
    original_u0, 
    all_models, 
    Model_no, 
    opt_positions_filtered, 
    method, 
    P, 
    d, 
    t_prediction_start::Float64, 
    end_time::Float64, 
    max_cost_threshold::Float64=0.25
)
    println("\n=======================================================")
    println("🔮 STARTING FORECAST HORIZON ANALYSIS (FIRST TEST: FROZEN) FROM DAY $t_prediction_start")
    println("=======================================================")

    if !isfile("param_matrix_backup.csv") || !isfile("midpoints_backup.csv")
        error("Error: Backups (param_matrix_backup.csv / midpoints_backup.csv) not found!")
    end

    param_matrix_loaded = Main.DelimitedFiles.readdlm("param_matrix_backup.csv", ',', Float64)
    midpoints_loaded = vec(Main.DelimitedFiles.readdlm("midpoints_backup.csv", ',', Float64))

    idx_freeze = findlast(m -> m <= t_prediction_start, midpoints_loaded)
    if idx_freeze === nothing
        error("No optimized parameter data found prior to the specified prediction start time!")
    end

    println("-> Using learned parameters from the first $idx_freeze intervals.")
    best_x_frozen = param_matrix_loaded[:, idx_freeze]

    println("--- Reconstructing u0 at prediction point step-by-step... ---")
    current_u0_pred = copy(original_u0)
    t_current = 0.0

    for j in 1:idx_freeze
        if j == 1
            t_start_local = 0.0
            t_end_local = midpoints_loaded[1] * 2.0
        else
            t_start_local = t_current
            t_end_local = midpoints_loaded[j] + (midpoints_loaded[j] - t_start_local)
        end
        
        x_all_local = copy(all_models[Model_no])
        x_all_local[opt_positions_filtered] = param_matrix_loaded[:, j]
        
        prob_local = Main.PDSProblem(P, d, current_u0_pred, (t_start_local, t_end_local), Tuple(x_all_local))
        sol_local = Main.solve(prob_local, method)
        
        current_u0_pred = copy(sol_local.u[end])
        t_current = t_end_local
    end

    println("-> System state u0 at day $t_current successfully prepared.")

    days_predicted = 0
    t_predict_limit = t_current

    while t_predict_limit < end_time
        t_next_target = min(t_predict_limit + 1.0, end_time)
        
        Main.eval(:(tspan = ($t_predict_limit, $t_next_target)))
        Main.eval(:(u0 = copy($current_u0_pred)))

        current_pred_cost = cost_func(best_x_frozen)
        println("🔮 Forecast Day $(days_predicted + 1) | Window: [$t_predict_limit -> $t_next_target] | Cost: $current_pred_cost")
        
        if current_pred_cost > max_cost_threshold
            println("🛑 Evaluation Stopped! Cost threshold exceeded ($current_pred_cost > $max_cost_threshold) at day $t_next_target.")
            break
        end

        x_all_pred = copy(all_models[Model_no])
        x_all_pred[opt_positions_filtered] = best_x_frozen
        
        prob_step = PDSProblem(P, d, current_u0_pred, (t_predict_limit, t_next_target), Tuple(x_all_pred))
        sol_step = solve(prob_step, method)
        
        current_u0_pred = copy(sol_step.u[end])
        t_predict_limit = t_next_target
        days_predicted += 1
    end

    println("\n=======================================================")
    println("🎯 FINAL RESULT (FIRST TEST: FROZEN):")
    println("   Starting from day $t_current, the model accurately forecast: 👉 $days_predicted days 👈")
    println("=======================================================")

    return days_predicted
end

"""
    test_prediction_horizon_weno_extrapolation(cost_func, original_u0, all_models, Model_no, opt_positions_filtered, method, P, d, weno3_interpolate, t_prediction_start::Float64, end_time::Float64, max_cost_threshold::Float64=0.25)

SECOND TEST: Predicts the system evolution day-by-day using active WENO extrapolation. 
At each step, future parameters are extrapolated via WENO, the non-autonomous PDS is solved,
and the history matrices are updated with the newly predicted coefficients to compute the next day's WENO curves.
"""
function test_prediction_horizon_weno_extrapolation(
    cost_func, 
    original_u0, 
    all_models, 
    Model_no, 
    opt_positions_filtered, 
    method, 
    P, 
    d,
    weno3_interpolate,
    t_prediction_start::Float64, 
    end_time::Float64, 
    max_cost_threshold::Float64=0.25
)
    println("\n=======================================================")
    println("🔮 STARTING FORECAST HORIZON ANALYSIS (SECOND TEST: WENO EXTRAPOLATION) FROM DAY $t_prediction_start")
    println("=======================================================")

    if !isfile("param_matrix_backup.csv") || !isfile("midpoints_backup.csv")
        error("Error: Backups (param_matrix_backup.csv / midpoints_backup.csv) not found!")
    end

    param_matrix_loaded = Main.DelimitedFiles.readdlm("param_matrix_backup.csv", ',', Float64)
    midpoints_loaded = vec(Main.DelimitedFiles.readdlm("midpoints_backup.csv", ',', Float64))

    idx_freeze = findlast(m -> m <= t_prediction_start, midpoints_loaded)
    if idx_freeze === nothing
        error("No optimized parameter data found prior to the prediction start time!")
    end

    current_midpoints = copy(midpoints_loaded[1:idx_freeze])
    current_param_matrix = copy(param_matrix_loaded[:, 1:idx_freeze])
    
    num_opt_vars = length(opt_positions_filtered)

    println("--- Reconstructing historical state u0 step-by-step... ---")
    current_u0_pred = copy(original_u0)
    t_current = 0.0

    for j in 1:idx_freeze
        if j == 1
            t_start_local = 0.0
            t_end_local = current_midpoints[1] * 2.0
        else
            t_start_local = t_current
            t_end_local = current_midpoints[j] + (current_midpoints[j] - t_start_local)
        end
        
        x_all_local = copy(all_models[Model_no])
        x_all_local[opt_positions_filtered] = current_param_matrix[:, j]
        
        prob_local = Main.PDSProblem(P, d, current_u0_pred, (t_start_local, t_end_local), Tuple(x_all_local))
        sol_local = Main.solve(prob_local, method)
        
        current_u0_pred = copy(sol_local.u[end])
        t_current = t_end_local
    end

    println("-> Base state u0 prepared at t = $t_current.")

    days_predicted = 0
    t_predict_limit = t_current

    while t_predict_limit < end_time
        t_next_target = min(t_predict_limit + 1.0, end_time)
        t_eval = (t_predict_limit + t_next_target) / 2.0
        
        # --- STEP A: Extrapolate parameters via WENO3 ---
        extrapolated_x = zeros(num_opt_vars)
        for j in 1:num_opt_vars
            y_vals = current_param_matrix[j, :]
            extrapolated_x[j] = weno3_interpolate(t_eval, current_midpoints, y_vals)
        end

        # --- STEP B: Validate using the cost function ---
        Main.eval(:(tspan = ($t_predict_limit, $t_next_target)))
        Main.eval(:(u0 = copy($current_u0_pred)))

        current_pred_cost = cost_func(extrapolated_x)
        println("🔮 Forecast Day $(days_predicted + 1) | Window: [$t_predict_limit -> $t_next_target] | Cost: $current_pred_cost")
        
        if current_pred_cost > max_cost_threshold
            println("🛑 Evaluation Stopped! Prediction cost exceeded the threshold ($current_pred_cost > $max_cost_threshold).")
            break
        end

        # --- STEP C: Solve the autonomous step ---
        x_all_pred = copy(all_models[Model_no])
        x_all_pred[opt_positions_filtered] = extrapolated_x
        
        prob_step = PDSProblem(P, d, current_u0_pred, (t_predict_limit, t_next_target), Tuple(x_all_pred))
        sol_step = solve(prob_step, method)
        
        # --- STEP D: Append predicted parameters to history for subsequent WENO updates ---
        push!(current_midpoints, t_eval)
        current_param_matrix = hcat(current_param_matrix, extrapolated_x)

        current_u0_pred = copy(sol_step.u[end])
        t_predict_limit = t_next_target
        days_predicted += 1
    end

    println("\n=======================================================")
    println("🎯 FINAL RESULT (SECOND TEST: WENO EXTRAPOLATION):")
    println("   Starting from day $t_current, the model accurately forecast: 👉 $days_predicted days 👈")
    println("=======================================================")

    return days_predicted, current_param_matrix, current_midpoints
end

end # module
module ParameterPredictor

using LinearAlgebra
using Statistics

export predict_until_time

"""
    predict_single_step(pred_matrix, pred_h, max_search_idx, L, min_duration, gamma)

Calculates the next parameters and step duration using a distance-weighted 
K-Nearest Neighbors (KNN) approach across historical states.
"""
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

"""
    predict_until_time(param_matrix_base, h_widths_base, X_interfaces_base, target_end_time; L, min_duration, gamma)

Runs iterative prediction forecasting step-by-step until the target time constraint is reached.
"""
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

end # module
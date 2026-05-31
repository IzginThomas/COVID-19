module ParameterPredictor

using LinearAlgebra
using Statistics

export predict_until_time

function predict_single_step(pred_matrix::Matrix{Float64}, pred_h::Vector{Float64}, max_search_idx::Int, L::Int, min_duration::Float64, gamma::Float64)
    num_vars, N_total = size(pred_matrix)
    
    aktueller_zustand = pred_matrix[:, end]
    aktueller_trend = pred_matrix[:, end] .- pred_matrix[:, end-1]
    
    abstaende = Float64[]
    valide_indizes = Int[]
    
    for t_idx in 2:max_search_idx
        historischer_zustand = pred_matrix[:, t_idx]
        historischer_trend = pred_matrix[:, t_idx] .- pred_matrix[:, t_idx-1]
        
        # 1. EUKLIDISCHER PARAMETER-ABSTAND (Summe der Quadrate)
        # Jeder Parameter liefert einen Beitrag, niemand wird ignoriert!
        sum_sq_dist = sum((aktueller_zustand .- historischer_zustand).^2)
        dist_abs = sqrt(sum_sq_dist)
        
        # 2. Euklidischer Trend-Abstand
        dist_trend = sqrt(sum((aktueller_trend .- historischer_trend).^2))
        
        # 3. Dynamischer Richtungsschutz (bleibt aktiv)
        richtungs_strafe = 0.0
        for v in 1:num_vars
            if aktueller_trend[v] * historischer_trend[v] < 0.0
                if abs(aktueller_trend[v]) > 1e-4
                    wichtigkeit = max(0.0, aktueller_zustand[v])
                    richtungs_strafe += 5.0 * wichtigkeit
                end
            end
        end
        
        # Gesamtabstand basiert nun auf der euklidischen Balance
        total_dist = dist_abs + gamma * dist_trend + richtungs_strafe
        
        push!(abstaende, total_dist)
        push!(valide_indizes, t_idx)
    end
    
    effektives_L = min(L, length(valide_indizes))
    if effektives_L == 0
        return pred_matrix[:, end], mean(pred_h[1:max_search_idx])
    end
    
    sortierte_indizes = sortperm(abstaende)
    L_naechste_indizes = valide_indizes[sortierte_indizes[1:effektives_L]]
    L_naechste_abstaende = abstaende[sortierte_indizes[1:effektives_L]]
    
    gewichte = zeros(effektives_L)
    anzahl_zwillinge = sum(L_naechste_abstaende .< 1e-12)
    
    if anzahl_zwillinge > 0
        for k in 1:effektives_L
            if L_naechste_abstaende[k] < 1e-12
                gewichte[k] = 1.0 / anzahl_zwillinge
            else
                gewichte[k] = 0.0
            end
        end
    else
        shift = mean(L_naechste_abstaende) * 0.1
        invers_dist = 1.0 ./ (L_naechste_abstaende .+ shift)
        gewichte = invers_dist / sum(invers_dist)
    end
    
    naechster_zustand = zeros(num_vars)
    naechste_dauer = 0.0
    
    for k in 1:effektives_L
        hist_idx = L_naechste_indizes[k]
        nachfolger_idx = hist_idx + 1
        w = gewichte[k]
        
        naechster_zustand .+= w .* pred_matrix[:, nachfolger_idx]
        
        for v in 1:num_vars
            if naechster_zustand[v] < 0.0
                naechster_zustand[v] = 0.0
            end
        end
        
        naechste_dauer += w * pred_h[nachfolger_idx]
    end
    
    if naechste_dauer < min_duration
        naechste_dauer = min_duration
    end
    
    return naechster_zustand, naechste_dauer
end

function predict_until_time(param_matrix_base::Matrix{Float64}, h_widths_base::Vector{Float64}, X_interfaces_base::Vector{Float64}, target_end_time::Float64; L::Int=3, min_duration::Float64=7.0, gamma::Float64=1.0)
    N_echt = size(param_matrix_base, 2)
    pred_matrix = copy(param_matrix_base)
    pred_h = copy(h_widths_base)
    pred_X = copy(X_interfaces_base)
    current_time = pred_X[end]
    
    if current_time >= target_end_time
        return pred_matrix, pred_h, pred_X
    end
    
    println("\n=== Starte iterative Prognose (Euklidisch balanciert) ===")
    step_counter = 0
    
    while current_time < target_end_time
        step_counter += 1
        max_search_idx = (N_echt - 1) - (step_counter - 1)
        
        if max_search_idx < 3
            rest_dauer = target_end_time - current_time
            current_time = target_end_time
            pred_h[end] += rest_dauer
            pred_X[end] = target_end_time
            break
        end
        
        next_params, next_duration = predict_single_step(pred_matrix, pred_h, max_search_idx, L, min_duration, gamma)
        
        if current_time + next_duration > target_end_time
            rest = target_end_time - current_time
            if rest < min_duration && step_counter > 1
                pred_h[end] += rest
                pred_X[end] = target_end_time
                current_time = target_end_time
                break
            else
                next_duration = rest
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
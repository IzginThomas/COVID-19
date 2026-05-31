module WenoPlotter

using PyCall

export plot_weno_parameters

# 5-point Gauss-Legendre quadrature to calculate cell average
function compute_cell_average(itp, x_l::Float64, x_r::Float64)
    nodes = [-0.906179845938664, -0.538469310105683, 0.0, 0.538469310105683, 0.906179845938664]
    weights = [0.236926885056189, 0.478628670499366, 0.568888888888889, 0.478628670499366, 0.236926885056189]
    mid = 0.5 * (x_l + x_r)
    half_width = 0.5 * (x_r - x_l)
    integral = 0.0
    for idx in 1:5
        t_val = mid + half_width * nodes[idx]
        integral += weights[idx] * itp(t_val)
    end
    return 0.5 * integral
end

"""
    plot_weno_parameters(midpoints::Vector{Float64}, param_matrix::Matrix{Float64}, interpolated_funcs::AbstractVector; opt_positions_filtered=nothing, original_values=nothing)

Plottet die WENO3-Interpolation für jeden optimierten Parameter.
Nutzt PyCall, um direkt auf das pyplot-Backend zuzugreifen.
Fügt den konstanten Originalwert als rote gestrichelte Linie hinzu.
"""
function plot_weno_parameters(midpoints::Vector{Float64}, param_matrix::Matrix{Float64}, interpolated_funcs::AbstractVector, weno_k::Int; opt_positions_filtered=nothing, original_values=nothing, backup_id::String="backupI", end_time::Float64=200.0)
    # PyPlot-Backend über PyCall importieren
    plt = pyimport("matplotlib.pyplot")
    
    num_opt_vars = size(param_matrix, 1) # Zeilen = Anzahl der Parameter
    N_cells = length(midpoints)
    
    # Aus den physikalischen Zellmitten (midpoints) die Zellgrenzen (Faces) konsistent berechnen
    X_f = Vector{Float64}(undef, N_cells + 1)
    X_f[1] = 0.0  # Annahme des Koordinatenursprungs
    for i in 1:N_cells
        X_f[i+1] = 2.0 * midpoints[i] - X_f[i]
    end
    
    max_per_fig = 4
    fig = nothing
    axs = nothing
    
    for j in 1:num_opt_vars
        # Index innerhalb der aktuellen Figure bestimmen (1 bis 4)
        local_idx = ((j - 1) % max_per_fig) + 1
        
        # Falls wir am Anfang einer neuen Figure stehen: Erstellen!
        if local_idx == 1
            fig, axs = plt.subplots(2, 2, figsize=(12, 8))
            axs = reshape(axs, 4) # Flachdrücken auf 1D-Array mit 4 Elementen
        end
        
        # 1. Diskrete Stützstellen aus der param_matrix holen (aktuelle Zeile j)
        discrete_y = param_matrix[j, :]
        
        # 2. Kontinuierliche Kurve abrufen
        itp = interpolated_funcs[j]
        
        # Aktuellen Achsen-Slot auswählen
        ax = axs[local_idx]
        
        # 3. Step-Koordinaten für originale Zellmittelwerte konstruieren
        step_x = Float64[]
        step_y = Float64[]
        for i in 1:N_cells
            push!(step_x, X_f[i])
            push!(step_y, discrete_y[i])
            push!(step_x, X_f[i+1])
            push!(step_y, discrete_y[i])
        end
        
        # 4. Zellweise Auswertung & Berechnung der rekonstruierten Mittelwerte
        rec_avgs = Vector{Float64}(undef, N_cells)
        for i in 1:N_cells
            # Gitterpunkte innerhalb der Zelle i (leicht nach innen verschoben)
            t_cell = range(X_f[i] + 1e-9, X_f[i+1] - 1e-9, length=50)
            interp_y_cell = [itp(t) for t in t_cell]
            
            # Segment der WENO-Kurve plotten
            ax.plot(t_cell, interp_y_cell, color="#2e7d32", linewidth=2.0, 
                    label=(i == 1 ? "WENO Reconstruction (k=$(weno_k))" : ""), zorder=4)
            
            # Berechne den Zellmittelwert der rekonstruierten Funktion
            rec_avgs[i] = compute_cell_average(itp, X_f[i], X_f[i+1])
        end
        
        # 5. Rekonstruierte Zellmittelwerte als gestrichelte Linie drüberlegen
        step_rec_x = Float64[]
        step_rec_y = Float64[]
        for i in 1:N_cells
            push!(step_rec_x, X_f[i])
            push!(step_rec_y, rec_avgs[i])
            push!(step_rec_x, X_f[i+1])
            push!(step_rec_y, rec_avgs[i])
        end
        
        max_err = maximum(abs.(rec_avgs .- discrete_y))
        err_str = max_err < 1e-11 ? "< 10⁻¹¹" : "$(round(max_err, sigdigits=3))"
        
        # --- Plot-Befehle via PyCall ---
        ax.plot(step_x, step_y, color="#1976d2", linewidth=3.0, zorder=5, label="Target Averages")
        ax.plot(step_rec_x, step_rec_y, color="#f57c00", linestyle="--", linewidth=1.5, zorder=6, 
                label="Rec. Averages (Err: $(err_str))")
        
        # 6. Sprünge an Zellgrenzen visualisieren
        first_jump = true
        for i in 2:N_cells
            y_left = itp(X_f[i] - 1e-9)
            y_right = itp(X_f[i] + 1e-9)
            jump_size = abs(y_left - y_right)
            
            if jump_size > 1e-5
                # Rote gepunktete Linie an der Grenze zeigt den Sprung
                ax.plot([X_f[i], X_f[i]], [y_left, y_right], color="#d32f2f", linestyle=":", linewidth=1.5, zorder=3,
                        label=(first_jump ? "Discontinuity (Jump)" : ""))
                first_jump = false
            end
            
            # Grenzwerte an der Schnittstelle markieren
            ax.scatter([X_f[i], X_f[i]], [y_left, y_right], color="#2e7d32", s=15, zorder=7, facecolors="none", edgecolors="#2e7d32")
        end
        
        # Vertikale Zellgrenzen einzeichnen (dezent gepunktet)
        for face in X_f
            ax.axvline(x=face, color="gray", linestyle=":", linewidth=0.8, alpha=0.3, zorder=1)
        end
        
        # Konstanten Originalwert aus x_model als horizontale Linie einzeichnen
        if original_values !== nothing
            ax.axhline(y=original_values[j], color="red", linestyle="--", linewidth=1.5, label="Initial Value")
        end
        
        # Titel setzen (Das "Var:" wurde hier komplett entfernt)
        name = opt_positions_filtered !== nothing ? "$(opt_positions_filtered[j])" : "Variable \$j\$"
        ax.set_title(name, fontsize=14, fontweight="bold")
        ax.set_xlabel("Days t", fontsize=12)
        ax.set_ylabel("Parameter", fontsize=12)
        ax.grid(true, linestyle=":", alpha=0.5)
        ax.legend(fontsize=8, loc="best")
        
        # Layout anpassen, wenn Figure voll ist oder am Ende angekommen ist
        if local_idx == max_per_fig || j == num_opt_vars
            if j == num_opt_vars && local_idx < max_per_fig
                for empty_idx in (local_idx + 1):max_per_fig
                    axs[empty_idx].set_visible(false)
                end
            end
            
            plt.tight_layout()
            
            # Construct a dynamic filename based on backup_id and end_time
            fig_num = div(j-1, max_per_fig) + 1
            fig_name = "weno_parameters_$(backup_id)_$(end_time)_fig$(fig_num).png"
            plt.savefig(fig_name, dpi=300)
            println("Saved parameter plot to $fig_name")
            
            if !isdefined(Main, :AntigravityHeadless)
                plt.show() 
            end
        end
    end
end

end # module

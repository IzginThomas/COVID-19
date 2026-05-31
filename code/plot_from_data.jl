using DelimitedFiles
using LinearAlgebra # Wichtig für pinv() in der neuen WENO-Funktion

if !@isdefined(BACKUP_ID)
    BACKUP_ID = "backupI"
end

base_p = all_models[Model_no];
println("--- Lade Optimierungsdaten aus den CSV-Backups (BACKUP_ID: $BACKUP_ID)... ---")
midpoints = vec(DelimitedFiles.readdlm("midpoints_$(BACKUP_ID).csv", ',', Float64))
param_matrix = DelimitedFiles.readdlm("param_matrix_$(BACKUP_ID).csv", ',', Float64)
println("Daten erfolgreich geladen!")

println("--- Rekonstruiere Zeit-Interpolationsfunktionen ---")
interpolated_funcs = []
num_opt_vars = length(opt_positions_filtered)

# --------------------------------===================================
# NEU: Gitter-Metrik aus den Mittelpunkten für das WENO-Modul ableiten
# --------------------------------===================================
# Da wir nur die Zellmittelpunkte gespeichert haben, rekonstruieren wir die Zellgrenzen (Interfaces).
# Für ein gleichmäßiges oder fast gleichmäßiges Optimierungs-Zeitgitter:
N_cells = length(midpoints)

if N_cells > 1
    # Berechne die Schrittweite zwischen den Mitten
    dt_mid = diff(midpoints)
    # Verwende die erste Differenz als Näherung für die Rand-Zellweiten
    dt_start = dt_mid[1]
    dt_end = dt_mid[end]
    
    # Rekonstruktion der Interfaces (Zellgrenzen)
    X_interfaces = zeros(N_cells + 1)
    X_interfaces[1] = 0.0
     for i in 1:N_cells
        X_interfaces[i+1] = 2.0 * midpoints[i] - X_interfaces[i]
    end
    
    # Exakte Zellweiten delta_x berechnen
    h_widths = diff(X_interfaces)
else
    X_interfaces = [midpoints[1] - 0.5, midpoints[1] + 0.5]
    h_widths = [1.0]
end

if !@isdefined(end_time)
    end_time = length(X_interfaces) > 1 ? X_interfaces[end] : 200.0
end





using Interpolations # Falls für spline-basierte Ansätze gewünscht

function evaluate_piecewise(t::Float64, midpoints::Vector{Float64}, h_widths::Vector{Float64}, 
                            X_interfaces::Vector{Float64}, y_values::Vector{Float64}, k::Int; 
                            method_type::Symbol=:weno, apply_zhang_shu::Bool=true, 
                            lower_bound::Float64=0.0)
    
    N = length(y_values)
    
    # 1. ZENTRALE INDEX-BESTIMMUNG
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
            # --- MATHEMATISCH KORREKTE KONSERVATIVE REKONSTRUKTION ---
            # Wir berechnen die Interface-Werte V_face so, dass Stetigkeit 
            # UND Mittelwerterhaltung über das Gleichungssystem gesichert sind.
            V_face = zeros(N + 1)
            
            # Da wir ein ungleichmäßiges Gitter haben, nutzen wir die Standard-
            # FVD-Herleitung (Finite Volume Reconstruction) für stetige lineare Profile:
            for j in 2:N
                # Abstand von der Mitte der linken Zelle zur Mitte der rechten Zelle
                dx_m2m = midpoints[j] - midpoints[j-1]
                # Gewichtung basierend auf den Abständen der Mittelpunkte zu den Grenzflächen
                V_face[j] = ((midpoints[j] - X_interfaces[j]) * y_values[j-1] + 
                             (X_interfaces[j] - midpoints[j-1]) * y_values[j]) / dx_m2m
            end
            
            # Jetzt zwingen wir die Randzellen zur exakten Mittelwerterhaltung
            V_face[1]   = 2.0 * y_values[1] - V_face[2]
            V_face[end] = 2.0 * y_values[end] - V_face[end-1]
            
            # Korrektur-Schritt für das Innere:
            # Damit in Zelle i der Mittelwert EXACT stimmt, muss der Verlauf korrigiert werden.
            # Ein rein interpoliertes V_face sichert Stetigkeit, verfälscht aber den Mittelwert leicht.
            # Wir verschieben die Linie parallel, um den exakten Mittelwert zu garantieren:
            v_L_temp = V_face[i]
            v_R_temp = V_face[i+1]
            temp_avg = 0.5 * (v_L_temp + v_R_temp)
            
            # Der Fehler zum gewünschten Zellmittelwert
            delta_avg = y_values[i] - temp_avg
            
            # Die korrigierten, exakt mittelwerterhaltenden Endpunkte für diese Zelle
            v_L = v_L_temp + delta_avg
            v_R = v_R_temp + delta_avg
            
            # Lineare Funktion im Intervall
            t_L = X_interfaces[i]
            v_raw = v_L + (v_R - v_L) * (t - t_L) / h_widths[i]
        end
        
   # ... (Zentrale Indexbestimmung bleibt gleich) ...

    elseif method_type == :constant
        if t <= X_interfaces[1]
            v_raw = y_values[1]
        elseif t >= X_interfaces[end]
            # ----------------------------------------------------------------
            # ZUKUNFTS-EXTRAPOLATION (t > end_time)
            # ----------------------------------------------------------------
            # Option A: Status Quo halten (Wert der letzten optimierten Zelle)
            v_raw = y_values[end]
            
            # Option B (Alternativ): Rückkehr zum statischen Basiswert aus base_p
            # Dafür müsstest du den Basiswert als 'default_value' an die Funktion übergeben:
            # v_raw = default_value 
        else
            # Normaler Bereich innerhalb des Daten-Gitters
            v_raw = y_values[i]
        end
    else
        error("Unbekannter Interpolationstyp: :$(method_type)")
    end

    # 3. ZENTRALER ZHANG-SHU LIMITER
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




println("--- Rekonstruiere Zeit-Interpolationsfunktionen ---")
interpolated_funcs = []
num_opt_vars = length(opt_positions_filtered)

# --- CONFIG-BEREICH FÜR DIE REKONSTRUKTION ---
gewaehlte_methode = :weno  # Optionen: :weno, :linear, :constant
weno_k           = 2     # k=2 für WENO-3, k=3 für WENO-5
limit_positivity = true   # Zhang-Shu-Limiter aktivieren

for j in 1:num_opt_vars
    y_values = param_matrix[j, :]
    
    if length(midpoints) > 1
        # Die anonyme Funktion greift jetzt dynamisch auf unseren Allround-Wrapper zu
        itp = t -> evaluate_piecewise(
            t, 
            midpoints, 
            h_widths, 
            X_interfaces, 
            y_values, 
            weno_k; 
            method_type = gewaehlte_methode,
            apply_zhang_shu = limit_positivity, 
            lower_bound = 0.0
        )
    else
        itp = t -> y_values[1]
    end
    push!(interpolated_funcs, itp)
end
println("Interpolationsfunktionen erfolgreich mit Methode [:", gewaehlte_methode, "] erstellt!")

# Lokale Wrapper für die nicht-autonomen Dynamiken
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

println("--- Löse die finale, nicht-autonome Differentialgleichung aus Backup ---")
if !@isdefined(original_u0)
    original_u0 = [30416000.0, 0.0, 5.0, 5.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
end
if !@isdefined(method)
    method = MPRK22(1.0)
end

final_prob_load = PDSProblem(P_non_autonomous_load, d_non_autonomous_load, original_u0, (0.0, end_time), base_p)
final_sol_load = solve(final_prob_load, method)

# Plot aufrufen (nutzt das in test_modul_bayesian_tools.jl definierte plot_res)
local_sols = @isdefined(local_solutions_history) ? local_solutions_history : nothing
plot_res(final_sol_load, base_p, interpolated_funcs, local_sols=local_sols)

# Save the final non-autonomous plot
fig_name = "final_non_autonomous_$(BACKUP_ID)_$(end_time).png"
plt.savefig(fig_name, dpi=300)
println("Saved final non-autonomous plot to $fig_name")
if !isdefined(Main, :AntigravityHeadless)
    plt.show()
end


# # 1. Hilfsfunktion zur Erzeugung valider LaTeX-Strings für Matplotlib
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

    # Wenn der Name extrem kurz ist (z. B. "μ", "γ", "x")
    if len < 3
        return "\$$(name)\$"
    end

    name_str = String(chars)

    # Regel für Variablen mit "Ve" am Ende (z. B. "rIVe", "rLVe", "μVe")
    if endswith(name_str, "Ve")
        if len >= 4
            base = string(chars[1])            # z. B. 'r'
            sup = String(chars[2:(end-2)])     # z. B. 'I' oder 'L'
            return "\$$(base)^$(sup)_{V_E}\$"
        else
            # Für 3-Zeichen-Parameter wie "μVe" -> \mu_{V_E}
            return "\$$(String(chars[1:(end-2)]))_{V_E}\$"
        end
    end

    # Regel für standardmäßige 3+ Zeichen-Variablen (z. B. "aEL" -> a^E_L)
    m = match(r"^([a-z])([A-Z])([A-Z])$", name_str)
    if m !== nothing
        return "\$$(m.captures[1])^$(m.captures[2])_$(m.captures[3])\$"
    end

    return "\$$(name_str)\$"
end

# 2. Originalliste aller Parameternamen definieren
parameters = [
    "k1", "k2", "cS", "cV", "cE", "cL", "cI", "cH", "cR", "cQ",
    "cVe", "δV", "δE", "δL", "δH", "δQ", "pv", "Λ", "μ", "φ",
    "ψ", "γ", "μVe", "αI", "αH", "aSV", "aVS", "aVE", "aEI", "aES",
    "aEQ", "aEL", "aLI", "aLR", "aLQ", "rLVe", "aIR", "aIQ", "aIH", "rIVe",
    "aHR", "aRS", "aQR", "b01","b02","b03","b04","b11","b12","b13", "b14","ω1", "ω2"
]

# 3. Die numerischen Originalwerte aus x_model für die optimierten Indizes ziehen
original_vals = x_model

# 4. Text-Labels aus dem Namens-Array holen und in LaTeX transformieren
raw_opt_names = parameters[opt_positions_filtered]
formatted_opt_names = format_param_latex.(raw_opt_names)

# # 5. Modul einbinden und den Plotter aufrufen
 include("WenoPlotter.jl")
 using .WenoPlotter: plot_weno_parameters

WenoPlotter.plot_weno_parameters(
    midpoints, 
    param_matrix, 
    interpolated_funcs,
    weno_k; 
    opt_positions_filtered = formatted_opt_names, # Übergabe der reinen LaTeX-Namen für den Titel
    original_values = original_vals,               # Übergabe der Konstanten für die axhline
    backup_id = BACKUP_ID,
    end_time = end_time
)

##

# ... (Laden der Backups wie gehabt) ...

# Bindest das neue Modul ein
include("ParameterPredictor.jl")
using .ParameterPredictor: predict_until_time

# --- DEINE WUNSCH-ENDZEIT FÜR DIE VORHERSAGE ---
wunsch_prognose_ende = 223.0  # Wie weit soll die Simulation INSGESAMT laufen? (z.B. 120 Tage)

# Führe die iterative Prognose aus
#wunsch_prognose_ende = 120.0  # Ziel-Zeitpunkt in Tagen (Gesamtlaufzeit)
L_nachbarn           = 1     # Wie viele historische Muster gemischt werden sollen
mindest_zeitfenster  = 1.0    # Jedes Zukunfts-Intervall ist mindestens 7 Tage lang
trend_gewicht        = 1.5

expanded_param_matrix, expanded_h_widths, expanded_X_interfaces = predict_until_time(
    param_matrix, 
    h_widths, 
    X_interfaces, 
    wunsch_prognose_ende; # <-- Semikolon beachten!
    L = L_nachbarn, 
    min_duration = mindest_zeitfenster, 
    gamma = trend_gewicht
)

# Da wir nun ein verlängertes, echtes Gitter haben, berechnen wir die erweiterten Mittelpunkte 
# für die interne Auswertung (falls benötigt, z.B. für Plots)
expanded_midpoints = [0.5 * (expanded_X_interfaces[i] + expanded_X_interfaces[i+1]) for i in 1:length(expanded_h_widths)]


# --- ERSTELLUNG DER FUNKTIONEN (Nutzt jetzt die erweiterten Daten!) ---
println("--- Erstelle zeit-kontinuierliche Wrapper für die ODE ---")
interpolated_funcs = []
num_opt_vars = length(opt_positions_filtered)

for j in 1:num_opt_vars
    # Wir übergeben der evaluate_piecewise die ERWEITERTEN Arrays aus dem Predictor
    y_values_expanded = expanded_param_matrix[j, :]
    
    itp = t -> evaluate_piecewise(
        t, 
        expanded_midpoints, 
        expanded_h_widths, 
        expanded_X_interfaces, 
        y_values_expanded, 
        weno_k; 
        method_type = :constant, # Da stückweise konstant am besten ist
        apply_zhang_shu = limit_positivity, 
        lower_bound = 0.0
    )
    push!(interpolated_funcs, itp)
end


# --- LÖSE DIE ODE BIS ZUR PROGNOSE-ENDZEIT ---
# println("--- Löse finale ODE mit prognostiziertem Parameter-Verlauf ---")
# final_prob_load = PDSProblem(
#     P_non_autonomous_load, 
#     d_non_autonomous_load, 
#     original_u0, 
#     (0.0, wunsch_prognose_ende), # <-- Hier deine neue Zielzeit
#     base_p
# )
# final_sol_load = solve(final_prob_load, method)

# Plotten (Der Plotter visualisiert automatisch das gesamte verlängerte Fenster)
# plot_res(final_sol_load, base_p, interpolated_funcs, local_sols=local_solutions_history)


# WenoPlotter.plot_weno_parameters(
#     expanded_midpoints, 
#     expanded_param_matrix, 
#     interpolated_funcs; 
#     opt_positions_filtered = formatted_opt_names, # Übergabe der reinen LaTeX-Namen für den Titel
#     original_values = original_vals                # Übergabe der Konstanten für die axhline
# )
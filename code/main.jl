include("opt_vars.jl");
include("bayesian_tools.jl");
include("get_data_new.jl");
include("weno.jl")
# include("weno_test.jl")

using .WenoInterpolation: weno_evaluate_non_uniform
using DelimitedFiles
using Revise
using PositiveIntegrators
using DataFrames
using Parameters
using PyPlot
using LinearAlgebra
using Interpolations
using Dates
using Trapz
# u0 erweitert um Index 11 (D_covid startet bei 0.0)
u0 = [30416000.0, 0.0, 5.0, 5.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]; # will be updated in optimization 
original_u0 = u0; # used for non-autonomous case to keep the same initial conditions
comp_labels = ["S", "V", "E", "L", "I", "H", "R", "Q", "D", "Ve", "D_covid"];
maxiteration = 150;
const BACKUP_ID = "backupI"; # identifier for backup files to distinguish from other runs

#tspan = (0.0, tref_IRHD[end]); # replace with real data time span
method = MPRK22(1.0);

# NEW
include("setup_model.jl"); # if included, check cost function (e.g. cummulative data?)

function getplotvar(t, u, p_all, idx)
    # u is of type Vector{Vector{Float64}} where each inner vector represents the compartment's values over time
    # p_all is a vector of parameters
    # idx is the index of the compartment to plot
    S = getindex.(u, 1)
    V = getindex.(u, 2)
    E = getindex.(u, 3)
    L = getindex.(u, 4)
    I = getindex.(u, 5)
    H = getindex.(u, 6)
    R = getindex.(u, 7)
    Q = getindex.(u, 8)
    D = getindex.(u, 9)
    Ve = getindex.(u, 10)
    D_covid = getindex.(u, 11)
    k1, k2, cS, cV, cE, cL, cI, cH, cR, cQ,
    cVe, δV, δE, δL, δH, δQ, pv, Λ, μ, φ,
    ψ, γ, μVe, αI, αH, aSV, aVS, aVE, aEI, aES,
    aEQ, aEL, aLI, aLR, aLQ, rLVe, aIR, aIQ, aIH, rIVe,
    aHR, aRS, aQR, b01, b02, b03, b11, b12, b13, ω = p_all


    if idx == 1 && flag_N # flags are in setup_model.jl
        y = (sum(u, dims=1))' - u[10, :] # sum over the compartments
    elseif (idx == 2 && flag_V && flag_cum_V) || (idx == 5 && flag_I && flag_cum_I)  # cumulative vaccinated/infected population using trapezoidal rule
        if idx == 5
            rate = aEI .* φ .* E + aLI .* γ .* L # production terms for I
        else
            rate = aSV .* S .+ pv .* Λ
        end
        y = [k == 1 ? rate[1] : trapz(t[1:k], rate[1:k]) for k in 1:length(t)]
        #elseif idx == 9 && flag_D
        #   rate = αI * I + αH * H
        #  y = [k == 1 ? rate[1] : trapz(t[1:k], rate[1:k]) for k in 1:length(t)]
    elseif idx == 5 && flag_I
        y = I # Infected (I) includes both I and L compartments
    else
        y = getindex.(u, idx)
    end
    return y
end


# END NEW
function P(u, p, t)
    S, V, E, L, I, H, R, Q, D, Ve, D_covid = u  # Erweitert um D_covid
    k1, k2, cS, cV, cE, cL, cI, cH, cR, cQ,
    cVe, δV, δE, δL, δH, δQ, pv, Λ, μ, φ,
    ψ, γ, μVe, αI, αH, aSV, aVS, aVE, aEI, aES,
    aEQ, aEL, aLI, aLR, aLQ, rLVe, aIR, aIQ, aIH, rIVe,
    aHR, aRS, aQR, b01, b02, b03, b04, b11, b12, b13, b14, ω1, ω2 = p

    pp = [cS, cV, cE, cL, cI, cH, cR, cQ, 0.0, cVe, 0.0] # Erweitert auf 11
    β0 = b01 * (1 + (b02 - 1)^2 * cos(pi * t / (365.0 * b03 + 1) + (b04 - pi)^2))
    β1 = (b11 - 1)^2 * (1 + (b12 - 1)^2 * cos(pi * t / (365.0 * b13 + 1) + (b14 - pi)^2))
    aSE = sum(u .* pp) * (sum(u))^round(k1) * β0 * (1 .+ β1 * cos(pi * t / (365.0 * ω1 + 1) + (ω2 - pi)^2))
    aVE = aVE * aSE^round(k2)

    Pmat = zeros(eltype(u), 11, 11) # Matrix auf 11x11 vergrößert

    # ==================== Production-rest terms ==================== #
    r1p = (1 - pv) * Λ
    r2p = pv * Λ
    r10p = rLVe * L + rIVe * I
    Pmat[1, 1] = r1p
    Pmat[2, 2] = r2p
    Pmat[10, 10] = r10p

    # ==================== Production terms ==================== #
    Pmat[1, 2] = aVS * V
    Pmat[1, 3] = aES * E
    Pmat[1, 7] = aRS * R
    Pmat[2, 1] = aSV * S
    Pmat[3, 1] = aSE * S
    Pmat[3, 2] = aVE * V
    Pmat[4, 3] = aEL * (1 - φ) * E
    Pmat[5, 3] = aEI * φ * E
    Pmat[5, 4] = aLI * γ * L
    Pmat[6, 5] = aIH * (1 - ψ) * I
    Pmat[7, 4] = aLR * (1 - γ) * L
    Pmat[7, 5] = aIR * ψ * I
    Pmat[7, 6] = aHR * H
    Pmat[7, 8] = aQR * Q
    Pmat[8, 3] = aEQ * E
    Pmat[8, 4] = aLQ * L
    Pmat[8, 5] = aIQ * I

    # Zustand 9 (D) behält die allgemeine/natürliche Mortalität
    Pmat[9, 1] = μ * S
    Pmat[9, 2] = μ * round(δV) * V
    Pmat[9, 3] = μ * round(δE) * E
    Pmat[9, 4] = μ * round(δL) * L
    Pmat[9, 5] = μ * I
    Pmat[9, 6] = μ * round(δH) * H
    Pmat[9, 7] = μ * R
    Pmat[9, 8] = μ * round(δQ) * Q

    # 🔥 NEU: Kompartiment 11 (D_covid) sammelt NUR die COVID-Todesfälle aus I und H
    Pmat[11, 5] = αI * I
    Pmat[11, 6] = αH * H

    return Pmat
end;
#=
function P(u, p, t)
    S, V, E, L, I, H, R, Q, D, Ve = u
    k1, k2, cS, cV, cE, cL, cI, cH, cR, cQ,
    cVe, δV, δE, δL, δH, δQ, pv, Λ, μ, φ,
    ψ, γ, μVe, αI, αH, aSV, aVS, aVE, aEI, aES,
    aEQ, aEL, aLI, aLR, aLQ, rLVe, aIR, aIQ, aIH, rIVe,
    aHR, aRS, aQR, b01,b02,b03,b04,b11,b12,b13,b14,ω1,ω2 = p

    pp = [cS, cV, cE, cL, cI, cH, cR, cQ, 0.0, cVe] 
    β0 = b01 * (1 + (b02-1)^2*cos(pi*t/(365.0*b03 + 1) + (b04-pi)^2)); # default: b01 = 1, b02 = 1, b03 = 1, b04 = 1
    β1 = (b11-1)^2 * (1 + (b12-1)^2*cos(pi*t/(365.0*b13 + 1) + (b14-pi)^2)); # default: b11 = 1, b12 = 1, b13 = 1, b14 = 1
    aSE = sum(u .* pp) * (sum(u))^round(k1) * β0 *(1 .+ β1*cos(pi*t/(365.0*ω1 + 1) + (ω2-pi)^2)); # default ω1= 1, ω2 = 1
    aVE = aVE * aSE^round(k2)

    Pmat = zeros(eltype(u), 10, 10) 

    # ==================== Production-rest terms ==================== #
    r1p = (1 - pv) * Λ
    r2p = pv * Λ
    r10p = rLVe * L + rIVe * I
    Pmat[1, 1] = r1p
    Pmat[2, 2] = r2p
    Pmat[10, 10] = r10p

    # ==================== Production terms ==================== #
    Pmat[1, 2] = aVS * V
    Pmat[1, 3] = aES * E
    Pmat[1, 7] = aRS * R
    Pmat[2, 1] = aSV * S
    Pmat[3, 1] = aSE * S
    Pmat[3, 2] = aVE * V
    Pmat[4, 3] = aEL * (1 - φ) * E
    Pmat[5, 3] = aEI * φ * E
    Pmat[5, 4] = aLI * γ * L
    Pmat[6, 5] = aIH * (1 - ψ) * I
    Pmat[7, 4] = aLR * (1 - γ) * L
    Pmat[7, 5] = aIR * ψ * I
    Pmat[7, 6] = aHR * H
    Pmat[7, 8] = aQR * Q
    Pmat[8, 3] = aEQ * E
    Pmat[8, 4] = aLQ * L
    Pmat[8, 5] = aIQ * I
    Pmat[9, 1] = μ * S
    Pmat[9, 2] = μ * round(δV) * V
    Pmat[9, 3] = μ * round(δE) * E
    Pmat[9, 4] = μ * round(δL) * L
    Pmat[9, 5] = (μ + αI) * I
    Pmat[9, 6] = (μ * round(δH) + αH) * H
    Pmat[9, 7] = μ * R
    Pmat[9, 8] = μ * round(δQ) * Q

    return Pmat
end;
=#

function d(u, p, t)
    S, V, E, L, I, H, R, Q, D, Ve, D_covid = u
    k1, k2, cS, cV, cE, cL, cI, cH, cR, cQ,
    cVe, δV, δE, δL, δH, δQ, pv, Λ, μ, φ,
    ψ, γ, μVe, αI, αH, aSV, aVS, aVE, aEI, aES,
    aEQ, aEL, aLI, aLR, aLQ, rLVe, aIR, aIQ, aIH, rIVe,
    aHR, aRS, aQR, b01, b02, b03, b04, b11, b12, b13, b14, ω1, ω2 = p
    res = zero(u)
    res[10] = μVe * Ve
    return res
end;


function cost(x)
    weights = ones(11) # weights for every compartment in the error calculation
    #weights[7] = 1.0; # weight for compartment R (recovered) set to one
    #weights[10] = 1.0; # weight for compartment S (susceptible) set to one

    params = ["k1", "k2", "cS", "cV", "cE", "cL", "cI", "cH", "cR", "cQ",
        "cVe", "δV", "δE", "δL", "δH", "δQ", "pv", "Λ", "μ", "φ",
        "ψ", "γ", "μVe", "αI", "αH", "aSV", "aVS", "aVE", "aEI", "aES",
        "aEQ", "aEL", "aLI", "aLR", "aLQ", "rLVe", "aIR", "aIQ", "aIH", "rIVe",
        "aHR", "aRS", "aQR", "b01", "b02", "b03", "b04", "b11", "b12", "b13", "b14", "ω1", "ω2"]
    # NEW
    x_all = all_models[Model_no] # get the full parameter vector for the chosen model
    #x_all = zeros(length(parameters)) # initialize full parameter vector with zeros
    #x_all[z_positions] .= 0.0; # set zero parameters to zero
    #x_all[nz_delta] .= 1.0
    x_all[opt_positions_filtered] = x # update only the parameters from the optimization

    # END NEW
    tab = DataFrame(
        Parameter=params,
        Model=x_all,
    )
    paras = Dict(Symbol.(tab.Parameter) .=> tab[!, Symbol("Model")]) # store all values

    @unpack k1, k2, cS, cV, cE, cL, cI, cH, cR, cQ,
    cVe, δV, δE, δL, δH, δQ, pv, Λ, μ, φ,
    ψ, γ, μVe, αI, αH, aSV, aVS, aVE, aEI, aES,
    aEQ, aEL, aLI, aLR, aLQ, rLVe, aIR, aIQ, aIH, rIVe,
    aHR, aRS, aQR, b01, b02, b03, b04, b11, b12, b13, b14, ω1, ω2 = paras


    p = (k1, k2, cS, cV, cE, cL, cI, cH, cR, cQ,
        cVe, δV, δE, δL, δH, δQ, pv, Λ, μ, φ,
        ψ, γ, μVe, αI, αH, aSV, aVS, aVE, aEI, aES,
        aEQ, aEL, aLI, aLR, aLQ, rLVe, aIR, aIQ, aIH, rIVe,
        aHR, aRS, aQR, b01, b02, b03, b04, b11, b12, b13, b14, ω1, ω2)
    prob = PDSProblem(P, d, u0, tspan, p)
    #if isdefined(Main, :dt0)
    sol = solve(prob, method)
    #else
    #    sol = solve(prob, method)
    #end;

    # u =  S, V, E, L, I, H, R, Q, D, Ve 
    # tref_V = 2 , tref_S = 1, tref_IRHD = 5, 6,7,9
    # error calculation
    err = zeros(eltype(weights), 11)
    for i in 1:11
        # default values for error computation
        tref = sol.t
        solref = getindex.(sol.u, i) # default to obtain zero error if flag_"Compartment" is false.

        if i == 1 && flag_N
            tref = tref_N
            solref = refN
        elseif i == 2 && flag_V
            tref = tref_V
            solref = refV
        elseif i == 5 && flag_I
            #tref = tref_IRHD
            tref = tref_I
            solref = refI
        elseif i == 6 && flag_H
            #tref = tref_IRHD
            tref = tref_H
            solref = refH
        elseif i == 7 && flag_R
            #tref = tref_IRHD
            tref = tref_R
            solref = refR
        elseif i == 11 && flag_D
            #tref = tref_IRHD
            tref = tref_D
            solref = refD
        end
        valid_indices = findall(x -> x >= tspan[1] && x <= tspan[end], tref)
        if !isempty(valid_indices)
            tref = tref[valid_indices]
            solref = solref[valid_indices]
        end
        dense_numsol = hcat(sol.(tref)) # approximations at real data time points using linear interpolation
        num_sol = getplotvar(tref, dense_numsol, p, i)
        ## Note: solref = num_sol if no data is available # maybe add min.(round.(num_sol) .- solref,0) +
        aux = num_sol .- solref
        denominator = max.(solref, 1.0) # add 1 to avoid division by zero
        err[i] = weights[i] .* (norm(aux) / norm(denominator)) # compare dense_numsol[:,j] with real data solref[j]
        if isnan(err[i]) || isinf(err[i])
            err[i] = 1e16 # if no data is available for this compartment, set error to zero
        end
    end
    #trans_numsol = dense_numsol # some transformation if needed to compare with data (e.g. summation of a compartment over time)

    return sum(err)
end
println("Cost of Model 1: ", cost(x_model))

#=
D = 43; # number of parameters
# set bounds dependend on dimension

bounds = fill(1.0, (2, D));
bounds[1, :] = bounds[1, :] * 0.0;
#bounds[2, 23:D] = bounds[2, 23:D] * 1e-1; # TRY TO SET REALISTIC BOUNDS FOR THE PROBLEM, THIS GREATLY AFFECTS THE OPTIMIZATION RESULTS
bounds[1, 1] = -1.0; # lower bound for k1
#bounds[2,1] = 0.0; #  k1 = 0 for model 1
#bounds[2,2] = 0.0; #  k2 = 0 for model 1
bounds[2, 18] = 2e3;     # upper bound for Λ  
bounds[2, 19] = 1e-3;     # upper bound for μ 
#bounds[1,18] = 1319.294;
bounds = bounds[:, nz_positions_delta]
=#

# NEW
bounds = fill(1.0, (2, length(x_model)));
bounds[1, :] = x_model' * 0.0; # lower bounds set to zero
bounds[2, :] = min.(x_model' * 2.0, 1.0); # upper bounds set to twice the initial values, but not exceeding 1.0
# END NEW


function plot_res(sol, base_p, parameter_kurven=nothing; local_sols=nothing)
    t = sol.t
    t_min = minimum(t)
    t_max = maximum(t)

    plt.figure()
    colors = []

    # 1. Schleife: Simulationsergebnisse (WENO3 Kontinuierlich) plotten
    for idx in 1:length(sol.u[1])
        if !(idx == 5 && flag_I || idx == 11 && flag_D) # Nur plotten, wenn Daten vorhanden
            continue
        end

        # HIER REKONSTRUIEREN WIR DIE DYNAMISCHEN PARAMETER FÜR GETPLOTVAR:
        if parameter_kurven !== nothing
            num_opt_vars = length(opt_positions_filtered)
            y_num = Vector{Float64}(undef, length(t))
            for k in eachindex(t)
                tk = t[k]
                uk = sol.u[k]

                p_dynamic = copy(base_p)
                for j in 1:num_opt_vars
                    p_dynamic[opt_positions_filtered[j]] = parameter_kurven[j](tk)
                end

                full_res = getplotvar([tk], [uk], p_dynamic, idx)
                y_num[k] = full_res[1]
            end
        else
            # Klassischer, statischer Aufruf (falls keine Kurven übergeben wurden)
            y_num = getplotvar(t, sol.u, base_p, idx)
        end

        # Kontinuierliche WENO3-Lösung plotten
        line, = plt.plot(t, y_num, label=comp_labels[idx], linewidth=2)
        current_color = line.get_color()
        push!(colors, (idx, current_color))

        # --- NEU: STÜCKWEISE LOKALE TRAJEKTORIEN HIER ZUFÜGEN ---
        if local_sols !== nothing
            for (b_idx, sol_block) in enumerate(local_sols)
                # Da die Blöcke statische Parameter nutzen, werten wir sie direkt aus
                # Wir holen uns die Parameter aus dem Block (falls im Problem gespeichert)
                # oder nutzen das statische base_p, da die lokalen Lösungen das p bereits intern tragen.
                y_block = getplotvar(sol_block.t, sol_block.u, sol_block.prob.p, idx)

                # Nur beim ersten Block ein Label setzen, um die Legende nicht zu überladen
                block_label = (b_idx == 1) ? (comp_labels[idx] * " (Piecewise)") : ""

                # Gepunktet (":") in der exakt gleichen Farbe wie die kontinuierliche Kurve
                plt.plot(sol_block.t, y_block, linestyle="-", color="black", alpha=0.9, linewidth=1.5, label=block_label)
            end
        end
    end

    # 2. Schleife: Reale Referenzdaten filtern und plotten
    for (i, col) in colors
        if i == 5 && flag_I
            tref = tref_I
            solref = refI
        elseif i == 2 && flag_V
            tref = tref_V
            solref = refV
        elseif i == 6 && flag_H
            tref = tref_H
            solref = refH
        elseif i == 7 && flag_R
            tref = tref_R
            solref = refR
        elseif i == 11 && flag_D
            tref = tref_D
            solref = refD
        else
            continue
        end

        valid_indices = findall(x -> x >= t_min && x <= t_max, tref)

        if !isempty(valid_indices)
            tref_filtered = tref[valid_indices]
            solref_filtered = solref[valid_indices]
            data_label = comp_labels[i] * " (Data)"
            plt.plot(tref_filtered, solref_filtered, linestyle="--", color=col, alpha=0.5, label=data_label)
        end
    end

    xlabel("Days")
    ylabel("Population")
    if parameter_kurven !== nothing
        title("Measurable Compartments")
    else
        title("Measurable Compartments")
    end
    legend()
    PyPlot.grid(true)
    tight_layout()
end

##

"""
Führt die fensterbasierte Parameteroptimierung in variablen Schritten aus.
Wenn der minimale Fehlerwert eines Fensters über 0.25 liegt, wird das Intervall 
halbiert und die Berechnung wiederholt –- es sei denn, die Schrittweite ist bereits 3.0.
"""
function run_rolling_optimization_and_solve(cost_func, end_time::Float64, bounds_matrix, max_iter::Int, initial_step_size::Float64; resume_from_backup::Bool=false)
    # Speicher für optimierte Parametervektoren, zeitliche Mittelpunkte und Kosten
    optimized_params_history = Vector{Vector{Float64}}()
    midpoints = Float64[]
    costs_history = Float64[] # Historie für die akzeptierten Kosten
    
    # --- Speicher für die stückweisen lokalen Lösungen ---
    local_solutions_history = []

    # Startwert u0 wird initial vom globalen u0-Vektor kopiert
    current_u0 = copy(u0)

    t_start = 0.0
    current_step_size = initial_step_size
    loop_counter = 1

    # --- RESTART / RESUME LOGIC ---
    if resume_from_backup && isfile("param_matrix_$(BACKUP_ID).csv") && isfile("midpoints_$(BACKUP_ID).csv")
        println("--- Reactivating rolling optimization from backup ---")

        # Load previously saved data
        param_matrix_loaded = DelimitedFiles.readdlm("param_matrix_$(BACKUP_ID).csv", ',', Float64)
        midpoints = vec(DelimitedFiles.readdlm("midpoints_$(BACKUP_ID).csv", ',', Float64))
        
        # Falls ein Kosten-Backup existiert, laden, ansonsten mit NaN füllen
        if isfile("costs_$(BACKUP_ID).csv")
            costs_history = vec(DelimitedFiles.readdlm("costs_$(BACKUP_ID).csv", ',', Float64))
        else
            costs_history = fill(NaN, length(midpoints))
        end

        # Convert matrix columns back into vector of vectors
        for col in 1:size(param_matrix_loaded, 2)
            push!(optimized_params_history, param_matrix_loaded[:, col])
        end

        println("-> $(length(midpoints)) intervals successfully recovered.")

        # Reconstruct local solutions piece by piece to update current_u0 to current state
        println("--- Calculating current u0 state over loaded intervals... ---")
        for j in 1:length(midpoints)
            t_end_local = midpoints[j] + (j == 1 ? midpoints[1] : midpoints[j] - midpoints[j-1]) / 2.0
            if j == 1
                t_start_local = 0.0
                t_end_local = midpoints[1] * 2.0
            else
                t_start_local = t_start
                t_end_local = midpoints[j] + (midpoints[j] - t_start_local)
            end

            best_x_local = optimized_params_history[j]
            x_all_local = copy(all_models[Model_no])
            x_all_local[opt_positions_filtered] = best_x_local
            p_tuple_local = Tuple(x_all_local)

            prob_local = PDSProblem(P, d, current_u0, (t_start_local, t_end_local), p_tuple_local)
            sol_local = solve(prob_local, method)

            push!(local_solutions_history, sol_local)
            current_u0 = copy(sol_local.u[end])
            t_start = t_end_local
            loop_counter += 1
        end
        println("-> u0 successfully updated to state at day $t_start. Continuing.")
    else
        println("--- Starting adaptive rolling optimization (Initial step size: $initial_step_size days) ---")
    end

    ## Variables for the fallback strategy during step size reduction
    lowest_cost_at_this_window = Inf
    best_x_at_this_window = nothing
    t_end_for_best_x = 0.0

    while t_start < end_time
        t_end = min(t_start + current_step_size, end_time)

        println("\n[Interval $loop_counter] Checking window: [$t_start, $t_end] (Step size: $current_step_size)")
        global tspan = (t_start, t_end)
        global u0 = copy(current_u0)

        result, cnt = bayesian_tools.run_opt(bounds_matrix, max_iter, cost_func)
        best_x = result.observed_optimizer
        current_cost = cost_func(best_x)
        println("-> Optimization complete. Best cost in this window: $current_cost")
        if current_cost < lowest_cost_at_this_window
            lowest_cost_at_this_window = current_cost
            best_x_at_this_window = copy(best_x)
            t_end_for_best_x = t_end
        end
        min_step_size = 1.0
        max_cost_threshold = 0.25
        
        # --- ADAPTIVE SCHRITTWEITENSTEUERUNG + FALLBACK ---
        if current_cost > max_cost_threshold
            if current_step_size <= min_step_size
                println("⚠️ Cost is high ($current_cost > $max_cost_threshold), but minimum step size ($current_step_size) is reached. Forcing step execution!")

                # Fallback: Use the absolute best attempt at this position
                if lowest_cost_at_this_window < current_cost
                    println("🔄 Best attempt at this window had a cost of $lowest_cost_at_this_window. Restoring this state!")
                    best_x = best_x_at_this_window
                    t_end = t_end_for_best_x
                    current_cost = lowest_cost_at_this_window
                else
                    println("🎯 Current attempt at step size $min_step_size is the best at this window (despite high cost). Continuing.")
                end
            else
                new_step_size = max(min_step_size, current_step_size / 2.0)
                current_step_size = new_step_size
                println("⚠️ Cost too high! Reducing step size to $current_step_size days and repeating interval.")
                continue # Repeat the same t_start with a smaller window
            end
        end

        # --- SCHRITT AKZEPTIERT ---
        push!(midpoints, (t_start + t_end) / 2.0)
        push!(optimized_params_history, best_x)
        push!(costs_history, current_cost)

        # === INTERNE ABSICHERUNG: BACKUP DIREKT NACH JEDEM SCHRITT SCHREIBEN ===
        param_matrix_tmp = hcat(optimized_params_history...)
        DelimitedFiles.writedlm("param_matrix_$(BACKUP_ID).csv", param_matrix_tmp, ',')
        DelimitedFiles.writedlm("midpoints_$(BACKUP_ID).csv", midpoints, ',')
        DelimitedFiles.writedlm("costs_$(BACKUP_ID).csv", costs_history, ',')

        x_all = copy(all_models[Model_no])
        x_all[opt_positions_filtered] = best_x
        p_tuple = Tuple(x_all)

        prob = PDSProblem(P, d, current_u0, (t_start, t_end), p_tuple)
        sol = solve(prob, method)

        push!(local_solutions_history, sol)

        #plot_res(sol, x_all)
        weno_k = 2; # ordnung 2k-1 
        #order = 2*weno_k - 1     
        # --- ZWISCHENPLOT DES NICHT-AUTONOMEN PROBLEMS (AB 4 INTERVALLEN) ---
        if length(midpoints) >= 5
            println("📊 [Intermediate Update] At least 5 intervals accepted. Simulating non-autonomous problem from t = 0.0 to t = $t_end...")
            num_opt_vars = length(opt_positions_filtered)
            current_param_matrix = hcat(optimized_params_history...)

            # ----------------------------------------------------------------
            # REKONSTRUKTION DER GITTERMETRIK MIT FIXEM STARTPUNKT X_interfaces[1] = 0.0
            # ----------------------------------------------------------------
            N_curr = length(midpoints)
            dt_mid = diff(midpoints)
            
            X_interfaces = zeros(N_curr + 1)
            X_interfaces[1] = 0.0  # Strikt bei Null fixiert
            
            for k in 1:N_curr-1
                X_interfaces[k+1] = midpoints[k] + 0.5 * dt_mid[k]
            end
            
            dt_end = X_interfaces[end-1] - (midpoints[end-1] - 0.5 * (dt_mid[end-1]))
            X_interfaces[end] = midpoints[end] + 0.5 * dt_end
            
            h_widths = diff(X_interfaces)
            

            limit_positivity = true  

            temp_funcs = []
            for j in 1:num_opt_vars
                y_vals = current_param_matrix[j, :]
                itp = t -> weno_evaluate_non_uniform(
                    t, 
                    y_vals, 
                    h_widths, 
                    X_interfaces, 
                    weno_k; 
                    apply_zhang_shu = limit_positivity, 
                    lower_bound = 0.0
                )
                push!(temp_funcs, itp)
            end

            function P_temp(u, p_static, t)
                p_dynamic = copy(p_static)
                for j in 1:num_opt_vars
                    p_dynamic[opt_positions_filtered[j]] = temp_funcs[j](t)
                end
                return P(u, Tuple(p_dynamic), t)
            end

            function d_temp(u, p_static, t)
                p_dynamic = copy(p_static)
                for j in 1:num_opt_vars
                    p_dynamic[opt_positions_filtered[j]] = temp_funcs[j](t)
                end
                return d(u, Tuple(p_dynamic), t)
            end

            base_p = all_models[Model_no]
            temp_prob = PDSProblem(P_temp, d_temp, original_u0, (0.0, t_end), base_p)
            temp_sol = solve(temp_prob, method)

            plot_res(temp_sol, base_p, temp_funcs, local_sols=local_solutions_history)
        end
        # ------------------------------------------------------------------------
        current_u0 = copy(sol.u[end])
        t_start = t_end
        loop_counter += 1

        lowest_cost_at_this_window = Inf
        best_x_at_this_window = nothing
        current_step_size = initial_step_size
    end

    println("\n--- Rolling optimization complete. Returning data. ---")
    param_matrix = hcat(optimized_params_history...)

    return param_matrix, midpoints, local_solutions_history, costs_history
end
#tspan = (0.0, maximum(tref_I)); # WILL BE GLOBAL VARIABLE!
cost_func = cost; # cost function needs to be defined in the global scope for the optimization loop
# end_time = maximum(tref_I)
end_time = 200.0 # for testing purposes, set to 6 months
bounds_matrix = bounds
max_iter = maxiteration
initial_step_size = 1.0
# 1. Optimierung ausführen
param_matrix, midpoints, local_solutions_history, cost_history = run_rolling_optimization_and_solve(cost_func, end_time, bounds_matrix, max_iter, initial_step_size; resume_from_backup= true)
# 2. Definition der nicht-autonomen P- und d-Funktionen für das finale PDSProblem

include("plot_from_data.jl")  # <plots parameter curves from csv files
include("run_rolling_forecast.jl") # <runs forecasts with different extrapolation methods and prediction horizons
#include("run_forecasts.jl")


## predictive checks
# include("ForecastingTools.jl")
# using .ForecastingTools

# # --- RUN FIRST TEST: Frozen Parameters ---
# days_test1 = ForecastingTools.test_prediction_horizon(
#     cost, original_u0, all_models, Model_no, opt_positions_filtered, method, P, d,
#     180.0, end_time, 0.25
# )

# # --- RUN SECOND TEST: Active WENO Extrapolation ---
# days_test2, weno_matrix, weno_midpoints = ForecastingTools.test_prediction_horizon_weno_extrapolation(
#     cost, original_u0, all_models, Model_no, opt_positions_filtered, method, P, d,
#     weno3_nonuniform,
#     180.0, end_time, 0.25
# )


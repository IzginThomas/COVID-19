module WenoInterpolation

export weno_evaluate_non_uniform

using LinearAlgebra
using ForwardDiff

# ==============================================================================
# Conservative stencil polynomial reconstruction
# ==============================================================================

function evaluate_stencil_p(
    x,
    v_bar,
    delta_x,
    x_nodes,
    i::Int,
    r::Int,
    k::Int
)

    T = typeof(x)

    p_val = zero(T)

    for m in 0:k

        # ----------------------------------------------------------------------
        # Conservative coefficient
        # ----------------------------------------------------------------------

        inner_sum = zero(T)

        for j in 0:(m-1)

            idx_j = i - r + j
            inner_sum += v_bar[idx_j] * delta_x[idx_j]

        end

        # ----------------------------------------------------------------------
        # Numerator
        # ----------------------------------------------------------------------

        numerator_sum = zero(T)

        for l in 0:k

            if l == m
                continue
            end

            prod_top = one(T)

            for q in 0:k

                if q == m || q == l
                    continue
                end

                idx_q = i - r + q
                prod_top *= (x - x_nodes[idx_q])

            end

            numerator_sum += prod_top

        end

        # ----------------------------------------------------------------------
        # Denominator
        # ----------------------------------------------------------------------

        denominator_prod = one(T)

        for l in 0:k

            if l == m
                continue
            end

            idx_m = i - r + m
            idx_l = i - r + l

            denominator_prod *= (
                x_nodes[idx_m] - x_nodes[idx_l]
            )

        end

        p_val += inner_sum * (
            numerator_sum / denominator_prod
        )

    end

    return p_val

end

# ==============================================================================
# Nonuniform smoothness indicators via quadrature + AD
# ==============================================================================

function compute_smoothness_indicator(
    v_bar,
    delta_x,
    x_nodes,
    i::Int,
    r::Int,
    k::Int
)

    # --------------------------------------------------------------------------
    # Reconstruction polynomial
    # --------------------------------------------------------------------------

    p(x) = evaluate_stencil_p(
        x,
        v_bar,
        delta_x,
        x_nodes,
        i,
        r,
        k
    )

    # --------------------------------------------------------------------------
    # Automatic differentiation
    # --------------------------------------------------------------------------

    dp(x) = ForwardDiff.derivative(p, x)

    # k=2 only needs first derivative
    if k == 2

        # ----------------------------------------------------------------------
        # Gauss-Legendre quadrature
        # ----------------------------------------------------------------------

        nodes = [
            -0.906179845938664,
            -0.538469310105683,
             0.0,
             0.538469310105683,
             0.906179845938664
        ]

        weights = [
            0.236926885056189,
            0.478628670499366,
            0.568888888888889,
            0.478628670499366,
            0.236926885056189
        ]

        x_l = x_nodes[i]
        x_r = x_nodes[i+1]

        mid = 0.5 * (x_l + x_r)
        half = 0.5 * (x_r - x_l)

        β = 0.0

        for q in eachindex(nodes)

            xq = mid + half * nodes[q]

            β += weights[q] * (
                dp(xq)^2
            )

        end

        return half * β

    end

    # --------------------------------------------------------------------------
    # k=3 uses first and second derivative
    # --------------------------------------------------------------------------

    d2p(x) = ForwardDiff.derivative(dp, x)

    nodes = [
        -0.906179845938664,
        -0.538469310105683,
         0.0,
         0.538469310105683,
         0.906179845938664
    ]

    weights = [
        0.236926885056189,
        0.478628670499366,
        0.568888888888889,
        0.478628670499366,
        0.236926885056189
    ]

    x_l = x_nodes[i]
    x_r = x_nodes[i+1]

    mid = 0.5 * (x_l + x_r)
    half = 0.5 * (x_r - x_l)

    dx_i = x_r - x_l

    β = 0.0

    for q in eachindex(nodes)

        xq = mid + half * nodes[q]

        β += weights[q] * (
            dp(xq)^2 +
            dx_i^2 * d2p(xq)^2
        )

    end

    return half * β

end

# ==============================================================================
# Compute nonuniform optimal weights
# ==============================================================================

function compute_non_uniform_d(
    x,
    v_bar,
    delta_x,
    x_nodes,
    i::Int,
    k::Int,
    p_r
)

    total_k = 2 * k - 1
    global_r = k - 1

    P_global = evaluate_stencil_p(
        x,
        v_bar,
        delta_x,
        x_nodes,
        i,
        global_r,
        total_k
    )

    A = zeros(eltype(x), 2, k)

    for r in 0:(k-1)

        A[1, r+1] = p_r[r+1]
        A[2, r+1] = 1.0

    end

    b = [P_global, 1.0]

    d = pinv(A) * b

    # positivity clipping
    for r in 1:k

        if d[r] < 0.0
            d[r] = 0.0
        end

    end

    s = sum(d)

    if s > 1e-14
        d = d ./ s
    else
        d .= 1.0 / k
    end

    return d

end

# ==============================================================================
# Pure WENO evaluation
# ==============================================================================

function weno_pure_at_cell_i(
    x,
    v_bar,
    delta_x,
    x_nodes,
    i::Int,
    k::Int
)

    eps = 1e-6

    p_r = zeros(eltype(x), k)
    beta = zeros(eltype(x), k)
    alpha = zeros(eltype(x), k)

    # --------------------------------------------------------------------------
    # candidate reconstructions
    # --------------------------------------------------------------------------

    for r in 0:(k-1)

        p_r[r+1] = evaluate_stencil_p(
            x,
            v_bar,
            delta_x,
            x_nodes,
            i,
            r,
            k
        )

        beta[r+1] = compute_smoothness_indicator(
            v_bar,
            delta_x,
            x_nodes,
            i,
            r,
            k
        )

    end

    # --------------------------------------------------------------------------
    # optimal weights
    # --------------------------------------------------------------------------

    d = compute_non_uniform_d(
        x,
        v_bar,
        delta_x,
        x_nodes,
        i,
        k,
        p_r
    )

    # --------------------------------------------------------------------------
    # nonlinear weights
    # --------------------------------------------------------------------------

    for r in 1:k

        alpha[r] = d[r] / (beta[r] + eps)^2

    end

    omega = alpha ./ sum(alpha)

    return sum(omega .* p_r)

end

# ==============================================================================
# Main WENO evaluator
# ==============================================================================

function weno_evaluate_non_uniform(
    x,
    v_bar_phys,
    delta_x_phys,
    x_nodes_phys,
    k::Int;
    apply_zhang_shu::Bool=true,
    lower_bound::Float64=0.0,
    upper_bound::Float64=Inf,
    extrapolation::Symbol=:constant
)
    # Berechne die exakt benötigte Anzahl an Ghost-Cells
    nghost = k - 1

    # --------------------------------------------------------------------------
    # Ghost cells für Zellmittelwerte (v_bar)
    # --------------------------------------------------------------------------
    if extrapolation == :linear && length(v_bar_phys) >= 2
        slope_left = v_bar_phys[2] - v_bar_phys[1]
        v_bar_left = [
            v_bar_phys[1] - (nghost - i + 1) * slope_left
            for i in 1:nghost
        ]

        slope_right = v_bar_phys[end] - v_bar_phys[end-1]
        v_bar_right = [
            v_bar_phys[end] + i * slope_right
            for i in 1:nghost
        ]

        v_bar = [v_bar_left; v_bar_phys; v_bar_right]
    else
        v_bar = [
            fill(v_bar_phys[1], nghost)
            v_bar_phys
            fill(v_bar_phys[end], nghost)
        ]
    end

    # Delta_x analog mit nghost erweitern
    delta_x = [
        fill(delta_x_phys[1], nghost)
        delta_x_phys
        fill(delta_x_phys[end], nghost)
    ]

    # --------------------------------------------------------------------------
    # Knoten (x_nodes) exakt erweitern
    # --------------------------------------------------------------------------
    dx_left = delta_x_phys[1]
    dx_right = delta_x_phys[end]

    x_nodes_left = [
        x_nodes_phys[1] - r * dx_left
        for r in nghost:-1:1
    ]

    x_nodes_right = [
        x_nodes_phys[end] + r * dx_right
        for r in 1:nghost
    ]

    x_nodes = [
        x_nodes_left
        x_nodes_phys
        x_nodes_right
    ]

    # --------------------------------------------------------------------------
    # Indizes für die physikalischen Grenzen anpassen
    # --------------------------------------------------------------------------
    N_total = length(delta_x)

    i_phys_min = nghost + 1
    i_phys_max = N_total - nghost

    i = searchsortedlast(x_nodes, x)

    if i < i_phys_min
        i = i_phys_min
    elseif i >= i_phys_max + 1
        i = i_phys_max
    end

    # --------------------------------------------------------------------------
    # WENO reconstruction
    # --------------------------------------------------------------------------

    v_weno_x = weno_pure_at_cell_i(
        x,
        v_bar,
        delta_x,
        x_nodes,
        i,
        k
    )

    # --------------------------------------------------------------------------
    # Zhang-Shu limiter
    # --------------------------------------------------------------------------

    if apply_zhang_shu

        x_left = x_nodes[i]
        x_right = x_nodes[i+1]

        M_i = v_weno_x
        U_i = v_weno_x

        n_points = 100

        for scale in range(0.0, 1.0, length=n_points)

            x_test = x_left + scale * (x_right - x_left)

            v_test = weno_pure_at_cell_i(
                x_test,
                v_bar,
                delta_x,
                x_nodes,
                i,
                k
            )

            M_i = min(M_i, v_test)
            U_i = max(U_i, v_test)

        end

        theta_lower = 1.0

        if M_i < lower_bound

            v_avg = v_bar[i]

            denom = v_avg - M_i

            theta_lower = (
                denom > 1e-14
            ) ? min(
                1.0,
                (v_avg - lower_bound) / denom
            ) : 1.0

        end

        theta_upper = 1.0

        if U_i > upper_bound

            v_avg = v_bar[i]

            denom = U_i - v_avg

            theta_upper = (
                denom > 1e-14
            ) ? min(
                1.0,
                (upper_bound - v_avg) / denom
            ) : 1.0

        end

        theta = min(
            1.0,
            max(0.0, theta_lower),
            max(0.0, theta_upper)
        )

        if theta < 1.0

            v_avg = v_bar[i]

            return theta * (
                v_weno_x - v_avg
            ) + v_avg

        end

    end

    return v_weno_x

end

end






#  module WenoInterpolation

#  export weno_evaluate_non_uniform

# # """
# #     weno3_interpolate(t::Float64, X::Vector{Float64}, Y::Vector{Float64}; eps_pos=1e-10)

# # Continuous Piecewise Linear Conservative Reconstruction for Non-Uniform Grids.
# # Guarantees global continuity across cell interfaces while strictly conserving 
# # the cell averages Y. Drop-in replacement for your ODE parameter tracking.
# # """
# # export weno_evaluate_non_uniform
# using LinearAlgebra

# """
# Berechnet das fundamentale Polynom p_r(x) für einen spezifischen Stencil r
# gemäß Gleichung (2.19) an einem beliebigen Punkt x.
# """
# function evaluate_stencil_p(x::Float64, v_bar::Vector{Float64}, delta_x::Vector{Float64}, x_nodes::Vector{Float64}, i::Int, r::Int, k::Int)
#     p_val = 0.0
#     for m in 0:k
#         inner_sum = 0.0
#         for j in 0:(m-1)
#             idx_j = i - r + j
#             inner_sum += v_bar[idx_j] * delta_x[idx_j]
#         end
        
#         numerator_sum = 0.0
#         for l in 0:k
#             if l == m continue end
#             prod_top = 1.0
#             @inbounds for q in 0:k
#                 if q == m || q == l continue end
#                 idx_q = i - r + q
#                 prod_top *= (x - x_nodes[idx_q])
#             end
#             numerator_sum += prod_top
#         end
        
#         denominator_prod = 1.0
#         for l in 0:k
#             if l == m continue end
#             idx_m = i - r + m
#             idx_l = i - r + l
#             denominator_prod *= (x_nodes[idx_m] - x_nodes[idx_l])
#         end
        
#         p_val += inner_sum * (numerator_sum / denominator_prod)
#     end
#     return p_val
# end

# """
# Berechnet den Glattheitsindikator β_r basierend auf den Zellmittelwerten.
# """
# function compute_smoothness_indicator(v_bar::Vector{Float64}, i::Int, r::Int, k::Int)
#     if k == 1
#         return 1.0
#     elseif k == 2 
#         v1 = v_bar[i - r]
#         v2 = v_bar[i - r + 1]
#         return (v2 - v1)^2
#     elseif k == 3
#         v1 = v_bar[i - r]
#         v2 = v_bar[i - r + 1]
#         v3 = v_bar[i - r + 2]
        
#         term1 = (13/12) * (v1 - 2*v2 + v3)^2
#         if r == 0
#             term2 = (1/4) * (3*v1 - 4*v2 + v3)^2
#         elseif r == 1
#             term2 = (1/4) * (v1 - v3)^2
#         else
#             term2 = (1/4) * (v1 - 4*v2 + 3*v3)^2
#         end
#         return term1 + term2
#     else
#         error("Ordnung k=$k ist nicht implementiert.")
#     end
# end

# """
# Bestimmt die gitter- und ortsabhängigen idealen Gewichte d_r(x) für ein 
# allgemeines, NICHT-ÄQUIDISTANTES Gitter.
# """
# function compute_non_uniform_d(x::Float64, v_bar::Vector{Float64}, delta_x::Vector{Float64}, x_nodes::Vector{Float64}, i::Int, k::Int, p_r::Vector{Float64})
#     total_k = 2 * k - 1
#     global_r = k - 1
    
#     P_global = evaluate_stencil_p(x, v_bar, delta_x, x_nodes, i, global_r, total_k)
    
#     A = zeros(2, k)
#     for r in 0:(k-1)
#         A[1, r+1] = p_r[r+1]
#         A[2, r+1] = 1.0
#     end
#     b = [P_global, 1.0]
    
#     d = pinv(A) * b
#     for r in 1:k
#         if d[r] < 0.0
#             d[r] = 0.0
#         end
#     end
#     s = sum(d)
#     if s > 1e-14
#         d = d ./ s
#     else
#         d = fill(1.0 / k, k)
#     end
#     return d
# end

# """
# Hilfsfunktion: Berechnet den reinen unlimitierten WENO-Wert an Stelle x für Zelle i.
# (Wird für die Minimumssuche des Zhang-Shu-Limiters benötigt).
# """
# function weno_pure_at_cell_i(x::Float64, v_bar::Vector{Float64}, delta_x::Vector{Float64}, x_nodes::Vector{Float64}, i::Int, k::Int)
#     eps = 1e-6
#     p_r = zeros(k)
#     beta = zeros(k)
#     alpha = zeros(k)
    
#     for r in 0:(k-1)
#         p_r[r+1] = evaluate_stencil_p(x, v_bar, delta_x, x_nodes, i, r, k)
#         beta[r+1] = compute_smoothness_indicator(v_bar, i, r, k)
#     end
    
#     d = compute_non_uniform_d(x, v_bar, delta_x, x_nodes, i, k, p_r)
    
#     for r in 1:k
#         alpha[r] = d[r] / (beta[r] + eps)^2
#     end
#     omega = alpha ./ sum(alpha)
#     return sum(omega .* p_r)
# end

# """
# HAUPTFUNKTION: Evaluiert die WENO-Approximation an einem beliebigen Punkt x
# für ein vollständig NICHT-ÄQUIDISTANTES Gitter. 
# Inklusive unbegrenzter Extrapolation und optionalem ZHANG-SHU LIMITER (z.B. Positivitätsschutz).
# """
# function weno_evaluate_non_uniform(
#     x::Float64, 
#     v_bar_phys::Vector{Float64}, 
#     delta_x_phys::Vector{Float64}, 
#     x_nodes_phys::Vector{Float64}, 
#     k::Int; 
#     apply_zhang_shu::Bool=true, 
#     lower_bound::Float64=0.0,
#     upper_bound::Float64=Inf,
#     extrapolation::Symbol=:constant
# )
#     N_phys = length(delta_x_phys)
    
#     # 1. Gitter- und Daten-Vektoren intern um k Ghost-Cells erweitern
#     # Wir nutzen Outflow/Neumann-Randbedingungen (konstant) oder lineare Extrapolation
#     if extrapolation == :linear && length(v_bar_phys) >= 2
#         slope_left = v_bar_phys[2] - v_bar_phys[1]
#         v_bar_left = [v_bar_phys[1] - (k - i + 1) * slope_left for i in 1:k]
        
#         slope_right = v_bar_phys[end] - v_bar_phys[end-1]
#         v_bar_right = [v_bar_phys[end] + i * slope_right for i in 1:k]
        
#         v_bar = [v_bar_left; v_bar_phys; v_bar_right]
#     else
#         v_bar = [fill(v_bar_phys[1], k); v_bar_phys; fill(v_bar_phys[end], k)]
#     end
    
#     delta_x = [fill(delta_x_phys[1], k); delta_x_phys; fill(delta_x_phys[end], k)]
    
#     # Für die Knoten (x_nodes) müssen wir die Koordinaten mathematisch 
#     # nach links und rechts mit der jeweiligen Rand-Schrittweite dx fortsetzen.
#     dx_left = delta_x_phys[1]
#     dx_right = delta_x_phys[end]
    
#     x_nodes_left = [x_nodes_phys[1] - r * dx_left for r in k:-1:1]
#     x_nodes_right = [x_nodes_phys[end] + r * dx_right for r in 1:k]
#     x_nodes = [x_nodes_left; x_nodes_phys; x_nodes_right]

#     # Total-Größe des erweiterten Gitters
#     N_total = length(delta_x)
    
#     # Die echten physikalischen Zellen rutschen durch das Anfügen nach innen
#     i_phys_min = k + 1
#     i_phys_max = N_total - k

#     # 2. Zelle 'i' mittels Binärsuche im erweiterten Gitter bestimmen
#     i = searchsortedlast(x_nodes, x)

#     # 3. Abfang-Logik (Clipping) für echte Extrapolationspunkte
#     # Wenn x außerhalb des echten Gitters liegt, zwingen wir die Auswertung 
#     # auf die jeweilige Randzelle, die nun intern dank Ghost-Cells voll funktionsfähig ist.
#     if i < i_phys_min
#         i = i_phys_min
#     elseif i >= i_phys_max + 1
#         i = i_phys_max
#     end

#     # 4. Berechne den unlimitierten WENO-Wert an Stelle x für Zelle i
#     v_weno_x = weno_pure_at_cell_i(x, v_bar, delta_x, x_nodes, i, k)

#     # 5. Zhang-Shu Limiter Logik (Doppelseitig für Unter- und Überschwinger)
#     if apply_zhang_shu
#         x_left = x_nodes[i]
#         x_right = x_nodes[i+1]
        
#         M_i = v_weno_x 
#         U_i = v_weno_x
#         n_points = 100
        
#         for scale in range(0.0, 1.0, length=n_points)
#             x_test = x_left + scale * (x_right - x_left)
#             v_test = weno_pure_at_cell_i(x_test, v_bar, delta_x, x_nodes, i, k)
#             if v_test < M_i
#                 M_i = v_test
#             end
#             if v_test > U_i
#                 U_i = v_test
#             end
#         end
        
#         theta_lower = 1.0
#         if M_i < lower_bound
#             v_avg = v_bar[i]
#             denom = v_avg - M_i
#             theta_lower = (denom > 1e-14) ? min(1.0, (v_avg - lower_bound) / denom) : 1.0
#         end
        
#         theta_upper = 1.0
#         if U_i > upper_bound
#             v_avg = v_bar[i]
#             denom = U_i - v_avg
#             theta_upper = (denom > 1e-14) ? min(1.0, (upper_bound - v_avg) / denom) : 1.0
#         end
        
#         theta = min(1.0, max(0.0, theta_lower), max(0.0, theta_upper))
#         if theta < 1.0
#             v_avg = v_bar[i]
#             return theta * (v_weno_x - v_avg) + v_avg
#         end
#     end
    
#     return v_weno_x
# end



# end # Ende des Moduls WenoInterpolation

module opt_vars

export get_opt_indices, matrix_map, parameters

matrix_map = Dict(
    # Compartment 1 (S):
    # In: r1p (pv, Λ), aVS, aES, aRS
    # Out: aSV, aSE (c-constants, k1), μ, Pmat[3,1], Pmat[9,1]
    "aSE" => [:cS, :cV, :cE, :cL, :cI, :cH, :cR, :cQ, :cVe, :k1], # This is the only parameter directly related to aSE, but k1 is also relevant for the transition from S to 
    
    "S"  => [:pv, :Λ, :aVS, :aES, :aRS, :aSV, :k1, :cS, :cV, :cE, :cL, :cI, :cH, :cR, :cQ, :cVe, :μ],# :b01,:b02,:b03,:b04,:b11,:b12,:b13, :b14,:ω1, :ω2],

    # Compartment 2 (V):
    # In: r2p (pv, Λ), aSV
    # Out: aVS, aVE (aSE^k2), μ, δV, Pmat[1,2], Pmat[3,2], Pmat[9,2]
    "V"  => [:pv, :Λ, :aSV, :aVS, :aVE, :k1, :k2, :cS, :cV, :cE, :cL, :cI, :cH, :cR, :cQ, :cVe, :μ, :δV],#, :b01,:b02,:b03,:b04,:b11,:b12,:b13, :b14,:ω1, :ω2],

    # Compartment 3 (E):
    # In: aSE, aVE
    # Out: aES, aEL, aEI, aEQ, μ, δE, Pmat[1,3], Pmat[4,3], Pmat[5,3], Pmat[8,3], Pmat[9,3]
    "E"  => [:k1, :k2, :cS, :cV, :cE, :cL, :cI, :cH, :cR, :cQ, :cVe, :aES, :aEL, :aEI, :aEQ, :φ, :μ, :δE],#, :b01,:b02,:b03,:b04,:b11,:b12,:b13, :b14,:ω1, :ω2],

    # Compartment 4 (L):
    # In: aEL
    # Out: aLI, aLR, aLQ, μ, δL, φ, Pmat[5,4], Pmat[7,4], Pmat[8,4], Pmat[9,4], Pmat[10,4]
    "L"  => [:aEL, :aLI, :aLR, :aLQ, :γ, :μ, :δL, :φ],

    # Compartment 5 (I):
    # In: aEI, aLI
    # Out: aIH, aIR, aIQ, φ, μ, αI, Pmat[6,5], Pmat[7,5], Pmat[8,5], Pmat[9,5], Pmat[10,5]
    "I"  => [:aEI, :aLI, :aIH, :aIR, :aIQ, :φ, :ψ, :γ, :μ, :αI],#:b01,:b02,:b03,:b04,:b11,:b12,:b13, :b14,:ω1, :ω2], # Added parameters for time-varying transmission

    # Compartment 6 (H):
    # In: aIH, ψ
    # Out: aHR, μ, δH, αH, Pmat[7,6], Pmat[9,6]
    "H"  => [:aIH, :aHR, :μ, :δH, :αH, :ψ],

    # Compartment 7 (R):
    # In: aLR, aIR, aHR, aQR
    # Out: aRS, μ, Pmat[1,7], Pmat[9,7]
    "R"  => [:aLR, :aIR, :aHR, :aQR, :aRS, :μ, :γ, :ψ],

    # Compartment 8 (Q):
    # In: aEQ, aLQ, aIQ
    # Out: aQR, μ, δQ, Pmat[7,8], Pmat[9,8]
    "Q"  => [:aEQ, :aLQ, :aIQ, :aQR, :μ, :δQ],

    # Compartment 9 (Mortality Sink):
    # This row collects all outflows. It is dependent on the death coefficients of all states.
    "D"  => [:μ, :δV, :δE, :δL, :δH, :δQ, :αI, :αH],

    # Compartment 10 (Recovery/Vaccine Pool):
    # In: r10p (rLVe, rIVe)
    # Out: μVe (Assuming μVe is the death/decay rate of this pool)
    "Ve" => [:rLVe, :rIVe, :μVe]
)

"""
Returns the positions in the 'parameters' array for a given set of rows.
Example: get_indices([1, 3]) returns params for S and E compartments.
"""
parameters = ["k1", "k2", "cS", "cV", "cE", "cL", "cI", "cH", "cR", "cQ",
              "cVe", "δV", "δE", "δL", "δH", "δQ", "pv", "Λ", "μ", "φ",
              "ψ", "γ", "μVe", "αI", "αH", "aSV", "aVS", "aVE", "aEI", "aES",
              "aEQ", "aEL", "aLI", "aLR", "aLQ", "rLVe", "aIR", "aIQ", "aIH", "rIVe",
              "aHR", "aRS", "aQR", "b01","b02","b03","b04","b11","b12","b13", "b14","ω1", "ω2"]

function get_opt_indices(targets, param_list, m_map)
    # Ensure targets is a collection (even if a single string is passed)
    target_list = targets isa String ? [targets] : targets
    
    selected_syms = Set{Symbol}()
    for t in target_list
        if haskey(m_map, t)
            union!(selected_syms, m_map[t])
        else
            @warn "Compartment $t not found in map."
        end
    end
    
    # Return sorted indices of where these parameters live in your list
    return sort([findfirst(==(string(s)), param_list) for s in selected_syms])
end
end
#=
# Usage:
idx = get_opt_indices("S", parameters, matrix_map)
# Or multiple:
idx_multi = get_opt_indices(["S", "E", "D"], parameters, matrix_map)


println("Indices for Row 1: ", idx)
println("Indices for Rows 1, 3, and 4: ", idx_multi)
#
=#
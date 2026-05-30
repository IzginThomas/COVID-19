#include("get_data_new.jl");
#include("opt_vars.jl");

Model_no = 1; # choose model number (1-5, Inf)


opt_compartments = ["S","V", "E", "L", "I", "H", "R", "Q","D", "Ve"];
#opt_compartments = ["aSE", "I", "D"];
opt_vars_pos = opt_vars.get_opt_indices(opt_compartments, opt_vars.parameters, opt_vars.matrix_map);

if Model_no < Inf
    # filter out zero and one parameters (fixed parameters - not to be optimized)
    nz_positions = findall(!iszero, all_models[Model_no]); # find positions of non-zero parameters for the chosen model)
    opt_positions = intersect(nz_positions, opt_vars_pos);
    special_symbols = ["φ", "ψ", "γ"] # (if =0 (filtered out already)or =1: fixed parameters - not to be optimized)
    special_idx = findall(p -> p in special_symbols, opt_vars.parameters)
    current_vals = all_models[Model_no]
    excl_special_idx = filter(i -> abs(current_vals[i]) == 1.0, special_idx) # filters k1 if it is -1.0, too
    delta_idx = findall(p -> occursin("δ", p), opt_vars.parameters) # indices of δV, δE, δL, δH, δQ
    opt_delta = intersect(opt_positions, delta_idx) # find non-zero δ parameters
    mu_idx = findall(p -> p == "μ", opt_vars.parameters) # indices of μ
    Lambda_idx = findall(p -> p == "Λ", opt_vars.parameters) # indices of Λ
    exclude_idx = union(delta_idx, excl_special_idx, mu_idx, Lambda_idx) # indices of parameters to exclude from optimization
    opt_positions_filtered = setdiff(opt_positions, exclude_idx)
    nopt_positions = setdiff(1:length(all_models[Model_no]), opt_positions_filtered) # find positions of zero parameters for the chosen model
    x_model = all_models[Model_no][opt_positions_filtered] # extract the non-zero, non-delta parameters for the chosen model to optimize
    if Model_no == 1
        flag_N = false# flag to indicate if total population data is available
        flag_D = true# flag to indicate if deaths data is available
        flag_V = false# flag to indicate if vaccinated data is available
        flag_cum_V = false# flag to indicate if cumulative vaccinated data is available
        flag_I = true# flag to indicate if infected data is available
        flag_cum_I = false# flag to indicate if cumulative infected data is available
        flag_H = false# flag to indicate if hospitalized data is available
        flag_R = false# flag to indicate if recovered data is available
        flag_Q = false# flag to indicate if total population data is available
        flag_Ve = false# 
        flag_E = false# flag to indicate if exposed data is available
        tspan = (0.0, maximum(tref_I)); # replace with real data time span
    else # TODO OTHER MODELS
        flag_N = false# flag to indicate if total population data is available
        flag_D = true# flag to indicate if deaths data is available
        flag_V = true# flag to indicate if vaccinated data is available
        flag_cum_V = true# flag to indicate if cumulative vaccinated data is available
        flag_I = true# flag to indicate if infected data is available
        flag_cum_I = false# flag to indicate if cumulative infected data is available
        flag_H = false# flag to indicate if hospitalized data is available
        flag_R = false# flag to indicate if recovered data is available
        flag_Q = false# flag to indicate if total population data is available
        flag_Ve = false# 
        flag_E = false# flag to indicate if exposed data is available
        tspan = (0.0, maximum([tref_I; tref_D;tref_V])); # replace with real data time span
    end
else
    opt_positions_filtered = 1:length(all_models[1]) # all parameters are non-zero for the comprehensive model
    nopt_positions = Int[] # no zero parameters
    #delta_idx = findall(p -> occursin("δ", p), parameters) # indices of δV, δE, δL, δH, δQ
    #opt_delta = findall(p -> occursin("δ", p), all_models[1]) # indices of δV, δE, δL, δH, δQ
    #opt_positions_delta = setdiff(opt_positions, delta_idx) # find non-zero parameters that are not δ
    x_model = all_models[1][opt_positions_filtered] # extract the non-zero parameters for the comprehensive model
end
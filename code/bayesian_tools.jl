

module bayesian_tools

export run_opt

using BayesianOptimization, GaussianProcesses, Distributions




function set_model(D)
    model = ElasticGPE(D,
                mean = MeanConst(1.0),
                kernel = Mat52Ard(fill(-1.0, D), 0.0),  # etwas neutraler Start
                logNoise = -10.0,                       # stabiler als -12
                capacity = 3000)                        # 🔥 reduziert (vorher 5000)
    
    set_priors!(model.mean, [Normal(0, 1)]) 

    # 🔥 Bounds angepasst für skalierten Input [0,1]
    lower_kern = vcat(fill(-2.5, D), -1.0)
    upper_kern = vcat(fill(-0.5, D), 3.0)

    modeloptimizer = MAPGPOptimizer(
        every = 50,                          # seltener optimieren
        noisebounds = [-12.0, -8.0],         # stabiler
        kernbounds  = [lower_kern, upper_kern],
        maxeval = 30                         # etwas leichter
    )

    return model, modeloptimizer
end


function run_opt(bounds, maxiteration, f)

    lower = bounds[1, :]
    upper = bounds[2, :]
    D = length(lower)
    ranges = upper .- lower

    function scaled_f(x_unit)
        x_real = lower .+ (x_unit .* ranges)
        return f(x_real)
    end

    model, modeloptimizer = set_model(D)
   # x_unit = (x_model .- lower) ./ ranges
    #y0 = f(x_model) 
    #append!(model, reshape(x_unit, :, 1), [y0])

    result = nothing
    stop = false
    cnt = 0
    oldmin = Inf

    while !stop
        if cnt == 0
            opt_Exp_Imp = BOpt(
                scaled_f,                        # 🔥 hier!
                model,
                ExpectedImprovement(0.0),
                modeloptimizer,
                zeros(D), ones(D),               # 🔥 jetzt [0,1] Raum
                repetitions = 1,
                maxiterations = maxiteration,
                sense = Min,
                acquisitionoptions = (
                    method = :LD_LBFGS,
                    restarts = 10,
                    maxtime = 0.1,
                    maxeval = 500
                ),
                verbosity = Silent
            )
            boptimize!(opt_Exp_Imp)
        else
            opt_Exp_Imp = BOpt(
                scaled_f,
                model,
                ExpectedImprovement(0.0),
                modeloptimizer,
                zeros(D), ones(D),
                repetitions = 1,
                maxiterations = maxiteration,
                sense = Min,
                acquisitionoptions = (
                    method = :LD_LBFGS,
                    restarts = 10,
                    maxtime = 0.1,
                    maxeval = 2000
                ),
                verbosity = Silent,
                initializer_iterations = 2*D + 2
            )
            boptimize!(opt_Exp_Imp)
        end

        opt_Prob_Imp = BOpt(
            scaled_f,
            model,
            ProbabilityOfImprovement(),
            modeloptimizer,
            zeros(D), ones(D),
            repetitions = 1,
            maxiterations = maxiteration,
            sense = Min,
            acquisitionoptions = (
                method = :LD_LBFGS,
                restarts = 1,
                maxtime = 0.1,
                maxeval = 2000
            ),
            verbosity = Silent,
            initializer_iterations = 0
        )

        result = boptimize!(opt_Prob_Imp)

        stop = result.observed_optimum + 1e-3 >= oldmin
        oldmin = result.observed_optimum
        cnt += 1
    end

    result = (
        observed_optimum = result.observed_optimum,
        observed_optimizer = lower .+ (result.observed_optimizer .* ranges),
        model_optimum = result.model_optimum,
        model_optimizer = lower .+ (result.model_optimizer .* ranges)
    )

    return result, cnt
end

end






#=

module bayesian_tools

export run_opt

using BayesianOptimization, GaussianProcesses, Distributions

#=
function set_model(D)
    model = ElasticGPE(D,                            # dimension = D input dimensions
                mean = MeanConst(0.),         
               kernel = SEArd(zeros(D), 5.), # size: input dimensions + 1
               logNoise = -10.,               # logNoise off
               capacity = 3000)              # the initial capacity of the GP is 3000 samples.
    set_priors!(model.mean, [Normal(1, 2)]) 
        #kernelbounds
    lower_kern = vcat(fill(-1.0, D), 0.0) #size input dimension + 1
    upper_kern = vcat(fill( 4.0, D), 10.0) 

    # Optimize the hyperparameters of the GP using maximum a posteriori (MAP) estimates every 50 steps
    modeloptimizer = MAPGPOptimizer(every = 10000, noisebounds = [-11, -10],       # bounds of the logNoise
                                kernbounds  = [lower_kern, upper_kern],  # bounds of the parameters GaussianProcesses.get_param_names(model.kernel)
                                maxeval = 40)
    return model, modeloptimizer
end
=#
function set_model(D)
    model = ElasticGPE(D,
                mean = MeanConst(0.),
                # Using 0.5 as a starting lengthscale (log(0.5) ≈ -0.7)
                kernel = Mat52Ard(fill(-0.7, D), 0.0), 
                logNoise = -12.0,                
                capacity = 5000)
    
    set_priors!(model.mean, [Normal(0, 1)]) 

    # Tighten bounds: 
    # Lower: e^-10 ≈ 4,5e-10 (allows some complexity)
    # Upper: e^-0.5 ≈ 0.9512 (prevents over-smoothing the 0-1 range)
    lower_kern = vcat(fill(-5.0, D), -1.0) # vorher -3 
    upper_kern = vcat(fill(-0.5, D), 3.0) # vorher -0.2

    modeloptimizer = MAPGPOptimizer(every = 25, 
                                    noisebounds = [-14.0, -10.0], 
                                    kernbounds  = [lower_kern, upper_kern],
                                    maxeval = 40)
    return model, modeloptimizer
end
# function run_opt
# to run the bayesian optimization, alternating the two acquisitionfunction ExpectedImprovement and ProbabilityOfImprovement
# input: bounds: 2 x D matrix containing the lower bound in first row and upper in second for every parameters
#        maxiteration: iterations of every call of BOpt
#        f(x):         costfunction with parameters stored in one vector 'x'
# output: result: optimized parameters, stored in: 
#                 @NamedTuple{observed_optimum::Float64, observed_optimizer::Vector{Float64}, model_optimum::Float64, model_optimizer::Vector{Float64}}
#                 cnt: counter of calls (ExpectedImprovement + ProbabilityOfImprovement)

function run_opt(bounds, maxiteration, f) #bounds 2 x D vector with lower bounds in first row, upper in second
lower = bounds[1, :]
    upper = bounds[2, :]
    D = length(lower)
    ranges = upper .- lower
#=
    # Internal wrapper: Always presents a [0, 1] space to the Optimizer
    function scaled_f(x_unit) # TODO use scaled_f for optimization. Dont forget to scale back the result
        x_real = lower .+ (x_unit .* ranges)
        return f(x_real)
    end
    =#
model, modeloptimizer = set_model(D)

result = nothing

stop = false
cnt = 0
oldmin = Inf

    while ~stop
        if cnt == 0
            opt_Exp_Imp = BOpt(f,
            model,
            ExpectedImprovement(0.01),                   # type of acquisition
            modeloptimizer,                        
            lower, upper,                     # lowerbounds, upperbounds         
            repetitions = 1,                          # evaluate the function for each input 5 times
            maxiterations = maxiteration,                      # evaluate at 100 input positions
            sense = Min,                              # minimize the function
            acquisitionoptions = (method = :LD_LBFGS, # run optimization of acquisition function with NLopts :LD_LBFGS method
                                     restarts = 15,       # run the NLopt method from 5 random initial conditions each time.
                                     maxtime = 1.0,      # run the NLopt method for at most 0.1 second each time
                                     maxeval = 1000),    # run the NLopt methods for at most 1000 iterations (for other options see https://github.com/JuliaOpt/NLopt.jl)
                                    verbosity = Progress)
        #result = boptimize!(opt_Exp_Imp)
        boptimize!(opt_Exp_Imp)
        else
            opt_Exp_Imp = BOpt(f, model, ExpectedImprovement(0.01), modeloptimizer,                        
                lower, upper,                     # lowerbounds, upperbounds         
                repetitions = 1,                          # evaluate the function for each input 5 times
                maxiterations = maxiteration,                      # evaluate at 100 input positions
                sense = Min,                              # minimize the function
                acquisitionoptions = (method = :LD_LBFGS,
                                        restarts = 15,
                                        maxtime = 1.0,
                                        maxeval = 1000),
                initializer_iterations = 2*D + 2) #important for any further iterations

            #result = boptimize!(opt_Exp_Imp)
            boptimize!(opt_Exp_Imp)
        end

    opt_Prob_Imp = BOpt(f, model, ProbabilityOfImprovement(),modeloptimizer,                        
            lower, upper,                     # lowerbounds, upperbounds         
            repetitions = 1,                          # evaluate the function for each input 5 times
            maxiterations = maxiteration,                      # evaluate at 100 input positions
            sense = Min,                              # minimize the function
            acquisitionoptions = (method = :LD_LBFGS,
                                    restarts = 1,
                                    maxtime = 0.1,
                                    maxeval = 10000),
            initializer_iterations = 0) #important for any further iterations

    result = boptimize!(opt_Prob_Imp)
    #boptimize!(opt_Prob_Imp)
    stop = opt_Prob_Imp.observed_optimum + 1e-3 >= oldmin
    oldmin = opt_Prob_Imp.observed_optimum
    cnt = cnt + 1
    end
    # Reconstruct the result with original scales
    #result.observed_optimizer = lower .+ (result.observed_optimizer .* ranges);
    
    return result, cnt

end


end

=#
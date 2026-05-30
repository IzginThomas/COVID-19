

module bayesian_tools

export run_opt

using BayesianOptimization, GaussianProcesses, Distributions




function set_model(D)
    model = ElasticGPE(D,
                mean = MeanConst(1.0),
                kernel = Mat52Ard(fill(-1.0, D), 0.0),  
                logNoise = -10.0,                       
                capacity = 3000)                        
    
    set_priors!(model.mean, [Normal(0, 1)]) 

    
    lower_kern = vcat(fill(-2.5, D), -1.0)
    upper_kern = vcat(fill(-0.5, D), 3.0)

    modeloptimizer = MAPGPOptimizer(
        every = 50,                          
        noisebounds = [-12.0, -8.0],         
        kernbounds  = [lower_kern, upper_kern],
        maxeval = 30                     
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


    result = nothing
    stop = false
    cnt = 0
    oldmin = Inf

    while !stop
        if cnt == 0
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



clc, clear, close all
warning('off','all')

maxiteration = 1000;

beta1 = optimizableVariable('beta1',[0.1 1]);
%absbeta2 = optimizableVariable('absbeta2',[0.05,0.4], Transform='log'); 
%beta1 = optimizableVariable('beta1',[0.1,1]);%, Transform='log');
beta2 = optimizableVariable('beta2',[-0.4 0.2]);%, Transform='log');
%beta3 = optimizableVariable('beta3',[0,0.1]);%, Transform='log');
%alpha2 = optimizableVariable('alpha2',[1/6, 0.5]);
%kappa2 = optimizableVariable('kappa2',[1,2],'Type','integer');
%vars = [beta1, beta2, beta3, alpha2, kappa2];
vars = [beta1, beta2];

fun = @fun2; % needs to be function handle
 %%% GRID SEARCCH YIELDS BETA1=0.1 AND BETA2=-0.09
results = bayesopt(fun, vars,'IsObjectiveDeterministic',true,...
    'MaxObjectiveEvaluations', maxiteration,...
    'AcquisitionFunctionName','probability-of-improvement',...
    'UseParallel',true);%,'MinWorkerUtilization',60);

delete(parcluster('Processes').Jobs)

save("Bayes_results.mat")

function cost = fun2(vars)
c = 50; %bound for tolerance in cost function
Alpha2 = 1/6;
Beta3 = 0;
Kappa2 = 1;
cost= 0;
solverstr = 'MPRK43IIadap';


%%%  Parameter of the method

switch solverstr
    case 'MPRK43Iadap'
        solvpara1 = 0.515;
        solvpara2 = 0.749;
    case 'MPRK43IIadap'
        solvpara1 = 0.563;
    case 'MPRK22adap'
        solvpara1 = 1;
end
    




%%% tolerances
%tolint = [1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1];
%tolint = [1e-4,1e-3,1e-2,1e-1];
%tolint = [1e-8];
%%% parameter of test problems PR3 and PR4
paraint = 0:0.2:0.5;
%paraint = [0.2];

%%% string of test problems
%tests = {{'PR3',paraint},{'PR4',paraint},'robertson','strat'};
tests ={'robertson','strat',{'PR4',paraint}};
%tests ={{'PR4',paraint}};
tests ={'robertson'};


for testpairs = tests
    testpair = testpairs{:};
    if length(testpair) == 2
        flag = true;
        Test = testpair{1};
        testpara =testpair{2};
        lengthpara = length(testpara);
    else
        flag = false;
        Test = testpair;
        lengthpara = 1;
    end

    if strcmp(Test, 'strat')
        tolint = [1e-6,1e-5,1e-4,1e-3,1e-2,1e-1];
    else
        tolint = [1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1];
    end
    
    lengthA = lengthpara*length(tolint);
    A = zeros(lengthA,2);
    k = 1;
    for tol = tolint
        if flag
            for para=testpara
                A(k,:) = [tol, para];
                k = k+1;
            end
        else
            A(k,:) = [tol, NaN];
            k = k +1;
        end
    end

    for i = 1:lengthA
        tol = A(i,1);
        para = A(i,2);
        %[succstep, cnt_rej, err] = ODESolverTestSuite(vars.beta1, vars.beta2, vars.beta3, vars.alpha2, tol, Test,para,solverstr,solvpara1,solvpara2);
        if exist('solvpara2','var')
            [succstep, cnt_rej, err] = ODESolverTestSuite(vars.beta1, vars.beta2, Beta3, Alpha2, tol, Test,para, Kappa2, solverstr,solvpara1,solvpara2);
        else
            [succstep, cnt_rej, err] = ODESolverTestSuite(vars.beta1, vars.beta2, Beta3, Alpha2, tol, Test,para, Kappa2, solverstr,solvpara1);
        end
        if isnan(cnt_rej)
            cnt_rej = 1e+6;
        elseif isnan(succstep)
            succstep = 1e+6;
        end

%cost = cost + succstep + cnt_rej + ( 1 - 1 / (1 + cnt_rej) ) + max(0, sign( log( err(1) / (c * tol) ) ) * 1e+6);
 cost = cost + succstep + cnt_rej + 1000*cnt_rej/(1+succstep)+ sign(max(0,1000 - succstep/(1+cnt_rej)))*1e+6 + max(0, sign( log( err(1) / (c * tol) ) ) * 1e+6);
    end
end

end
% Optimization ROBERTSON PI & ALPHA2=1/6 completed.
% MaxObjectiveEvaluations of 1000 reached.
% Total function evaluations: 1000
% Total elapsed time: 1926.9225 seconds
% Total objective function evaluation time: 4026.5182
% 
% Best observed feasible point:
%      beta1      beta2  
%     _______    ________
% 
%     0.19258    -0.18674
% 
% Observed objective function value = 25264
% Estimated objective function value = 48102.5938
% Function evaluation time = 3.5467
% 
% Best estimated feasible point (according to models):
%      beta1      beta2  
%     _______    ________
% 
%     0.19187    -0.18849
% 
% Estimated objective function value = 9967.2754
% Estimated function evaluation time = 3.3004
% 
% 


clc, clear, close all
%warning('off','all')
solverset = {'ROS2','MPRK43IIadap','MPRK43Iadap', 'MPRK22adap'};
maxiteration = 500;

beta1 = optimizableVariable('beta1',[-5 5]);
beta2 = optimizableVariable('beta2',[-3 3]);%, Transform='log');
beta3 = optimizableVariable('beta3',[-2 2]);%, Transform='log');
alpha2 = optimizableVariable('alpha2',[-3, 3]);
kappa2 = optimizableVariable('kappa2',[1,4],'Type','integer');
vars = [beta1, beta2, beta3, alpha2, kappa2];
%vars = [beta1, beta2];
for solver = solverset
    % solverstr = 'MPRK43Iadap';
    solverstr = solver{:};
    fun = @(vars) costfunction(vars,solverstr); % needs to be function handle

    stop = false;

    oldmin = inf;
    cnt = 0;
    while ~stop
        if cnt == 0
            results = bayesopt(fun, vars,'IsObjectiveDeterministic',true,...
                'MaxObjectiveEvaluations', maxiteration,... %'AcquisitionFunctionName','probability-of-improvement',...
                'UseParallel',true,'MinWorkerUtilization',8);
        else
            results = resume(results, 'IsObjectiveDeterministic',true,'MaxObjectiveEvaluations',maxiteration);
        end

        % delete(parcluster('Processes').Jobs)
        results = resume(results, 'IsObjectiveDeterministic',true,'MaxObjectiveEvaluations',maxiteration, 'AcquisitionFunctionName','probability-of-improvement');

        stop = results.MinObjective + 1e-3 >= oldmin; % no significant improvement after 1000 iterations
        oldmin = results.MinObjective;
        cnt = cnt + 1;
    end

    save("Bayes_results" + string(solverstr) + "_" + string(num2str(cnt)) + "000_nov2025.mat")
end

function cost = costfunction(vars, solverstr)
% solverstr = 'MPRK43IIadap';
c = 1;
%Alpha2 = 0;
%Beta3 = 0;
%kappa2 = 1;
cost= 0;



%%%  Parameter of the method

switch solverstr
    case 'MPRK43Iadap'
        solvpara1 = 0.5;
        solvpara2 = 0.75;
        p = 3;
    case 'MPRK43IIadap'
        solvpara1 = 0.563;
        p = 3;
    case 'MPRK22adap'
        solvpara1 = 1;
        p = 2;
    case 'ROS2'
        % solvpara1 = 1 + 0.5*sqrt(2);
        solvpara1 = 1/(2+sqrt(2));
        p = 2;
end





%%% tolerances
%tolint = [1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1];
%tolint = [1e-4,1e-3,1e-2,1e-1];
%tolint = [1e-8];
%%% parameter of test problems PR3 and PR4
paraint = 0:0.1:0.5;
%paraint = 0.4;

%%% string of test problems
%tests = {{'PR3',paraint},{'PR4',paraint},'robertson','strat'};
tests ={'hires8eq', 'npzd','robertson',{'PR4',paraint}};
%tests ={{'PR4',paraint}};
%tests ={'hires8eq'};

test_index = 1;
cost_old = 0;
cost_test = 0;
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

    %  if strcmp(Test, 'strat')
    %      tolint = [1e-6,1e-5,1e-4,1e-3,1e-2,1e-1];
    %  else
    tolint = [1e-1,1e-2,1e-3,1e-4,1e-5,1e-6,1e-7,1e-8];
    %  end

    lengthA = lengthpara*length(tolint);
    A = zeros(lengthA,2);
    k = 1;

    if flag
        for para=testpara
            for tol = tolint
                A(k,:) = [tol, para];
                k = k+1;
            end
        end
    else
        for tol = tolint
            A(k,:) = [tol, NaN];
            k = k +1;
        end
    end

    penalty = 1e6;
    flag_cont = true;
    flag_aux = true;

    for i = 1:lengthA
        if flag_cont
            if mod(i-1,length(tolint)) == 0
                err_old = penalty;
                steps_old = 1;
            end
            tol = A(i,1);
            para = A(i,2);
       	   % if vars.alpha2 < 1/6
           %	vars.alpha2 = 0;
           % end
           %[succstep, cnt_rej, err] = ODESolverTestSuite(vars.beta1, vars.beta2, vars.beta3, vars.alpha2, tol, Test,para,solverstr,solvpara1,solvpara2);
           if exist('solvpara2','var')
               [succstep, cnt_rej, err] = ODESolverTestSuite(vars.beta1, vars.beta2, vars.beta3, vars.alpha2, tol, Test, para, vars.kappa2, solverstr,solvpara1,solvpara2);
           else
               [succstep, cnt_rej, err] = ODESolverTestSuite(vars.beta1, vars.beta2, vars.beta3, vars.alpha2, tol, Test, para, vars.kappa2, solverstr,solvpara1);
           end
           if isnan(cnt_rej)
               cnt_rej = 10*penalty;
               err(1) = 1e10*c*tol; % this implies (in line 125) max(0, log( err(1) / (c * tol) ) / log(10) ) * penalty = 10*penalty
               flag_aux = false;
           elseif isnan(succstep)
               succstep = 10*penalty;
               err(1) = 1e10*c*tol; % this implies (in line 125) max(0, log( err(1) / (c * tol) ) / log(10) ) * penalty = 10*penalty
               flag_aux = false;
           end
           if isnan(err(1))
               cnt_rej = 10*penalty;
               succstep = 10*penalty;
               err(1) = 1e10*c*tol;
               flag_aux = false;
           end
           steps = succstep + cnt_rej;
           if flag_aux
               dy = log(err_old/err(1));
               dx = log(steps / steps_old);
               if  dy <= 0 ||   dx <= 0
                   flag_cont = false;
               elseif dy / dx <= 0.7 && tol < 1e-2
                   flag_cont = false;
               elseif dy / dx <= 0.35 && tol >= 1e-2
                   flag_cont = false;
               end
           end
           flag_aux = true;
           err_old = err(1);
           steps_old = steps;

           %cost = cost + succstep + cnt_rej + ( 1 - 1 / (1 + cnt_rej) ) + max(0, sign( log( err(1) / (c * tol) ) ) * 1e+6);
           %cost = cost + succstep + cnt_rej + sign(max(0,1000 - succstep/(1+cnt_rej)))*1e+6 + max(0, sign( log( err(1) / (c * tol) ) ) * 1e+6);
           %cost = cost + succstep + cnt_rej + max(0, log( err(1) / (c * tol) ) / log(10) ) * penalty;

           %cost = cost + (((log(succstep + cnt_rej))^2 + (log( err(1) / (c * tol) ) )^2) / log(10)^2 )^2 ;
           cost = cost + p * log(succstep + cnt_rej) +  log( err(1) / tol ) + max(0, log( err(1) / (c * tol) ))   ;
           if mod(i,length(tolint)) == 0
               cost_test(test_index) = atan((cost - cost_old)/100)^2;
               cost_old =  cost;
               if strcmp(Test,'PR4')
                   cost_test(test_index) = cost_test(test_index) / length(paraint);
               end
               test_index = test_index + 1;
           end
        else
            cost_test(end + 1) = atan((cost - cost_old)/100)^2;
            break
        end
    end
    if not(flag_cont)
        cost = sum(cost_test) + 10;
        break
    end
end
if flag_cont
    cost = sum(cost_test);
end

end

% function cost = fun2(vars)
% solverstr = 'MPRK43Iadap';
% c = 50; % An error of c*tol is tolerated in the cost function
% %Alpha2 = 0;
% %Beta3 = 0;
% %kappa2 = 1;
% cost= 0;
%
%
%
% %%%  Parameter of the method
%
% switch solverstr
%     case 'MPRK43Iadap'
%         solvpara1 = 0.5;
%         solvpara2 = 0.75;
%     case 'MPRK43IIadap'
%         solvpara1 = 0.563;
%     case 'MPRK22adap'
%         solvpara1 = 1;
% end
%
%
%
%
%
% %%% tolerances
% %tolint = [1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1];
% %tolint = [1e-4,1e-3,1e-2,1e-1];
% %tolint = [1e-8];
% %%% parameter of test problems PR3 and PR4
% paraint = 0:0.2:0.5;
% %paraint = [0.2];
%
% %%% string of test problems
% %tests = {{'PR3',paraint},{'PR4',paraint},'robertson','strat'};
% tests ={'robertson',{'PR4',paraint}, 'hires8eq'};
% %tests ={{'PR4',paraint}};
% %tests ={'robertson'};
%
%
% for testpairs = tests
%     testpair = testpairs{:};
%     if length(testpair) == 2
%         flag = true;
%         Test = testpair{1};
%         testpara =testpair{2};
%         lengthpara = length(testpara);
%     else
%         flag = false;
%         Test = testpair;
%         lengthpara = 1;
%     end
%
%   %  if strcmp(Test, 'strat')
%   %      tolint = [1e-6,1e-5,1e-4,1e-3,1e-2,1e-1];
%   %  else
%         tolint = [1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1];
%   %  end
%
%     lengthA = lengthpara*length(tolint);
%     A = zeros(lengthA,2);
%     k = 1;
%     for tol = tolint
%         if flag
%             for para=testpara
%                 A(k,:) = [tol, para];
%                 k = k+1;
%             end
%         else
%             A(k,:) = [tol, NaN];
%             k = k +1;
%         end
%     end
% penalty = 1e6;
%     for i = 1:lengthA
%         tol = A(i,1);
%         para = A(i,2);
%         %[succstep, cnt_rej, err] = ODESolverTestSuite(vars.beta1, vars.beta2, vars.beta3, vars.alpha2, tol, Test,para,solverstr,solvpara1,solvpara2);
%         if exist('solvpara2','var')
%             [succstep, cnt_rej, err] = ODESolverTestSuite(vars.beta1, vars.beta2, vars.beta3, vars.alpha2, tol, Test, para, vars.kappa2, solverstr,solvpara1,solvpara2);
%         else
%             [succstep, cnt_rej, err] = ODESolverTestSuite(vars.beta1, vars.beta2, vars.beta3, vars.alpha2, tol, Test, para, vars.kappa2, solverstr,solvpara1);
%         end
%         if isnan(cnt_rej)
%             cnt_rej = 10*penalty;
%             err(1) = 1e10*c*tol; % this implies (in line 125) max(0, log( err(1) / (c * tol) ) / log(10) ) * penalty = 10*penalty
%         elseif isnan(succstep)
%             succstep = 10*penalty;
%             err(1) = 1e10*c*tol; % this implies (in line 125) max(0, log( err(1) / (c * tol) ) / log(10) ) * penalty = 10*penalty
%         end
%         if isnan(err(1))
%             cnt_rej = 10*penalty;
%             succstep = 10*penalty;
%             err(1) = 1e10*c*tol;
%         end
%
% %cost = cost + succstep + cnt_rej + ( 1 - 1 / (1 + cnt_rej) ) + max(0, sign( log( err(1) / (c * tol) ) ) * 1e+6);
% %cost = cost + succstep + cnt_rej + sign(max(0,1000 - succstep/(1+cnt_rej)))*1e+6 + max(0, sign( log( err(1) / (c * tol) ) ) * 1e+6);
%  cost = cost + succstep + cnt_rej + max(0, log( err(1) / (c * tol) ) / log(10) ) * penalty;
%     ends
% end
%
% end



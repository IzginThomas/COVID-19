clc, clear, close all
warning('off','all')

Matlab_comp_bool = 1; % Work-precision Diagram comparing optimized MPRK and ode-solver
saveplot = true;
plotopt = {'LineWidth',4};
txtopt = {'FontSize',20};
txtoptlgnd = {'FontSize',18};
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%% Choose test problems %% If PR4 is used, put it always as the last test problem
paraint = [0.1 0.3 0.4 0.5];
paraint = 0.25;
%paraint = 0.4;
tests = {'PME', 'FokkerPlanck', 'robertson', 'hires8eq','npzd', 'brusselator',{'PR4',paraint}};
%tests = {'robertson',{'PR4',paraint}};
%tests = {'PME'};
%tests = {'robertson', 'hires8eq','npzd', {'PR4',paraint}};
%tests = {'brusselator',{'PR4',[0.1 0.3 0.5]}};

%% Choose method
solverstr = 'ROS2';
% The parameters will be set automatically:
switch solverstr
    case 'MPRK43Iadap'
        solvpara = [0.5,0.75];
        points = {[1.8476    -0.11068    -0.2863    -0.24606      2]; [1,0,0,0,1]; [0.6,-0.2,0,0,1];[0.7,-0.4,0,0,1]};
    case 'MPRK43IIadap'
        solvpara = 0.563;
        points = {[1.5193    -0.42576    -0.078535    -0.29465      2]; [1,0,0,0,1]; [0.6,-0.2,0,0,1];[0.7,-0.4,0,0,1]};
    case 'MPRK22adap' % [  0.98561    -0.2473    0.099065    0.19584      2],  c=5:[  0.99907    -0.068063    0.054175    0.17639      2]
        solvpara = 1;
        points = {[1.5193    -0.42576    -0.078535    -0.29465      2];[2,-1,0,-1,1]; [0.6,-0.2,0,0,1];[0.7,-0.4,0,0,1]};
    case 'ROS2'
        solvpara = 1/(2+sqrt(2));
        points = {[1.4783   -1.027    -0.27743    -0.23112     3];[1,0,0,0,1];[2,-1,0,-1,1];[0.6,-0.2,0,0,1];[0.7,-0.4,0,0,1]};
end



Markers= {'--d','--o','--*','--x','--+','--v','--^','--s','-->','--<'};
Markeropt ={'Markersize',10};

%% Choose which parameters to compare

%points={...
%  [0.6,-0.2,0,0,1];[0.7,-0.4,0,0,1];[2,-1,0,-1,1];[1/6,-1/3,0,0,1];[1/6,1/6,0,0,1];...
% [1,0,0,0,1];[1/2,1/2,0,1/2,1];[1/18,1/9,1/18,0,1];[1/4,1/4,1/4,0,1]}; %standard
points = {}; % To skip work-precision for comparing parameters

if ~isempty(points)
    % figure
    %% creating the legend

    lgnd = cell(1,length(points));

    for i=1:length(points)
        point = points{i};
        beta1 = point(1);
        beta2 = point(2);
        beta3 = point(3);
        alpha2 = point(4);
        kappa2 = point(5);
        lgnd{i} = sprintf('(%.3g, %.3g, %.3g, %.3g, %.3g)',beta1,beta2,beta3,alpha2,kappa2);
    end
    testset={};
    for i=1:length(tests)
        ind = 1;
        if length(tests{i}) == 2
            while ind <= length(tests{i}{2})
                teststr = tests{i}{1};
                testpara = tests{i}{2}(ind);
                testset{end + 1} = {teststr, num2str(testpara)};
                ind = ind + 1;
            end
        else
            testset{end + 1} = tests{i};
        end
    end
    %% for each point compute the costs, number of succesful steps, rejected steps and errors
    % format:
    % "costs" is an array, containing one element for each parameter combination in
    % "points". Each element in costs is a matrix with a new row for a different test
    % problem. Each row is of the form [cost_step, cost_tol, cost_step +
    % cost_tol] containing the different parts of the cost function sumed over
    % all test cases and tolerances and the cost function itself.
    %
    %"steps" is an array, each array element containing the information on the number
    % of successful steps for each parameter combination in "points".
    % Each array element is a matrix, again with a new row for a different test
    % problem. Each row contains the number of successful steps for each of the
    % eight tolerances [1e-8, ..., 1e-1] (Also for strat!)
    %
    %"rejected" and "error" are similar to "steps" but contain the information
    % of rejected steps and errors, respectively.

    for i=1:length(points)
        point = points{i};
        beta1 = point(1);
        beta2 = point(2);
        beta3 = point(3);
        alpha2 = point(4);
        kappa2 = point(5);
        [costs{i}, steps{i}, rejected{i}, errors{i},~] = workprecision(beta1,beta2,beta3,alpha2,kappa2,tests,solverstr,solvpara);
        for j = 1:size(errors{i},1)
            test_errors = errors{i}(j, :);
            total_steps = steps{i}(j, :) + rejected{i}(j, :);
            subplot(1,size(errors{i},1), j);
            loglog(total_steps,test_errors,Markers{i}, plotopt{:}, Markeropt{:})
            if length(string(testset{j})) == 2
                teststr = testset{j}{1};
                testpara = testset{j}{2};
                title([teststr,' p = ',testpara]);
            else
                switch tests{j}
                    case 'hires8eq'
                        title('HIRES');
                    case 'npzd'
                        title('NPZD');
                    case 'strat'
                        title('Strat');
                    case 'robertson'
                        title('Robertson');
                    case 'brusselator'
                        title('Brusselator');
                    case 'PME'
                        title('PME');
                    case 'FokkerPlanck'
                        title('Fokker-Planck');
                    otherwise
                        title('SetTitle!');
                end
            end
            xlabel("Number of total steps",txtopt{:},'Interpreter','latex');
            ylabel("Error",txtopt{:},'Interpreter','latex');
            grid on
            l = legend(lgnd,'location','sw','ItemHitFcn',@cb_legend);
            set(gca,txtoptlgnd{:})
            l.Visible = 'on';
            hold on
        end

    end
end
if saveplot
    Sub_Fig_Divider(solverstr)
end
hold off
%%
str = {'PME'}; % Choose test problem for comparison with build-in solvers
run(['test_' str{1}])
if length(str) == 2
    [~,~,f] = PD(0,y0,str{2});
else
    [~,~,f] = PD(0,y0);
end
tint = [tstart tend];
kend = 10;
for k = 1:kend
    if Matlab_comp_bool
        k
        lgnd = cell(1,5);
        tolint = [1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1];


        [~, steps2, rejected2, errors2, times2(:,k)] = workprecision(1.4783,   -1.027,    -0.27743,   -0.23112 ,    3,str,'ROS2',1/(2+sqrt(2)));
        total_nfevals2 = 2*(steps2 + rejected2);
        solves2 = steps2 - length(y0) - 1 + rejected2;
        if k == kend
            figure
            loglog(total_nfevals2, errors2, Markers{1}, plotopt{:}, Markeropt{:})
            lgnd{1} = sprintf('ROS2');
            hold on
        end
        [~, steps22, rejected22, errors22, times22(:,k)] = workprecision(1.5193 ,   -0.42576 ,   -0.078535   , -0.29465   ,   2,str,'MPRK22adap',1);
        total_nfevals22 = 2*(steps22 + rejected22);
        solves22 = 2*(steps22 + rejected22);
        if k == kend
            loglog(total_nfevals22, errors22, Markers{1}, plotopt{:}, Markeropt{:})
            lgnd{2} = sprintf('MPRK22');
        end

        [~, steps43I, rejected43I, errors43I, times43I(:,k)] = workprecision(1.8476, -0.11068 , -0.2863 ,   -0.24606,      2,str,'MPRK43Iadap',[0.5, 0.75]);
        total_nfevals43I =3*(steps43I + rejected43I);
        solves43I = 4*(steps43I + rejected43I);
        if k == kend
            loglog(total_nfevals43I, errors43I, Markers{2},plotopt{:}, Markeropt{:})
            lgnd{3} = sprintf('MPRK43I');
        end

        [~, steps43II, rejected43II, errors43II, times43II(:,k)] = workprecision(1.5193 ,   -0.42576 ,   -0.078535   , -0.29465   ,   2,str,'MPRK43IIadap',0.563);
        total_nfevals43II = 3*(steps43II + rejected43II);
        solves43II = 4*(steps43II + rejected43II);
        if k == kend
            loglog(total_nfevals43II, errors43II, Markers{3}, plotopt{:}, Markeropt{:})
            lgnd{4} = sprintf('MPRK43II');
        end

        errode15 = zeros(1,8);
        errode23 = errode15;
        total_nfevals_ode15 = errode15;
        total_nfevals_ode23 = errode15;
        for i =1:length(tolint)
            tol = tolint(i);
            options.AbsTol = tol;
            options.RelTol = tol;
            options.NonNegative = 1;
            tic
            solode15 = ode15s(f,tint,y0,options);
            total_etime_ode15(i,k) = toc;
            total_solves_ode15(i) = solode15.stats.nsolves;
            total_nfevals_ode15(i) = solode15.stats.nfevals;
            tic
            solode23 = ode23s(f,tint,y0,options);
            total_etime_ode23(i,k) = toc;
            total_solves_ode23(i) = solode15.stats.nsolves;
            total_nfevals_ode23(i) = solode23.stats.nfevals;

            opt.AbsTol = 1e-13;
            opt.RelTol = 1e-13;
            opt.NonNegative = [];
            [~, yrefode15] = ode15s(f, solode15.x, y0, opt);
            [~, yrefode23] = ode15s(f, solode23.x, y0, opt);
            if exist('pde_plot','var') && pde_plot
                errode15(i) = L2err_rel2d(solode15.x, x, solode15.y, yrefode15);
                errode23(i) = L2err_rel2d(solode23.x, x, solode23.y, yrefode23);    
            else
                errode15(i) = L2err_rel(solode15.x, solode15.y, yrefode15);
                errode23(i) = L2err_rel(solode23.x, solode23.y, yrefode23);
            end
        end
        if k == kend
            loglog(total_nfevals_ode15, errode15, Markers{4}, plotopt{:}, Markeropt{:})
            lgnd{5} = sprintf('ode15s');
            loglog(total_nfevals_ode23, errode23, Markers{5}, plotopt{:}, Markeropt{:})
            lgnd{6} = sprintf('ode23s');
        end
        switch str{1}
            case 'hires8eq'
                titlestr = 'WP Diagram: HIRES';
            case 'npzd'
                titlestr = 'WP Diagram: NPZD';
            case 'strat'
                titlestr = 'WP Diagram: Strat';
            case 'robertson'
                titlestr = 'WP Diagram: Robertson';
            case 'brusselator'
                titlestr = 'WP Diagram: Brusselator';
            case 'PR4'
                titlestr = strcat('WP Diagram: PR4 p=', num2str(str{2}));
            case 'FokkerPlanck'
                titlestr = 'WP Diagram: Fokker-Planck';
            case 'PME'
                 titlestr = 'WP Diagram: PME';
            otherwise
                titlestr = 'SetTitle!';
        end
        if k == kend
            title(titlestr)
            xlabel("Number of RHS evaluations",txtopt{:},'Interpreter','latex');
            ylabel("Error",txtopt{:},'Interpreter','latex');
            grid on
            l = legend(lgnd,'location','ne','ItemHitFcn',@cb_legend);
            set(gca,txtoptlgnd{:})
            l.Visible = 'on';
            hold off
            if length(str) == 2
                print('-depsc2',['./save/' str{1} num2str(str{2}) '_WP_RHS.eps']);
            else
                print('-depsc2',['./save/' str{1} '_WP_RHS.eps']);
            end
            figure
            loglog(solves2, errors2, Markers{1}, plotopt{:}, Markeropt{:})
            lgnd{1} = sprintf('ROS2');
            hold on
            loglog(solves22, errors22, Markers{1}, plotopt{:}, Markeropt{:})
            lgnd{2} = sprintf('MPRK22');
            loglog(solves22, errors43I, Markers{2},plotopt{:}, Markeropt{:})
            lgnd{3} = sprintf('MPRK43I');
            loglog(solves22, errors43II, Markers{3}, plotopt{:}, Markeropt{:})
            lgnd{4} = sprintf('MPRK43II');
            loglog(total_solves_ode15, errode15, Markers{4}, plotopt{:}, Markeropt{:})
            lgnd{5} = sprintf('ode15s');
            loglog(total_solves_ode23, errode23, Markers{5}, plotopt{:}, Markeropt{:})
            lgnd{6} = sprintf('ode23s');


            title(titlestr)
            xlabel("Total number of solved linear systems",txtopt{:},'Interpreter','latex');
            ylabel("Error",txtopt{:},'Interpreter','latex');
            grid on
            l = legend(lgnd,'location','ne','ItemHitFcn',@cb_legend);
            set(gca,txtoptlgnd{:})
            l.Visible = 'on';
            hold off
            
            if length(str) == 2
                print('-depsc2',['./save/' str{1} num2str(str{2}) '_WP_solves.eps']);
            else
                print('-depsc2',['./save/' str{1} '_WP_solves.eps']);
            end
        end
    end
end
if Matlab_comp_bool
    times2(:,1) = sum(times2,2) / size(times2,1);
    times22(:,1) = sum(times22,2) / size(times22,1);
    times43I(:,1) = sum(times43I,2) / size(times43I,1);
    times43II(:,1) = sum(times43II,2) / size(times43II,1);
    total_etime_ode15(:,1) = sum(total_etime_ode15,2);
    total_etime_ode23(:,1) = sum(total_etime_ode23,2);
    figure
    loglog(times2(:,1), errors2, Markers{1}, plotopt{:}, Markeropt{:})
    lgnd{1} = sprintf('ROS2');
    hold on
    loglog(times22(:,1), errors22, Markers{1}, plotopt{:}, Markeropt{:})
    lgnd{2} = sprintf('MPRK22');
    loglog(times43I(:,1), errors43I, Markers{2},plotopt{:}, Markeropt{:})
    lgnd{3} = sprintf('MPRK43I');
    loglog(times43II(:,1), errors43II, Markers{3}, plotopt{:}, Markeropt{:})
    lgnd{4} = sprintf('MPRK43II');
    loglog(total_etime_ode15(:,1), errode15, Markers{4}, plotopt{:}, Markeropt{:})
    lgnd{5} = sprintf('ode15s');
    loglog(total_etime_ode23(:,1), errode23, Markers{5}, plotopt{:}, Markeropt{:})
    lgnd{6} = sprintf('ode23s');

    title(titlestr)
    xlabel("Averaged elapsed time [s]",txtopt{:},'Interpreter','latex');
    ylabel("Error",txtopt{:},'Interpreter','latex');
    grid on
    l = legend(lgnd,'location','sw','ItemHitFcn',@cb_legend);
    set(gca,txtoptlgnd{:})
    l.Visible = 'on';
    hold off
    if length(str) == 2
        print('-depsc2',['./save/' str{1} num2str(str{2}) '_WP_time.eps']);
    else
        print('-depsc2',['./save/' str{1} '_WP_time.eps']);
    end
end
function cb_legend(~,evt)
if strcmp(evt.Peer.Visible,'on')
    evt.Peer.Visible = 'off';
else
    evt.Peer.Visible = 'on';
end
end




function [costs, steps, rejected, errors, times] = workprecision(beta1,beta2,beta3,alpha2,kappa2,tests,solverstr,solvpara)
c = 50; %bound for tolerance in cost function

if length(solvpara) == 1
    solvpara1 = solvpara;
elseif length(solvpara) == 2
    solvpara1 = solvpara(1);
    solvpara2 = solvpara(2);
else
    error('Too many parameter inputs for method')
end
cnt_loop = 1;
%if strcmp(Test, 'strat')
%   tolint = [1e-6,1e-5,1e-4,1e-3,1e-2,1e-1];
% else
tolint = [1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1];
%tolint = [1e-2,1e-1];
% end
lengthtolint = length(tolint);

cnt_para = 0;
total_paras = 0;
for testpairs = tests
    testpair = testpairs{:};
    if length(testpair) == 2
        cnt_para = cnt_para + 1;
        testpara = testpair{2};
        total_paras = total_paras + length(testpara);
    end
end
steps = zeros(length(tests) - cnt_para + total_paras, lengthtolint);
costs = zeros(length(tests) - cnt_para + total_paras, 3);
errors = zeros(length(tests) - cnt_para + total_paras, lengthtolint);




for testpairs = tests
    testpair = testpairs{:};
    if length(testpair) == 2
        flag = true;
        Test = testpair{1};
        testpara = testpair{2};
        lengthpara = length(testpara);
    else
        flag = false;
        Test = testpair;
        lengthpara = 1;
    end


    cost_step = 0;
    %cost_rej = 0;
    cost_tol = 0;

    rejected = steps;



    lengthA = lengthpara * lengthtolint;


    A = zeros(lengthA,2);
    k = 1;


    if flag
        for para=testpara
            for tol = tolint
                A(k,:) = [tol, para];
                k = k + 1;
            end
        end
    else
        for tol = tolint
            A(k,:) = [tol, NaN];
            k = k + 1;
        end
    end

    for i = 1:lengthA
        tol = A(i,1);
        para = A(i,2);
        k = mod(i - 1, lengthtolint) + 1;
        %[succstep, cnt_rej, err] = ODESolverTestSuite(vars.beta1, vars.beta2, vars.beta3, vars.alpha2, tol, Test,para,solverstr,solvpara1,solvpara2);
        if exist('solvpara2','var')
            [succstep, cnt_rej, err] = ODESolverTestSuite(beta1, beta2, beta3, alpha2, tol, Test, para, kappa2, solverstr,solvpara1,solvpara2);
        else
            [succstep, cnt_rej, err] = ODESolverTestSuite(beta1, beta2, beta3, alpha2, tol, Test, para, kappa2, solverstr,solvpara1);
        end
        times(i) = err(2);
        if isnan(cnt_rej)
            cnt_rej = 1e+6;
            err(1) = 1;  %NaN %10^10 * c * tol;
        elseif isnan(succstep)
            succstep = 1e+6;
            err(1) = 1; %NaN %10^10 * c * tol;
        end
        steps(cnt_loop, k) = steps(cnt_loop, k) + succstep;
        rejected(cnt_loop,k) = rejected(cnt_loop, k) + cnt_rej;
        errors(cnt_loop, k) = err(1);

        cost_step = cost_step + succstep + cnt_rej;
        % cost_rej = cost_rej; %+ max(0, log(cnt_rej / p*succstep)) * 1e+6;%sign( max( 0, 1000 - succstep / ( 1 + cnt_rej ) ) ) * 1e+6;
        cost_tol = cost_tol + max(0,  log( err(1) / (c * tol) ) / log(10) ) * 1e+5;

        if flag && mod(i, lengthtolint) == 0
            %             costs(cnt_loop, :) = [cost_step, cost_rej, cost_tol, cost_step + cost_rej + cost_tol];
            costs(cnt_loop, :) = [cost_step, cost_tol, cost_step + cost_tol];
            cnt_loop = cnt_loop + 1;
        end

    end
    if ~flag
        %costs(cnt_loop, :) = [cost_step, cost_rej, cost_tol, cost_step + cost_rej + cost_tol];
        costs(cnt_loop, :) = [cost_step, cost_tol, cost_step + cost_tol];
        cnt_loop = cnt_loop + 1;
    end
end

end

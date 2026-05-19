function [varargout] = ODESolverTestSuite(varargin)
warning('off','all')
%clc, clear, close all
if nargin > 7
    Beta1 = varargin{1};
    Beta2 = varargin{2};
    Beta3 = varargin{3};
    Alpha2 = varargin{4};
    tol = varargin{5};
    test = varargin{6};
    para = varargin{7};
    kappa2 = varargin{8};

    assert(~and( or(strcmp(test,'PR3'),strcmp(test,'PR4')), isnan(para)),'invalid range for para')

    solverstr = varargin{9};
    solverpara1 = varargin{10};
    if nargin >10
        solverpara2 = varargin{11};
    end
    task = 'data';
else
    clc, clear, close all
    solverstr ='ROS2';
    switch solverstr
        case 'MPRK43Iadap'
            solverpara1 = 0.5;
            solverpara2 = 0.75;
            Beta1 =1.8476;
            Beta2 =      -0.11068;
            Beta3 =  -0.2863;
            Alpha2 = -0.24606;
            kappa2 = 2;
        case 'MPRK43IIadap'
            solverpara1 = 0.563;

            Beta1 =1.5193 ;
            Beta2 =    -0.42576    ;
            Beta3 = -0.078535  ;
            Alpha2 =  -0.29465 ;
            kappa2 = 2;
        case 'MPRK22adap'
            solverpara1 = 1;
            Beta1 = 1.5193;
            Beta2 =      -0.42576;
            Beta3 =  -0.078535;
            Alpha2 = -0.29465;
            kappa2 = 2;
        case 'ROS2'
            % solverpara1 = 1 + 0.5*sqrt(2);
            solverpara1 = 1/(2+sqrt(2));
            Beta1 = 1.4783;
            Beta2 =    -1.027;
            Beta3 =  -0.27743;
            Alpha2 = -0.23112;
            kappa2 = 3;
        otherwise
            % solverpara1 = 1 + 0.5*sqrt(2);
            Beta1 = 0.7;
            Beta2 = -0.4;
            Beta3 =  0;
            Alpha2 = 0;
            kappa2 = 1;
    end
    tol = 1e-6;
    test ='npzd';
    if strcmp(test,'PR3') || strcmp(test,'PR4')
        para = 0.4;
    else
        para = NaN;
    end
    task = 'plot';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ODE Solver Test Suite                                               %%%
%%%                                                                     %%%
%%% Author: Stefan Kopecz (kopecz@mathematik.uni-kassel.de)             %%%
%%% Date: 05/07/2019                                                    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% References %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% [ 1]  Burchardt et. al.; A high-order conservative Patankar-type
%%%       discretisation for stiff systems of production–destruction
%%%       equations; APNUM; 2003
%%% [ 2]  Kopecz, Meister; On order conditions for modified
%%%       Patankar-Runge-Kutta schemes; APNUM; 2018
%%% [ 3]  Kopecz, Meister; Unconditionally positive and conservative third
%%%       order modified Patankar-Runge-Kutta discretizations of
%%%       production-destruction systems; BIT; 2018
%%% [ 4]  Meister, Butcher; Sensitivity of modified Patankar-type schemes
%%%       for systems of conservative productiondestruction equations,
%%%       AIP Conference Proceedings 1863; 2017
%%% [ 5]  Martin-Vaquero et. al; Higher-order nonstandard finite difference
%%%       schemes for a MSEIR model for a malware propagation; Journal of
%%%       computational and applied mathematics; 2017
%%% [ 6]  Mickens, Washington; NSFD discretizations of interacting
%%%       population models satisfying conservation laws; Computers and
%%%       Mathematics with Applications; 2013
%%% [ 7]  Formaggia, Scotti; POSITIVITY AND CONSERVATION PROPERTIES OF SOME
%%%       INTEGRATION SCHEMES FOR MASS ACTION KINETICS; SIAM J. Numer.
%%%       Anal.; 2011
%%% [ 8]  Friedmann et. al; Well-Posedness of a Linear Spatio-Temporal
%%%       Model of the JAK2/STAT5 Signaling Pathway; Comm. Math. Anal.;
%%%       2013
%%% [ 9]  Leveque; Finite Difference Methods for Ordinary and Partial
%%%       Differential Equations; SIAM; 2007
%%% [10]  Shu, Osher: Efficient implementation of essentially
%%%       non-oscillatory shock-capturing schemes; J. Comput. Phys.; 1988
%%% [11]  Hundsdorfer, Verwer; Numerical solution of time-dependent
%%%       Advection-Diffusion-Reaction equations
%%% [12]  Bruggeman et. al.; A second-order, unconditionally positive,
%%%       mass-conserving integration scheme for biochemical systems;
%%%       APNUM; 2007
%%% [13]  Broekhuizen et. al.; An improved and generalized second order,
%%%       unconditionally positive , mass conserving integration scheme for
%%%       biochemical systems; APNUM; 2008
%%% [14]  Schippmann, Burchard; Rosenbrock methods in biogeochemical
%%%       modelling - A comparison to Runge-Kutta methods and modified
%%%       Patankar schemes; Ocean Modelling; 2011
%%% [15]  Bertolazzi; Positive and Conservative Schemes for Mass Action
%%%       Kinetics; Computers Math. Applic.; 1996
%%% [16]  Burchard et. al; Application of modified Patankar schemes to
%%%       stiff biogeochemical models for the water column; Ocean Dynamics;
%%%       2005
%%% [17]  Huang, Shu: Positivity-Preserving Time Discretizations for
%%%       Production–Destruction Equations with Applications to
%%%       Non-equilibrium Flows; J. Sci. Comput; 2018
%%% [18]  https://en.wikipedia.org/wiki/Lotka%E2%80%93Volterra_equations#A_simple_example
%%% [19]  Kopecz: Ein gekoppeltes Finite-Elemente/Discontinuous-Galerkin-
%%%       Verfahren für Strömungs-Transport-Probleme; 2012
%%% [20]  Cresson, Pierret: Non standard finite difference scheme preserving
%%%       dynamical properties; J. Comput. Appl. Math.; 2016
%%% [21]  Luke, N.S., DeVito, M.J., Shah, I. et al. Bull. Math. Biol. (2010)
%%%       72: 1799. https://doi.org/10.1007/s11538-010-9508-5
%%% [22]  Voth: Modellierung und Simulation mit gewöhnlichen Differential-
%%%       gleichungen für die Therapie retinaler Erkrankungen, BA FH
%%%       Bielefeld, 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% Tasks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% plot:       Plots approximations computed with the different solvers
%%%             specified below
%%% plotsave:   Same as plot + EPS files are saved in the subdirectory
%%%             'save'
%%% error:      Generates error plots for the different solvers
%%%             specified below
%%% efficiency: Plots errors over time
%%% adaptivity: Plots errors over CPU time for different Tolerances
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% Select Test Case %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Autonomous test cases
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% linmod:         Linear model from [1]
%%% nonlinmod:      Nonlinear model from [1]
%%% robertson:      Robertson test problem from [1]
%%% brusselator:    Brusselator test problem from [2]
%%% sir:            SIR model from [6]
%%% mseir:          MSEIR model from [5]
%%% leveque:        Stiff linear kinetics model from [9], Example 8.2
%%% stifflin:       Stiff linear test
%%% stiffnonlin:    Stiff nonlinear test
%%% upwind:         1st order upwind advection
%%% bertolazzi:     Nonlinear problem from [15]
%%% kavitation:     Nonlinear problem from [19]
%%% npzd:           Nonlinear model from [16]
%%% npd0d:          Nonlinear model from [14]
%%% nonlinBBKS:     Nonlinear test from [12]
%%% conservation:   Nonlinear test with dim ker S^T > 1
%%% rifampicin:     Problem from [21]
%%% vegf:           Problem from [22]
%%% babooncheetah:  Lotka-Volterra system from [19]
%%% CressonPierret: Problem from [20]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Non-autonomous test cases
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% lin_nonauto:   Non-autonomous linear test
%%% jakstat:       JAKSTAT model from [8], see equations (3.10 - 3.14)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% TESTCASES CAN BE ADAPTED IN THE FILES test_<testcase_name>.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%% Select Methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Runge-Kutta and Rosenbrock schemes:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FE:            Forward Euler method
%%% RK22:          Explict 2nd order RK
%%%                {'RK22', alpha} with alpha ~= 0
%%% RK3I:          Explicit 3rd order RK, see [3, Lemma 2.3]
%%%                {'RK3I', alpha, beta} with alpha ~= beta, alpha ~= 0,
%%%                beta ~= 0, alpha ~= 2/3
%%% RK3II:         Explicit 3rd order RK, see [3, Lemma 2.3]
%%%                {'RK3II, alpha} with alpha ~= 0
%%% SSP33:         Optimal third order SSP scheme, see [10]
%%% SSP33adap:     Time adaptive variant of SSP33 using Heun's 2nd order
%%%                scheme to estimate the error. Entspricht RKF2(3).
%%%                {'SSP33adap', [RelTol AbsTol]}
%%% RK4:           Explicit Runge-Kutta 4th order
%%% BE:            Backward Euler method
%%% SDIRK2:        SDIRK 2nd order
%%% SDIRK3:        SDIRK 3rd order
%%% ROS1:          Linearized theta method, see [11].
%%% ROS2 :
%%% ode23wrapper:  Wrapper for MATLAB's ode23 method
%%% ode23swrapper: Wrapper for MATLAB's ode23s method
%%% ode15swrapper: Wrapper for MATLAB's ode15s method
%%% ode45wrapper:  Wrapper for MATLAB's ode45 method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Multistep methods:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% BDF3:         BDF3 scheme
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Modified Patankar-Runge-Kutta schemes:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MPE:          Modified Patankar-Euler, see [1]
%%% MPRK:         Modified Patankar-Runge-Kutta, see [1], is equivalent to
%%%               {'MPRK22',1}
%%% MPRK22:       Second order MPRK with delta = 1, see [2]
%%%               {'MPRK22', alpha} with alpha >= 0.5
%%% MPRK22adap:   Second order MPRK with delta = 1, see [2]
%%%               {'MPRK22', alpha, [RelTol AbsTol]} with positive alpha
%%%               when problem has remainder
%%% MPRK22ncs:    Second order MPRK with delta = 0, see [2]
%%%               {'MPRK22', alpha} with alpha >= 0.5
%%% MPRK32:       Three stage, 2nd order scheme, similar to MPRK43,
%%%               unpublished
%%% MPRK32ncs:    Three stage, 2nd order scheme, similar to MPRK43ncs,
%%%               unpublished
%%% NSMPRK22:     Non standard 2nd order MPRK scheme, unpublished
%%%               {'NSMPRK22', alpha} with alpha >= 0.5
%%% MPRK43I:      Third order MPRK with delta = 1, see [3]
%%%               {'MPRK43I, alpha, beta}, see Lemma 2.4 in [3] for
%%%               restrictions on alpha and beta
%%% MPRK43Incs:   Third order MPRK with delta = 0, see [3]
%%%               {'MPRK43Incs, alpha, beta}, see Lemma 2.4 in [3] for
%%%               restrictions on alpha and beta
%%% MPRK43II:      Third order MPRK with delta = 1, see [3]
%%%               {'MPRK43II', alpha} with 0.375 <= alpha <= 0.75
%%% MPRK43IIncs:   Third order MPRK with delta = 0, see [3]
%%%               {'MPRK43II', alpha} with 0.375 <= alpha <= 0.75
%%% MPRK43Iadap:  Adaptive MPRK43I scheme, unpublished
%%%               {'MPRK43Iadap', alpha, beta, [RelTol AbsTol]}
%%% MPRK43IIadap: Adaptive MPRK43II scheme, unpublished
%%%               {'MPRK43Iadap', alpha, [RelTol AbsTol]}
%%% SSPMPRK22:    SSP variant of MPRK22
%%%               {'SSPMPRK22', alpha, beta}
%%% mpdec:        Modified Patankar Deferred Correction schemes
%%%               {'mpdec',M}, M = 0, 1, 2, ...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MPRK schemes for specific test problems:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MPElin:       1-stage, 2nd order scheme for linmod, see [2]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Other schemes:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PE:           Patankar-Euler scheme from [1]
%%% CMPE:         CMPE scheme from [4]
%%%               {'CMPE', alpha} with 0 <= alpha <= 1
%%% CMPRK:        CMPRK scheme similiar to CMPE, unpublished
%%%               {'CMPRK', alpha} with 0 <= alpha <= 1
%%% PCMP:         2nd to 3rd order multistep Patankar-type scheme, see [7]
%%% BBKS1:        1st order scheme from [12]
%%% BBKS2:        2nd order scheme from [12]
%%% mBBKS1:       1st order scheme from [13]
%%% mBBKS2:       2nd order scheme from [13]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Other schemes for specific test problems:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SDIRK2stifflin:
%%%               SDIRK2 for the stifflin test problem. L-stable method.
%%% ImplTrapstifflin:
%%%               Implicit trapezoidal rule for stifflin test problem.
%%%               A-stable but not L-stable.
%%% NSFDSIR:      Nonstandard FD scheme for the SIR test, see [6] equations
%%%               (50-52)
%%% NSFDMSEIR:    Four nonstandard FD schemes to solve the MSEIR test,
%%%               see [5]
%%%               {'NSFDMSEIR', n} with n = 1,2,3,4
%%% NSFDMSEIRextp2:
%%%               2nd order extrapolated versions of NSFDMSEIR, see (10) in
%%%               [5]. Not unconditionally positive!!!
%%%               {'NSFDMSEIRextp2', n} with n = 1,2,3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SELECT SOLVERS HERE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%
if exist('solverpara1','var')
    if exist('solverpara2','var')
        solvers = {{solverstr,solverpara1,solverpara2,[tol tol]}};
    else
        solvers = {{solverstr,solverpara1,[tol tol]}};
    end
else
    solvers = {solverstr};
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MAIN PROGRAM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Read Test Description
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


run(['test_' test])
tint = [tstart tend];

if ~exist('xscale','var')
    xscale = 'linear';
end
if ~exist('yscale','var')
    yscale = 'linear';
end
if ~exist('scaleode','var')
    scaleode = 1;
end
if ~exist('pde_plot','var')
    pde_plot = false;
end
if ~exist('one_dim_prob','var')
    one_dim_prob = false;
end
defaultMarkers = {'d','o','*','x','+','v','-^','-s','->','-<'};
smoothMarkers = {'-','-','-','-','-','-','-','-s','->','-<'};
errMarkers = {'-d','-o','-*','-x','-+','-v','-^','-s','->','-<'};
% Markersref = {'--d','--o','--*','--x','--+','--v','--^','--s','-->','--<'};
Markersref = {'--','--','--','--','--','--','--','--s','-->','--<'};
%Markersref = {'--d','--o','--*','--x','--+','--v','--^','--s','-->','--<'};
Markers= errMarkers;
Markeropt = {'Markersize',1};
default = 1;
Markerint = default;
plotopt = {'LineWidth',2};
txtopt = {'FontSize',14};
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0; %check if octave
switch task
    %%% Plot approximations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case {'plot','plotsave', 'data'}
        if ~strcmp(task,'data')
            if ~exist('exsol', 'var')
                if ~isOctave
                    options.AbsTol = 1e-13;
                    options.RelTol = 1e-13;
                    if isnan(para)
                        [~,~,f] = PD(0,y0/scaleode,scaleode);
                    else
                        [~,~,f] = PD(0,y0/scaleode,scaleode,para);
                    end
                    options.NonNegative = 1:length(y0);
                    options.NonNegative = [];
                    [tex,yex] = ode15s(f,tint,y0,options);
                    %[tex,yex] = ode23s(f,tint,y0,options);
                else
                    warning('Octave is currently not supported.')
                    return
                end
            else
                if ~ pde_plot
                    tex = linspace(tstart,tend,(tend-tstart)*1000);
                    yex =exsol(tex');
                else
                    tex = linspace(tstart,tend,(tend-tstart)*1000)';
                    yex = exsol(tex,x);
                end
            end
            yex = yex/scaleode;

            %%% Plot options
            plotopt = {'LineWidth',2};
            txtopt = {'FontSize',14};

            %%% Reference solution


            if ~exist('E','var')
                E = ones(1,size(yex,2));
            end

            if isempty(E)
                sumyex = [];
            else
                sumyex = E*yex';
            end

            yex = [yex, sumyex'];

            if exist('visfac','var')
                yex = yex.*visfac;
            end

            %%% Plotting
            lgnd = cell(1,2*(length(y0) + size(E,1)));
            lgnderr = lgnd;
            for i = 1:length(y0)
                if exist('visfac','var') && visfac(i) ~= 1
                    lgnd{i} = sprintf('%i y_%i', visfac(i), i);
                    lgnd2{i} = sprintf('%i err_%i', visfac(i), i);
                else
                    lgnd{i} = sprintf('y_%i', i);
                    lgnd2{i} = sprintf('err_%i', i);
                end
                if exist('visfac','var') && visfac(i) ~= 1
                    lgnd{length(y0) + size(E,1) + i} = sprintf('%i y_%i num', visfac(i), i);
                else
                    lgnd{length(y0) + size(E,1) + i} = sprintf('y_%i num', i);
                end
            end
            if exist('visfac','var') && visfac(end) ~= 1
                for i = 1:size(E,1)
                    lgnd{length(y0) + i} = [num2str(visfac(end)) '\Sigma y_i'];
                    lgnd{2*length(y0) + size(E,1) + i} = [num2str(visfac(end)) '\Sigma y_i'];
                end
            else
                if all(all(E))
                    for i = 1:size(E,1)
                        lgnd{length(y0) + i} = '\Sigma y_i';
                        lgnd{2*length(y0) + size(E,1) + i} = '\Sigma y_i num';
                    end
                else
                    for i = 1:size(E,1)
                        lgnd{length(y0) + i} = sprintf('\\Sigma e_{%ii}y_i',i);
                        lgnd{2*length(y0) + size(E,1) + i} = sprintf('\\Sigma e_{%ii}y_i num',i);
                    end
                end
            end

            if ~exist('loc','var')
                loc = 'ne';
            end

            % figure
            if strcmp(test,'hires')
                hdl1 = semilogx(tex,yex,plotopt{:});
            elseif strcmp(test,'strat')
                tstart = tstart/3600;
                tend = tend/3600 +15;
                %tend = tend +15;
                %tex = tex + 12;
                tex_temp = tex/3600 + 12;

                %yex = (diag([9.906e1; 6.624e8; 5.326e11; 1.697e16; 4.000e16; 1.093e9])*yex')';
                subplot(3,2,5)
                plot(tex_temp,yex(:,1),plotopt{:})
                title('O^{1D}')
                xlim([tstart tend])
                xlabel('t',txtopt{:})

                subplot(3,2,2)
                plot(tex_temp,yex(:,2),plotopt{:})
                title('O')
                xlim([tstart tend])
                xlabel('t',txtopt{:})

                subplot(3,2,3)
                plot(tex_temp,yex(:,3),plotopt{:})
                title('O_3')
                xlim([tstart tend])
                xlabel('t',txtopt{:})

                subplot(3,2,1)
                plot(tex_temp,yex(:,4) - 2*1.697e+16,plotopt{:})
                title('O_2 - 2*1.697e+16')
                xlim([tstart tend])
                xlabel('t',txtopt{:})

                subplot(3,2,4)
                plot(tex_temp,yex(:,5),plotopt{:})
                title('NO')
                xlim([tstart tend])
                xlabel('t',txtopt{:})

                subplot(3,2,6)
                plot(tex_temp,yex(:,6),plotopt{:})
                title('NO_2')
                xlim([tstart tend])
                xlabel('t',txtopt{:})


            elseif strcmp(test,'reducednucreact')

                subplot(2,2,1)
                plot(tex,yex(:,1),plotopt{:})
                title('Neutron Density')
                grid on;

                subplot(2,2,2)
                plot(tex,yex(:,2),plotopt{:})
                title('Precursor Density')
                grid on;

                subplot(2,2,3)
                plot(tex,yex(:,3),plotopt{:})
                title('Reactivity')
                grid on;

                subplot(2,2,4)
                plot(tex,T0 + (rho0 - yex(:,3))/alphatilde,plotopt{:})
                title('Temperature')
                grid on;
                sgtitle('Reference solution',txtopt{:})
            else
                if ~pde_plot
                    hdl1 = plot(tex,yex,plotopt{:});
                    ax = gca;
                    ax.XScale = xscale;
                    ax.YScale = yscale;
                end
            end



            if ~pde_plot && ~strcmp(test,'strat')
                xlim([tstart tend])
                xlabel('t',txtopt{:})
                if ~strcmp(test,'reducednucreact')
                    title('Reference solution',txtopt{:})
                    legend(lgnd(1:length(y0)+size(E,1)),'location',loc)
                    set(gca,txtopt{:})
                    if exist('yrange','var')
                        ylim(yrange);
                    end
                end
            end
            if ~pde_plot && strcmp(task,'plotsave')
                print('-depsc2',['./save/' test '_ref.eps'])
            end
        end

        if ~exist('dtfunc','var')
            dtfunc = @(dt) dt;
        end

        for i = 1:length(solvers)
            meth = solvers{i};
            try
                if ~iscell(meth)
                    str = meth;
                    fh = str2func(str);
                    [t,Y] = fh(PD, y0, dt0, tint, dtfunc);
                else
                    switch length(meth)
                        case 2
                            str = meth{1};
                            theta = meth{2};

                            fh = str2func(str);
                            [t,Y] = fh(PD, y0, dt0, tint, theta, dtfunc,Beta1, Beta2, Beta3, Alpha2, kappa2, para, scaleode);
                            if strcmp(str,'SSP33adap')
                                str = sprintf('%s',str);
                            elseif strcmp(str,'ode23wrapper')
                                str = 'ODE23';
                            elseif strcmp(str,'ode45wrapper')
                                str = 'ODE45';
                            elseif strcmp(str,'ode15swrapper')
                                str = 'ODE15s';
                            elseif strcmp(str,'BE')
                                str = sprintf('%s (%s)',str,upper(theta));
                            elseif strcmp(str,'SDIRK2')
                                str = sprintf('%s (%s)',str,upper(theta));
                            else
                                str = sprintf('%s(%.3f)',str,theta);
                            end
                        case 3
                            str = meth{1};
                            theta = meth{2};
                            omega = meth{3};
                            fh = str2func(str);
                            tic
                            if isnan(para)
                                [t,Y,cnt_rej,lengtht] = fh(PD, y0, dt0, tint, theta, omega, dtfunc, Beta1, Beta2, Beta3, Alpha2, kappa2, para, scaleode);
                            else
                                [t,Y,cnt_rej,lengtht] = fh(PD, y0, dt0, tint, theta, omega, dtfunc, Beta1, Beta2, Beta3, Alpha2, kappa2, para);
                            end
                            times = toc;
                            if strcmp(str,'MPRK43IIadap') || strcmp(str,'MPRK22adap') || strcmp(str,'mBBKS2adap') || strcmp(str,'ROS2')
                                str = sprintf('%s(%.3f)',str,theta);
                            elseif strcmp(str,'ode23swrapper')
                                str = 'ODE23s';
                            elseif strcmp(str,'ode45wrapper')
                                str = 'ODE45';
                            elseif strcmp(str,'ode15swrapper')
                                str = 'ODE15s';
                            else
                                str = sprintf('%s(%.3f, %.3f)',str,theta,omega);
                            end
                        case 4
                            str = meth{1};
                            para1 = meth{2};
                            para2 = meth{3};
                            para3 = meth{4};
                            fh = str2func(str);
                            tic
                            if isnan(para)
                                [t,Y,cnt_rej,lengtht] = fh(PD, y0, dt0, tint, para1, para2, para3, Beta1, Beta2, Beta3, Alpha2, kappa2, para, scaleode);
                            else
                                [t,Y,cnt_rej,lengtht] = fh(PD, y0, dt0, tint, para1, para2, para3, Beta1, Beta2, Beta3, Alpha2, kappa2, para);
                            end
                            times = toc;
                            str = sprintf('%s(%.3f,%.3f)',str,para1,para2);
                    end
                end

                if strcmp(task,'plot')
                    figure
                    loglog(t,sum(Y)./sum(y0))
                    title('Conservativity')
                    hold off
                    if strcmp(test, 'strat')
                        figure
                        ttemp = t/3600 + 12;
                        semilogy(ttemp,(Y(5,:) + 0.5 * Y(6,:))./(y0(5) + 0.5*y0(6)))
                        title('NO + 0.5 NO2')
                        hold off
                    end
                end
                %         if strcmp(test,'strat')
                %             figure
                %             loglog(t,(2*Y(:,5)+Y(:,6))./(2*y0(5)+y0(6)))
                %             title('Second Invariant')
                %             hold off
                %         end
                %sumY = sum(Y);
                Y=Y/scaleode;
                if ~isempty(E)
                    sumY = E*Y;
                else
                    sumY = [];
                end

                % if ~exist('cnt_rej', 'var')
                %      cnt_rej = 0;
                %  end


                % fprintf('%.3f_%.3f_%.3f\n',Beta1,Beta2,Beta3)
                Y = [Y; sumY];
                if exist('visfac','var')
                    Y = Y.*visfac';
                end


                if ~strcmp(task,'data')
                    if  ~pde_plot
                        figure
                        %hdl1 = plot(tex,yex,'--',tex,sumyex,'--',plotopt{:});
                        if strcmp(test,'hires')
                            hdl1 = semilogx(tex,yex,'--',plotopt{:});
                            hold on
                        elseif strcmp(test,'strat')
                            linespec = 'o-';
                            %t = t/3600 + 12;
                            %t = t+12;
                            %yex = (diag([9.906e1; 6.624e8; 5.326e11; 1.697e16; 4.000e16; 1.093e9])*yex')';
                            subplot(3,2,5)
                            plot(tex_temp,yex(:,1),plotopt{:})
                            title('O^{1D}')
                            ax = gca;
                            ax.ColorOrderIndex = 1;
                            grid on;
                            hold on
                            plot(ttemp,(Y(1,:))',linespec,plotopt{:})
                            hold off
                            xlim([tstart tend])
                            xlabel('t',txtopt{:})

                            subplot(3,2,2)
                            plot(tex_temp,yex(:,2),plotopt{:})
                            title('O')
                            ax = gca;
                            ax.ColorOrderIndex = 1;
                            grid on;
                            hold on
                            plot(ttemp,(Y(2,:))',linespec,plotopt{:})
                            hold off
                            xlim([tstart tend])
                            xlabel('t',txtopt{:})

                            subplot(3,2,3)
                            plot(tex_temp,yex(:,3),plotopt{:})
                            title('O_3')
                            ax = gca;
                            ax.ColorOrderIndex = 1;
                            grid on;
                            hold on
                            plot(ttemp,(Y(3,:))',linespec,plotopt{:})
                            hold off
                            xlim([tstart tend])
                            xlabel('t',txtopt{:})

                            subplot(3,2,1)
                            plot(tex_temp,yex(:,4) - 2*1.697e+16,plotopt{:})
                            title('O_2 - 2*1.697e+16')
                            ax = gca;
                            ax.ColorOrderIndex = 1;
                            grid on;
                            hold on
                            plot(ttemp,(Y(4,:))' - 2*1.697e+16,linespec,plotopt{:})
                            hold off
                            xlim([tstart tend])
                            xlabel('t',txtopt{:})

                            subplot(3,2,4)
                            plot(tex_temp,yex(:,5),plotopt{:})
                            title('NO')
                            ax = gca;
                            ax.ColorOrderIndex = 1;
                            grid on;
                            hold on
                            plot(ttemp,(Y(5,:))',linespec,plotopt{:})
                            hold off
                            xlim([tstart tend])
                            xlabel('t',txtopt{:})

                            subplot(3,2,6)
                            plot(tex_temp,yex(:,6),plotopt{:})
                            title('NO_2')
                            ax = gca;
                            ax.ColorOrderIndex = 1;
                            grid on;
                            hold on
                            plot(ttemp,(Y(6,:))',linespec,plotopt{:})
                            hold off
                            xlim([tstart tend])
                            xlabel('t',txtopt{:})
                            sgtitle(str,txtopt{:})
                        elseif strcmp(test,'reducednucreact')
                            linespec = 'o-';
                            subplot(2,2,1)
                            plot(tex,yex(:,1),'--',plotopt{:})
                            title('Neutron Density')
                            ax = gca;
                            ax.ColorOrderIndex = 1;
                            grid on;
                            hold on
                            plot(t,(Y(1,:))',linespec,plotopt{:})
                            hold off

                            subplot(2,2,2)
                            plot(tex,yex(:,2),'--',plotopt{:})
                            title('Precursor Density')
                            ax = gca;
                            ax.ColorOrderIndex = 1;
                            grid on;
                            hold on
                            plot(t,(Y(2,:))',linespec,plotopt{:})
                            hold off

                            subplot(2,2,3)
                            plot(tex,yex(:,3),'--',plotopt{:})
                            title('Reactivity')
                            ax = gca;
                            ax.ColorOrderIndex = 1;
                            grid on;
                            hold on
                            plot(t,(Y(3,:))',linespec,plotopt{:})
                            hold off

                            subplot(2,2,4)
                            plot(tex,T0 + (rho0 - yex(:,3))/alphatilde,'--',plotopt{:})
                            title('Temperature')
                            ax = gca;
                            ax.ColorOrderIndex = 1;
                            grid on;
                            hold on
                            plot(t,T0 + (rho0 - (Y(3,:))')/alphatilde,linespec,plotopt{:})
                            hold off
                            sgtitle(str,txtopt{:})


                        else
                            hold off
                            hdl1 = plot(tex,yex,'--',plotopt{:});
                            hold on
                        end


                        %hdl2 = plot(t,Y',linespec,t,sumY,linespec,plotopt{:});
                        ax = gca;
                        ax.ColorOrderIndex = 1;

                        if ~exist('linespec','var')
                            linespec = 'o-';
                        end
                        if strcmp(test,'hires')
                            hdl2 = semilogx(t,Y',linespec,plotopt{:});
                        else

                            if ~strcmp(test,'reducednucreact') && ~strcmp(test,'strat')
                                hdl2 = plot(t,Y',linespec,plotopt{:});
                            end
                        end

                        xlim([tstart tend])
                        xlabel('t',txtopt{:})
                        if ~strcmp(test,'reducednucreact') && ~strcmp(test,'strat')
                            title(str,txtopt{:})
                            l = legend(lgnd,'location',loc);
                            l.Visible = 'on';
                            if exist('auxlines','var')
                                ax = gca;
                                ax.ColorOrderIndex = 1;
                                plot([tstart tend],auxlines([1 1],:),'-');
                            end
                            ax.XScale = xscale;
                            ax.YScale = yscale;
                            hold off
                            if exist('yrange','var')
                                ylim(yrange);
                            end
                            set(gca,txtopt{:})
                        end
                        if strcmp(task,'plotsave')
                            print('-depsc2',['./save/' test '_' str '.eps'])
                        end

                        if exist('exsol', 'var')
                            figure
                            if ~isempty(E)
                                Y(end,:)=[];
                            end
                            hdl3 = semilogy(t,abs(Y'-exsol(t'))./abs(exsol(t')),'-',plotopt{:});
                            title('Relative Error',txtopt{:})
                            l2 = legend(lgnd2,'location',loc);
                            l2.Visible = 'on';
                        end
                        if strcmp(task,'plotsave')
                            print('-depsc2',['./save/' test '_err_' str '.eps'])
                        end
                    else
                        for idx = 1: size(inds,1)
                            if exist('exsol','var')
                                lgnd{1}='exact';
                            else
                                lgnd{1}='ref';
                            end
                            lgnd{2} = 'num';
                            loc = 'ne';
                            plot(x,yex(end,inds(idx,:)),Markers{1},plotopt{:},'MarkerIndices',1:20*Markerint:length(x), Markeropt{:})
                            xlim([xstart xend])
                            xlabel('x',txtopt{:},'Interpreter','latex')
                            if exist('ylabels','var')
                                ylabel(ylabels{idx},txtopt{:},'Interpreter','latex')
                            end
                            title('Reference solution',txtopt{:})
                            if strcmp(task,'plotsave')
                                print('-depsc2',['./save/' test '_ref.eps'])
                            end
                            figure
                            plot(x,yex(end,inds(idx,:)),Markersref{1},plotopt{:},'MarkerIndices',1:20*Markerint:length(x), Markeropt{:})
                            hold on
                            plot(x,Y(inds(idx,:),end),Markers{1},plotopt{:},'MarkerIndices',1:Markerint:length(x),'Markersize',4);

                            hold off
                            xlim([xstart xend])
                            xlabel('x',txtopt{:},'Interpreter','latex')
                            if exist('ylabels','var')
                                ylabel(ylabels{idx},txtopt{:},'Interpreter','latex')
                            end
                            print('-depsc2',['./save/' test '_' str '.eps'])
                            if one_dim_prob
                                figure
                                [T,X]=meshgrid(t,x);
                                [Tex,Xex] = meshgrid(tex,x);
                                if ~strcmp(task,'plotsave')
                                    subplot(1, 2, 2)
                                    meshc(T,X,Y);
                                    subplot(1, 2, 1)
                                end
                                meshc(Tex,Xex,yex');
                                %title('Exact solution', 'FontSize', 12', 'FontWeight', 'bold')
                                xlabel('t', 'FontSize', 12')
                                ylabel('x', 'FontSize', 12')
                                if strcmp(task,'plotsave')
                                    print('-depsc2',['./save/' test '_ref3d.eps'])
                                end
                            end

                        end
                    end
                end

                %         if ~and(strcmp(test,'strat'), strcmp(task,'plot'))
                if ~exist('exsol', 'var')
                    if ~isOctave
                        options.AbsTol = 1e-13;
                        options.RelTol = 1e-13;
                        if isnan(para)
                            [~,~,f] = PD(0,y0/scaleode,scaleode);
                        else
                            [~,~,f] = PD(0,y0,para);
                        end
                        %                     if strcmp(test,'mseir')
                        %                         %[tex,yex] = ode23s(f,tint,y0,options);
                        %                         [t,yex] = ode15s(f,t,y0,options);
                        %                     elseif strcmp(test,'bertolazzi')
                        %                         opt.RelTol = 1e-5;
                        %                         opt.AbsTol = opt.RelTol;
                        %                         opt.NonNegtaive = [0 1 1];
                        %                         [t,yex] = ode23s(f,t,y0,opt);
                        %                     else
                        %options.NonNegative = 1:length(y0);
                        %options.NonNegative = [];

                        if any(diff(t)<=0)
                            [tex,yy] = ode15s(f,tint,y0,options);
                            yy = yy/scaleode;
                            y = zeros(size(Y'));
                            for i = 1:length(t)
                                tarray = repmat(t(i), [length(tex) 1]);
                                [~,closestIndex] = min(abs(tarray-tex));
                                y(i,:) = yy(closestIndex,:);
                            end
                            %varargout{3} = NaN;
                        else
                            if length(t) > 1
                                [tex,y] = ode15s(f,t,y0,options);
                                if length(t) == 2
                                    y = [y(1,:); y(end,:)];
                                end
                                y = y/scaleode;
                            end
                        end
                        if exist('visfac','var')
                            Y= Y./visfac';
                        end
                        %err2exsol = mean( vecnorm(y-Y')./vecnorm(y) )
                        if length(t) > 1
                            if pde_plot
                                err2exsol = L2err_rel2d(t,x,Y,y); %L2err_rel(x,(Y(:,end))',(y(end,:))');% norm(Y(:,end) - (y(end,:))',Inf) /  norm(y(end,:),Inf);
                            else
                                err2exsol =  L2err_rel(t,Y,y);
                            end
                        else
                            err2exsol = NaN;
                        end
                        %                         err2lininv = NaN;
                        if strcmp(task, 'data') && strcmp(test, 'strat') %% relative error to second linear invariant
                            % err2lininv =  L2err_rel(t,(Y(5,:) + 0.5 * Y(6,:))./(y0(5) + 0.5*y0(6)), (y(:,5) + 0.5 * y(:,6))./(y0(5) + 0.5*y0(6)));
                            %err2lininv =  L2err_rel(t,Y(5,:) + 0.5 * Y(6,:), y(:,5) + 0.5 * y(:,6));
                        end

                        %                         [tex,y] = ode15s(f,t,y0,options);
                        %                         if length(t) == 2
                        %                             y = [y(1,:); y(end,:)];
                        %                         end
                        %[tex,yex] = ode23s(f,tint,y0,options);
                        % end
                    else
                        warning('Octave is currently not supported.')
                        return
                    end
                else
                    %err2exsol = mean( vecnorm(exsol(t')-Y')./vecnorm(exsol(t')) );
                    if length(t) > 1
                        if pde_plot
                            err2exsol = L2err_rel2d(t,x,Y,exsol(t',x));%L2err_rel(x,(Y(:,end))',exsol(tend,x')) %norm(Y(:,end)-exsol(tend,x'),Inf) / norm(exsol(tend,x'));
                        else
                            err2exsol = L2err_rel(t,Y,exsol(t'));
                        end
                    else
                        err2exsol = NaN;
                    end
                    %                     err2lininv = NaN;
                end
                %         else
                %             err2exsol = -inf; %%% dummy output if strat is plotted
                %         end

                if nargout == 3
                    if isnan(lengtht)
                        varargout{1} = NaN;
                        %disp(lengtht)
                    else
                        varargout{1} = length(t) - 1;
                        if strcmp(solverstr,'ROS2')
                            varargout{1} = varargout{1} + length(y0) + 1; % will be multiplied by 2 in workprecision.m to obtain the rhs evaluations
                        end
                    end
                    varargout{2} = cnt_rej;
                    varargout{3} = [err2exsol, times];

                elseif strcmp(task,'plot')
                    fprintf('%s\n',str);
                    fprintf('%i successful steps\n',length(t)-1);
                    fprintf('%i failed attempts\n\n',cnt_rej);
                    % if pde_plot
                    fprintf('%1.2e relative L2 error\n\n',err2exsol);
                    % else
                    %     fprintf('%1.2e relative L2 error\n\n',err2exsol);
                    % end
                end



            catch ME
                fprintf('\n')
                fprintf('Method: %s\n',str)
                warning(ME.message)
                fprintf('\n')
            end
        end




        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Compute error %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'error'
        %%% Time step sizes
        dtvec = (tend - tstart)./2.^(4:12);
        dtvec = (tend - tstart)./2.^(4:10);
        %dtvec = 1./2.^(0:5);
        %dtvec = [2 1 .5 .25 .125];
        %dtvec = (tend - tstart)./2.^(10:17);
        %dtvec = 0.1./2.^(0:5);

        %%% Plot options
        plotopt = {'LineWidth',2};
        txtopt = {'FontSize',14};

        if ~exist('E','var')
            E = ones(1,length(y0));
        end

        tint = [tstart tend];

        legc = cell(1,length(solvers));
        err = inf*ones(length(solvers),length(dtvec));
        cons = err;
        for k = 1:length(dtvec)
            dt = dtvec(k);

            for i = 1:length(solvers)
                meth = solvers{i};
                if ~iscell(meth)
                    str = meth;
                    fh = str2func(str);
                    [t,C] = fh(PD, y0, dt, tint);
                else
                    switch length(meth)
                        case 2
                            str = meth{1};
                            theta = meth{2};
                            fh = str2func(str);
                            [t,C] = fh(PD, y0, dt, tint, theta);

                            str = sprintf('%s(%.2f)',str,theta);
                        case 3
                            str = meth{1};
                            theta = meth{2};
                            omega = meth{3};
                            fh = str2func(str);
                            [t,C] = fh(PD, y0, dt, tint, theta, omega);

                            str = sprintf('%s(%.2f, %.2f)',str,theta,omega);
                    end
                end

                opt.RelTol = 1e-10;
                opt.AbsTol = 1e-12;
                %opt.NonNegative = 1:length(y0);
                [~,~,f] = PD(0,y0);
                [~,cex] = ode15s(f, t, y0, opt);

                if isreal(C)
                    err(i,k) = L2err_rel(t,C,cex);
                    %err(i,k) = comp_err(t,C,cex,2,3,false);
                    %err(i,k) = norm(C(:,end) - cex(end,:)');
                else
                    err(i,k) = NaN;
                end

                %N = size(C,2);
                %cons(i,k) = 1/sum(y0)*sqrt(1/N* sum((sum(y0)-sum(C)).^2));
                if ~isempty(E)
                    cons(i,k) = L2err_rel(t,E*C,(E*y0(:,ones(1,size(C,2))))');
                end
                legc{i} = str;
            end

        end

        %%

        loglog(dtvec,err(1,:),'o-',plotopt{:})
        hold on
        for i = 2:length(solvers)
            loglog(dtvec,err(i,:),'o-',plotopt{:})
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %     ax = gca;
        %     ax.ColorOrderIndex = 1;
        %     Caux = err(:,1)./(dtvec(1).^(2:6))';
        %     for i = 1:length(solvers)
        %       loglog(dtvec,Caux(i)*dtvec.^(i+1),'--',plotopt{:})
        %     end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        hold off
        l = legend(legc,'location','se');
        l.ItemHitFcn = @legend_hitcallback;
        xlabel('\Delta t',txtopt{:})
        title('Click items in the legend to show or hide the associated chart')
        set(gca,txtopt{:})
        %ylim([-Inf 1])

        figure(3)
        semilogx(dtvec,err(1,:),'o-',plotopt{:})
        hold on
        for i = 2:length(solvers)
            semilogx(dtvec,err(i,:),'o-',plotopt{:})
        end
        hold off
        l = legend(legc,'location','se');
        l.ItemHitFcn = @legend_hitcallback;
        xlabel('\Delta t',txtopt{:})
        title('Click items in the legend to show or hide the associated chart')
        set(gca,txtopt{:})
        %ylim([-Inf 1])

        figure(2)
        loglog(dtvec,cons(1,:),plotopt{:})
        hold on
        for i = 2:length(solvers)
            loglog(dtvec,cons(i,:),plotopt{:})
        end
        hold off
        l = legend(legc,'location','ne');
        l.ItemHitFcn = @legend_hitcallback;
        xlabel('\Delta t',txtopt{:})
        title('Click items in the legend to show or hide the associated chart')
        set(gca,txtopt{:})
        figure(1)


        P = [];
        for i = 1:length(solvers)
            p = log(err(i,2:end)./err(i,1:end-1))./log(dtvec(2:end)./dtvec(1:end-1));
            P = [P p'];
        end
        P

        %%% Efficiency %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'efficiency',
        %%% Time step sizes
        dtvec = (tend - tstart)./2.^(4:12);

        %%% Plot options
        plotopt = {'LineWidth',2};
        txtopt = {'FontSize',14};

        legc = cell(1,length(solvers));
        err = inf*ones(length(solvers),length(dtvec));
        CPUtime = zeros(length(solvers),length(dtvec));

        jmax = 10;
        for j = 1:jmax
            for k = 1:length(dtvec)
                dt = dtvec(k);

                opt.RelTol = 1e-10;
                opt.AbsTol = 1e-10;
                [~,~,f] = PD(0,y0);
                tint = [tstart tend];
                [~,cex] = ode15s(f, tstart:dt:tend, y0, opt);

                for i = 1:length(solvers)
                    meth = solvers{i};
                    if ~iscell(meth)
                        str = meth;
                        fh = str2func(str);
                        time0 = tic;
                        [t,C] = fh(PD, y0, dt, tint);
                        telapsed = toc(time0);
                    else
                        switch length(meth)
                            case 2,
                                str = meth{1};
                                theta = meth{2};
                                fh = str2func(str);
                                time0 = tic;
                                [t,C] = fh(PD, y0, dt, tint, theta);
                                telapsed = toc(time0);
                                str = sprintf('%s(%.2f)',str,theta);
                            case 3,
                                str = meth{1};
                                theta = meth{2};
                                omega = meth{3};
                                fh = str2func(str);
                                time0 = tic;
                                [t,C] = fh(PD, y0, dt, tint, theta, omega);
                                telapsed = toc(time0);
                                str = sprintf('%s(%.2f, %.2f)',str,theta,omega);
                        end
                    end

                    N = length(C(1,:)) - 1;
                    %err(i,k) = sqrt(1/N*sum((cex(:,1)' - C(1,:)).^2))/(1/N*sum(cex(:,1)));
                    err(i,k) = L2err_rel(t,C,cex);
                    CPUtime(i,k) = CPUtime(i,k) + telapsed;
                    legc{i} = str;
                end
            end
        end
        CPUtime = CPUtime/jmax;

        loglog(CPUtime(1,:),err(1,:),'o-',plotopt{:})
        hold on
        for i = 2:length(solvers)
            loglog(CPUtime(i,:),err(i,:),'o-',plotopt{:})
        end
        hold off
        l = legend(legc,'location','sw');
        l.ItemHitFcn = @legend_hitcallback;
        xlabel('CPU Time',txtopt{:})
        ylabel('Error',txtopt{:})
        title('Click items in the legend to show or hide the associated chart')
        set(gca,txtopt{:})


        %%% Adaptivity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'adaptivity'
        tint = [tstart tend];

        %%% Plot options
        plotopt = {'LineWidth',2};
        txtopt = {'FontSize',14};

        m = 2:4;
        for k = 1:length(m)
            RelTol = 10^(-m(k));
            AbsTol = 10^(-m(k)-2);
            dt0 = RelTol;

            Tol = [RelTol AbsTol];
            stiffsolvers = {{'MPRK43Iadap',1,.5,Tol},{'MPRK43Iadap',.5,.75,Tol}, ... %{'MPRK43IIadap',.5,Tol},
                {'ode23swrapper',Tol},{'ode15swrapper',Tol}};
            nonstiffsolvers = {{'SSP33adap',Tol},{'ode45wrapper',Tol}};
            solvers = [stiffsolvers nonstiffsolvers];
            solvers = stiffsolvers;

            for i = 1:length(solvers)
                meth = solvers{i};
                try
                    if ~iscell(meth)
                        str = meth;
                        fh = str2func(str);
                        t0 = tic;
                        [t,Y] = fh(PD, y0, dt0, tint);
                        t1 = toc(t0);
                    else
                        switch length(meth)
                            case 2,
                                str = meth{1};
                                theta = meth{2};
                                fh = str2func(str);
                                t0 = tic;
                                [t,Y] = fh(PD, y0, dt0, tint, theta);
                                t1 = toc(t0);
                                if strcmp(str,'SSP33adap')
                                    str = sprintf('%s',str);
                                elseif strcmp(str,'ode23wrapper')
                                    str = 'ODE23';
                                elseif strcmp(str,'ode23swrapper')
                                    str = 'ODE23s';
                                else
                                    str = sprintf('%s(%.3f)',str,theta);
                                end
                            case 3,
                                str = meth{1};
                                theta = meth{2};
                                omega = meth{3};
                                fh = str2func(str);
                                t0 = tic;
                                [t,Y] = fh(PD, y0, dt0, tint, theta, omega);
                                t1 = toc(t0);
                                if ~strcmp(str,'MPRK43IIadap')
                                    str = sprintf('%s(%.3f, %.3f)',str,theta,omega);
                                else
                                    str = sprintf('%s(%.3f)',str,theta);
                                end
                            case 4,
                                str = meth{1};
                                para1 = meth{2};
                                para2 = meth{3};
                                para3 = meth{4};
                                fh = str2func(str);
                                t0 = tic;
                                [t,Y] = fh(PD, y0, dt0, tint, para1, para2, para3);
                                t1 = toc(t0);
                                str = sprintf('%s(%.3f,%.3f)',str,para1,para2);
                        end
                    end

                    time(i,k) = t1;

                    [~,~,f] = PD(0,y0);
                    opt.RelTol = 1e-12;
                    opt.AbsTol = 1e-12;
                    if ~strcmp(test,'mseir')
                        [~,Yex] = ode15s(f, t, y0, opt);
                    else
                        [~,Yex] = ode23s(f, t, y0, opt);
                    end
                    err(i,k) = L2err_rel(t,Y,Yex);

                    lgnd{i} = str;

                catch ME
                    fprintf('\n')
                    warning(ME.message)
                    fprintf('\n')
                end
            end
        end
        %%
        for i = 1:length(solvers)
            semilogx(time(i,:),log10(err(i,:)),'o-',plotopt{:})
            hold on
        end
        hold off
        xlabel('CPU time')
        ylabel('log_{10}(error)')
        l = legend(lgnd,'location','ne');
        l.ItemHitFcn = @legend_hitcallback;
        set(gca,txtopt{:})
        grid
    otherwise
        error('Task %s is unknown.\n',task)
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %            _
% % Triangle  |/
%
% x1 = 1e-2;
% x2 = 1e-1;
% y2 = 1;
% p = 1;
% C = exp(log(y2)-p*log(x2));
% plot([x1 x2 x2 x1 x1],[y2 y2 C*x2^p C*x1^p y2],'Color',[.25 .25 .25])
% text(sqrt(x1*x2),y2,'1','VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',12)
% text(x1,sqrt(y2*C*x1^p),num2str(p),'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',12)
% %
% % Triangle   /|
% %
% x1 = 2*1e-2;
% x2 = 2*1e-1;
% y1 = 1e-7;
% p = 3;
% C = exp(log(y1)-p*log(x1));
% plot([x2 x1 x1 x2 x2],[C*x2^p C*x1^p C*x1^p C*x1^p C*x2^p],'Color',[.25 .25 .25])
% text(sqrt(x1*x2),C*x1^p,'1','VerticalAlignment','top','HorizontalAlignment','center','FontSize',12)
% text(x2,sqrt(C*x2^p*C*x1^p),num2str(p),'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',12)
end
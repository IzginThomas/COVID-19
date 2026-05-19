%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ODE Solver Test Suite                                               %%%
%%%                                                                     %%%
%%% Author: Stefan Kopecz (kopecz@mathematik.uni-kassel.de)             %%%
%%% Date: 10/26/2017                                                    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Robertson problem                                                   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Right hand side
PD = @PD_robertson;

% Time interval
tstart = 0;
tend = 1e8;
deltta = 1e-1;
% Initial values
y0 = [deltta; deltta; 1-2*deltta];
y0 = [1; 0;0];

% Time step size
dt0 = 1e-6;
%dt0 = 1e-3;

% Adaptive time stepping
dtfunc = @(dt) dt*2;
%dtfunc = @(dt) dt;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot function
%plot = @semilogx;

%Linear invariant
E =[];

% Visualization amplification factors
visfac = [1 1e4 1];

% Plot range
yrange =  [-.1 1.1];

xscale = 'log';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tend = 1e2;
% dt0 = 5;
% y0 = [1; 0.37*1e-4; eps]

% tend = 1e5;
% dt0 = 3e+03;
% y0 = [0.9; 2e-5;0.1]; 




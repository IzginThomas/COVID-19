%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ODE Solver Test Suite                                               %%%
%%%                                                                     %%%
%%% Author: Stefan Kopecz (kopecz@mathematik.uni-kassel.de)             %%%
%%% Date: 10/26/2017                                                    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Brusselator problem                                                 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% Right hand side
PD = @PD_brusselator;

% Time interval
tstart = 0;
tend = 10;

% Initial values
y0 = [10; 10; realmin; realmin; 0.1; 0.1];
%y0 = [10; 10; 0; 0; 0.1; 0.1];

% Time step size
dt0 = .1;

% Visualization amplification factors
%visfac = [1 1 1 1 1 1 0.5];
%visfac = [1 1 1 1 1 1 1 1 0.5];

% Plot range
%yrange =  [-1 11];

% Compound composition matrix
%E = [1 0 0 1 1 1; 0 1 1 0 0 0; 1 1 1 1 1 1];
E = [];
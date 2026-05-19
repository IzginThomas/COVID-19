

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ODE Solver Test Suite                                               %%%
%%%                                                                     %%%
%%% Author: Stefan Kopecz (kopecz@mathematik.uni-kassel.de)             %%%
%%% Date: 09/02/2018                                                    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% NPZD problem                                                        %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Right hand side
PD = @PD_npzd;

% Time interval
tstart = 0;
tend = 5;
%tend = 50;

% Initial values
y0 = [8; 2; 1; 4];
%y0 = [7;eps;eps;1];

% Time step size
dt0 = 1;
%dt0 = 0.1;
%dt0 = 5;

% Plot function
%plot = @semilogx;

% Plot range
yrange =  [-1 11];

% Line specification
%linespec = 'o-';
E = [];

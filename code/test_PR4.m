%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ODE Solver Test Suite                                               %%%
%%%                                                                     %%%
%%% Author: Thomas Izgin (izgin@mathematik.uni-kassel.de)               %%%
%%% Date: 09/29/2022                                                    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Nonlinear problem 4D                                                %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% Right hand side
PD = @PD_PR4;

% Time interval
tstart = 0; % in seconds
tend = 20*pi;

% Initial values
y0 = [2; 2; 1; 1]; % =exsol(0)


% Time step size
dt0 = 1;

% Visualization amplification factors
%visfac = [1 1 1 1];

% Plot range
yrange =  [-1 6.1];

% Compound composition matrix
%E = [1 0 0 1 1 1; 0 1 1 0 0 0; 1 1 1 1 1 1];
%E = [1 1 1 1];
E=[];

% exact solution
oscifac =@(t) 0.1*cos(0.5.*t);
oscidiff =@(t) -0.05*sin(.5.*t); % derivative of oscifac(t)
proddiff =@(t) oscidiff(t).*t + oscifac(t); % derivative of oscifac(t)*t
k1 = 0.3; %k1<=0.5 for positivity

%A = [-1 1-para para 0; para -1 0 1-para; 1-para 0 -1 para;0 para 1-para -1];

exsol=@(t) [2+k1.*sin(oscifac(t).*t) 2+sin(oscifac(t).*t) 1-sin(oscifac(t).*t) 1-k1.*sin(oscifac(t).*t)];
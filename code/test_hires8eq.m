% Right hand side
PD = @PD_hires8eq;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Time interval
 tstart = 0;
 tend = 321.8122;
 %tend =5;

 % (Initial) Time step size
 dt0 = .5e-3;

 % Initial values
 y0=[1;0;0;0;0;0;0;0.0057];
 % y0(1) = 1;
 % y0(8) = 0.0057;
 % y0 = y0' + realmin*ones(8,1);
% y0(9)=0;
% y0=y0' + realmin*ones(9,1);
 

% Visualization amplification factors
visfac = [1 5 5 1 5 1 1e2 1e2];

xscale = 'log';
% yscale = 'log';


% Compound composition matrix
 E=[];
 

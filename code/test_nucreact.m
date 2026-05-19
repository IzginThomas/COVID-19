% Right hand side
PD = @PD_nucreact;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Time interval
 tstart = 0;
 tend = 150;


 % (Initial) Time step size
 dt0 = .5e-2;

 % Initial values

lambda1 = 0.0124;
lambda2 = 0.0305;
lambda3 = 0.111;
lambda4 = 0.301;
lambda5 = 1.13;
lambda6 = 3.0;
lambdavector = [lambda1;lambda2;lambda3;lambda4;lambda5;lambda6];

beta1 = 0.00021;
beta2 = 0.00141;
beta3 = 0.00127;
beta4 = 0.00255;
beta5 = 0.00074;
beta6 = 0.00027;
betavector = [beta1;beta2;beta3;beta4;beta5;beta6];
beta = sum(betavector);

LAMBDA = 0.0005;
Phi0 = 10;
cvector = Phi0/LAMBDA*(betavector./lambdavector);
T0 = 300;
rho0 = 0.5*beta;

y0 = [Phi0;cvector;rho0]; % length(y0)=1 + 6 + 1 = 8

% Visualization amplification factors
%visfac = [1 5 5 1 5 1 1e2 1e2 1];

%yscale = 'log';
% Compound composition matrix
 E=[];
 

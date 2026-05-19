

PD = @PD_FokkerPlanck;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Time interval
tstart = 0;
tend = 10;

% Spatial domain
Nx = 200;
xstart = -1;
xend = 1;
xspan = xend - xstart;

% x= -5.0:0.0525:5.4475;
% dx =0.0525;
% xstart = -5;
% xend = 5.4475;
% Nx = length(x);
bc = 'noflux';
% dx = (xend-xstart)/Nx;
% x = linspace(xstart,xend-dx,Nx);
x_interfaces =  linspace(xstart,xend,Nx + 1);
dx = x_interfaces(2) - x_interfaces(1);
x = x_interfaces(1:end-1) + dx/2;
para = 0.2; % sigma^2



% (Initial) Time step size
CFL = 0.8;
dt0 = CFL*dx; % CFL

% Initial values

y01_func = @(x) exp(-30*(x+0.5).^2)+exp(-30*(x-0.5).^2);
beta = 1 / integral(y01_func,xstart,xend,'RelTol',1e-14,'AbsTol',1e-14);
y0_func = @(x) beta * y01_func(x);
y0 = y0_func(x');
% exsol = @(t,x) 0.*y0;

% y0x = @(x) x .* y0_func(x);
% uu = integral(y0x,xstart,xend,'RelTol',1e-14,'AbsTol',1e-14);
% y_stat_1 = @(x) 1./(1-x.^2).^2 .* ((1+x)./(1-x)).^(uu ./ (2.*para)) .*exp(-(1-uu.*x)./(para .* (1-x.^2)));
% K = 1./integral(y_stat_1,xstart,xend,'RelTol',1e-14,'AbsTol',1e-14);
% y_stat = @(x) K .* y_stat_1(x);
% exsol = @(x) y_stat(x);

% Compound composition matrix
E=[];
inds = 1:Nx;
pde_plot = true;
one_dim_prob = true;
% % Relaxation
entropy_flag = true;
if entropy_flag
    % Entropy
    eta = str2func(['@(Y) sum(Y .* log(Y) - Y)*' num2str(dx)]);
    % Entropy prime
    eta_prime =str2func(['@(y) (log(y)) *' num2str(dx)]);

    entropy_cons = false;

end
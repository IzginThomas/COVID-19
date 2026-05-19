


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ODE Solver Test Suite                                               %%%
%%%                                                                     %%%
%%% Author: Stefan Kopecz (kopecz@mathematik.uni-kassel.de)             %%%
%%% Date: 04/20/2018                                                    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dt,flag,eps1] = compute_dt(y,yn,sigma,dt,RTOL,ATOL,errctrl,p,varargin)

if nargin >8
    epsminus = varargin{1};
    eps0 = varargin{2};
    Beta1 = varargin{3};
    Beta2 = varargin{4};
    Beta3 = varargin{5};
    Alpha2 = varargin{6};
    kappa2 = varargin{7};
    dtold = varargin{8};
end

flag = false;

fmin = 0.1;
fmax = 10;
fsafety = 0.81;
%fsafety = 0.3;
eps1 = [];
kappa = @(x) 1+kappa2*atan((x-1)/kappa2); % Continuously bounding dt_factor

%Beta1 = 0.7;
%Beta2 = -0.4;
%Beta3 = 0;



switch errctrl
  case 'EPUS'
    %r = norm((y - sigma)./y,inf);
    r = norm(y - sigma);
    TOL = RTOL*norm(y) + ATOL;
    if r > TOL
      dt = dt*max([fmin, fsafety*(TOL/r)^(1/(p+1))]);
      flag = true;
    else
      dt = dt*min([fmax, fsafety*(TOL/r)^(1/(p+1))]);
    end
  case 'SWP'
    sk = ATOL + max(abs(y),abs(sigma))*RTOL;
    err = norm((y - sigma)./sk)/sqrt(length(sigma));
    dt = dt*min(fmax,max(fmin, fsafety*(1/err)^(1/(p+1))));
    if err > 1
      flag = true;
    end
  case 'PID'
    sk = ATOL + max(abs(y),abs(sigma))*RTOL;
    w = norm((y - sigma)./sk)./sqrt(length(y)); 
    %w = max(abs(y - sigma)./sk);
    w = max(w,eps); % We want w>=eps=2.2204e-16
    eps1 = 1 / w;
    dt_factor = kappa( ( eps1^Beta1* eps0^Beta2 * epsminus^Beta3 * (dt/dtold)^(- Alpha2) )^(1/p) );
    %dt_factor = (eps1^Beta1*eps0^Beta2*epsminus^Beta3)^(1/p);
    dt = dt*dt_factor;
    if dt_factor < fsafety
        flag = true;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ODE Solver Test Suite                                               %%%
%%%                                                                     %%%
%%% Author: Stefan Kopecz (kopecz@mathematik.uni-kassel.de)             %%%
%%% Date: 10/26/2017                                                    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t, Y, varargout] = MPRK43IIadap(PD, y0, dt0, tint, varargin)
%warning('off','all')

if nargin > 4, gamma = varargin{1}; else gamma = 0; end
if nargin > 5
  RTOL = varargin{2}(1); ATOL = varargin{2}(2);
else
  RTOL = 1e-2; ATOL = 1e-2;
end
if nargin > 7
    Beta1 =  varargin{4};
    Beta2 =  varargin{5};
    Beta3 =  varargin{6};
    Alpha2 = varargin{7};
    kappa2 = varargin{8};
end
if length(varargin) >8
    para = varargin{9};
end
if length(varargin) == 10
    scaleode = varargin{10};
end
errctrl = 'PID'; % PID, EPUS, SWP, off


assert(3/8 <= gamma && gamma <= 3/4)

a21 = 2/3;
a31 = 2/3 - 1/(4*gamma);
a32 = 1/(4*gamma);
b1 = 1/4;
b2 = 3/4 - gamma;
b3 = gamma;

beta2 = 1/(2*a21);
beta1 = 1-beta2;

p = 1/(3*a21*(a31+a32)*b3);
q = 1/a21;

assert(abs(b2*a21+b3*(a31+a32)-1/2)<1e-14)
assert(abs(b2*a21^2+b3*(a31+a32)^2-1/3)<1e-14)
assert(abs(a21*a32*b3-1/6)<1e-14)
assert(abs(beta1+beta2-1)<1e-14)
assert(abs(a21*beta2-1/2)<1e-14)

N = length(y0);
Idmat = speye(N);
tstart = tint(1);
tend = tint(2);

t = tstart;
dt = dt0;
dtold = dt;

if strcmp(errctrl,'PID') && t(end) == tstart
    epsminus = 1; %eps_{-1}
    eps0 = 1; %eps_0
end

Y = y0;
y = y0 + realmin;

%global cnt_rej;
cnt_rej = 0;
lengtht = 0;

 %deltat = 0;
while t(end) < tend
  if t(end) + dt > tend
    dt = tend - t(end);
  end
 y = y + realmin;
  yn = y;
  
  if isnan(para)
      [P,D,~,~,~,Rrem, Rminus] = PD(t(end),y,scaleode);
  else
      [P,D,~,~,~,Rrem, Rminus] = PD(t(end),y,para);
  end
  % M2 = zeros(N,N);
  % for i = 1:N
  %   M2(i,i) = 1 + a21*dt*(sum(D(i,:)) + Rminus(i))/y(i);
  %   for j = 1:N
  %     if j ~= i
  %       M2(i,j) = -a21*dt*P(i,j)/y(j);
  %     end
  %   end
  % end
  f1 = diag( sum(D,2) + Rminus ) - P;
  M2 = Idmat + a21*dt*f1*diag(1./y);
  y2 = M2\(y + a21*dt*Rrem) + realmin;
  

  if isnan(para)
      [P2, D2,~,~,~,R2rem,R2minus] = PD(t(end)+2/3*dt,y2,scaleode);
  else
      [P2, D2,~,~,~,R2rem,R2minus] = PD(t(end)+2/3*dt,y2,para);
  end
  % M3 = zeros(N,N);
  rho = y.*(y2./y).^p+ realmin;
  % for i = 1:N
  %   M3(i,i) = 1 + dt*((a31*sum(D(i,:)) + Rminus(i)) + a32*(sum(D2(i,:)) + R2minus(i)))/rho(i);
  %   for j = 1:N
  %     if j ~= i
  %       M3(i,j) = -dt*(a31*P(i,j) + a32*P2(i,j))/rho(j);
  %     end
  %   end
  % end
  f2 = diag( sum(D2,2) + R2minus ) - P2;
  M3 = Idmat + (a31*dt*f1 + a32*dt*f2)*diag(1./rho);
  y3 = M3\(y + dt*(a31*Rrem + a32*R2rem))+realmin;
 
  

  % M4= zeros(N,N);
  mu = y.*(y2./y).^q+realmin;
  % for i = 1:N
  %   M4(i,i) = 1 + dt*(beta1*(sum(D(i,:)) + Rminus(i)) + beta2*(sum(D2(i,:)) + R2minus(i)) )/mu(i);
  %   for j = 1:N
  %     if j ~= i
  %       M4(i,j) = -dt*(beta1*P(i,j) + beta2*P2(i,j))/mu(j);
  %     end
  %   end
  % end
  M4 = Idmat + (beta1*dt*f1 + beta2*dt*f2)*diag(1./mu);
  sigma = M4\(y+dt*(beta1*Rrem+beta2*R2rem))+realmin;

  if isnan(para)
      [P3, D3,~,~,~,R3rem,R3minus] = PD(t(end)+2/3*dt,y3,scaleode);
  else
      [P3, D3,~,~,~,R3rem,R3minus] = PD(t(end)+2/3*dt,y3,para);
  end
  % M = zeros(N,N);
  % for i = 1:N
  %   M(i,i) = 1 + dt*(b1*(sum(D(i,:)) + Rminus(i)) + b2*(sum(D2(i,:)) + R2minus(i)) + b3*(sum(D3(i,:)) + R3minus(i)))/sigma(i);
  %   for j = 1:N
  %     if j ~= i
  %       M(i,j) = -dt*(b1*P(i,j) + b2*P2(i,j) + b3*P3(i,j))/sigma(j);
  %     end
  %   end
  % end
  f3 = diag( sum(D3,2) + R3minus ) - P3;
  M = Idmat + (b1*dt*f1 + b2*dt*f2 + b3*dt*f3)*diag(1./sigma) ;
  y = M\(y+ dt*(b1*Rrem + b2*R2rem + b3*R3rem))+realmin;

  %assert(isreal(P2),'test')

  Y(:,end+1) = y;
  t(end+1) = t(end) + dt;

%   if t(end)-t(end-1) < 1e-15
%   %    t(end) = [];
%   %    Y(:,end) = [];
%      % disp(t(end)-t(end-1))
%       deltat = NaN;
%     %  break
%   else
%       deltat = 0;
%   end

  pp = 3; % order of main method  
  if strcmp(errctrl,'PID')
    dt_temp = dt;
    [dt,flag,eps1] = compute_dt(y, yn, sigma, dt, RTOL, ATOL, errctrl, pp, epsminus, eps0, Beta1, Beta2, Beta3, Alpha2, kappa2,dtold);
    if ~flag
        dtold = dt_temp;
    end
  elseif strcmp(errctrl,'off')
      flag = false;
  else
    [dt,flag] = compute_dt(y,yn,sigma,dt,RTOL,ATOL,errctrl,p);
  end
  if flag
    Y(:,end) = [];
    t(end) = [];
    y = Y(:,end);
    cnt_rej = cnt_rej + 1;
  else
      if strcmp(errctrl,'PID')
          epsminus = eps0;
          eps0 = eps1;
      end
  end
  %   %r = norm((y - sigma)./y,inf);
  %   TOL = RTOL*norm(y,'inf') + ATOL;
  %   r = norm(y - sigma,inf);
  %   if r > TOL
  %     Y(:,end) = [];
  %     t(end) = [];
  %     y = Y(:,end);
  %     dt = dt*max([fmin, fsafety*(TOL/r)^(1/4)]);
  %     cnt_rej = cnt_rej + 1;
  %   else
  %     dt = dt*min([fmax, fsafety*(TOL/r)^(1/4)]);
  %   end
 % fprintf('%i failed attempts\n\n',cnt_rej);
 
 if length(t) == 1e+6  || dt < 1e-100
     lengtht = NaN;
     break
 else
     lengtht = 0;
 end
 if cnt_rej == 1e+4 || cnt_rej > 1e2 * length(t)
     cnt_rej = NaN;
     break
 end

  
end
if nargout > 1
      varargout{1} = cnt_rej;
      varargout{2} = lengtht;
   %   varargout{3} = deltat;
end
%fprintf('%s(%.3f)\n',mfilename,gamma)
%fprintf('%i successful steps\n',length(t)-1)
% fprintf('%i failed attempts\n\n',cnt_rej);
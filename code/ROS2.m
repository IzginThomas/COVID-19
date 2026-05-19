%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ODE Solver Test Suite                                               %%%
%%%                                                                     %%%
%%% Author: Stefan Kopecz (kopecz@mathematik.uni-kassel.de)             %%%
%%% Date: 29/08/2018                                                    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t, Y, varargout] = ROS2(PD, y0, dt0, tint,varargin)
%
% P     Matrix der Produktionsterme mit p_{i,j} = P(i,j),
% D     Matrix der Produktionsterme mit d_{i,j} = D(i,j),

if nargin > 4, gamma = varargin{1}; else gamma = 1/(2+sqrt(2)); end
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

tstart = tint(1);
tend = tint(2);

t = tstart;
dt = dt0;
dtold = dt;


if strcmp(errctrl,'PID') && t(end) == tstart
    epsminus = 1; %eps_{-1}
    eps0 = 1; %eps_0
end
N = length(y0);

Y = y0;
y = y0;

if isnan(para)
    [~,~,f] = PD(0,y,scaleode);
else
    [~,~,f] = PD(0,y,para);
end

h = 1e-8;
% Ay = zeros(N);

cnt_rej = 0;
lengtht = 0;

while t(end) < tend
    if t(end) + dt > tend
        dt = tend - t(end);
    end
    
    yn = y;
%    tic;
    % for i=1:N
    %   e = zeros(N,1);
    %   e(i) = 1;
    %   Ay(:,i)=(f(t(end),y + h*e) - f(t(end),y - h*e))/(2*h);
    % end
    E = eye(N);
    Ay = (f(t(end), y + h*E) - f(t(end), y - h*E)) / (2*h);
    % Ay = cell2mat(arrayfun(@(i) ...
    %     (f(t(end), y + h*E(:,i)) - f(t(end), y - h*E(:,i))) / (2*h), ...
    %     1:N, 'UniformOutput', false));
  %  toc;
    At = (f(t(end) + h, y) - f(t(end) - h, y))/(2*h);

    % k1 = (eye(N) - gamma*dt*Ay)\(dt*f(t(end),y) + gamma*dt^2*At);
    % k2 = (eye(N) - gamma*dt*Ay)\(dt*f(t(end)+dt,y + k1) - 2*k1);
    W = (E - gamma*dt*Ay);
    k1 = W\(f(t(end),y) + gamma*dt*At);
    k2 = W\(f(t(end)+0.5*dt,y + 0.5*dt*k1) - k1) + k1;
    sigma = y + dt*k1;
    % y = y + 1.5*k1 + 0.5*k2;
    y = y + dt*k2;

    Y(:,end+1) = y;
    t(end+1) = t(end) + dt;
    if any(y<0)
        warning('ODESolverTestSuite:ROS2 violated positivity constraint\n')
        %break
    end

    p = 2; % order of main method
    if strcmp(errctrl,'PID')
        dt_temp = dt;
        [dt,flag,eps1] = compute_dt(y,yn,sigma,dt,RTOL,ATOL,errctrl,p,epsminus,eps0,Beta1,Beta2,Beta3, Alpha2, kappa2, dtold);
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
        % disp(cnt_rej)
    else
        if strcmp(errctrl,'PID')
            epsminus = eps0;
            eps0 = eps1;
        end
    end
    if length(t) == 1e+6  || dt < 1e-100
        lengtht = NaN;
        break
    else
        lengtht = 0;
    end
    if cnt_rej == 1e+4 || cnt_rej > 1e2 * length(t)
    	cnt_rej =NaN;
        break
    end


end
if nargout > 1
    varargout{1} = cnt_rej;
    varargout{2} = lengtht;
    %varargout{3} = deltat;
end
end

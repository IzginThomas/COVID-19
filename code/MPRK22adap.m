%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ODE Solver Test Suite                                               %%%
%%%                                                                     %%%
%%% Author: Stefan Kopecz (kopecz@mathematik.uni-kassel.de)             %%%
%%% Date: 10/26/2017                                                    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t, Y,varargout] = MPRK22adap(PD, y0, dt0, tint, varargin)
warning('off','all')

if nargin > 4, a21 = varargin{1}; else a21=1; end
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


b2 = 1/(2*a21);
b1 = 1 - b2;

assert(a21 ~= 0)

I = length(y0);
Idmat = speye(I);
tstart = tint(1);
tend = tint(2);

t = tstart;
dt = dt0;
dtold = dt;


if strcmp(errctrl,'PID') && t(end) == tstart
    epsminus = 1; %eps_{-1}
    eps0 = 1; %eps_0
end


y = y0 + realmin;
Y = y;


cnt_rej = 0;
lengtht = 0;
% deltat = 0;

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

    % M2 = zeros(I);
    % for i = 1:I
    %     if a21 > 0
    %         M2(i,i) = 1 + a21*dt*(sum(D(i,:)) + Rminus(i))/y(i);
    %     else
    %         M2(i,i) = 1 - a21*dt*(sum(P(i,:)))/y(i);
    %     end
    %     for j = 1:I
    %         if j ~= i
    %             if a21 > 0
    %                 M2(i,j) = -a21*dt*P(i,j)/y(j);
    %             else
    %                 M2(i,j) =  a21*dt*D(i,j)/y(j);
    %             end
    %         end
    %     end
    % end
    f1 = diag( sum(D,2) + Rminus ) - P;
    M2 = Idmat + a21*dt*f1*diag(1./y);
    y2 = M2\(y + a21*dt*Rrem + realmin);

    if isnan(para)
        [P2, D2,~,~,~,R2rem,R2minus] = PD(t(end)+a21*dt,y2,scaleode);
    else
        [P2, D2,~,~,~,R2rem,R2minus] = PD(t(end)+a21*dt,y2,para);
    end

    % M = zeros(I);
    sigma = y.*(y2./y).^(1/a21) + realmin;
    %sigma = y2.^(1/a21).*y.^(1-1/a21);
    % for i = 1:I
    %     if b1 >= 0
    %         M(i,i) = 1 + dt*b1*(sum(D(i,:)) + Rminus(i))/(sigma(i));
    %     else
    %         M(i,i) = 1 - dt*b1*(sum(P(i,:)))/(sigma(i));
    %     end
    %     if b2 >= 0
    %         M(i,i) = M(i,i) + dt*b2*(sum(D2(i,:)) + R2minus(i))/(sigma(i));
    %     else
    %         M(i,i) = M(i,i) - dt*b2*(sum(P2(i,:)))/(sigma(i));
    %     end
    %     for j=1:I
    %         if j ~= i
    %             if b1 >= 0
    %                 M(i,j) = -dt*b1*P(i,j)/(sigma(j));
    %             else
    %                 M(i,j) = dt*b1*D(i,j)/(sigma(j));
    %             end
    %             if b2 >= 0
    %                 M(i,j) = M(i,j) - dt*b2*P2(i,j)/(sigma(j));
    %             else
    %                 M(i,j) = M(i,j) + dt*b2*D2(i,j)/(sigma(j));
    %             end
    %         end
    %     end
    % end
    f2 = diag( sum(D2,2) + R2minus ) - P2;
    M = Idmat + dt*(b1*f1 + b2*f2)*diag(1./sigma) ;
    rem = b1*dt*Rrem + b2*dt*R2rem;

%     if t(end) >= 800
%         xx=1;
%     end
    y = M\(y + rem + realmin) + realmin;



    Y(:,end+1) = y;
    t(end+1) = t(end) + dt;

%     if t(end)-t(end-1) < 1e-14
%         %t(end) = [];
%         %Y(:,end) = [];
%         % disp(t(end)-t(end-1))
%         deltat = NaN;
%         %break
%     else 
%         deltat = 0;
%     end

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
%fprintf('%s(%.3f)\n',mfilename,a21)
%fprintf('%i successful steps\n',length(t)-1)
%fprintf('%i failed attempts\n\n',cnt_rej);
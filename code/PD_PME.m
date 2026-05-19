%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ODE Solver Test Suite                                               %%%
%%%                                                                     %%%
%%% Author: Thomas Izgin (izgin@mathematik.uni-kassel.de)             %%%
%%% Date: 09/29/2022                                                    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% HIRES problem                                                   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [P, D, f, S, r, Rrem,Rminus] = PD_PME(t,y,varargin)
% getpara();
test_PME(); % getting Nx and dx here
m = para;
a = @(y) m * y.^(m-1);
cs = @(y,k) circshift(y,k);

yplus = @(y) cs([zeros(1,size(y,2));y(2:end,:)],-1);
yminus = @(y) cs([y(1:end-1,:);zeros(1,size(y,2))],1);
interior= @(y) y(2:end-1,:);

P = sparse(Nx,Nx);
y_vecup = (yplus(y).*(a(yplus(y)) + a(y)));
y_veclow = (yminus(y).*(a(yminus(y)) + a(y)));
P_aux = spdiags(cs(y_vecup,1),1, Nx, Nx) + spdiags(cs(y_veclow,-1),-1, Nx, Nx);
P(2:end-1,2:end-1) = P_aux(2:end-1,2:end-1) / (2*dx^2);
P(1,2) = a(y(2))*y(2)/(2*dx^2);
P(end,end-1) = a(y(end-1))*y(end-1)/(2*dx^2);

D = P';



if nargout > 2
    f = @(t,y) [a(y(2,:)).*y(2,:)./(2*dx^2); ...
                interior( ( ( a(y) + a(yplus(y)) ) .* yplus(y) - ...
                ( a(yminus(y)) + 2 .* a(y) + a(yplus(y)) ) .* y + ...
                ( a(yminus(y)) + a(y)) .* yminus(y) ) ./ (2*dx^2) );... 
                a(y(end-1,:)).*y(end-1,:)/(2*dx^2)];
end

if nargout > 3 %% Need to be done
 S=zeros(length(y));
  
  
   r = @(y) zeros(length(y),1);

end

if nargout > 5 
Rrem = zeros(length(y),1); % remainder
Rminus = Rrem;
% Rminus(2) = a(y(2))*y(2)/(2*dx^2);
% Rminus(end-1) = a(y(end-1))*y(end-1)/(2*dx^2);
end


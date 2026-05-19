%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ODE Solver Test Suite                                               %%%
%%%                                                                     %%%
%%% Author: Thomas Izgin (izgin@mathematik.uni-kassel.de)             %%%
%%% Date: 09/29/2022                                                    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% HIRES problem                                                   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [P, D, f, S, r, Rrem,Rminus] = PD_FokkerPlanck(t,y,varargin)

test_FokkerPlanck(); % getting x,x_interfaces, para, Nx,dx

P = production_FP_ranocha(t,y,x,x_interfaces, para, Nx,dx);

D = P';


if nargout > 2
    % f = @(t,y) - (num_flux(y,phys_flux) - circshift(num_flux(y,phys_flux),1)) / dx;
    f = @(t,y) rhs_FP(t,y,x,x_interfaces, para,Nx,dx); %sum(production_FP_ranocha(t,y,x,x_interfaces, para, Nx,dx) - (production_FP_ranocha(t,y,x,x_interfaces, para, Nx,dx))', 2);
end
%if max(abs(sum(P-D,2)-f(t,y))) > 1e-14
    % disp(max(abs(sum(P-D,2)-f(t,y))));
%end
if nargout > 3 %% Need to be done
    S=zeros(length(y));


    r = @(y) zeros(length(y),1);

end

% lmean_y = @(y) log_mean(y,@(y) y);
% y_aux = @(y) (reshape(y,Nx,d))';
% pvec = lmean_y(y).*(reshape((A*y_aux(y))', ...
%     d*Nx,1)) /(dx^2);
%
% P = sparse(d*Nx,d*Nx);
% for k=1:d
%     P((k-1)*Nx + 1: k*Nx, (k-1)*Nx + 1: k*Nx) = spdiags(pvec((k-1)*Nx + 1: k*Nx),1, Nx, Nx);
% end
% no-flux boundary condition for d=3 (u^k_1,u^k_Nx=const, k=1,2,3)
% P([1,Nx-1,Nx+1,2*Nx-1,2*Nx+1,3*Nx-1],:) = sparse(6,3*Nx);
%
% D = P';
%
% v1=ones(d*Nx,1);
% v1([1,Nx,Nx+1,2*Nx,2*Nx+1,3*Nx])=0;
%
%
% if nargout > 2
%
%     f = @(t,y) v1.*lmean_y(y).*reshape((A*(circshift(y_aux(y),-1) - circshift(y_aux(y),1) ))', d*Nx,1) /dx^2;
% end
% [a,b]=max(f(t,y)-sum(P - P', 2))
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


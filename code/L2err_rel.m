%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ODE Solver Test Suite                                               %%%
%%%                                                                     %%%
%%% Author: Stefan Kopecz (kopecz@mathematik.uni-kassel.de)             %%%
%%% Date: 04/23/2018                                                    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function err = L2err_rel(t,Y,Yref)

% N = size(Y,1);
%
% err = 0;
% for j = 1:N
%   err = err + sqrt(sum(diff(t).*abs(Yref(2:end,j)' - Y(j,2:end)).^2))/sqrt(sum(diff(t)'.*abs(Yref(2:end,j)).^2));
%   %err = err + sqrt(sum(diff(t).*abs(Yref(2:end,j)' - Y(j,2:end)).^2));
% end
% err = err/N;

p = 2;
aux = sum(abs(Yref').^p,1);
l2norm_Yex = (sum(diff(t).*(0.5*(aux(1:end-1)+aux(2:end)))))^(1/p);

aux = sum(abs(Y - Yref').^p,1);
l2norm_Ydiff = (sum(diff(t).*(0.5*(aux(1:end-1)+aux(2:end)))))^(1/p);
err= l2norm_Ydiff/l2norm_Yex;







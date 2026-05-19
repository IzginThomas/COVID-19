%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ODE Solver Test Suite                                               %%%
%%%                                                                     %%%
%%% Author: Stefan Kopecz (kopecz@mathematik.uni-kassel.de)             %%%
%%% Date: 04/23/2018                                                    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function err = L2err_rel2d(t,x,Y,Yref)

% Y: space x time
% Yref: time x space

Y = Y';
dts = diff(t);
dxs = diff(x);
aux = abs(Yref).^2;
l2norm_Yex = sum(diag(dts) * (0.25 * (aux(1:end-1,1:end-1) + aux(2:end,1:end-1) + aux(1:end-1,2:end) + aux(2:end,2:end)) ),1);
l2norm_Yex = (sum(dxs .* l2norm_Yex))^(1/2);

aux = abs(Y - Yref).^2;
l2norm_Ydiff = sum(diag(dts) * (0.25 * (aux(1:end-1,1:end-1) + aux(2:end,1:end-1) + aux(1:end-1,2:end) + aux(2:end,2:end)) ),1);
l2norm_Ydiff = (sum(dxs .* l2norm_Ydiff))^(1/2);
err= l2norm_Ydiff/l2norm_Yex;
end







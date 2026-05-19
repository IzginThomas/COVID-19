%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ODE Solver Test Suite                                               %%%
%%%                                                                     %%%
%%% Author: Stefan Kopecz (kopecz@mathematik.uni-kassel.de)             %%%
%%% Date: 09/02/2018                                                    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% NPZD problem                                                        %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [P, D, f, S, r, Rrem, Rminus,varargout] = PD_npzd(t,y,varargin)
% N = y1, P = y2, Z = y3, D = y4

alpha = 0.01;
rmax = 1.0;
gmax = 0.5;
Ipar = 25;
Imin = 25;
Iopt = max(0.25*Ipar,Imin);
Iv = 1.1;
rpn = 0.01;
rzn = 0.01;
rdn = 0.003;
rpd = 0.05;
rzd = 0.02;

dnp = rmax*Ipar/Iopt*exp(1 - Ipar/Iopt)*y(1)/(alpha + y(1))*y(2);
dpz = gmax*(1 - exp(-Iv^2*y(2)^2))*y(3);
dpn = rpn*y(2);
dzn = rzn*y(3);
ddn = rdn*y(4);
dpd = rpd*y(2);
dzd = rzd*y(3);

% D = [0 dnp 0 0;
%      dpn 0 dpz dpd;
%      dzn 0 0 dzd;
%      ddn 0 0 0];     
% 
% P = D';

P = [ 0  dpn dzn ddn;
     dnp  0   0   0 ;
      0  dpz  0   0;
      0  dpd dzd  0]; 
% 
%      [0 dnp 0 0;
%      dpn 0 dpz dpd;
%      dzn 0 0 dzd;
%      ddn 0 0 0]; 

%D = [dnp;dpn + dpz + dpd; dzn + dzd; ddn];
D = P';

if nargout > 2
dpn = @(t,y) rpn*y(2,:);
dzn = @(t,y) rzn*y(3,:);
ddn = @(t,y) rdn*y(4,:);
dnp = @(t,y) rmax*Ipar/Iopt.*exp(1 - Ipar/Iopt).*y(1,:)./(alpha + y(1,:)).*y(2,:);
dpz = @(t,y) gmax.*(1 - exp(-Iv^2.*y(2,:).^2)).*y(3,:);
dzd = @(t,y) rzd*y(3,:);
dpd = @(t,y) rpd*y(2,:);

f = @(t,y) [dpn(t,y) + dzn(t,y) + ddn(t,y) - dnp(t,y);
            dnp(t,y) - (dpn(t,y) + dpz(t,y) + dpd(t,y));
            dpz(t,y) - (dzn(t,y) + dzd(t,y));
            dpd(t,y) + dzd(t,y) - ddn(t,y)];
end

if nargout > 3
  S = [ 1  1  1 -1  0  0  0; 
       -1  0  0  1 -1  0 -1;
        0 -1  0  0  1 -1  0;
        0  0 -1  0  0  1  1];
  r = @(y) [dpn(t,y); dzn(t,y); ddn(t,y); dnp(t,y); dpz(t,y); dzd(t,y); dpd(t,y)];
end


if nargout > 5 

Rrem = zeros(length(y),1); % remainder
Rminus = Rrem;

end
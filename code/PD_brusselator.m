%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ODE Solver Test Suite                                               %%%
%%%                                                                     %%%
%%% Author: Stefan Kopecz (kopecz@mathematik.uni-kassel.de)             %%%
%%% Date: 10/26/2017                                                    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Brusselator problem                                                 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [P, D, f, S, r, Rrem,Rminus] = PD_brusselator(t,y,varargin)

P = [  0      0       0    0        0           0;
       0      0       0    0        0           0;
       0  y(2)*y(5)   0    0        0           0;
       0      0       0    0      y(5)          0;
     y(1)     0       0    0        0      y(5)^2*y(6);
       0      0       0    0    y(2)*y(5)       0];

D = P';

if nargout > 2
  f = @(t,y) [-y(1,:);
              -y(2,:).*y(5,:);
               y(2,:).*y(5,:);
               y(5,:);
               y(1,:) - y(2,:).*y(5,:) + y(5,:).^2.*y(6,:) - y(5,:);
               y(2,:).*y(5,:) - y(5,:).^2.*y(6,:)]; 
end
if nargout > 3 %% Need to be done
 S=zeros(length(y));
   r = @(y) zeros(length(y),1);
end

if nargout > 5 
Rrem = zeros(length(y),1); % remainder
Rminus = Rrem;
end

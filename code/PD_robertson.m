%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ODE Solver Test Suite                                               %%%
%%%                                                                     %%%
%%% Author: Stefan Kopecz (kopecz@mathematik.uni-kassel.de)             %%%
%%% Date: 10/26/2017                                                    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Robertson problem                                                   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [P, D, f, S, r, Rrem, Rminus] = PD_robertson(t,y,varargin)

A = 1e4;
B = 4e-2;
C = 3e7;

p12 = A*y(2)*y(3);
p21 = B*y(1);
p32 = C*y(2)*y(2);

% P = [  0    A*y(2)*y(3)   0;
%      B*y(1)      0        0;
%        0    C*y(2)*y(2)   0];        

P = [0 p12 0; p21 0 0; 0 p32 0];

% D = P';
D = [0 p21 0;p12 0 p32; 0 0 0];

if nargout > 2
f = @(t,y) [A*y(2,:).*y(3,:) - B*y(1,:);
            B*y(1,:) - A*y(2,:).*y(3,:) - C*y(2,:).*y(2,:);
            C*y(2,:).*y(2,:)];
end

if nargout > 3
  S = [1 -1 0;-1 1 -1; 0 0 1];
  r = @(y) [A*y(2)*y(3); B*y(1); C*y(2)*y(2)];
end
 if nargout > 5
  Rrem = zeros(length(y),1);
  Rminus = Rrem;
end
% [msg, id] = lastwarn;          
% if (A~=1e4 || B~=4e-2 || C~=3e7) && ...
%     ~strcmp(id,'ODESolverTestSuite:PARAMETER:robertson')
%   warning('ODESolverTestSuite:PARAMETER:robertson', ...
%     'Parameters of Robertson problem changed:\n  A=%e, B=%e, C=%e\n', ...
%     A, B, C);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ODE Solver Test Suite                                               %%%
%%%                                                                     %%%
%%% Author: Thomas Izgin (izgin@mathematik.uni-kassel.de)               %%%
%%% Date: 09/29/2022                                                    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Prothero & Robinson Problem with 4 Equations                        %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [P, D, f, S, r, Rrem, Rminus,varargout] = PD_PR4(t,y,varargin)

if nargin >2 && ~isempty(varargin{1})
    para = varargin{1};
else
    error('para is not specified in PD_nonlin4')
end
oscifac =@(t) 0.1*cos(0.5*t);
oscidiff =@(t) -0.05*sin(.5*t); % derivative of oscifac(t)
proddiff =@(t) oscidiff(t)*t + oscifac(t); % derivative of oscifac(t)*t
k1 = 0.3; %k1<=0.5 for positivity

A = [-1 1-para para 0; para -1 0 1-para; 1-para 0 -1 para;0 para 1-para -1];

exsol=@(t) [2+k1*sin(oscifac(t)*t);2+sin(oscifac(t)*t);1-sin(oscifac(t)*t);1-k1*sin(oscifac(t)*t)];
g = exsol(t);
gprime=[k1*proddiff(t)*cos(oscifac(t)*t); proddiff(t)*cos(oscifac(t)*t);-proddiff(t)*cos(oscifac(t)*t); -k1*proddiff(t)*cos(oscifac(t)*t)];


P(1,2) = y(2);
P(1,3) = g(1);
P(1,4) = para*y(3) + para*g(2) + min(0,gprime(1));
P(2,1) = g(2);
P(2,3) = para*g(4)+ para*y(1) + min(0,gprime(2));   
P(2,4) = y(4);
P(3,1) = y(1);%  + k1*min(0,gprime(3));
P(3,2) = para*(g(1) + y(4)) + min(0,gprime(3));
P(3,4) = g(3);
P(4,1) = para*(y(2) + g(3)) + min(0,gprime(4)); 
P(4,2) = g(4);
P(4,3) = y(3);

P=real(P);

D = P';

if nargout > 2
f = @(t,y) A*(y - [2+k1*sin(oscifac(t)*t);2+sin(oscifac(t)*t);1-sin(oscifac(t)*t);1-k1*sin(oscifac(t)*t)]) + [k1*proddiff(t)*cos(oscifac(t)*t);proddiff(t)*cos(oscifac(t)*t);-proddiff(t)*cos(oscifac(t)*t);-k1*proddiff(t)*cos(oscifac(t)*t)];
end

if nargout > 3 %% Need to be done
 S=zeros(length(y));
  
  
   r = @(y) zeros(length(y),1);

end

if nargout > 5 

Rrem = zeros(length(y),1); % remainder
Rminus = Rrem;

end
if nargout > 7 

varargout{1} = exsol;

end

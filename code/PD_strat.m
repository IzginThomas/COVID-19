%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ODE Solver Test Suite                                               %%%
%%%                                                                     %%%
%%% Author: Thomas Izgin (izgin@mathematik.uni-kassel.de)               %%%
%%% Date: 09/29/2022                                                    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Stratospheric reaction problem    RESCALED to be CONSERVATIVE       %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [P, D, f, S, r, Rrem,Rminus] = PD_strat(t,y,scaleode)

T = @(t) mod(t/3600 + 12,24); % t in seconds, t=0 --> start at noon
Tr = 4.5;
Ts = 19.5;

function res = sig(t)
    if (Tr <= T(t)) && (T(t) <= Ts)
        %sig =@(t) ( 0.5*(T(t) - Tr)./abs(T(t) - Tr) + 0.5*(Ts - T(t))./abs(Ts - T(t))) .* (0.5 + 0.5*cos( pi * abs((2*T(t) - Tr - Ts)./(Ts - Tr)) .* (2*T(t) - Tr - Ts)./(Ts - Tr) ));
        res = 0.5 + 0.5*cos( pi * abs((2*T(t) - Tr - Ts)./(Ts - Tr)) .* (2*T(t) - Tr - Ts)./(Ts - Tr) );
    else
        res = 0;
    end
end


k1 = @(t)  (sig(t)).^3*2.643e-10;
r1 = @(t,y) k1(t)*y(4);

k2 = 8.018e-17;
r2 = @(y) k2*y(4)*y(2)/scaleode;

k3 = @(t) (sig(t))*6.120e-4;
r3 = @(t,y) k3(t)*y(3);

k4 = 1.576e-15;
r4 = @(y) k4*y(2)*y(3)/scaleode;

k5 = @(t) (sig(t)).^2*1.070e-3;
r5 = @(t,y) k5(t)*y(3);

M = 8.120e+16;
k6 = 7.110e-11;
r6 = @(y) k6*M*y(1);

k7 = 1.200e-10;
r7 = @(y) k7*y(1)*y(3)/scaleode;

k8 = 6.062e-15;
r8 = @(y) k8*y(3)*y(5)/scaleode;

k9 = 1.069e-11;
r9 = @(y) k9*y(2)*y(6)/scaleode;

k10 = @(t) (sig(t))*1.289e-2;
r10 = @(t,y) k10(t)*y(6);

k11 = 1.0e-8;
r11 = @(y) k11*y(5)*y(2)/scaleode;



D(1,2) = r6(y);
D(1,4) = 1/3*r7(y);
D(2,3) = 1/2*r2(y);
D(2,4) = 1/3*r4(y);
D(2,5) = 1/2*r9(y);
D(2,6) = r11(y);
D(3,1) = 1/3*r5(t,y);
D(3,2) = 1/3*r3(t,y);
D(3,4) = 2/3*r3(t,y) + r4(y)  + r7(y) + 2/3*r8(y) + 2/3*r5(t,y);
D(3,6) = 1/3*r8(y);
D(4,2) = r1(t,y);
D(4,3) = r2(y);
D(5,6) = r11(y) + 1/3*r8(y);
D(6,2) = 1/2*r10(t,y);
D(6,4) = r9(y);
D(6,5) = 1/2*r10(t,y);


     
P = D';


if nargout > 2
    
f = @(t,y) [1/3*r5(t,y) - (r6(y) + 1/3*r7(y));
           1/2*r10(t,y) + r1(t,y) + 1/3*r3(t,y) + r6(y) - (r11(y) + 1/2*r2(y) + 1/3*r4(y) + 1/2*r9(y));
           3/2*r2(y) - (r3(t,y) + r4(y)  + r7(y) + r8(y) + r5(t,y));
           r9(y) + 2/3*r3(t,y) + 4/3*r4(y)  + 4/3*r7(y) + 2/3*r8(y) + 2/3*r5(t,y) - (r2(y) + r1(t,y));
           1/2*r10(t,y) + 1/2*r9(y) - (r11(y) + 1/3*r8(y));
           2/3*r8(y) + 2*r11(y) - (r9(y) + r10(t,y))];
%f = @(t,y) [r5(t,y) - r6(y) - r7(y);
%          2*r1(t,y) - r2(y) + r3(t,y) - r4(y) + r6(y) - r9(y) + r10(t,y) - r11(y);
 %          r2(y) - r3(t,y) - r4(y) - r5(t,y) - r7(y) - r8(y);
  %         -r1(t,y) - r2(y) + r3(t,y) + 2*r4(y) + r5(t,y) + 2*r7(y) + r8(y) + r9(y);
   %        -r8(y) + r9(y) + r10(t,y) - r11(y);
    %       r8(y) - r9(y) - r10(t,y) + r11(y)];
%f = @(t,y) diag([9.906e1; 6.624e8; 5.326e11; 1.697e16; 4.000e16; 1.093e9])\g(t,diag([9.906e1; 6.624e8; 5.326e11; 1.697e16; 4.000e16; 1.093e9])*y);
end

if nargout > 3 %% Need to be done
 S=zeros(length(y));
  
  
   r = @(y) zeros(length(y),1);

end

if nargout > 5 

Rrem = zeros(length(y),1); % remainder
Rminus = Rrem;

end
end

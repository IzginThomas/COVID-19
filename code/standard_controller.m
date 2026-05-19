function [varargout] = standard_controller()

s1.beta1=0.6;
s1.beta2=-0.2;
s1.beta3=0;
s1.alpha2=0;
s1.kappa2=1;

s2.beta1=0.7;
s2.beta2=-0.4;
s2.beta3=0;
s2.alpha2=0;
s2.kappa2=1;

s3.beta1=1/6;
s3.beta2=-1/3;
s3.beta3=0;
s3.alpha2=0;
s3.kappa2=1;

s4.beta1=1/6;
s4.beta2=1/6;
s4.beta3=0;
s4.alpha2=0;
s4.kappa2=1;

s5.beta1=1;
s5.beta2=0;
s5.beta3=0;
s5.alpha2=0;
s5.kappa2=1;

s6.beta1=2;
s6.beta2=-1;
s6.beta3=0;
s6.alpha2=-1;
s6.kappa2=1;

s7.beta1=0.5;
s7.beta2=0.5;
s7.beta3=0;
s7.alpha2=0.5;
s7.kappa2=1;

s8.beta1=0.056;
s8.beta2=0.111;
s8.beta3=0.056;
s8.alpha2=0;
s8.kappa2=1;

s9.beta1=0.25;
s9.beta2=0.25;
s9.beta3=0.25;
s9.alpha2=0;
s9.kappa2=1;

varargout{1} = s1;
varargout{2} = s2;
varargout{3} = s3;
varargout{4} = s4;
varargout{5} = s5;
varargout{6} = s6;
varargout{7} = s7;
varargout{8} = s8;
varargout{9} = s9;
end
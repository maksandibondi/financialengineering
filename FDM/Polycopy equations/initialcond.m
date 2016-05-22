function [out] =  initialcond(x)
out = exp(-x^2)*x*(1-x); 
return


% cos(2*pi*x)+x^2-2
% 3/2*(x^2)
% exp(-x^2)*x*(1-x)
% x^2
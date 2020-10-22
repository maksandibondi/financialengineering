function [out] =  initialcond_hyper(x)
%sigma = 0.005;

if (x>=0 && x<(1/3))
    
out = 3*x;

elseif (x>=(1/3) && (x<=1))
        
out = (3-3*x)/2;

else
    
out = 0;

end;

 

return

% out = exp((-(x-1)^2)/(2*sigma));  WITH  sigma = 0.005 ; for physical
% if x>0 & x<1/3
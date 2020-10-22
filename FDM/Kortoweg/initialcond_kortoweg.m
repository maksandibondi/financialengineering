function [out] =  initialcond_kortoweg(x)

V(1) = 400; V(2) = 200; V(3) = 150;

out = V(1)/(2*(cosh((V(1)^(1/2))*(x-1)/2))^2) + V(2)/(2*(cosh((V(2)^(1/2))*(x-1.8)/2))^2) + V(3)/(2*(cosh((V(3)^(1/2))*(x-2.5)/2))^2);


return
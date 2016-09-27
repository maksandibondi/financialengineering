function [path] = EconomicVarSimul(T,N,Ro,sigma2,initial_value,bmpath)

path(1) = initial_value;
delta_t = T/N;

for i = 2:N
path(i) = path(i-1) - Ro*path(i-1)*delta_t + sigma2*(bmpath(i)-bmpath(i-1));
end;
 
return;



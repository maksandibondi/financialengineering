function [path] = RateSimulator(T,N,r0,nu,gama,sigma)

delta_t = T/N;
path(1) = r0;


bmpath = BMsimulator(T,N,0);

for i = 2:N
path(i) = path(i-1) + (nu-gama*path(i-1))*delta_t+sigma*(bmpath(i)-bmpath(i-1));
end;
 
return;


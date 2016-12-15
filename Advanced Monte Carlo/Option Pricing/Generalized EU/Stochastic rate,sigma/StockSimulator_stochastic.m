function [path] = StockSimulator_stochastic(T,N,S0,r,sigma)

delta_t = T/N;

path(1) = S0;

bmpath = BMsimulator(T,N,0);

for i = 2:N
    
    path(i) = path(i-1) + path(i-1)*(r(i-1)*delta_t+sigma(i-1)*(bmpath(i)-bmpath(i-1)));

end;
 
return;
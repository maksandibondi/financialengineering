function [path] = WealthSimul(T,N,r,nu,lambda,sigma1,initial_value_M,Ro,sigma2,initial_value_X)

teta = (1:N+1)*0.001;
path(1) = initial_value_M;
delta_t = T/N;


bmpath = BMsimulator(T,N,0);
xpath = EconomicVarSimul(T,N,Ro,sigma2,initial_value_X,bmpath);

for i = 2:N
path(i) = path(i-1) + path(i-1)*(r*(1-teta(i))+teta(i)*(nu+lambda*xpath(i)))*delta_t + teta(i)*sigma1*(bmpath(i)-bmpath(i-1));
end;


return;


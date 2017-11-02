function [path] = FundSimul(T,N,nu,lambda,sigma1,Ro,sigma2,initial_value_X,initial_value_P)

path(1) = initial_value_P;
delta_t = T/N;


bmpath = BMsimulator(T,N,0);
xpath = EconomicVarSimul(T,N,Ro,sigma2,initial_value_X,bmpath);

for i = 2:N
path(i) = path(i-1) + path(i-1)*(nu+lambda*xpath(i))*delta_t + sigma1*path(i-1)*(bmpath(i)-bmpath(i-1));
end;


return;



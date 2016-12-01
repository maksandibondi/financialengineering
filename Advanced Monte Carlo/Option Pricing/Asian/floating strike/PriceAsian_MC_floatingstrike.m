function [AsianCall, AsianPut] = PriceAsian_MC_floatingstrike(S0,T,N,r,sigma)

num_of_iter = 1000;
AsianCall = 0;
AsianPut = 0;

for i = 1:num_of_iter
    [call, put] = payoff_asian_floatingstrike(T,N,S0,r,sigma);
    AsianCall = AsianCall + exp(-r*T)*call/num_of_iter;
    AsianPut = AsianPut + exp(-r*T)*put/num_of_iter;
end; 

return;
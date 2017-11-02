function [AsianCall, AsianPut] = PriceAsian_MC(S0,At,t,T,N,r,sigma,K)

num_of_iter = 100;
AsianCall = 0;
AsianPut = 0;

for i = 1:num_of_iter
    [call, put] = payoff_asian(At,T,t,N,S0,r,sigma,K);
    AsianCall = AsianCall + exp(-r*T)*call/num_of_iter;
    AsianPut = AsianPut + exp(-r*T)*put/num_of_iter;
end; 

return;
function [LookBack] = Price_lookback_MC(S0,T,N,r,sigma,K)

num_of_iter = 100;
LookBack = 0;

for i = 1:num_of_iter
    payoff = payoff_lookback(T,N,S0,r,sigma,K);
    LookBack = LookBack+ exp(-r*T)*payoff/num_of_iter;
end; 

return;
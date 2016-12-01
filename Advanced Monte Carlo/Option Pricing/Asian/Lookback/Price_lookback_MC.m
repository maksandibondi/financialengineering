function [LookBackCall, LookBackPut] = Price_lookback_MC(S0,T,N,r,sigma,K)

num_of_iter = 100;
LookBackCall = 0;
LookBackPut = 0;

for i = 1:num_of_iter
    [call, put] = payoff_lookback(T,N,S0,r,sigma,K);
    LookBackCall = LookBackCall + exp(-r*T)*call/num_of_iter;
    LookBackPut = LookBackPut + exp(-r*T)*put/num_of_iter;
end; 

return;
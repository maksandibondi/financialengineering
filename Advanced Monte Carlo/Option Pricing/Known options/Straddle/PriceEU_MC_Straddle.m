function [Straddle] = PriceEU_MC_Straddle(S0,K,r,sigma,T)

num_of_iter = 500;
sum_straddle = 0;

for k = 1:num_of_iter
    S_end(k) = S0*exp((r-(sigma^2)/2)*T+sigma*sqrt(T)*randn());
    sum_straddle = sum_straddle + payoff_straddle(S_end(k),K)*exp(-r*T);
end;

Straddle = sum_straddle/num_of_iter;

return;


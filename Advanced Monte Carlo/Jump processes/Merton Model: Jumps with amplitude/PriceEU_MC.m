function [price] = PriceEU_MC(S0,K,r,sigma,T)

num_of_iter = 1000;
summ = 0;

for k = 1:num_of_iter
    S_end(k) = S0*exp((r-(sigma^2)/2)*T+sigma*sqrt(T)*randn());
    summ = summ + payoff(S_end(k),K)*exp(-r*T);
end;

price = summ/num_of_iter;

return;


function [Butterfly] = PriceEU_Butterfly_MC(S0,K,r,sigma,T)

num_of_iter = 100;
sum = 0;

for k = 1:num_of_iter
    S_end(k) = S0*exp((r-(sigma^2)/2)*T+sigma*sqrt(T)*randn());
    sum = sum + payoff_butterfly(S_end(k),K)*exp(-r*T);
end;

Butterfly = sum/num_of_iter;

return;


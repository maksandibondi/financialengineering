function [price] = PriceEU_forward(S0,r,sigma,T,T1,alpha)

num_of_iter = 1000;
summ = 0;

for k = 1:num_of_iter
    S_T1 = S0*exp((r-(sigma^2)/2)*T1+sigma*sqrt(T1)*randn());
    K = alpha*S_T1;
    S_end = S_T1*exp((r-(sigma^2)/2)*(T-T1)+sigma*sqrt(T-T1)*randn());
    summ = summ + payoff_forwardstart(S_end,K)*exp(-r*T);
end;

price = summ/num_of_iter;

return;


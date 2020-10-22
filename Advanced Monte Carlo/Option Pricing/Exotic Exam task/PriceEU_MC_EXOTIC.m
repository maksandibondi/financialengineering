function [price] = PriceEU_MC_EXOTIC(S1_0,S2_0,K,r1,r2,sigma1,sigma2,T)

num_of_iter = 10000;
summ = 0;

for k = 1:num_of_iter
    S1_middle = S1_0*exp((r1-(sigma1^2)/2)*(T/2)+sigma1*sqrt(T/2)*randn());
    S2_middle = S2_0*exp((r2-(sigma2^2)/2)*(T/2)+sigma2*sqrt(T/2)*randn());
    
    if S2_middle<=exp(-1/4)
        sigma1_new = 1/2;
    else
        sigma1_new = 2;
    end;
    
    S1_end(k) = S1_middle*exp((r1-(sigma1_new^2)/2)*(T/2)+sigma1_new*sqrt(T/2)*randn());
    summ = summ + payoff(S1_end(k),K)*exp(-r1*T);
end;

price = summ/num_of_iter;

return;


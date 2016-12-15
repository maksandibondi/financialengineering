function [price] = PriceEU_MC_stochasticRate(S0,N,K,sigma,r0,nu_r,gama_r,sigma_r,T)

num_of_iter = 100;
summ = 0;

delta_t = T/N;

for k = 1:num_of_iter
    
    sum_r = 0; % will be used for actualization
    r(k,:) = RateSimulator(T,N,r0,nu_r,gama_r,sigma_r);
    S(k,:) = StockSimulator_stochastic(T,N,S0,r(k,:),sigma);
    
    sum_r = sum(r(k,:).*delta_t);
    actualization_rate = exp(-sum_r);
    summ = summ + payoff(S(k,end),K)*actualization_rate;
    
end;

price = summ/num_of_iter;

return;


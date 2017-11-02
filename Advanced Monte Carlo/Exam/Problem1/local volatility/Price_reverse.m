function [price] = Price_reverse(S0,K,r,sigma,T,N)

num_of_iter = 1000;
summ = 0;

for k = 1:num_of_iter
    S = StockSimulator_localvol(T,N,S0,K,r,sigma);
    S_end(k) = S(end);
    summ = summ + payoff(S_end(k),K)*exp(-r*T);
end;

price = summ/num_of_iter;

return;


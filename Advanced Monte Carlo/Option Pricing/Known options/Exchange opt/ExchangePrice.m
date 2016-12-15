function [price] = ExchangePrice(T,N,S10,S20,r,sigma1,sigma2,ro)

num_of_iter = 100;
sum_opt = 0;

for n = 1:num_of_iter
    path(n,:,:) = CorrStockSimulator(T,N,S10,S20,r,sigma1,sigma2,ro);
    sum_opt = sum_opt+payoff_exchange(path(n,1,end),path(n,2,end))*exp(-r*T);
    % where path(n,i,end) is a price of i-th stock on n-th iteration at time end = T
    
end;

price = sum_opt/num_of_iter;


return;
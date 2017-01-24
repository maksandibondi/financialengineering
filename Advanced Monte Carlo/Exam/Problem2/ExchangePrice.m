function [price] = ExchangePrice(T,N,S10,S20,r,sigma1,sigma2,ro)

num_of_iter = 1000;
summ = 0;

for n = 1:num_of_iter
    [path1(n,:),path2(n,:)] = CorrStockSimulator(T,N,S10,S20,r,sigma1,sigma2,ro);
    summ = summ + payoff_exchange(path1(n,end),min(path2(n,:)))*exp(-r*T);    
end;

price = summ/num_of_iter;


return;
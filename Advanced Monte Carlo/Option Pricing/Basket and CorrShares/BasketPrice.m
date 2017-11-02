function [price] = BasketPrice(T,N,K,S10,S20,r,sigma1,sigma2,lamb1,lamb2,ro)

num_of_iter = 500;
summ = 0;

for n = 1:num_of_iter
    path(n,:,:) = CorrStockSimulator(T,N,S10,S20,r,sigma1,sigma2,ro);
    summ = summ + payoff_basket(path(n,1,end),path(n,2,end),lamb1,lamb2,K)*exp(-r*T);
    % where path(n,i,end) is a price of i-th stock on n-th iteration at time end = T
    
end;

price = summ/num_of_iter;


return;
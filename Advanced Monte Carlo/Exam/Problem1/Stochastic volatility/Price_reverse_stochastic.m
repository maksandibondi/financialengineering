function [price,S,path_sigma] = Price_reverse_stochastic(S0,K,r,sigma,T,N,ro)

num_of_iter = 1000;
summ = 0;

for k = 1:num_of_iter
    [path_sigma,S] = CorrProcessSimulator(T,N,S0,r,sigma,ro);
    S_end(k) = S(end);
    summ = summ + payoff(S_end(k),K)*exp(-r*T);
end;

price = summ/num_of_iter;

return;


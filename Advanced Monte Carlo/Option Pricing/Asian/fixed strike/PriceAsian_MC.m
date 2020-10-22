function [Asian] = PriceAsian_MC(S0,T,N,r,sigma,K)

num_of_iter = 100;
Asian = 0;

for i = 1:num_of_iter
    payoff = payoff_asian(T,N,S0,r,sigma,K);
    Asian = Asian + exp(-r*T)*payoff/num_of_iter;
 
end; 

return;
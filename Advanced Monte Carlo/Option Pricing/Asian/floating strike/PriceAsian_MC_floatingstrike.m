function [Asian] = PriceAsian_MC_floatingstrike(S0,T,N,r,sigma)

num_of_iter = 500;
Asian = 0;

for i = 1:num_of_iter
    payoff = payoff_asian_floatingstrike(T,N,S0,r,sigma);
    Asian = Asian + exp(-r*T)*payoff/num_of_iter;
end; 

return;
function [C,P] = PriceEU_MC(S0,K,r,sigma,T)

num_of_iter = 2000;
sum_call = 0; sum_put = 0;


for k = 1:num_of_iter
    S_end(k) = S0*exp((r-(sigma^2)/2)*T+sigma*sqrt(T)*randn());
    sum_call = sum_call + max(S_end(k)-K,0)*exp(-r*T);
    sum_put = sum_put + max(K-S_end(k),0)*exp(-r*T);
end;

C = sum_call/num_of_iter;
P = sum_put/num_of_iter;

return;


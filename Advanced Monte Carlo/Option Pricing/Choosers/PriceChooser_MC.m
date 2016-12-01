function [C,P] = PriceChooser_MC(S0,K,r,sigma,T_deс,T)

num_of_iter = 2000;
sum_call = 0; sum_put = 0;


for k = 1:num_of_iter
    S_dec_end(k) = S0*exp((r-(sigma^2)/2)*T_des+sigma*sqrt(T_des)*randn());
    S_end(k) = S_des_end(k)*exp((r-(sigma^2)/2)*(T-T_des)+sigma*sqrt(T-T_des)*randn());
    
    sum_call = sum_call + max(S_end(k)-K,0)*exp(-r*T);
    sum_put = sum_put + max(K-S_end(k),0)*exp(-r*T);
end;

C = sum_call/num_of_iter;
P = sum_put/num_of_iter;

return;


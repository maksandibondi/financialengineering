function [C,P] = PriceEU_MC_discrete(S0,K,r,sigma,jump_freq,t0,T,discretization_num_t,v1,p1,v2,p2)

num_of_iter = 100;
sum_call = 0; sum_put = 0;

EV_Y = v1*p1+v2*p2;
drift = r - jump_freq*EV_Y;

delta_t = (T-t0)/discretization_num_t;
% definition of the t-values on axis
t(1) = t0;
for q = 2:1:discretization_num_t
    t(q) = t(q-1) + delta_t; 
end;


for k = 1:num_of_iter
    
    S(k,1) = S0;
    
    for i = 2:discretization_num_t
        
        S(k,i) = S(k,i-1)*(Amplitude_discrete(jump_freq,t(i)-t(i-1),v1,p1,v2,p2)*exp((drift-(sigma^2)/2)*delta_t+sigma*sqrt(delta_t)*randn()));
        
    end;
    
        sum_call = sum_call + max(S(k,end)-K,0)*exp(-r*T);
        sum_put = sum_put + max(K-S(k,end),0)*exp(-r*T);
    
end;

C = sum_call/num_of_iter;
P = sum_put/num_of_iter;

return;


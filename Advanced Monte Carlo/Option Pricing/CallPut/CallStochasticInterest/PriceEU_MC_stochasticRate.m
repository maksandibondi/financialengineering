function [C,P] = PriceEU_MC_stochasticRate(S0,N,K,sigma,r0,nu_r,gama_r,sigma_r,T)

num_of_iter = 2000;
sum_call = 0; sum_put = 0; 

delta_t = T/N;

for k = 1:num_of_iter
    sum_r = 0; % will be used for actualization
    path_r(k,:) = RateSimulator(T,N,r0,nu_r,gama_r,sigma_r);
    
    S(k,1) = S0;
    for i = 2:N
        S(k,i) = S(k,i-1)+S(k,i-1)*(path_r(k,i-1)*delta_t+sigma*sqrt(delta_t)*randn());
    end;
        sum_r = sum(path_r(k,:).*delta_t);
    
    actualization_rate = exp(-sum_r);
    sum_call = sum_call + max(S(k,end)-K,0)*actualization_rate;
    sum_put = sum_put + max(K-S(k,end),0)*actualization_rate;
    
end;

C = sum_call/num_of_iter;
P = sum_put/num_of_iter;

return;


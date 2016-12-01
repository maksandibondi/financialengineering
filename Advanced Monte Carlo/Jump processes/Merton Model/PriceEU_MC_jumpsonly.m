function [C,P] = PriceEU_MC_jumpsonly(S0,K,r,jump_height,jump_freq,T)

num_of_iter = 500;
sum_call = 0; sum_put = 0;
drift = r + jump_freq*(exp(jump_height)-1);

for k = 1:num_of_iter
    S_end(k) = S0*exp(drift*T+jump_height*Jumps_t(jump_freq,T));
    sum_call = sum_call + max(S_end(k)-K,0)*exp(-r*T);
    sum_put = sum_put + max(K-S_end(k),0)*exp(-r*T);
end;

C = sum_call/num_of_iter;
P = sum_put/num_of_iter;

return;


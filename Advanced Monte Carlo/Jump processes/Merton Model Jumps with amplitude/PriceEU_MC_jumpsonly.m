function [price] = PriceEU_MC_jumpsonly(S0,K,r,jump_height,jump_freq,T)

num_of_iter = 500;
summ = 0;
drift = r + jump_freq*(exp(jump_height)-1);

for k = 1:num_of_iter
    S_end(k) = S0*exp(drift*T+jump_height*Jumps_t(jump_freq,T));
    summ = summ + payoff(S_end(k),K)*exp(-r*T);
end;

price = summ/num_of_iter;

return;


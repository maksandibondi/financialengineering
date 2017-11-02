function [price] = PriceEU_MC_Merton(S0,K,r,sigma,m_y,sigma_y,sigma_j,jump_freq,T)

num_of_iter = 500;
summ = 0;
EV_Y = exp(m_y+sigma_y^2/2)-1;
drift = r - jump_freq*EV_Y;

for k = 1:num_of_iter
    % Instead of amplitude function we use just a power of EV as we need
    % only the last value of S
    S_end(k) = S0*(Amplitude(m_y,sigma_j,jump_freq,T)*exp((drift-sigma^2/2)*T+sigma*sqrt(T)*randn()));
    summ = summ + payoff(S_end(k),K)*exp(-r*T);
end;

price = summ/num_of_iter;


return;


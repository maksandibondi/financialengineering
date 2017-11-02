function [payoff] = payoff_exchange(S1T,M2T)

payoff = max(S1T-M2T,0);

return;
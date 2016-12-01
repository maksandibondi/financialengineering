function [payoff] = payoff_exchange(S1T,S2T)

payoff = max(S1T-S2T,0);

return;
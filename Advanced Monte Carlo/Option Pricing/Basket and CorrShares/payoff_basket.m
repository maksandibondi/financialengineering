function [payoff] = payoff_basket(S1T,S2T,lamb1,lamb2,K)

payoff = max(S1T*lamb1+S2T*lamb2-K,0);

return;
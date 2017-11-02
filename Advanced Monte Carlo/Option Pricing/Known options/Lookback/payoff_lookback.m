function [payoff] = payoff_lookback(T,N,S0,r,sigma,K)

    delta_t = T/N;
    Spath = StockSimulator(T,N,S0,r,sigma);

%% Payoff calculation
payoff = max(max(Spath)-K,0);

return;
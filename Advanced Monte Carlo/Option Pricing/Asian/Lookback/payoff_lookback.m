function [Call,Put] = payoff_lookback(T,N,S0,r,sigma,K)

    delta_t = T/N;
    Spath = StockSimulator(T,N,S0,r,sigma);

%% Payoff calculation
Call = max(max(Spath)-K,0);
Put = max(K-min(Spath),0);


return;
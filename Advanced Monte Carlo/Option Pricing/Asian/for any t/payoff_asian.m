function [Call,Put] = payoff_asian(At,T,t,N,S0,r,sigma,K)

    delta_t = (T-t)/N;
    Spath = StockSimulator(T,N,S0,r,sigma);

    A = 0;
   
%% mean value of the share price path
for i = 1:N
    A = A + Spath(i)*delta_t/T;  
end;

%% Payoff calculation

f = At*t+A;
Call = max(f-K,0);
Put = max(K-f,0);


return;
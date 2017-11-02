function [payoff] = payoff_asian(T,N,S0,r,sigma,K)

    delta_t = T/N;
    Spath = StockSimulator(T,N,S0,r,sigma);

    A = 0;
   
%% mean value of the share price path
for i = 1:N
    A = A + Spath(i)*delta_t/T;  
end;

%% Payoff calculation
payoff = max(A-K,0);

return;
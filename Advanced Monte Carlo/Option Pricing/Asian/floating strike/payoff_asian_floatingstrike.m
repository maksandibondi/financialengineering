function [Call,Put] = payoff_asian_floatingstrike(T,N,S0,r,sigma)

    delta_t = T/N;
    Spath = StockSimulator(T,N,S0,r,sigma);

    A = 0;
   
%% mean value of the share price path
for i = 1:N
    A = A + Spath(i)*delta_t/T;  
end;

%% Payoff calculation
Call = max(A-Spath(end),0);
Put = max(Spath(end)-A,0);


return;
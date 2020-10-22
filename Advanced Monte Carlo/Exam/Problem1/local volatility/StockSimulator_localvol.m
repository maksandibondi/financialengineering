function [S] = StockSimulator_localvol(T,N,S0,K,r,sigma0)

delta_t = T/N;

S(1) = S0;

bmpath = BMsimulator(T,N,0);

for i = 2:N
    
    if S(i-1)<K
       sigma(i-1) = sigma0*(1+S(i-1)/K);
    else
       sigma(i-1) = sigma0*2;
    end;
    
    S(i) = S(i-1) + S(i-1)*(r*delta_t+sigma(i-1)*(bmpath(i)-bmpath(i-1)));

end;
 
return;
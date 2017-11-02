function[out]=d1(D,S,t,T,K,r,sigma)
    out=(log(S/K)+(r-D+0.5*sigma^2)*(T-t))/(sigma*sqrt(T-t));
end
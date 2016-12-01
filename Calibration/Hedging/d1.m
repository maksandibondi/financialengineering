function[out]=d1(S,t,T,K,r,sigma)
    out=(log(S/K)+(r+0.5*sigma^2)*(T-t))/(sigma*sqrt(T-t));
end
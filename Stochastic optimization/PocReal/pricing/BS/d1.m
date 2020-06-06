function[out]=d1(sigma,S0,r,T,K)
    out=(log(S0/K)+(r+0.5*sigma^2)*T)/(sigma*sqrt(T));
end
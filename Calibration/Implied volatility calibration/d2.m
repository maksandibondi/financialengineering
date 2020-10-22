function[out2]=d2(sigma,S0,r,T,K)
    out2=(log(S0/K)+(r-0.5*sigma^2)*T)/(sigma*sqrt(T));
end
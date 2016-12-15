function[out2]=d2(S,t,T,K,r,sigma)
    out2=(log(S/K)+(r-0.5*sigma^2)*(T-t))/(sigma*sqrt(T-t));
end
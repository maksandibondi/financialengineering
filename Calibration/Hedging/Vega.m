function[v]=Vega(sigma,S0,r,T,K)
    v=S0*sqrt(0.5*T/pi)*exp(-0.5*d1(sigma,S0,r,T,K)^2);
end
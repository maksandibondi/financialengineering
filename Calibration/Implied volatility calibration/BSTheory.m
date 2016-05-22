function[P]=BSTheory(sigma,S0,r,T,K)
    P=S0*normal(d1(sigma,S0,r,T,K))-K*exp(-r*T)*normal(d2(sigma,S0,r,T,K));
end
function[P] = BSTheory_put(S0,t,T,K,r,sigma)
    P = S0*normal(d1(S0,t,T,K,r,sigma))-K*exp(-r*T)*normal(d2(S0,t,T,K,r,sigma)); % for call
    P = P + K*exp(-r*T) - S0; % for put
end
function[P] = BSTheory_call_div(D,S0,t,T,K,r,sigma)
    P = exp(-D*(T-t))*S0*normal(d1(D,S0,t,T,K,r,sigma))-K*exp(-r*(T-t))*normal(d2(D,S0,t,T,K,r,sigma)); % for call
    % P = P + K*exp(-r*T) - S0; % for put
end
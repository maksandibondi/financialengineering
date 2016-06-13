function[P] = BSTheory(sigma,S0,r,T,K)
    P = S0*normal(d1(sigma,S0,r,T,K))-K*exp(-r*T)*normal(d2(sigma,S0,r,T,K)); % for call
    % P = P + K*exp(-r*T) - S0; % for put
end
function [S] = GBM_ITOSimulator(S0,mu,sigma,T,delta_t)

discretization_num_t = T/delta_t;

t(1) = 0;
for i = 2:discretization_num_t
    t(i) = t(i-1)+delta_t;
end;

W = BMsimulator(T,discretization_num_t,'Polar');
S(1) = S0;
for i = 2:discretization_num_t
    S(i) = S(1)*exp((mu-sigma^2/2)*t(i)+sigma*W(i));
end;

return;
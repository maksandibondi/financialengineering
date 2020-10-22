function [R] = Ornstein_FDMSimulator(R0,a,b,sigma,T,delta_t)

discretization_num_t = T/delta_t;

t(1) = 0;
for i = 2:discretization_num_t
    t(i) = t(i-1)+delta_t;
end;

W = BMsimulator(T,discretization_num_t,'Polar');
R(1) = R0;
for i = 2:discretization_num_t
    R(i) = R(i-1) + a*(b-R(i-1))*delta_t + sigma*(W(i)-W(i-1));
end;

return;
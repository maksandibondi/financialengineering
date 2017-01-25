clear; clc;

R0 = 0.1;
mu = 0.1;
sigma = 0.5;
T = 0.5;
delta_t = 0.005;
a = 0.5;
b = 0.15;

discretization_num_t = T/delta_t;

t(1) = 0;
for i = 2:discretization_num_t
    t(i) = t(i-1)+delta_t;
end;

W = BMsimulator(T,discretization_num_t,'Reject');
R(1) = R0;
for i = 2:discretization_num_t
    
    for s = 2:i
        integral(s) = exp(a*(t(s)-t(i)))*(W(s)-W(s-1));
    end;
    R(i) = R(1)*exp(-a*t(i))+b*(1-exp(-a*t(i)))+sigma*sum(integral);
end;

plot(t,R);
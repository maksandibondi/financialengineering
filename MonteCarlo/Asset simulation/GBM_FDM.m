clear; clc;

S0 = 10;
mu = 0.1;
sigma = 0.5;
T = 0.5;
delta_t = 0.005;


discretization_num_t = T/delta_t;

t(1) = 0;
for i = 2:discretization_num_t
    t(i) = t(i-1)+delta_t;
end;

W = BMsimulator(T,discretization_num_t,'TCL');
S(1) = S0;
for i = 2:discretization_num_t
    S(i) = S(i-1) + S(i-1)*((mu-sigma^2/2)*delta_t+sigma*(W(i)-W(i-1)));
end;

plot(t,S);
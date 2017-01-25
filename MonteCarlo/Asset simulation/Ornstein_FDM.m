clear; clc;

R0 = 0.2;
a = 0.3;
b = 0.8;
sigma = 0.2;
T = 5;
delta_t = 0.005;


discretization_num_t = T/delta_t;

t(1) = 0;
for i = 2:discretization_num_t
    t(i) = t(i-1)+delta_t;
end;

W = BMsimulator(T,discretization_num_t,'Reject');
R(1) = R0;
for i = 2:discretization_num_t
    R(i) = R(i-1) + a*(b-R(i-1))*delta_t + sigma*(W(i)-W(i-1));
end;

plot(t,R);
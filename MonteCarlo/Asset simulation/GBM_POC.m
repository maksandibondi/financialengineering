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


S = GBM_ITOSimulator(S0,mu,sigma,T,delta_t);
S2 = GBM_FDMSimulator(S0,mu,sigma,T,delta_t);

plot(t,S,'r');
hold on;
plot(t,S2,'b');
xlabel('t'); ylabel('S');
title('Geometric brownian motion')
legend('Trajectory using Ito solution', 'Trajectory using FDM');


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

R = Ornstein_ITOSimulator(R0,a,b,sigma,T,delta_t);
R2 = Ornstein_FDMSimulator(R0,a,b,sigma,T,delta_t);

plot(t,R,'r');
hold on;
plot(t,R2,'b');
xlabel('t'); ylabel('R');
title('Ornstein-Ulenbeck')
legend('Trajectory using Ito solution', 'Trajectory using FDM');

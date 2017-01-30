clear; clc;

T = 10;
discretization_num_t = 500;

delta_t = T/discretization_num_t;
t(1) = 0;
for q = 2:discretization_num_t
    t(q) = t(q-1) + delta_t; 
end;

hold on;

methods = {'Random walk','Polar','Reject','TCL'};

for k = 1:size(methods,2)
  W(:,k) = BMsimulator(T,discretization_num_t,methods(k));
  plot(t,W(:,k));
end;
xlabel('t'); ylabel('W(t)');
title('Brownian trajectories')
legend('Random Walk', 'Polar', 'Reject', 'TCL');

%% We use this mean and variance caclul in case of creation of multiple trajectories for one method
% EV = mean(W(end,:))
% Var = var(W(end,:))


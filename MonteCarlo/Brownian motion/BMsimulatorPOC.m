clear; clc;

T = 0.2;
discretization_num_t = 500;

delta_t = T/discretization_num_t;
t(1) = 0;
for q = 2:discretization_num_t
    t(q) = t(q-1) + delta_t; 
end;

hold on;

for k = 1:1000
  W(:,k) = BMsimulator(T,discretization_num_t,'Reject');
  
  plot(t,W(:,k));
end;

EV = mean(W(end,:))
Var = var(W(end,:))


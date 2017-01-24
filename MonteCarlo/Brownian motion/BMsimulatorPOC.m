clear; clc;

T = 0.8;
discretization_num_t = 300;

delta_t = T/discretization_num_t;
t(1) = 0;
for q = 2:discretization_num_t
    t(q) = t(q-1) + delta_t; 
end;

for k = 1:1000
  W(:,k) = BMsimulator(T,discretization_num_t,'Polar');
end;

EV = mean(W(end,:))
Var = var(W(end,:))

hold on;
plot(t,W);

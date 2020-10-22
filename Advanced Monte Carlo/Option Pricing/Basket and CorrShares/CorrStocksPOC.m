%% Correlated stock simulation POC
clear; clc;

%% Variables intervals and discretization
T = 0.5;
discretization_num_t = 100;
S10 = 10; % initial value of stock 1
S20 = 10; % initial value of stock 1
r = 0.02;
sigma1 = 0.5;
sigma2 = 0.5;
ro = [-1, 0, 1];

%% t-vector
delta_t = T/discretization_num_t;
t(1) = 0;
for q = 2:1:discretization_num_t
    t(q) = t(q-1) + delta_t; 
end;


%% Main
    path1 = CorrStockSimulator(T,discretization_num_t,S10,S20,r,sigma1,sigma2,ro(1));
    path2 = CorrStockSimulator(T,discretization_num_t,S10,S20,r,sigma1,sigma2,ro(2));
    path3 = CorrStockSimulator(T,discretization_num_t,S10,S20,r,sigma1,sigma2,ro(3));


%% Plotting
hold on;
plot(t,path1(1,:),'b');
plot(t,path1(2,:),'b');
% plot(t,path2(1,:),'r');
% plot(t,path2(2,:),'r');
% plot(t,path3(1,:),'black');
% plot(t,path3(2,:),'black');
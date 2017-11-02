%% Variables intervals and discretization
clear; clc;

Div = 0.1;
B0 = 1;
t0 = 0;
T = 5;
r = 0.05;
sigma = 0.5;
K = 1.5;
S0 = 1;
hedge_freq = 1;

discretization_num_t = 100; % number of discretization

delta_t = (T - t0)/(discretization_num_t+1); % delta y

%% X and T arrays creation

% definition of the t-values on axis
t(1) = t0;
for q = 2:1:discretization_num_t
    t(q) = t(q-1) + delta_t; 
end;

%% Main   

[V, P, A, B, PL] = hedging_div(Div, B0, t0, T, discretization_num_t, r, sigma, K, S0,hedge_freq);

%% Calculation of D

for i = 2:discretization_num_t
D(i) = A(i) - A(i-1);
end;

%% plotting

% Plotting deltas and cash deposit
figure;
hold on;
plot(t,A,'black'); 
plot(t,B,'blue'); 
legend('A','B');

% Plotting the difference
figure;
plot(t,D,'-o');
legend('D');


% Ploting prices portfolio
figure;
hold on;
plot(t(2:end),P(2:end),'blue');
plot(t,V,'black');
legend('P','V');

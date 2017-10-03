%% Variables intervals and discretization
clear; clc;

B_call0 = 1;
B_put0 = 1;
t0 = 0;
T = 5;
r = 0.05;
sigma = 0.5;
K = 1;
S0 = 1;
hedge_freq = 10;

discretization_num_t = 100; % number of discretization

delta_t = (T - t0)/(discretization_num_t+1); % delta y

%% X and T arrays creation

% definition of the t-values on axis
t(1) = t0;
for q = 2:1:discretization_num_t
    t(q) = t(q-1) + delta_t; 
end;

%% Main   
%[V_call, V_put, P_call, P_put, A_call, A_put, B_call, B_put, PLcall, PLput] = hedging(B_call0, B_put0, t0, T, discretization_num_t, r, sigma, K, S0,hedge_freq);
[V_call, V_put, P_call, P_put, A_call, A_put, B_call, B_put, PLcall, PLput] = hedging_exact_repl(t0, T, discretization_num_t, r, sigma, K, S0,hedge_freq);

%% Calculation of D

for i = 2:discretization_num_t
D_call(i) = A_call(i) - A_call(i-1);
D_put(i) = A_put(i) - A_put(i-1);
end;

%% plotting

% Plotting deltas and cash deposit
figure;
hold on;
plot(t,A_call,'black'); % plot(t,A_put,'red');
plot(t,B_call,'blue'); % plot(t,B_put,'yellow'); % Should be the same?
legend('A call','B call');

% Plotting the difference
figure;
plot(t,D_call,'-o');
% plot(t,D_put,'black');
legend('D call');


% Ploting prices portfolio
figure;
hold on;
plot(t(2:end),P_call(2:end),'blue');
% plot(t(2:end),P_put(2:end),'yellow');
plot(t,V_call,'black');
% plot(t,V_put,'red');
legend('P call','V call');

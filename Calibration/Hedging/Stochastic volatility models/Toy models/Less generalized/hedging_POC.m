%% Variables intervals and discretization
clear; clc;

B0 = 1;

t0 = 0;
T = 5;
r = 0.05;
K = 1;
S0 = 1;
hedge_freq = 1;

discretization_num_t = 100; % number of discretization


delta_t = (T - t0)/(discretization_num_t); % delta y

%% sigma and T arrays creation

% definition of the t-values on axis
t(1) = t0;

sigma1(1) = volatility1();
% sigma1(1) = 0.3;
% sigma(1) = 0.3;

for q = 2:1:discretization_num_t
    t(q) = t(q-1) + delta_t; 
    sigma1(q) = volatility1();
    % sigma1(q) = volatility2(sigma1(q-1));
    % sigma(q) = sigma(q-1)+sigma(q-1)*sqrt(delta_t)*randn();
    % sigma(q) = sigma(q-1);
end;

figure;
plot(t,sigma1);


%% Main   
[V, P, A, B, PL] = hedging_exact_repl(t0, T, discretization_num_t, r, sigma1, K, S0,hedge_freq);

%% plotting

% Plotting deltas and cash deposit
figure;
hold on;
plot(t,A,'black'); % plot(t,A_put,'red');
plot(t,B,'blue'); % plot(t,B_put,'yellow'); % Should be the same?

legend('A','B');

% Ploting prices portfolio
figure;
hold on;
plot(t(2:end),P(2:end),'blue');
% plot(t(2:end),P_put(2:end),'yellow');
plot(t,V,'black');
% plot(t,V_put,'red');

legend('P','V');


difference = V - P;
figure;
plot(t,difference,'black');

%% Variables intervals and discretization
clear; clc;

B_call0 = 1;
B_put0 = 1;
t0 = 0;
T = 5;
r = 0.05;
sigma = 0.5;
K = 1.5;
S0 = 1;
hedge_freq = [1 2 5 10];
num_of_iter = 5000;

discretization_num_t = 100; % number of discretization

delta_t = (T - t0)/(discretization_num_t+1); % delta y

%% X and T arrays creation

% definition of the t-values on axis
t(1) = t0;
for q = 2:1:discretization_num_t
    t(q) = t(q-1) + delta_t; 
end;

%% Main  

s = size(hedge_freq,2);

for k = 1:s

    for i = 1:num_of_iter

        [V_call, V_put, P_call, P_put, A_call, A_put, B_call, B_put, PLcall(i,k), PLput(i,k)] = hedging(B_call0, B_put0, t0, T, discretization_num_t, r, sigma, K, S0, hedge_freq(k));

    end;

end;

for k = 1:s
    EV_PLcall(k) = sum(PLcall(:,k))/num_of_iter;
    Var_PLcall(k) = sum(PLcall(:,k).^2)/num_of_iter - EV_PLcall(k)^2;
end;

display(EV_PLcall);
display(Var_PLcall);

plot(hedge_freq, EV_PLcall);

figure;
plot(hedge_freq, Var_PLcall);
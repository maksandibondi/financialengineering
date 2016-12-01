%% Variables intervals and discretization
clear;clc;
S_0 = 0;  S_f = 20; %interval for x
t_0 = 0;  t_f = 0.5;  % interval for t 
r = 0.1;
sigma = 0.5;
K = 10;
discretization_num_S = 20; % number of discretization 
discretization_num_t = 70; % number of discretization

delta_S = (S_f - S_0)/(discretization_num_S); % delta x
delta_t = (t_f - t_0)/(discretization_num_t); % delta y

%% X and T arrays creation
% definition of the x-values on axis
S(1) = S_0;
for q = 2:1:discretization_num_S
    S(q) = S(q-1) + delta_S;   
end;

% definition of the t-values on axis
t(1) = t_0;
for q = 2:1:discretization_num_t
    t(q) = t(q-1) + delta_t; 
end;


N = discretization_num_t;

%% Main
for k = 1:discretization_num_t
    for i = 1:discretization_num_S
        [LookBackCall(k,i),LookBackPut(k,i)] = Price_lookback_MC(S(i),t(k),N,r,sigma,K);
        [EUCall(k,i), EUPut(k,i)] = PriceEU_MC(S(i),K,r,sigma,t(k));
    end;
end;

%% Final condition
hold on;
plot(S,LookBackCall(end,:));
plot(S,LookBackPut(end,:));
plot(S,EUCall(end,:));
plot(S,EUPut(end,:));

%% Surface
figure;
surf(S,t,LookBackCall);

%% Variables intervals and discretization
clear;clc;
S_0 = 0;  S_f = 20; %interval for x
t_0 = 0;  t_f = 5;  % interval for t 

% rate params for vasicek model
r0 = 0.1;
K = 10;
nu_r = 0.016;
gama_r = 0.2;
sigma_r = 0.02;

discretization_num_S = 20; % number of discretization 
discretization_num_t = 70; % number of discretization
N = discretization_num_t;
delta_S = (S_f - S_0)/discretization_num_S; % delta x
delta_t = (t_f - t_0)/discretization_num_t; % delta y

%% Stochastic volatility sigma array creation
sigma(1) = 0.3;
for q = 2:discretization_num_t
    sigma(q) = volatility_model(sigma(q-1)); 
end;

%% X and T arrays creation
% definition of the x-values on axis
S(1) = S_0;
for q = 2:discretization_num_S
    S(q) = S(q-1) + delta_S;   
end;

% definition of the t-values on axis
t(1) = t_0;
for q = 2:discretization_num_t
    t(q) = t(q-1) + delta_t; 
end;

%% Main
for k = 1:discretization_num_t
    for i = 1:discretization_num_S
        price(i,k) = PriceEU_MC_stochasticRate(S(i),N,K,sigma,r0,nu_r,gama_r,sigma_r,t(k));
        price_det(i,k) = PriceEU_MC(S(i),K,r0,sigma(1),t(k));
    end;
end;

%% Final condition
hold on;
plot(S,price(:,end),'red');
plot(S,price_det(:,end),'black');
legend('stochastic','deterministic');

%% Surface
figure;
surf(t,S,price);
figure;
surf(t,S,price_det);


%% Variables intervals and discretization
clear; clc;

S1_0 = 0;  S1_f = 20; %interval for x1
S2_0 = 0;  S2_f = 20; %interval for x1
t_0 = 0;  t_f = 0.5;  % interval for t
T = t_f - t_0;
r = 0.1;
sigma1 = 0.5;
sigma2 = 0.4;
ro = -1;
discretization_num_S = 40; % number of discretization 
discretization_num_t = 25; % number of discretization

delta_S1 = (S1_f - S1_0)/(discretization_num_S); % delta x
delta_S2 = (S2_f - S2_0)/(discretization_num_S); % delta x
delta_t = (t_f - t_0)/(discretization_num_t); % delta y

%% X and T arrays creation
% definition of the x-values on axis
S1(1) = S1_0;
S2(1) = S2_0;
for q = 2:discretization_num_S
    S1(q) = S1(q-1) + delta_S1;  
    S2(q) = S2(q-1) + delta_S2;
end;

% definition of the t-values on axis
t(1) = t_0;
for q = 2:discretization_num_t
    t(q) = t(q-1) + delta_t; 
end;

%% Main
for k = 1:discretization_num_S
    for i = 1:discretization_num_S
        price(k,i) = ExchangePrice(T,discretization_num_t,S1(k),S2(i),r,sigma1,sigma2,ro);
    end;
end;

%% Surface
figure;
surf(S1,S2,price);


%% Variables intervals and discretization
clear; clc;
S_0 = 0;  S_f = 20; %interval for x
t_0 = 0;  t_f = 0.5;  % interval for t 
r = 0.4;
sigma = 0.5;

T = t_f;
jump_freq = 0.5;
K = 10;
discretization_num_S = 20; % number of discretization 
discretization_num_t = 100; % number of discretization
v1 = 0.8; p1 = 0.8; v2 = -0.7; p2 = 0.2;


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

%% Main
for k = 1:discretization_num_t
    for i = 1:discretization_num_S
        price_Discrete(k,i) = PriceEU_MC_discrete(S(i),K,r,sigma,jump_freq,t(k),T,discretization_num_t,v1,p1,v2,p2);
        price(k,i) = PriceEU_MC(S(i),K,r,sigma,t(k));
    end;
end;

%% Final condition
plot(S,price_Discrete(end,:),'black');
hold on;
plot(S,price(end,:),'blue');

%% Surface
figure;
surf(S,t,price_Discrete);

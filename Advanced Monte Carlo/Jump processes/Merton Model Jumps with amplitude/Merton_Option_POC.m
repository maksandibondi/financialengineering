%% Variables intervals and discretization
clear; clc;
S_0 = 0;  S_f = 20; %interval for x
t_0 = 0;  t_f = 1;  % interval for t 
r = 0.1;
sigma = 0.5;

jump_freq = 0.5;
m_y = 0.1; % expected value of variable Y s.t 1+Y = exp(j(i))
sigma_y = 0.3;
sigma_j = 0.5; % sigma of j s.t j(i) = m_j + sigma_j*randn()
K = 10;
discretization_num_S = 99; % number of discretization 
discretization_num_t = 100; % number of discretization



delta_S = (S_f - S_0)/(discretization_num_S); % delta x
delta_t = (t_f - t_0)/(discretization_num_t); % delta y

%% For other models (simple BS, jumps of constant height)
jump_height = 0.05;


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
        price_Merton(k,i) = PriceEU_MC_Merton(S(i),K,r,sigma,m_y,sigma_y,sigma_j,jump_freq,t(k));
        price_jumpsonly(k,i) = PriceEU_MC_jumpsonly(S(i),K,r,jump_height,jump_freq,t(k));
        price(k,i) = PriceEU_MC(S(i),K,r,sigma,t(k));
    end;
end;

%% Final condition
plot(S,price_Merton(end,:),'black');
hold on;
plot(S,price_jumpsonly(end,:),'r');
plot(S,price(end,:),'blue');

%% Surface
figure;
surf(S,t,price_Merton);

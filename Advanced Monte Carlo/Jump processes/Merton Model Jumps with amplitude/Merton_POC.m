%% Variables intervals and discretization
clear; clc;
t_0 = 0;  t_f = 1;  % interval for t

r = 0.1;
sigma = 0.5;
jump_freq = 0.5;
m_y = 0.1;
sigma_j = 0.5;
K = 10;
S0 = 10;

discretization_num_t = 100; % number of discretization

delta_t = (t_f - t_0)/(discretization_num_t); % delta y

%% X and T arrays creation

% definition of the t-values on axis
t(1) = t_0;
for q = 2:1:discretization_num_t
    t(q) = t(q-1) + delta_t; 
end;

%% Main

num_of_iter = 10;
EV_Y = exp(m_y+sigma_j^2/2);
drift = r - jump_freq*EV_Y;


for k = 1:num_of_iter
    S(k,1) = S0;
    
    for i = 2:discretization_num_t
    
    S(k,i) = S(k,i-1)*Amplitude(m_y,sigma_j,jump_freq,delta_t)*exp((drift-sigma^2/2)*delta_t+sigma*sqrt(delta_t)*randn());
    
    end;
    
    hold on;
    plot(t,S(k,:));
    
end;









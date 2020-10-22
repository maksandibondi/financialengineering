%% Variables intervals and discretization
clear; clc;
S0 = 10;
t_0 = 0;  t_f = 0.5;  % interval for t 
r = 0.4;
sigma = 0.5;

T = t_f;
jump_freq = 1;

K = 10;
discretization_num_S = 99; % number of discretization 
discretization_num_t = 100; % number of discretization
v1 = 0.5; p1 = 0.8; v2 = -0.7; p2 = 0.2;
EV_Y = v1*p1+v2*p2;
drift = r - jump_freq*EV_Y;

delta_t = (t_f - t_0)/(discretization_num_t); % delta y

num_of_iter = 100;

%% X and T arrays creation

% definition of the t-values on axis
t(1) = t_0;
for q = 2:1:discretization_num_t
    t(q) = t(q-1) + delta_t; 
end;

%% Main
for k = 1:num_of_iter
    
    S(k,1) = S0;
    
    for i = 2:discretization_num_t
        
        
        S(k,i) = S(k,i-1)*(Amplitude_discrete(jump_freq,t(i)-t(i-1),v1,p1,v2,p2)*exp((drift-(sigma^2)/2)*delta_t+sigma*sqrt(delta_t)*randn()));
        
    end;
    
    plot (t,S(k,:));
    hold on;
    
end;

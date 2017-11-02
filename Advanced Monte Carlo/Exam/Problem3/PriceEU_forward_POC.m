%% Variables intervals and discretization
clear; clc;
S_0 = 0;  S_f = 40; %interval for x
t_0 = 0;  t_f = 0.5;  % interval for t 
T = t_f;
T1 = 0.25;
alpha = 0.8;
r = 0.1;
sigma = 0.3;
discretization_num_S = 20; % number of discretization 
discretization_num_t = 100; % number of discretization

delta_S = (S_f - S_0)/discretization_num_S; % delta x
delta_t = (t_f - t_0)/discretization_num_t; % delta y

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
% we go only to discretization_num_t/2 because t(k) should be smaller than T1
for k = 1:discretization_num_t/2
    for i = 1:discretization_num_S
        price(k,i) = PriceEU_forward(S(i),r,sigma,T-t(k),T1,alpha);
    end;
end;

%% Task1: price at t=0 for S0 = 10
% we take price(1,6) because S0 = 10 is 6-th element of S array
option_price = price(1,6)

%% Task2
hold on;
plot(S,price(1,:),'b');
plot(S,price(end,:),'r');
xlabel('S0'); ylabel('price');
legend('price at t = 0','price at t = T1');

%% Surface
figure;
surf(S,t(1:50),price);
xlabel('S'); ylabel('t'); zlabel('price');

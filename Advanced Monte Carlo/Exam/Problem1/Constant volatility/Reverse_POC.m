%% Variables intervals and discretization
clear; clc;
K = 10;
S_0 = 0;  S_f = 3*K; %interval for x
t_0 = 0;  t_f = 0.5;  % interval for t 
T = t_f;
r = 0.4;
sigma = 0.5;
discretization_num_S = 30; % number of discretization 
discretization_num_t = 100; % number of discretization

delta_S = (S_f - S_0)/discretization_num_S; % delta x
delta_t = (t_f - t_0)/discretization_num_t; % delta y

%% S and T arrays creation
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
        price(k,i) = Price_reverse(S(i),K,r,sigma,T-t(k));
    end;
end;

%% Task1
Option_price = price(1,21)

%% Task2
hold on;
plot(S,price(1,:));
plot(S,price(end,:));
xlabel('S'); ylabel('price at t = 0');

%% Task3
figure;
surf(S,t,price);
xlabel('S'); ylabel('t'); zlabel('price');

%% Task4
figure;
plot(S,price(51,:));
xlabel('S'); ylabel('price at t = T/2');

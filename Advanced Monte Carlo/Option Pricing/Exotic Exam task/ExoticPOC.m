%% Variables intervals and discretization
clear; clc;
S1_0 = 0;  S1_f = 2; %interval for x
S2_0 = 0;  S2_f = 2; %interval for x2
t_0 = 0;  t_f = 1;  % interval for t 
r1 = 0;
r2 = 0;
sigma2 = 1;
sigma1 = 1;
K = 1;
discretization_num_S = 40; % number of discretization 
discretization_num_t = 100; % number of discretization

delta_S1 = (S1_f - S1_0)/discretization_num_S; % delta x
delta_S2 = (S2_f - S2_0)/discretization_num_S; % delta x
delta_t = (t_f - t_0)/discretization_num_t; % delta y

%% X and T arrays creation
% definition of the x-values on axis
S1(1) = S1_0;
for q = 2:1:discretization_num_S
    S1(q) = S1(q-1) + delta_S1;   
end;
S2(1) = S2_0;
for q = 2:1:discretization_num_S
    S2(q) = S2(q-1) + delta_S2;   
end;

% definition of the t-values on axis
t(1) = t_0;
for q = 2:1:discretization_num_t
    t(q) = t(q-1) + delta_t; 
end;

%% Main
for k = 1:discretization_num_S
    for i = 1:discretization_num_S
        price(k,i) = PriceEU_MC_EXOTIC(S1(i),S2(k),K,r1,r2,sigma1,sigma2,t_f);
    end;
end;

%% Surface
figure;
surf(S1,S2,price);
xlabel('S1'); ylabel('S2'); zlabel('price');
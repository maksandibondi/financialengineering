%% Variables intervals and discretization
clear; clc;

S1_0 = 40;  % S1_f = 20; %interval for x1
S2_0 = 50;  % S2_f = 20; %interval for x1
t_0 = 0;  t_f = 1.5;  % interval for t
T = t_f - t_0;
r = 0.1;
sigma1 = 0.2;
sigma2 = 0.1;
K = 80;
lamb1 = 1;
lamb2 = 1;
ro_0 = -1; ro_f = 1; 

discretization_num_ro = 100;
discretization_num_t = 100; % number of discretization

delta_t = (t_f - t_0)/discretization_num_t; % delta y
delta_ro = (ro_f - ro_0)/discretization_num_ro;

%% X and T arrays creation

ro(1) = ro_0;
for q = 2:1:discretization_num_ro
    ro(q) = ro(q-1) + delta_ro; 
end;


% definition of the t-values on axis
t(1) = t_0;
for q = 2:1:discretization_num_t
    t(q) = t(q-1) + delta_t; 
end;




%% Main

for j = 1:discretization_num_ro
        price(j) = BasketPrice(T,discretization_num_t,K,S1_0,S2_0,r,sigma1,sigma2,lamb1,lamb2,ro(j));
end;

%% Final condition
figure;
plot(ro,price);

% hold on;
% plot(S2,price(1,:));
% 
% %% Surface
% figure;
% surf(S1,S2,price);
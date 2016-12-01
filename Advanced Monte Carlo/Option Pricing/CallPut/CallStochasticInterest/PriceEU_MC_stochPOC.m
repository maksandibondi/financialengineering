%% Variables intervals and discretization
clear;clc;
S_0 = 0;  S_f = 20; %interval for x
t_0 = 0;  t_f = 5;  % interval for t 
r0 = 0.1;
r = 0.1;
sigma = 0.5;
K = 10;
nu_r = 0.016;
gama_r = 0.2;
sigma_r = 0.02;
discretization_num_S = 100; % number of discretization 
discretization_num_t = 100; % number of discretization

delta_S = (S_f - S_0)/discretization_num_S; % delta x
delta_t = (t_f - t_0)/discretization_num_t; % delta y


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
% for k = 1:discretization_num_t
    for i = 1:discretization_num_S
        
        % [priceC(k,i),priceP(k,i)] = PriceEU_MC_stochasticRate(S(i),discretization_num_t,K,sigma,r0,nu_r,gama_r,sigma_r,t(k));
        % [priceC_det(k,i),priceP_det(k,i)] = PriceEU_MC(S(i),K,r,sigma,t(k));
        [priceC(i),priceP(i)] = PriceEU_MC_stochasticRate(S(i),discretization_num_t,K,sigma,r0,nu_r,gama_r,sigma_r,t_f);
        [priceC_det(i),priceP_det(i)] = PriceEU_MC(S(i),K,r,sigma,t_f);
    end;
    
% end;

%% Final condition
hold on;
plot(S,priceC,'red');
plot(S,priceC_det,'black');
hold on;
plot(S,priceP,'red');
plot(S,priceP_det,'black');
legend('stochastic','deterministic');

%% Surface
% figure;
% surf(S,t,priceC);
% surf(S,t,priceC_det);
% figure;
% surf(S,t,priceP);
% surf(S,t,priceP_det);
% 

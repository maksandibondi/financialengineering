%% Variables intervals and discretization
S_0 = 0;  S_f = 20; %interval for x
t_0 = 0;  t_f = 0.5;  % interval for t 
r = 0.1;
sigma = 0.5;
K = 10;
discretization_num_S = 99; % number of discretization 
discretization_num_t = 100; % number of discretization

delta_S = (S_f - S_0)/(discretization_num_S+1); % delta x
delta_t = (t_f - t_0)/(discretization_num_t+1); % delta y

%% X and T arrays creation
% definition of the x-values on axis
S(1) = S_0;
for q = 2:1:discretization_num_S + 2
    S(q) = S(q-1) + delta_S;   
end;

% definition of the t-values on axis
t(1) = t_0;
for q = 2:1:discretization_num_t + 2
    t(q) = t(q-1) + delta_t; 
end;

%% Main
for k = 1:discretization_num_t+2
    for i = 1:discretization_num_S+2
        [priceC(k,i),priceP(k,i)] = PriceEU_MC(S(i),K,r,sigma,t(k));
    end;
end;

%% Final condition
plot(S,priceC(end,:));
hold on;
plot(S,priceP(end,:));

%% Surface
figure;
surf(S,t,priceC);
figure;
surf(S,t,priceP);


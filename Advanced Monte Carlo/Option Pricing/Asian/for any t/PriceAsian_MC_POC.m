%% Variables intervals and discretization
S_0 = 0;  S_f = 20; %interval for x
At_0 = 0;  At_f = 30;  % interval for t 
T = 0.5;
t = 0.2;

r = 0.1;
sigma = 0.5;
K = 10;
discretization_num_S = 20; % number of discretization 
discretization_num_At = 70; % number of discretization

delta_S = (S_f - S_0)/discretization_num_S; % delta x
delta_At = (At_f - At_0)/discretization_num_At; % delta y

%% X and T arrays creation
% definition of the x-values on axis
S(1) = S_0;
for q = 2:1:discretization_num_S
    S(q) = S(q-1) + delta_S;   
end;

% definition of the t-values on axis
At(1) = At_0;
for q = 2:1:discretization_num_t
    At(q) = At(q-1) + delta_At; 
end;


N = discretization_num_At;

%% Main
for k = 1:discretization_num_At
    for i = 1:discretization_num_S
        [AsianCall(k,i),AsianPut(k,i)] = PriceAsian_MC(S(i),At(k),t,T,N,r,sigma,K);
        [EUCall(k,i), EUPut(k,i)] = PriceEU_MC(S(i),K,r,sigma,t);
    end;
end;

%% Final condition
hold on;
plot(S,AsianCall(end,:));
plot(S,AsianPut(end,:));
plot(S,EUCall(end,:));
plot(S,EUPut(end,:));

%% Surface
figure;
surf(S,At,AsianCall);

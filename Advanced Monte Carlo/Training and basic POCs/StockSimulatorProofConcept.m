clc; clear;

%% Initial data
N = 100;
T = 0.5;
r = 0.1;
sigma = 0.2;
S0 = 10;
num_of_iter = 5000;

%% t-vector creation
delta_t = T/N;
t_vector(1) = 0;
for k = 2:N
    t_vector(k) = t_vector(k-1)+delta_t;
end;

%% Multiple path creation
for i = 1:num_of_iter
    path(:,i) = StockSimulator(T,N,S0,r,sigma);
end;

%% EV and Var calculation
EV = sum(path(N,:))/num_of_iter;
display(EV);
EVtheory = S0*exp(r*T);
display(EVtheory);
Var = sum(path(N,:).^2)/num_of_iter - EV^2;
display(Var);

%% Plotting
for i = 1:num_of_iter 
    plot(t_vector,path(:,i));
    hold on;
end;

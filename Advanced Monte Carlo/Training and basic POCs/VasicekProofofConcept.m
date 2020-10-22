clear; clc;

%% Initial data
r0 = 0.1;
nu = 0.016;
gama = 0.2;
omega = 0.02;
T = 5;
N = 100;
num_of_iter = 5000;

%% t-vector creation
delta_t = T/N;
t_vector(1) = 0;
for k = 2:N
    t_vector(k) = t_vector(k-1)+delta_t;
end;

%% Multiple path creation
for i = 1:num_of_iter
    path(:,i) = RateSimulator(T,N,r0,nu,gama,omega);
end;

%% EV and Var calculation
EV = sum(path(N,:))/num_of_iter;
display(EV);
Var = sum(path(N,:).^2)/num_of_iter - EV^2;
display(Var);

%% Plotting
for i = 1:num_of_iter 
    plot(t_vector,path(:,i));
    hold on;
end;

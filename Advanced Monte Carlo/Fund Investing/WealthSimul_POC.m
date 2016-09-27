clc; clear;

%% Initial data
N = 100;
T = 5;
Ro = 0.6;
sigma2 = 0.2;
X0 = 5;
r = 0.1;
nu = 0.0001;
lambda = 0.0005;
M0 = 1;
sigma1 = 0.2;
num_of_iter = 6000;

%% t-vector creation
delta_t = T/N;
t_vector(1) = 0;
for k = 2:N
    t_vector(k) = t_vector(k-1)+delta_t;
end;

%% Multiple path creation
for i = 1:num_of_iter
    path(:,i) = WealthSimul(T,N,r,nu,lambda,sigma1,M0,Ro,sigma2,X0);
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
clear; clc;

%% Initial data
T = 5;
N = 100;
num_of_iter = 1000;
initial_point = 0;

%% t-vector creation
delta_t = T/N;
t_vector(1) = 0;
for k = 2:N
    t_vector(k) = t_vector(k-1)+delta_t;
end;

%% Multiple path creation
for i = 1:num_of_iter
      path(:,i) = BMsimulator(T,N,initial_point);
end;

%% EV and Var calculation
EV = sum(path(N,:))/num_of_iter;
Var = sum(path(N,:).^2)/num_of_iter - EV^2;
display(EV);
display(Var);

%% Plotting
for i = 1:num_of_iter 
    plot(t_vector,path(:,i));
    hold on;
end;

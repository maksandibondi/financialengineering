clear; clc;
%% Variables intervals and discretization
x_0 = 0;  x_f = 20; %interval for x
t_0 = 0;  t_f = 0.5;  % interval for t 
r = 0.1;
sigma = 0.5;
Strike = 10;
discretization_num_x = 99; % number of discretization 
discretization_num_t = 4999; % number of discretization
delta_x = (x_f - x_0)/(discretization_num_x+1); % delta x
delta_t = (t_f - t_0)/(discretization_num_t+1); % delta y

%% X and T arrays creation
% definition of the x-values on axis
x(1) = x_0;
for q = 2:1:discretization_num_x + 2
    x(q) = x(q-1) + delta_x;   
end;

% definition of the t-values on axis
t(1) = t_0;
for q = 2:1:discretization_num_t + 2
    t(q) = t(q-1) + delta_t; 
end;

%% Final condition realization and plotting
%u(discretization_num_t + 2,1) = initialcond_bs(x(1),Strike); % by initial condition and knowing first valuse of x we can find the first value of u
for i = 1:discretization_num_x + 2;
    u(discretization_num_t + 2,i) = initialcond_bs(x(i), Strike); % g(x) function. Go through all values of u for given t = t_0
end;
plot(x, u(discretization_num_t + 2,:), 'r');

%% Dirichlet boundary condition realization
%u(1,discretization_num_x + 2) = boundaryfinal_bs(t(1), Strike, r, t_f, x_f); % by boundary condition and knowing first value of t we can find the last value of u
for n = 1:discretization_num_t + 1;
   u(n,1) = boundarystart_bs(t(n), Strike, r, t_f, x_f);  % Go through all values of u for given x = x_0
   u(n,discretization_num_x + 2) = boundaryfinal_bs(t(n), Strike, r, t_f, x_f); % Go through all values of u for given x = x_f 
end;

figure; plot(t,u(:,1),'r'); 
figure; plot(t,u(:,discretization_num_x + 2), 'b');

%% Main program
% till this moment we have defined u(1,:), u(:,1), u(:,end)
for n = discretization_num_t + 2: -1: 2
    for i = 2:discretization_num_x + 1
        u(n-1,i) = u(n,i+1)*(delta_t/2)*((sigma^2)*(x(i)^2)/(delta_x)^2 + r*x(i)/delta_x)+u(n,i)*(1-delta_t*((sigma^2)*x(i)^2/(delta_x)^2+r))+u(n,i-1)*(delta_t/2)*((sigma^2)*(x(i))^2/(delta_x)^2-r*x(i)/delta_x);
    end;
end;
 
%% Grid resizing and Plotting the solution
figure;
plot(x,u(1,:));
display(u(discretization_num_t+2,:));
display(u(1,:));

figure;
surf(x,t,u);

%plotting and normalizing vectors
for j=1:101                      
    W(j,:) = u(50*(j-1)+1,:);
    tnew(j) = t(50*(j-1)+1);
end
%tnew = (1:101)*t/(101);
figure;
surfl(x,tnew,W,'light'); xlabel('variable x'); ylabel('variable t'); zlabel('function u');
title('solution eqation 3 part4');   

 %display(showsolutionInPoint(0.5, 0.5, u, x_0, delta_x, t_0, delta_t));
 
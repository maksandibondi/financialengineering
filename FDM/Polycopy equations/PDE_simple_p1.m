clear; clc;
%% Variables intervals and discretization
x_0 = 0;  x_f = 1; %interval for x
t_0 = 0;  t_f = 1;  % interval for t 
discretization_num_x = 20; % number of discretization 
discretization_num_t = 1000; % number of discretization
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

%% Initial condition realization and plotting
u(1,1) = initialcond(x(1)); % by initial condition and knowing first valuse of x we can find the first value of u
for i= 2:discretization_num_x + 2;
    u(1,i) = initialcond(x(i-1)); % g(x) function. Go through all values of u for given t = t_0
end;
plot(x, u(1,:), 'r');

%% Dirichlet boundary conditions realization and plotting
u(1,discretization_num_x + 2) = boundaryfinal(t(1)); % by boundary condition and knowing first value of t we can find the last value of u
for n = 2:discretization_num_t + 2;
   u(n,1) = boundarystart(t(n));  % Go through all values of u for given x = x_0
   u(n,discretization_num_x + 2) = boundaryfinal(t(n)); % Go through all values of u for given x = x_f 
end;

figure; plot(t,u(:,1),'r'); 
figure; plot(t,u(:,discretization_num_x + 2), 'b');
display (u(:,discretization_num_x + 2));

%% Main program
% till this moment we have defined u(1,:), u(:,1), u(:,end)
for n = 1:discretization_num_t + 1
    for i = 2:discretization_num_x + 1
        u(n+1,i) = u(n,i) + (delta_t/(delta_x)^2)*(u(n,i+1) + u(n,i-1) - 2*u(n,i)) - 2*delta_t;
    end;
end;
 
%% Grid resizing and Plotting the solution
for j=1:101                            
    W(j,:)=u(10*(j-1)+1,:);
    tnew(j) = t(10*(j-1)+1);
end
figure;
surfl(x,tnew,W,'light'); xlabel('variable x'); ylabel('variable t'); zlabel('function u'); 
title('solution eqation 1 part1');       
        
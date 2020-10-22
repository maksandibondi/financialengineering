clear; clc;
%% Variables intervals and discretization
x_0 = 0;  x_f = 20; %interval for x
t_0 = 0;  t_f = 0.5;  % interval for t 
r = 0.2;
sigma = 0.9;
Strike = 10;
discretization_num_x = 100; % number of discretization 
discretization_num_t = 5000; % number of discretization
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
% Final bs-condition definition
%u(discretization_num_t + 2,1) = initialcond_bs(x(1),Strike); % by initial condition and knowing first valuse of x we can find the first value of u
for i = 1:discretization_num_x + 2;
    u(discretization_num_t + 2,i) = initialcond_bs(x(i), Strike); % g(x) function. Go through all values of u for given t = t_0
end;
plot(x, u(discretization_num_t + 2,:), 'r');

display (u(:,discretization_num_x + 2));

%% Main program
% till this moment we have defined u(1,:), u(:,1), u(:,end) (defining all
% the surface including boundary conditions implementation - discretized)
% To change to put we need just to change x(i)-Strike to the Strike-x(i)
for n = discretization_num_t + 2: -1: 2
    for i = 2:discretization_num_x + 1
        u(n-1,i) = max(u(n,i+1)*(delta_t/2)*((sigma^2)*(x(i))^2/(delta_x)^2+r*x(i)/delta_x) + u(n,i)*(1-delta_t*((sigma^2)*(x(i))^2/(delta_x)^2+r)) + u(n,i-1)*(delta_t/2)*((sigma^2)*(x(i))^2/(delta_x)^2-r*x(i)/delta_x), max(x(i)-Strike,0));
    end;
%% Neumann boundary conditions realization
        u(n-1,1) = u(n-1,2) - boundarystart_bs(t(n), Strike, r, t_f, x_f)*delta_x; % Go through all values of u for given x = x_0
        u(n-1,discretization_num_x + 2) =  u(n-1,discretization_num_x + 1) + boundaryfinal_bs(t(n), Strike, r, t_f, x_f)*delta_x; % Go through all values of u for given x = x_f 
end;
  
%% Grid resizing and Plotting the solution and boundary conditions
figure; plot(t,u(:,1),'r');  xlabel('variable t'); ylabel('variable u');
title('solution eqation 3 part 4 T = 0.5 boundary condition x = 1'); 

figure; plot(t,u(:,discretization_num_x + 2), 'b'); xlabel('variable t'); ylabel('variable u');
title('solution eqation 3 part 4 T = 0.5 boundary condition x = 22');
display (u(:,discretization_num_x + 2));

%plotting and normalizing vectors
for j=1:101                      
    W(j,:) = u(50*(j-1)+1,:);
    tnew(j) = t(50*(j-1)+1);
end;


%tnew = (1:101)*t/(101);
figure;
plot(x,W(101,:));
hold on;
plot (x, u(discretization_num_t + 2,:), 'r');

figure;
surfl(x,tnew,W,'light'); xlabel('variable x'); ylabel('variable t'); zlabel('function u');
title('solution eqation 3 part4');  


 %display(showsolutionInPoint(0.5, 0.5, u, x_0, delta_x, t_0, delta_t));
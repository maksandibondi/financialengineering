clear; clc;
%% Variables intervals and discretization
x_0 = 0;  x_f = 1; %interval for x
t_0 = 0;  t_f = 10;  % interval for t 
discretization_num_x = 20; % number of discretization 
discretization_num_t = 10000; % number of discretization
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
figure;
plot(x, u(1,:), 'r'); xlabel('variable x'); ylabel('variable u');
title('solution eqation 2 part 5 T = 10 initial condition t = 1'); 

%% Main program
% till this moment we have defined u(1,:), u(:,1), u(:,end)
for n = 1:discretization_num_t + 1
    for i = 2:discretization_num_x + 1
        u(n+1,i) = u(n,i) + (delta_t/(delta_x)^2)*(u(n,i+1) + u(n,i-1) - 2*u(n,i)) - delta_t*(1+t(n)*u(n,i)) + (delta_t/delta_x)*x(i)*(u(n,i+1)-u(n,i)) - (x(i)^2)*delta_t/2;
    end;
%% Neumann boundary conditions realization and plotting
        u(n+1, 1) = u(n+1,2) - boundarystart(t(n))*delta_x;   % newmann condition
        u(n+1, discretization_num_x+2) = u(n+1, discretization_num_x+1) + boundaryfinal(t(n))*delta_x; % newmann condition
    
end;

%% Grid resizing and Plotting the solution and boundary conditions
figure; plot(t,u(:,1),'r');  xlabel('variable t'); ylabel('variable u');
title('solution eqation 2 part 5 T = 10 boundary condition x = 1'); 

figure; plot(t,u(:,discretization_num_x + 2), 'b'); xlabel('variable t'); ylabel('variable u');
title('solution eqation 2 part 5 T = 10 boundary condition x = 22');
display (u(:,discretization_num_x + 2));
   
for j=1:101                      
    W(j,:)=u(10*(j-1)+1,:);
    tnew(j) = t(10*(j-1)+1);
end
%tnew = (1:101)*t/(101);
figure;
surfl(x,tnew,W,'light'); xlabel('variable x'); ylabel('variable t'); zlabel('function u');
title('solution eqation 2 part 5 T = 10');     

display(showsolutionInPoint(0.5, 0.5, u, x_0, delta_x, t_0, delta_t));
 
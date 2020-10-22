clear; clc;
%% Variables intervals and discretization
x_0 = 0;  x_f = 4; %interval for x
t_0 = 0;  t_f = 10^(-2);  % interval for t 
discretization_num_x = 100; % number of discretization 
discretization_num_t = 9999; % number of discretization
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
for i = 1:discretization_num_x + 2; % Initial condition definition
    u(1,i) = initialcond_kortoweg(x(i)); 
end;
% Plotting the initial condition
plot(x, u(1,:), 'b'); xlabel('variable x'); ylabel('variable u');
title('solution eqation kortoweg initial condition t = 1'); 

%% Definition of the second level
for i = 3:discretization_num_x % Second level definition
    u(2,i) = u(1,i)+delta_t*(-(u(1,i+1)+u(1,i)+u(1,i-1))*(u(1,i+1)-u(1,i-1))/delta_x-(u(1,i+2)-2*u(1,i+1)+2*u(1,i-1)-u(1,i-2))/(2*(delta_x)^3));
end;

%% Periodic boundary conditions and plotting
u(2,1) = u(2,discretization_num_x-1); 
u(2,2) = u(2,discretization_num_x);
u(2,discretization_num_x+1) = u(2,3);
u(2,discretization_num_x+2) = u(2,4);
% Plotting second level 
hold on;
plot(x, u(2,:),'r'); xlabel('variable x'); ylabel('variable u');
title('solution eqation kortoweg initial condition t = 2'); 

%% Main program
for n = 2:discretization_num_t + 1 % main program to calculate all levels valking through all t and x
    for i = 3:discretization_num_x
       u(n+1,i) = u(n-1,i)+2*delta_t*(-(u(n,i+1)+u(n,i)+u(n,i-1))*(u(n,i+1)-u(n,i-1))/delta_x-(u(n,i+2)-2*u(n,i+1)+2*u(n,i-1)-u(n,i-2))/(2*(delta_x)^3));
    end;
    
%% Periodic boundary conditions realization
    u(n+1,1) = u(n+1,discretization_num_x-1);
    u(n+1,2) = u(n+1,discretization_num_x);
    u(n+1,discretization_num_x+1) = u(n+1,3);
    u(n+1,discretization_num_x+2) = u(n+1,4);
end;

%% Grid resizing and Plotting the solution and boundary conditions
% Plotting boundary values (x = 1, x = discretization_num_t+2)
figure; plot(t,u(:,1),'r');  xlabel('variable t'); ylabel('variable u');
title('solution eqation kortoweg boundary condition x = 1'); 

figure; plot(t,u(:,discretization_num_x + 2), 'b'); xlabel('variable t'); ylabel('variable u');
title('solution eqation kortoweg boundary condition x = end');

%% Movie
figure;
for n = 1:50:discretization_num_t+2
    plot(x,u(n,:));
   axis([0 4 0 300]);
    MD(n) = getframe;
end;
figure;
movie(MD);
title('solution eqation kortoweg movie');
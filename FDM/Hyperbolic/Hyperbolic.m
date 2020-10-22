clear; clc;
%% Variables intervals and discretization
x_0 = 0;  x_f = 1; %interval for x
t_0 = 0;  t_f = 1;  % interval for t 
discretization_num_x = 100; % number of discretization 
discretization_num_t = 9000; % number of discretization
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
% Final bs-condition definition
%u(discretization_num_t + 2,1) = initialcond_bs(x(1),Strike); % by initial condition and knowing first valuse of x we can find the first value of u
for i = 1:discretization_num_x + 2;
    u(1,i) = initialcond_hyper(x(i)); % g(x) function. Go through all values of u for given t = t_0
    u(2,i) = u(1,i)+delta_t*Velocity(x(i));
end;
plot(x, u(1,:), 'r');

%% Non-periodic boundary conditions realization
for n = 1:discretization_num_t + 2;
   u(n,1) = boundarystart_hyper(t(n));  % Go through all values of u for given x = x_0
   u(n,discretization_num_x + 2) = boundaryfinal_hyper(t(n)); % Go through all values of u for given x = x_f 
end;

%% Main program
% till this moment we have defined u(1,:), u(:,1), u(:,end) (defining all
% the surface including boundary conditions implementation - discretized)
% To change to put we need just to change x(i)-Strike to the Strike-x(i)
for n = 2:1:discretization_num_t + 1
    for i = 2:discretization_num_x + 1
       u(n+1,i) = 2*u(n,i)-u(n-1,i)+((delta_t/delta_x)^2)*(u(n,i+1)+u(n,i-1)-2*u(n,i));
    end;
%% periodic boundary conditions realization
%     u(n,1) = u(n,discretization_num_x+1); % periodic boundarystart condition
%     u(n,discretization_num_x+2) = u(n,2); % periodic boundaryfinal condition
end;
  
%% Grid resizing and Plotting the solution and boundary conditions
figure; plot(t,u(:,1),'r');  xlabel('variable t'); ylabel('variable u');
title('solution eqation 3 part 4 T = 0.5 boundary condition x = 1'); 

figure; plot(t,u(:,discretization_num_x + 2), 'b'); xlabel('variable t'); ylabel('variable u');
title('solution eqation 3 part 4 T = 0.5 boundary condition x = 22');

%% Movie
for n = 1:50:discretization_num_t+2
    plot(x,u(n,:));
    axis([0 1 -1 1]);
    MD(n) = getframe;
end;
movie(MD);








% 
% 
% %plotting and normalizing vectors
% for j=1:101                      
%     W(j,:) = u(90*(j-1)+1,:);
%     tnew(j) = t(90*(j-1)+1);
% end;
% 
% 
% 
% 
% tnew = (1:101)*t/(101);
% % figure;
% % plot(x,W(101,:));
% % % hold on;
% % plot (x, u(discretization_num_t + 2,:), 'r');
% % 
% figure;
% surfl(x,tnew,W,'light');
% xlabel('variable x');
% ylabel('variable t');
% zlabel('function u');
% % title('solution eqation 3 part4');  
% 
% 
%  %display(showsolutionInPoint(0.5, 0.5, u, x_0, delta_x, t_0, delta_t));
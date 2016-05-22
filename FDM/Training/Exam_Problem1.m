x_0 = 0;  x_f = 1; %interval for x
t_0 = 0;  t_f = 0.1;  % interval for t 

discretization_num_x = 49; % number of discretization 
discretization_num_t = 999; % number of discretization
delta_x = (x_f - x_0)/(discretization_num_x+1); % delta x
delta_t = (t_f - t_0)/(discretization_num_t+1); % delta y

% definition of the x-values on axis
x(1) = x_0;
for q = 2:discretization_num_x + 2
    x(q) = x(q-1) + delta_x;   
end;

% definition of the t-values on axis
t(1) = t_0;
for q = 2:discretization_num_t + 2
    t(q) = t(q-1) + delta_t; 
end;

% Initial condition
for i = 1:discretization_num_x + 2;
    u(1,i) = sin(2*pi*x(i)); % g(x) function. Go through all values of u for given t = t_0
end;
figure; plot(x, u(1,:), 'r'); xlabel('variable t'); ylabel('variable u'); title('solution problem 1 T = 0.1 initial condition t = MIN'); 


for n = 1:discretization_num_t + 2;
   u(n,1) = -8*t(n);  % Go through all values of u for given x = x_0
   u(n,discretization_num_x+2) = 10*t(n);
end;

% till this moment we have defined u(1,:), u(:,1), u(:,end) (defining all
% the surface including boundary conditions implementation - discretized)
for n = 1:discretization_num_t+1
    for i = 2:discretization_num_x + 1
        u(n+1,i) = u(n,i)+(delta_t/(delta_x)^2)*(u(n,i+1)+u(n,i-1)-2*u(n,i))+(delta_t/delta_x)*x(i)*(u(n,i+1)-u(n,i-1))+3*delta_t*t(n)*u(n,i)+delta_t*(t(n)^2)*(x(i)-1);
    end;
end;
  
figure; plot(t,u(:,1),'r');  xlabel('variable t'); ylabel('variable u');
title('solution problem 1 T = 0.1 boundary condition x = MIN'); 
figure; plot(t,u(:,discretization_num_x + 2), 'b'); xlabel('variable t'); ylabel('variable u'); title('solution problem 1 T = 0.1 boundary condition x = MAX');

%plotting and normalizing vectors
for j=1:101                      
    W(j,:) = u(10*(j-1)+1,:);
    tnew(j) = t(10*(j-1)+1);
end;

%tnew = (1:101)*t/(101);
figure;
surfl(x,tnew,W,'light'); xlabel('variable x'); ylabel('variable t'); zlabel('function u'); title('solution problem 1');   


figure; plot(x,u(1,:)); xlabel('variable x'); ylabel('variable t'); zlabel('function u'); title('u(t,x) when t = 0'); 
figure; plot(x,u((discretization_num_t+1)/2,:)); xlabel('variable x'); ylabel('variable t'); zlabel('function u'); title('u(t,x) when t = T/2'); 
figure; plot(x,u(discretization_num_t+2,:)); xlabel('variable x'); ylabel('variable t'); zlabel('function u'); title('u(t,x) when t = T'); 

answer_6 = showsolutionInPoint(1/2, t_f/4, u, x_0, delta_x, t_0, delta_t);
display(answer_6);
x_0 = 0;  x_f = 1; %interval for x
t_0 = 0;  t_f = 0.01;  % interval for t 

discretization_num_x = 100; % number of discretization 
discretization_num_t = 9999; % number of discretization
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

% Final bs-condition definition
%u(discretization_num_t + 2,1) = initialcond_bs(x(1),Strike); % by initial condition and knowing first valuse of x we can find the first value of u
for i = 1:discretization_num_x + 2;
    u(1,i) = initialcond_training(x(i)); % g(x) function. Go through all values of u for given t = t_0
end;
plot(x, u(1,:), 'r');

display (u(:,discretization_num_x + 2));

for n = 1:discretization_num_t + 2;
   u(n,1) = 100*t(n);  % Go through all values of u for given x = x_0
end;

% till this moment we have defined u(1,:), u(:,1), u(:,end) (defining all
% the surface including boundary conditions implementation - discretized)
for n = 1:discretization_num_t+1
    for i = 2:discretization_num_x + 1
        u(n+1,i) = u(n,i)+(delta_t/(delta_x)^2)*(u(n,i+1)+u(n,i-1)-2*u(n,i))+(delta_t/(2*delta_x))*(u(n,i+1)-u(n,i-1))+delta_t*exp(2*t(n))*x(i)*(x(i)-1);
    end;
        u(n+1,discretization_num_x + 2) = u(n+1,discretization_num_x) - 4*t(n)*delta_x; 
end;
  
figure; plot(t,u(:,1),'r'); 
xlabel('variable t');
ylabel('variable u');
title('solution eqation 3 part 4 T = 0.5 boundary condition x = 1'); 
figure; plot(t,u(:,discretization_num_x + 2), 'b');
xlabel('variable t');
ylabel('variable u');
title('solution eqation 3 part 4 T = 0.5 boundary condition x = 22');
display (u(:,discretization_num_x + 2));

%plotting and normalizing vectors
for j=1:101                      
    W(j,:) = u(100*(j-1)+1,:);
    tnew(j) = t(100*(j-1)+1);
end;


%tnew = (1:101)*t/(101);
figure;
surfl(x,tnew,W,'light');
xlabel('variable x');
ylabel('variable t');
zlabel('function u');
title('solution eqation 3 part4');   

 %display(showsolutionInPoint(0.5, 0.5, u, x_0, delta_x, t_0, delta_t));
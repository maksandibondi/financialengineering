x_0 = 0;  x_f = 40; %interval for x
t_0 = 0;  t_f = 0.5;  % interval for t 
sigma = 0.5;
r = 0.1; 
K = 10;

discretization_num_x = 100; % number of discretization 
discretization_num_t = 1999; % number of discretization
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
for i = 1:discretization_num_x + 2;
    if (x(i)<=8 || x(i)>=12)
    u(discretization_num_t+2,i) = 0;
    else
    u(discretization_num_t+2,i) = 1;
    end;
end;
figure; plot(x, u(discretization_num_t+2,:), 'r'); xlabel('variable t'); ylabel('variable u'); title('solution Problem 3 part 2 T = 0.5 final condition t = MAX');


% dirichlet boundaries
for n = 1:discretization_num_t + 2;
   u(n,1) = 0;  % Go through all values of u for given x = x_0
   u(n,discretization_num_x + 2) = 0;
end;



% till this moment we have defined u(1,:), u(:,1), u(:,end) (defining all
% the surface including boundary conditions implementation - discretized)
for n = discretization_num_t+2:-1:2
    for i = 2:discretization_num_x + 1
        u(n-1,i) = u(n,i)-delta_t*r*u(n,i)+(delta_t/(2*delta_x))*r*x(i)*(u(n,i+1)-u(n,i-1))+(delta_t/(2*delta_x^2))*(sigma^2)*(x(i)^2)*(u(n,i+1)-2*u(n,i)+u(n,i-1));
    end; 
        u(n-1,1) = u(n-1,2);
        u(n-1,discretization_num_x+2) = u(n-1,discretization_num_x+1)+delta_x;
end;
  

figure; plot(t,u(:,1),'r');  xlabel('variable t'); ylabel('variable u'); title('solution Problem 3 part 2 T = 0.5 boundary condition x = MIN'); 
figure; plot(t,u(:,discretization_num_x + 2), 'b'); xlabel('variable t'); ylabel('variable u'); title('solution Problem 3 part 2 T = 0.5 boundary condition x = MAX');

%plotting and normalizing vectors
for j=1:101                      
    W(j,:) = u(20*(j-1)+1,:);
    tnew(j) = t(20*(j-1)+1);
end;


%tnew = (1:101)*t/(101);
figure;
surfl(x,tnew,W,'light'); xlabel('variable x'); ylabel('variable t'); zlabel('function u'); title('solution Problem 3 part 2');   

figure; plot(x,u(1,:)); xlabel('variable x'); ylabel('variable t'); zlabel('function u'); title('u(t,x) when t = 0');
figure; plot(x,u(discretization_num_t+2,:)); xlabel('variable x'); ylabel('variable t'); zlabel('function u'); title('u(t,x) when t = T'); 
display(showsolutionInPoint(10, t_0, u, x_0, delta_x, t_0, delta_t));
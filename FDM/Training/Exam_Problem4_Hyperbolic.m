x_0 = 4;  x_f = 8; %interval for x
t_0 = 0;  t_f = 5;  % interval for t 
A = 5; sigma = 0.5; Q = 4*pi;


discretization_num_x = 100; % number of discretization 
discretization_num_t = 4999; % number of discretization
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
    u(1,i) = A*exp(-((x(i)-x_0)^2)/(2*sigma^2))*cos(Q*(x(i)-x_0)); % g(x) function. Go through all values of u for given t = t_0
    u(2,i) = u(1,i)+delta_t*velocity(x(i),x_0,A,sigma,Q); % velocity
end;
figure; plot(x, u(1,:), 'r'); xlabel('variable t'); ylabel('variable u'); title('solution problem 4 initial condition t = MIN');


for n = 1:discretization_num_t + 2;
   u(n,1) = 0;  % Go through all values of u for given x = x_0 
   u(n,discretization_num_x+2) = 0;
end;

% Main
for n = 2:discretization_num_t+1
    for i = 2:discretization_num_x + 1
        u(n+1,i) = 2*u(n,i)-u(n-1,i)+((delta_t/delta_x)^2)*(u(n,i+1)+u(n,i-1)-2*u(n,i));
    end;
    u(n,1) = u(n,discretization_num_x+1); % periodic boundarystart condition
    u(n,discretization_num_x+2) = u(n,2); % periodic boundaryfinal condition
end;
 

%plotting and normalizing vectors
for j=1:101                      
    W(j,:) = u(50*(j-1)+1,:);
    tnew(j) = t(50*(j-1)+1);
end;


figure;
for n = 1:50:discretization_num_t+2
    plot(x,u(n,:));
    axis([0 10 -5 5]);
    MD(n) = getframe;
end;
movie(MD);

figure; plot(x,u(1,:)); xlabel('variable x'); ylabel('variable t'); zlabel('function u'); title('u(t,x) when t = 0'); 
figure; plot(x,u((discretization_num_t+1)/2,:)); xlabel('variable x'); ylabel('variable t'); zlabel('function u'); title('u(t,x) when t = T/2'); 
figure; plot(x,u(discretization_num_t+2,:)); xlabel('variable x'); ylabel('variable t'); zlabel('function u'); title('u(t,x) when t = T'); 
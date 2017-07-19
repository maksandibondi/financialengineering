clear; clc;

%% Initial data
t_f = 0.5; t_0 = 0;
discretization_num_t = 20;
T = t_f-t_0;

%% Time and BM arrays creation
delta_t = T/discretization_num_t;
t(1) = t_0;
for i = 2:discretization_num_t
    t(i) = t(i-1)+delta_t;
end;

%% Market data 

Rm = [0.01 0.108 0.039 0.038 0.05 0.1291 0.1915 0.1914 0.1238 0.1587 0.3521 0.4047 0.4650 0.3700 0.4116 0.5142 0.4010 0.5389 0.8104 0.6796];

%% Initial values definition
d = [1 1];
beta = [1 1]; % params to calibrate
A = beta(1); B = beta(2);
precision = 10^(-9); lambda = 0.01;
identity = eye(2,2);

%% Levenberg Marquart algorithm

while norm(d) > precision

    % Giving concrete values to symbolic parameters
    for p = 2:discretization_num_t
           
           Rmarket = Rm(p); r = Rm(p-1);
        
           residual(p) = Rmarket - (A * r + B);
                       
           J(p,1) = -Rm(p);  %Derivative in respect of beta(1)
           J(p,2) = -1;     %Derivative in respect to beta(2) 
           
    end;
    
    % Calculating d using the algorithm

    M = J' * J + lambda * identity;
    d = -inv(M)*J'*residual';
    
    beta = beta + d';
    A = beta(1); 
    B = beta(2); 
 
end;
display(beta);


%% Determining theoretical curve and plotting
Rth(1) = Rm(1);
for p = 2:discretization_num_t
    Rth(p) = Rm(p-1)*beta(1)+beta(2);
end;

figure;
plot(Rm*beta(1)+beta(2), Rm); % line
hold on;
plot(Rm(2:end),Rm(1:end-1),'*'); % line




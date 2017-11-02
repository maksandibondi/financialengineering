clear; clc;

%% Initial data
t_f = 30; t_0 = 0;
r0 = 0.1; 
discretization_num_t = 500;
T = t_f-t_0;

%% Time and BM arrays creation
delta_t = T/discretization_num_t;
t(1) = t_0;
for i = 2:discretization_num_t
    t(i) = t(i-1)+delta_t;
end;

%% Market data simulation (using BM and arbitary parameters to be calibrated later)
Rm(1) = r0;
eta_m = 0.25; % long term mean
gamma_m = 0.25; % speed of reversion
sigma_m = 0.03; % volatility
A = exp(-gamma_m*delta_t);
B = (eta_m/gamma_m)*(1-A);
Dev = sigma_m*sqrt((1-exp(-2*gamma_m*delta_t))/(2*gamma_m));

for k = 2:discretization_num_t 
    Rm(k) = Rm(k-1)*A + B + Dev*randn();
end;
plot(t,Rm,'*');
        
%% Initial values definition
d = [1 1];
beta = [A B]; % params to calibrate
A_cal = beta(1); B_cal = beta(2);
precision = 10^(-9); lambda = 0.01;
identity = eye(2,2);

%% Levenberg Marquart algorithm

while norm(d) > precision

    % Giving concrete values to symbolic parameters
    for p = 2:discretization_num_t
           
           Rmarket = Rm(p); r = Rm(p-1);
        
           residual(p) = Rmarket - (A_cal * r + B_cal);
                       
           J(p,1) = -Rm(p);  %Derivative in respect of beta(1) [Ac] ??????????????????????
           J(p,2) = -1;     %Derivative in respect to beta(2) [Bc]
           
    end;
    
    % Calculating d using the algorithm

    M = J' * J + lambda * identity;
    d = -inv(M)*J'*residual';
    
    beta = beta + d';
    A_cal = beta(1); 
    B_cal = beta(2); 
 
end;
display(beta);

%% Comparison of market and calibrated parameters
gamma_cal = -(log(beta(1))/delta_t)
gamma_m
eta_cal = (beta(2) * gamma_cal) / (1 - beta(1))
eta_m

variance = sum(residual.^2)/discretization_num_t;
sigma_cal = sqrt(variance)*sqrt(2*gamma_cal/(1 - exp(-2 * gamma_cal * delta_t)))
sigma_m

%% Determining theoretical curve and plotting
Rth(1) = Rm(1);
for p = 2:discretization_num_t
    Rth(p) = Rm(p-1)*beta(1)+beta(2);
end;

figure;
plot(Rm*beta(1)+beta(2), Rm); % line
hold on;
plot(Rm(2:end),Rm(1:end-1),'*'); % line




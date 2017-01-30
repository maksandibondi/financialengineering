clear; clc;

%% Initial data
t2 = 5; t1 = 0;
r0 = 0.3; 
discretization_num_t = 100;
T = t2-t1;

%% Time and BM arrays creation
W = BMsimulator(T,discretization_num_t,'TCL');
Ws2_(1) = 0;
Ws1_(1) = 0;
delta_t = T/discretization_num_t;
t(1) = t1;
for i = 2:discretization_num_t
    t(i) = t(i-1)+delta_t;
    Ws2_(i) = W(i);
    Ws1_(i) = W(i-1);
end;

%% Market data simulation (using BM and arbitary parameters to be calibrated later)
Rm(1) = r0;
eta_m = 0.4; % long term mean
gamma_m = 0.6; % speed of reversion
sigma_m = 0.25; % volatility
for k = 2:discretization_num_t 
    Rm(k) = Rm(k-1)*exp(-gamma_m*(t(k)-t(k-1)))+(eta_m/gamma_m)*(1-exp(-gamma_m*(t(k)-t(k-1))))+sigma_m*exp(-gamma_m*(t(k)-t(k-1)))*(Ws2_(k)-Ws1_(k));
end;
plot(t,Rm,'*');
% display(Rm);

%% Symbolic calculations
        syms eta gamma sigma t2 t1 Ws2 Ws1 res Rth Rmarket;
        
        R_th_term_1 = r0*exp(-gamma*(t2-t1));
        R_th_term_2 = (eta/gamma)*(1-exp(-gamma*(t2-t1)));
        % R_th_term_3 = sigma*exp(-gamma*(t2-t1))*(Ws2-Ws1);
        
        Rth = R_th_term_1 + R_th_term_2;
        res = Rmarket - Rth;
        % Symbolic derivatives of residual function to minimize
        drEta = diff(res,eta);
        drGamma = diff(res,gamma);
        drSigma = diff(res,sigma);
        
%% Initial values definition
d = [1 1];
beta = [2 2]; % params to calibrate

precision = 0.001; lambda = 1;
identity = eye(2,2);

%% Levenberg Marquart algorithm
    
while norm(d) > precision

    % Giving concrete values to symbolic parameters
    for p = 1:discretization_num_t
           Rmarket = Rm(p); eta = beta(1); gamma = beta(2); 
           t2 = t(p); t1 = 0;
           Ws2 = Ws2_(p);
           Ws1 = Ws1_(p);
           
    % jacobian matrix evaluation    
           J(p,1) = eval(drEta); 
           J(p,2) = eval(drGamma); 
           residual(p) = eval(res);  
    end;
    
    % Calculating d using the algorithm
    d = -inv((J'*J)+(lambda.*identity))*J.'*residual';
    display(beta);
    beta = beta+d';
end;
display(beta);

%% Determining theoretical curve and plotting
variance = 0;
for j = 1:discretization_num_t
    Rth(j) = r0*exp(-beta(2)*(t(j)-t1))+(beta(1)/beta(2))*(1-exp(-beta(2)*(t(j)-t1)));
    variance = variance + (Rm(j)-Rth(j))^2;
end;
figure;
plot(t,Rth,'*'); hold on; scatter(t,Rm,'r');

%% Comparison of market and calibrated parameters
beta(1)
eta_m
beta(2)
gamma_m
sigma_calibrated = eval(sqrt(variance))
sigma_m


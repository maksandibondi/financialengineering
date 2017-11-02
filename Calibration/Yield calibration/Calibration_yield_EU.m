clear; clc;

%% Market data

R = csvread('rate_R2.csv',0,0);
T = csvread('Maturity_T2.csv',0,0);

% plot(T,R,'*')

%% Symbolic calculations
        syms B1 B2 B3 B4 dT r_ R_h R_m;
        
        R_hat_term_1 = B1;
        R_hat_term_2 = ((B2 * B4)/dT)*(1-exp(-dT/B4));
        R_hat_term_3 = - B3 * exp(-dT/B4);
        R_hat_term_4 = ((B3*B4)/dT)*(1-exp(-dT/B4));
        
        R_h = R_hat_term_1 + R_hat_term_2 + R_hat_term_3 + R_hat_term_4;
        r_ = R_m - R_h;
        % Symbolic derivatives of residual function to minimize
        drB1 = diff(r_,B1);
        drB2 = diff(r_,B2);
        drB3 = diff(r_,B3);
        drB4 = diff(r_,B4);
        
%% Initial values definition

d = [1 1 1 1]; beta = [-1 1 1 1];
sz = size(R,1);

precision = 0.001; t = 2; lambda = 1;
identity = eye(4,4);

%% Newton algorithm for curve calibration

while norm(d) > precision

    % Jacobian matrix evaluation
    for p = 1:sz
           R_m = R(p); B1 = beta(1); B2 = beta(2); B3 = beta(3); B4 = beta(4);
           dT = T(p); 
        
           J(p,1) = eval(drB1); J(p,2) = eval(drB2); J(p,3) = eval(drB3); J(p,4) = eval(drB4);
           r(p) = eval(r_);  
    end

    d = -inv((J'*J)+(lambda.*identity))*J.'*r';
    display(beta);
    beta = beta+d';
end
display(beta);

%% Plotting
for j = 1:sz
    R_hat(j) = beta(1)+beta(2)*((beta(4)/(T(j)-t))*(1-exp(-(T(j)-t)/beta(4))))+beta(3)*((beta(4)/(T(j)-t))*(1-exp(-(T(j)-t)/beta(4)))-exp(-(T(j)-t)/beta(4)));
end;
plot(T,R_hat); hold on; scatter(T,R,'r');
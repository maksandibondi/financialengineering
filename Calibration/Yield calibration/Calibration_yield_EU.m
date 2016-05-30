clear; clc;

%% market data

R = csvread('rate_R.dat',0,0);
T = csvread('Maturity_T.dat',0,0);

% plot(T,R,'*')

%% initial values

d = [1 1 1 1];
beta = [-1 1 1 1];
sz = size(R,1);
precision = 0.0000001;

t = 0;
lambda = 1;
identity = eye(4,4);

%% Newton algorithm

while norm(d) > precision

    for p = 1:sz
        J(p,1) = -1;
        J(p,2) = -((beta(4)/T(p))*(1-exp(-T(p)/beta(4))));
        J(p,3) = -(beta(4)/T(p))*(1-exp(-T(p)/beta(4))) + exp(-T(p)/beta(4));
        J(p,4) = -((beta(2)+beta(3))/T(p))*(1-exp(-T(p)/beta(4))) + (beta(2)+beta(3)*(T(p)/beta(4))+beta(3))*((exp(-T(p)/beta(4)))/beta(4));
        
        R_hat = beta(1)+beta(2)*((beta(4)/(T(p)-t))*(1-exp(-(T(p)-t)/beta(4))))+beta(3)*((beta(4)/(T(p)-t))*(1-exp(-(T(p)-t)/beta(4)))-exp(-(T(p)-t)/beta(4)));
        
        r(p) = R(p)-(R_hat);
    end
    
    d = -inv((J'*J)+(lambda.*identity))*J.'*r';
    beta = beta+d';
end

%% values of Beta
display(beta);

%% Plotting
for j = 1:sz
    R_hat(j) = beta(1)+beta(2)*((beta(4)/(T(j)-t))*(1-exp(-(T(j)-t)/beta(4))))+beta(3)*((beta(4)/(T(j)-t))*(1-exp(-(T(j)-t)/beta(4)))-exp(-(T(j)-t)/beta(4)));
end;
plot(T,R_hat); hold on; scatter(T,R,'r');
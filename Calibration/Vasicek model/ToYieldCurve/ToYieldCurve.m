clear; clc;

%% Market data
r0 = 0.023;
R = csvread('rate_R2.csv',0,0);
T = csvread('Maturity_T2.csv',0,0);

%% Symbolic calculations

syms eta variance gamma A B r tao T_ t_ Ytheory Ymarket res 

tao = T_ - t_;
B = (1-exp(-gamma*tao))/gamma;
A = (B-tao)*(eta*gamma-variance/2)/gamma - variance*(B^2)/(4*gamma);
Ytheory = -(A-r*B)/T_;

res = Ymarket - Ytheory;

dEta = diff(res, eta);
dVar = diff(res, variance);
dGamma = diff(res, gamma);

%% Initial values

d = [1 1 1];
beta = [1 1 1];
precision = 0.001; lambda = 0.01;
identity = eye(3,3);
discretization_num_t = size(T,1);
t_ = 0; 

%% Levenberg-Marquart

while norm(d) > precision
    
    for i = 1:discretization_num_t
        
        Ymarket = R(i); T_ = T(i); r = r0;
        eta = beta(1); variance = beta(2); gamma = beta(3);
        
        J(i,1) = eval(dEta); 
        J(i,2) = eval(dVar); 
        J(i,3) = eval(dGamma);
        residual(i) = eval(res);
        
    end;
    
    d = -inv((J'*J)+(lambda*identity))*J'*residual';
    norm(d)
    beta = beta+d';
    
end;

Eta = beta(1)
Var = beta(2)
Gamma = beta(3)


for i = 1:discretization_num_t
    r = r0; T_ = T(i); 
    Rth(i) = eval(-(A-r*B)/T_);
end;

scatter([t_; T],[r0; R]);
hold on;
plot([t_; T],[r0 Rth]);




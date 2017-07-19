clear; clc;

%% Market data
r0 = 0.2;
R = csvread('rate_R2.csv',0,0);
T = csvread('Maturity_T2.csv',0,0);

%% Symbolic calculations

syms a variance A B r tao T_ t_ Ytheory Ymarket res 

tao = T_ - t_;
B = tao;
A = -(1/2)*a*tao^2+(1/6)*variance*tao^3;
Ytheory = -(A-r*B)/T_;

res = Ymarket - Ytheory;

da = diff(res, a);
dVar = diff(res, variance);

%% Initial values

d = [1 1];
beta = [1 1];
precision = 10^(-9); lambda = 0.01;
identity = eye(2,2);
discretization_num_t = size(T,1);
t_ = 0; 

%% Levenberg-Marquart

while norm(d) > precision
    
    for i = 1:discretization_num_t
        
        Ymarket = R(i); T_ = T(i); r = r0;
        a = beta(1); variance = beta(2);
        
        J(i,1) = eval(da); 
        J(i,2) = eval(dVar); 
        residual(i) = eval(res);
        
    end;
    
    d = -inv((J'*J)+(lambda*identity))*J'*residual';
    norm(d);
    beta = beta+d';
    
end;

a = beta(1)
Var = beta(2)


for i = 1:discretization_num_t
    r = r0; T_ = T(i); 
    Rth(i) = eval(-(A-r*B)/T_);
end;

scatter([t_; T],[r0; R]);
hold on;
plot([t_; T],[r0 Rth]);
xlabel('T'); ylabel('Y(t,T)');



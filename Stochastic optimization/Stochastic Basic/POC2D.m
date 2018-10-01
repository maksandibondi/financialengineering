clear; clc;
tic;
% Stochastic algoruthm of optimimztion basic

%% Input data
T0 = 500;
t = 0.2;
Nmc = 2^12;
M = 64;

x1 = (-1+2*rand(30,1))*5.12; % Initial population of 30 candidates x1,x2 Rastragin
x2 = (-1+2*rand(30,1))*5.12;

sz = size(x1,1);


%func = @(x1, x2) 10*2 + (x1.^2-10*cos(2*pi.*x1)) + (x2.^2-10*cos(2*pi.*x2));
func = @(x1,x2) x1.^2+x2.^2;
cost = @(x1, x2) abs(func(x1, x2));
gibbs = @(x1, x2, T) exp(-cost(x1, x2)/T);

T = T0;
for k = 1:Nmc
    if (min(cost(x1, x2))>0.001)
    T = T*(1-k*(mod(k,M)==0)/Nmc)^4;
    
    % Selection
    for i = 1:sz  
        if (rand() > gibbs(x1(i),x2(i),T))
            for l = 1:sz
                % Proba of chosing another l-th candidate. sum of these proba == 1
                p(l) = gibbs(x1(l),x2(l),T)/sum(gibbs(x1,x2,T)); 
            end;
            idx = ChoosenIdx(rand(),p);
            x1(i) = x1(idx);
            x2(i) =  x2(idx);
        end;
    end;
    
    
    % Mutation
    for i = 1:sz
        x1(i) = x1(i) + t*(2*rand()-1);
        x2(i) = x2(i) + t*(2*rand()-1);
    end;
    
    else 
        break;
    end;
    m = min(cost(x1,x2))
    %scatter(m,func(m)); xlabel('x'); ylabel('y'); hold on; 
end;
toc;
%figure; hold on;
res = min(cost(x1,x2))
%y = -20:0.1:20;
x1_ = -5.2:0.1:5.2;
x2_ = -5.2:0.1:5.2;

[X1,X2]=meshgrid(x1_,x2_);
Z = func(X1,X2);
surf(Z); hold on;
k

scatter3(x1,x2,func(x1,x2));

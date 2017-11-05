clear; clc;
% Stochastic algoruthm of optimimztion basic

%% Input data
T0 = 1000;
t = 0.2;
Nmc = 2^16;
M = 64;
x = (-1+2*rand(30,1))*5; % Initial population of 10 candidates
%y = (-1+2*rand(10,1))*5; % Initial population of 10 candidates
sz = size(x);

%beale = @(x,y) (1.5-x+x.*y).^2+(2.25-x+x.*(y.^2)).^2+(2.625-x+x.*(y.^3)).^2;
%gibbs = @(x,y,T) exp(-beale(x,y)/T);

func = @(x) x.^2-6.*x-3;
cost = @(x) abs(func(x)-(-12)); % min if -6 at x=3
gibbs = @(x,T) exp(-cost(x)/T);

T = T0;
for k = 1:Nmc
    if (sum(cost(x)>0.1)>0)
    T = T*(1-k*(mod(k,M)==0)/Nmc)^4;
    
    % Selection
    for i = 1:sz  
        if (rand() > gibbs(x(i),T))
            for l = 1:sz
                % Proba of chosing another candidate. sum of these proba == 1
                p(l) = gibbs(x(l),T)/sum(gibbs(x,T)); 
            end;
            x(i) = x(ChoosenIdx(rand(),p));
        end;
    end;
    
    
    % Mutation
    for i = 1:sz
        x(i) = x(i) + t*(2*rand()-1);
    end;
    
    else 
        break;
    end;
    
end;

figure; hold on;
y = -10:0.1:10;
plot(y,func(y));

k
x
cost(x)

scatter(sort(x),func(sort(x)));
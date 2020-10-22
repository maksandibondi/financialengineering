clear; clc;
tic;
% Stochastic algoruthm of optimimztion basic

%% Input data
T0 = 1000;
t = 0.05;
Nmc = 2^14;
M = 64;
%x = (-1+2*rand(30,1))*20; % Initial population of 30 candidates
x = (-1+2*rand(30,1))*5; % Initial population of 30 candidates Rastragin
%y = (-1+2*rand(10,1))*5; % Initial population of 10 candidates
sz = size(x);

%beale = @(x,y) (1.5-x+x.*y).^2+(2.25-x+x.*(y.^2)).^2+(2.625-x+x.*(y.^3)).^2;
%gibbs = @(x,y,T) exp(-beale(x,y)/T);

%func = @(x) 10 + (x.^2-10*cos(2*pi.*x));
%cost = @(x) abs(func(x));
func = @(x) x.^2-6.*x-3;
cost = @(x) abs(func(x)-(-12)); % min if -6 at x=3

gibbs = @(x,T) exp(-cost(x)/T);

T = T0;
for k = 1:Nmc
    if (min(cost(x))>0.001)
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
   % m = min(x);
    %scatter(m,func(m)); xlabel('x'); ylabel('y'); hold on; 
end;
toc;
%figure; hold on;
res = min(cost(x))
%y = -20:0.1:20;
y = -5:0.1:5;
plot(y,func(y)); hold on;

k
x;
func(x); 

scatter(sort(x),func(sort(x)));

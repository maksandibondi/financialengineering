clear; clc;

% BS parameters
S0 = 15;
mu = 0.08;
r = 0;
sigma = 0.2;
T = 1; 
delta_t = 0.02;
K = 18;

% Algo parameters
Nmc = 500; % total num of paths

% creating t-array
discretization_num_t = T/delta_t; 
 t(1) = 0; 
 for i = 2:discretization_num_t 
     t(i) = t(i-1)+delta_t; 
 end; 

 %initializing first betas (known)


PriceEU0 = 0;
PriceEU = 0;
payoff_call0(1) = 0;
payoff_call(1) = 0;

for i = 2:Nmc+1
  
    X0 = BSStockSimulator(S0,mu,sigma,T,delta_t);
    X = BSStockSimulator(S0,mu,sigma,T,delta_t);
    payoff_call0(i) = payoff_EU(X0(end), K, 'Call')*exp(-r*T);
    payoff_call(i) = (payoff_EU(X(end), K, 'Put'))*exp(-r*T) + S0 - exp(-r*T)*K; % Problem - it has a negative value sometimes
    
    PriceEU0 = PriceEU +  payoff_call0(i)/Nmc; % with var reduction       
    PriceEU = PriceEU +  payoff_call(i)/Nmc; % with var reduction
    
  
end;

payoff_call
payoff_call0

Pricemean0 = mean(payoff_call0(2:end))
Pricemean = mean(payoff_call(2:end))


Pricevar0 = var(payoff_call0(2:end))
Pricevar = var(payoff_call(2:end))

PriceEU0
PriceEU
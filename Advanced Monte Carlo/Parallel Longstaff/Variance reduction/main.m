clear; clc;

% BS parameters
S0 = 14;
r = 0.08;
sigma = 0.1;
T = 1; 
delta_t = 0.005;
K = 18;

% Algo parameters
Nmc = 40000; % total num of paths

% creating t-array
discretization_num_t = T/delta_t; 
 t(1) = 0; 
 for i = 2:discretization_num_t 
     t(i) = t(i-1)+delta_t; 
 end; 

 %initializing first betas (known)


PriceEUcall0 = 0;
PriceEUcall = 0;
PriceEUput0 = 0;
PriceEUput = 0;

payoff_call0(1) = 0;
payoff_call(1) = 0;
payoff_put0(1) = 0;
payoff_put(1) = 0;


for i = 2:Nmc+1
  
   % X0 = BSStockSimulator(S0,r,sigma,T,delta_t);
    X0 = S0*exp((r-(sigma^2)/2)*T+sigma*sqrt(T)*randn());
    
    %% Pour call
    payoff_call0(i) = payoff_EU(X0, K, 'Call')*exp(-r*T);
    payoff_call(i) = payoff_EU(X0, K, 'Put')*exp(-r*T) + S0 - exp(-r*T)*K; % Problem - it has a negative value sometimes 
    
    %% Pour put
    payoff_put0(i) = payoff_EU(X0, K, 'Put')*exp(-r*T);
    payoff_put(i) = payoff_EU(X0, K, 'Call')*exp(-r*T) - S0 + exp(-r*T)*K; % Problem - it has a negative value sometimes 
    
    
    PriceEUcall0 = PriceEUcall0 +  payoff_call0(i)/Nmc; % with var reduction       
    PriceEUcall = PriceEUcall +  payoff_call(i)/Nmc; % with var reduction
    PriceEUput0 = PriceEUput0 +  payoff_put0(i)/Nmc; % with var reduction       
    PriceEUput = PriceEUput +  payoff_put(i)/Nmc; % with var reduction
    
    
  
end;

PriceEUcall0
PriceEUcall
Pricevarcall0 = var(payoff_call0(2:end))
Pricevarcall = var(payoff_call(2:end))

PriceEUput0
PriceEUput
Pricevarput0 = var(payoff_put0(2:end))
Pricevarput = var(payoff_put(2:end))

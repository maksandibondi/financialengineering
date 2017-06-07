clear; clc;

% BS parameters
S0 = 15;
mu = 0.08;
r = 0;
sigma = 0.6;
T = 1; 
delta_t = 0.02;
K = 18;

% Algo parameters
n = 50; % number of iterations
Nmc = 500; % total num of paths
NpathIter = Nmc/n; % number of paths in one iteration

% creating t-array
discretization_num_t = T/delta_t; 
 t(1) = 0; 
 for i = 2:discretization_num_t 
     t(i) = t(i-1)+delta_t; 
 end; 

 %initializing first betas (known)

% Initial beta values
for k = 1:discretization_num_t
  beta(1,k) = 1; %  initial a(0,k) in the article
end;
% Weights of iterations
for i = 1:n
   IterWeight(i) = 1;
end;
PriceFinal = 0;
qfinal = 0;
P = zeros(NpathIter, discretization_num_t);
u = zeros(n,NpathIter,discretization_num_t);
v = zeros(n,NpathIter,discretization_num_t);
U = zeros(n, discretization_num_t);
V = zeros(n, discretization_num_t);


for i = 2:n
  
  for j = 1:NpathIter
    
    X = BSStockSimulator(S0,mu,sigma,T,delta_t);
    Pnext = max(payoff_AM(X(end), K, 'Put'),0);
    for k = discretization_num_t:-1:1
      F = max(payoff_AM(X(k), K, 'Put'),0);
      C = beta(i-1,k)*F;
      u(i,j,k) = PathWeight(X(k), K, 'Put')*(F^2);
      v(i,j,k) = PathWeight(X(k), K, 'Put')*F*exp(-r*delta_t)*Pnext;
      
      if (F>C) % comparison of exercise F against continuation C
        P(j,k) = F; % current payoff (exercise)
      else
        P(j,k) = Pnext*exp(-r*delta_t); % discounted payoff from the next period (continuation)
      end;
      
      Pnext = P(j,k);
      U(i,k) = U(i-1,k) + u(i,j,k);
      V(i,k) = V(i-1,k) + v(i,j,k);
   
   end;
    
  end;
  
  PriceContribution(i) = sum(P(:,1)); % Sum of all paths payoffs at time 1
  beta(i,:) = V(i,:)./U(i,:);
  
  PriceFinal = PriceFinal + IterWeight(i)*PriceContribution(i);
  qfinal = qfinal + IterWeight(i)*(NpathIter);

end;

AMprice = PriceFinal/qfinal
clear; clc;

%% Data reading
data=dlmread('sp-index.csv',';',1,0);
% B=csvread('sp-index.dat',1,0);  %the matrix with all the data
A = data(59:115,:); %the matrix we study  (1:58,:)  (59:115,:)  (116:149,:) (150:171,:) (172:214,:)  (215:237,:)  (238:265,:)  (266:280,:)
display(A);

%% Parameters' arrays creation
T = A(1,1); % Maturity value
q = 0.0217;           
S0 = 1260.36*exp(-q*T); % Initial price
r = A(1,7)/100; % Rate
K = A(:,2); % Strike
M = (A(:,3)+A(:,4))/2; % Market price data
e = 0.0001; % eps
sz = size(K);
maturity = T*ones(1,length(K));

%% Main program
for j=1:sz
    
%% Arbitrage condition
    if max(S0-K*exp(-r*T),0) < M < S0 
    
        sigma = sqrt(2*abs((log(S0/K(j))+r*T)/T));
        
%% Newton algorithm
        while abs(BSTheory(sigma,S0,r,T,K(j))-M(j)) > e
        
            sigma = sigma-(BSTheory(sigma,S0,r,T,K(j))-M(j))/Vega(sigma,S0,r,T,K(j));
        
        end
    
        sigma_implied(j) = sigma;
    
    else
        
    sigma_implied(j) = 0;   
    
    end;
    
    disp(sigma_implied(j));
   
end

%% Filtering the data before plotting (bad data was marked as 0)
k = 1;
for j = 1:sz
    if sigma_implied(j) ~= 0
        sigma_new(k) = sigma_implied(j);
        K_new(k) = K(j);
        maturity_new(k) = T;
        k = k + 1;
    end;
end;

%% Plotting

% plot3(maturity, K, sigma_implied)
plot3(maturity_new, K_new, sigma_new);

   

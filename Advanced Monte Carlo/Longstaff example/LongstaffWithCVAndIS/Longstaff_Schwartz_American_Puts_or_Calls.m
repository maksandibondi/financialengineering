% Longstaff-Schwartz for American puts and calls under Black SCholes

% Fabrice Douglas Rouah, FRouah.com and Volopta.com

clc; clear;

%% Inline function for the Black Scholes call and put
BSCall = @(s,K,r,q,v,T) s*exp(-q*T)*normcdf((log(s/K) + (r-q+v^2/2)*T)/v/sqrt(T)) - K*exp(-r*T)*normcdf((log(s/K) + (r-q+v^2/2)*T)/v/sqrt(T) - v*sqrt(T));
BSPut  = @(s,K,r,q,v,T) K*exp(-r*T)*normcdf(-(log(s/K) + (r-q+v^2/2)*T)/v/sqrt(T) + v*sqrt(T)) - s*exp(-q*T)*normcdf(-(log(s/K) + (r-q+v^2/2)*T)/v/sqrt(T));

%% Stock price features.
Spot = 36;               % Spot price.
r = 0.05;                 % Risk free rate.
q = 0.02;

% Option features.
K = 40;                  % Strike price.
v = 0.30;                % Volatility.
T = 1;                % Maturity in years.

% Simulation features.
NT = 100;                % Number of time increments.
dt = T/NT;               % Time increment.
%NS = 1e4;                % Number of simulated paths.
NS = 15000;

% Trinomial tree features.
n = 250;                 % Time steps for trinomial tree.
EuroAmer = 'A';          % Flavor indicator for trinomial tree.
PutCall = 'p';           % Select a put in the trinomial tree.

% IS and CV features
IS = 'false';
theta = 0;
CV = 'call';



%% Simulate the stock price process under Black Scholes
FirstCol = Spot.*ones(NS,1);

randomT = randn(NS,NT-1);
E1 = exp((r-q-0.5*v^2)*dt + sqrt(dt)*v.*randomT);
E2 = exp((r-q+v*theta-0.5*v^2)*dt + sqrt(dt)*v.*randomT);

All = [FirstCol E1];
S1 = cumprod(All,2);

All2 = [FirstCol E2];
S2 = cumprod(All2,2);
clear FirstCol E1 All1 E2 All2

%% Calculate the BlackScholes European price for comparison.
if strcmp(PutCall,'P')
	EuroPrice = BSPut(Spot,K,r,q,v,T);
else
	EuroPrice = BSCall(Spot,K,r,q,v,T);
end
clear BSPut BSCall

%% Longstaff-Schwartz price
% Design matrix for the LSM algorithm
XmatrixHandle = {@(y)ones(length(y),1), @(y)(1-y),@(y)1./2.*(2-4.*y-y.^2)};
% Run the LSM algorithm
[EuroPriceLSM, AmerPriceLSM, Variance] = BlackScholesLSM(S1,K,r,q,v,T,NS,NT,dt,PutCall,XmatrixHandle,randomT, 'null', 'false', theta);
[EuroPriceLSM2, AmerPriceLSM2, Variance2] = BlackScholesLSM(S2,K,r,q,v,T,NS,NT,dt,PutCall,XmatrixHandle,randomT, CV, IS, theta);
% Early exercise premium
PremiumLSM = AmerPriceLSM - EuroPriceLSM;
% Control variate price
ControlVariate = EuroPrice + PremiumLSM;

%% Trinomial Tree price
AmerPriceTRI = TrinomialTree(Spot,K,r,q,v,T,n,'A',PutCall);
EuroPriceTRI = TrinomialTree(Spot,K,r,q,v,T,n,'E',PutCall);
PremiumTRI   = AmerPriceTRI - EuroPriceTRI;

%% Display the results
fprintf('Method              EuroPrice   AmerPrice   Premium \n')
fprintf('----------------------------------------------------\n')
fprintf('Longstaff-Schwartz  %8.4f   %8.4f   %8.4f \n',EuroPriceLSM,AmerPriceLSM,PremiumLSM);
fprintf('Trinomial Tree      %8.4f   %8.4f   %8.4f \n',EuroPriceTRI,AmerPriceTRI,PremiumTRI);
fprintf('LSM Control Variate    n/a     %8.4f      \n',ControlVariate);
fprintf('----------------------------------------------------\n')



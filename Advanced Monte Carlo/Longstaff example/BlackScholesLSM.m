function [EuroPrice AmerPrice] = BlackScholesLSM(S,K,r,q,T,NS,NT,dt,PutCall,XmatrixHandle)

% Longstaff-Schwartz method for American options
% Fabrice Douglas Rouah, FRouah.com and Volopta.com

% INPUTS
%  S = Matrix of simulated stock price, size NSxNT (rows are paths)
%  K = Strike
%  r = Risk free rate
%  q = Dividend yield
%  T = Maturity
%  NS = Number of stock price paths
%  NT = Number of time steps per path
%  dt = Time increment (T/NT)
%  PutCall = 'P' or 'C'

% Number of columns in the X matrix
NX = length(XmatrixHandle);

% Initialize the Cash Flows
CF = zeros(NS,NT);

% Set the last cash flows to the intrinsic value
if strcmp(PutCall,'P')
	CF(:,NT) = max(K - S(:,NT), 0);
elseif strcmp(PutCall,'C')
	CF(:,NT) = max(S(:,NT) - K, 0);
end

% European option value
EuroPrice = exp(-r*T)*mean(CF(:,NT));

% Work backwards through the stock prices until time j=2.
% We could work through to time j=1 but the regression will not be
% of full rank at time 1, so this is cleaner.
for j=NT-1:-1:2
	if strcmp(PutCall,'P')
		I = find(S(:,j) < K);           % Indices for Puts
	else
		I = find(S(:,j) > K);           % Indices for Calls
	end
	X = S(I,j);                         % X-vector = time j stock prices.
	Y = CF(I,j+1)*exp(-r*dt);           % Y-vector = time j+1 discounted CF.
    Z = zeros(length(X),NX);            % Design matrix for regression to predict cash flows
    for k=1:NX
        Z(:, k) = feval(XmatrixHandle{k}, X);
    end
    display(Z);
	beta = Z\Y;                       	% Regression parameters.
	P = Z*beta;                         % Regression predicted CF.
	if strcmp(PutCall,'P')
		J = max(K - X, 0) > P;          % Early exercise for puts.
	else
		J = max(X - K, 0) > P;          % Early exercise for calls.
    end
 	E = I(J);                           % Stock price indices where immediate exercise is optimal.
	C = setdiff((1:NS),E)';             % Stock price indices where continuation is optimal.
	if strcmp(PutCall,'P')
		CF(E,j) = max(K - X(J), 0);     % Replace with early exercise for puts
	else
		CF(E,j) = max(X(J) - K, 0);
	end
	CF(C,j) = exp(-r*dt)*CF(C,j+1);     % Continued CF are discounted back one period.
end

% American option value
AmerPrice = exp(-r*dt)*mean(CF(:,2));


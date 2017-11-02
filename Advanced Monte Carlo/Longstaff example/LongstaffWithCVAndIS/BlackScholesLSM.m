function [EuroPrice, AmerPrice, Variance] = BlackScholesLSM(S,K,r,q,v,T,NS,NT,dt,PutCall,XmatrixHandle, randomT, CV, IS, theta)


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

if strcmp(CV,'call')
ContVar(:,NT) = max(S(:,NT) - K, 0);
else
ContVar(:,NT) = zeros(NS,1);
end;

% European option value
EuroPrice = exp(-r*T)*mean(CF(:,NT));

% Time vector creation
  t(1) = 0;
  for k=2:NT
  t(k) = dt*(k-1);
  end;

% Work backwards through the stock prices until time j=2.
% We could work through to time j=1 but the regression will not be
% of full rank at time 1, so this is cleaner.
for j=NT-1:-1:1

    if strcmp(PutCall,'P')
		I = find(S(:,j) < K); % Indices for Puts
    else
        I = find(S(:,j) > K); % Indices for Calls
    end;

	  X = S(I,j);                     % X-vector = time j stock prices.
	
    Y = CF(I,j+1)*exp(-r*dt);           % Y-vector = time j+1 discounted CF.
    Z = zeros(length(X),NX);            % Design matrix for regression to predict cash flows
    for k=1:NX
        Z(:, k) = feval(XmatrixHandle{k}, X);
    end;
	beta = Z\Y;                       	% Regression parameters.
	P = Z*beta;                         % Regression predicted CF.
  
  
  % Calculating IS weight

  if strcmp(IS,'true') 
      for k = 1:NS
        sumGaussian(k) = sum(randomT(k,1:j));
      end;   
    weight = exp(-theta * sqrt(dt)* sumGaussian' - 0.5 * (theta^2) * t(j)*ones(NS,1));
  else
    weight = ones(NS,1);
  end;
      
      if strcmp(PutCall,'P')
        J = max(K - X, 0).*(weight(I)) > P;   % Early exercise for puts.
      else
        J = max(X - K, 0).*(weight(I)) > P;   % Early exercise for calls.
      end;
      
 	  E = I(J);                           % Stock price indices where immediate exercise is optimal.
	  C = setdiff((1:NS),E)';             % Stock price indices where continuation is optimal.
  
       if strcmp(PutCall,'P')
           CF(E,j) = max(K - X(J), 0).*weight(J);     % Replace with early exercise for puts
       else
           CF(E,j) = max(X(J)-K, 0).*weight(J);
       end;
       
      if strcmp(CV,'call')
        ContVar(E,j) = exp(-q*(T-t(j)))*normcdf((log(X(J)/K) + (r-q+v^2/2)*(T-t(j)))/v/sqrt((T-t(j)))).*X(J) - K*exp(-r*(T-t(j)))*normcdf((log(X(J)/K) + (r-q+v^2/2)*(T-t(j)))/v/sqrt((T-t(j))) - v*sqrt((T-t(j))));
      end;
    
	  CF(C,j) = exp(-r*dt)*CF(C,j+1);     % Continued CF are discounted back one period.
      ContVar(C,j) = exp(-r*dt)*ContVar(C,j+1);
  
end

% Calculate the variance of CF's
covarCFCV = cov(CF(:,1),ContVar(:,1));
varCV = var(ContVar(:,1));
if (varCV == 0)
  lambda = 1;
else 
  lambda = covarCFCV(2,1)/varCV ;
end;

display(lambda);

Variance = var(CF(:,1)-lambda*(ContVar(:,1)-mean(ContVar(:,1)))); 

% American option value
AmerPrice = exp(-r*dt)*mean(CF(:,1));
display(Variance);
display(AmerPrice);


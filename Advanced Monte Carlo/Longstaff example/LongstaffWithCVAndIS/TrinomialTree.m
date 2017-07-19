function y = TrinomialTree(Spot,K,r,q,v,T,N,EuroAmer,PutCall)

% Trinomial tree for Black Scholes
% Fabrice Douglas Rouah, FRouah.com and Volopta.com

%  Spot = Spot Price
%  K = Strike price
%  r = Risk free rate
%  q = Dividend yield
%  v = Volatility
%  T = Maturity
%  N = Number of time steps
%  EuroAmer = 'E'uropean or 'A'merican
%  PutCall = 'P'ut or 'C'all

% Trinomial tree parameters and probabilities.
dt = T/N;
u = exp(v*sqrt(2*dt));
d = 1/u;
pu = (exp((r-q)*dt/2) - exp(-v*sqrt(dt/2)))^2/(exp(v*sqrt(dt/2)) - exp(-v*sqrt(dt/2)))^2;
pd = (exp(v*sqrt(dt/2)) - exp((r-q)*dt/2))^2/(exp(v*sqrt(dt/2)) - exp(-v*sqrt(dt/2)))^2;
pm = 1 - pu - pd;

% Initialize the stock prices
S = zeros(2*N+1,N+1);
S(1,1) = Spot;

% Calculate all the stock prices.
for j=2:N+1
    for i=1:2*j-1
        S(i,j) = S(1,1)*u^j*d^i;
	end
end

% Initialize the option prices.
V = zeros(2*N+1,N+1);

% Calculate terminal option prices.
switch PutCall
	case 'C'
		V(:,N+1) = max(S(:,N+1) - K, 0);
	case 'P'
		V(:,N+1) = max(K - S(:,N+1), 0);
end

% Calculate Remaining entries for Calls and Puts
for j=N:-1:1
	for i=1:2*j-1
		switch EuroAmer
			case 'A'
				if strcmp(PutCall,'C')
					V(i,j) = max(S(i,j) - K, exp(-r*dt)*(pu*V(i,j+1) + pm*V(i+1,j+1) + pd*V(i+2,j+1)));
				else
					V(i,j) = max(K - S(i,j), exp(-r*dt)*(pu*V(i,j+1) + pm*V(i+1,j+1) + pd*V(i+2,j+1)));
				end
			case 'E'
				V(i,j) = exp(-r*dt)*(pu*V(i,j+1) + pm*V(i+1,j+1) + pd*V(i+2,j+1));
		end
	end
end

% Option price is at the first node.
y = V(1,1);

	
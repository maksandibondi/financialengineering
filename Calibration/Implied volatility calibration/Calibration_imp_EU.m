function[] = Calibration_imp_EU(filename,matur,S0)

%% Data reading
% data = dlmread(filename,';',1,0);
data = dlmread(filename,',',1,0);
display(data);
sz = size(data,1);

% Filtering the data with the maturity needed
startindex = 0; endindex = 0;
for n = 1:sz
    if (data(n,1) == matur) && (startindex == 0)
        startindex = n; endindex = n;
    elseif (data(n,1) == matur) && (startindex ~= 0)
        endindex = n;
    elseif (data(n,1) ~= matur) && (startindex ~= 0)
        break;
    end;
end;
display(startindex); display(endindex);
A = data(startindex:endindex,:);
display(A);

%% Parameters' arrays creation
T = A(1,1); % Maturity value  
if strcmp(S0,'no')
S0 = 1260.36*exp(-0.0217*T); % Initial price
end
display(S0);
r = A(1,7)/100; % Rate
K = A(:,2); % Strike
M = (A(:,3)+A(:,4))/2; % Market price data for call
% M = (A(:,5)+A(:,6))/2; % Market price data for put
e = 0.0001; % eps
sz = size(K);
maturity = T*ones(1,length(K));

%% Main program
for j=1:sz
    
%% Arbitrage condition
   if max(S0-K*exp(-r*T),0) < M < S0 % for call
   % if max(K*exp(-r*T)-S0,0) < M < S0 % for put       
    
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

   

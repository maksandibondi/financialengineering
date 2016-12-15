function [V, P, A, B, PL] = hedging(B, t0, T,discretization_num_t, r, sigma, K, S0, hedge_freq)

%% Variables initialization

B(1) = B;

t_0 = t0;  t_f = T;  % interval for t 

T = t_f-t_0;

delta_t = (t_f - t_0)/(discretization_num_t+1); % delta y

%% X and T arrays creation

% definition of the t-values on axis
t(1) = t_0;
for q = 2:1:discretization_num_t
    t(q) = t(q-1) + delta_t; 
end;

%% Main
V0 = BSTheory(S0,t_0,T,K,r,sigma(1));

V(1) = V0;

S(1) = S0;
A(1) = normal(d1(S0,t_0,T,K,r,sigma(1))); % delta
P(1) = A(1)*S(1) + B(1);

% For put
% A(1) = normal(d1(S0,t_0,T,K,r,sigma(1)))-1; % delta
% P(1) = A(1)*S(1) + B(1);


j = 1;
for k = 2:hedge_freq:discretization_num_t
    hedge_array(j) = k;
    j = j+1;
end;


for i = 2:discretization_num_t
    
    % Stock price simulation
    S(i) = S(i-1)*exp((r-sigma(i)^2/2)*delta_t+sigma(i)*sqrt(delta_t)*randn());
    P(i) = S(i)*A(i-1)+B(i-1)*(1+r*delta_t); % recalculating hedging prtf value using known alpha
    
    
    % Recalculating hedging portfolio
    if ismember(i,hedge_array)
    
        A(i) = normal(d1(S(i),t(i),T,K,r,sigma(i))); % New alpha recalculation
        
        % For put
        % A(i) = normal(d1(S(i),t(i),T,K,r,sigma(i)))-1; % New alpha recalculation

    else

        A(i) = A(i-1);
     
    end;
    
    % Parameters recalculation
    B(i) = P(i)-S(i)*A(i); % New Cash Deposit value
    
    % Discounting
    P(i) = P(i) + (V(1) - P(1))*exp(r*t(i)); % Hedging Portfolio discounting (for final info)
    
    V(i) = BSTheory(S(i),t(i),T,K,r,sigma(i)); 
    
end;
    
    PL = P(end)-V(end);
     
return;
    
    
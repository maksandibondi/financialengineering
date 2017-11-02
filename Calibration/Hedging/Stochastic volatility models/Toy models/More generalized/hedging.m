function [V_call, V_put, P_call, P_put, A_call, A_put, B_call, B_put, PLcall, PLput] = hedging(B_call, B_put, t0, T,discretization_num_t, r, sigma, K, S0,hedge_freq)

%% Variables initialization

B_call(1) = B_call;
B_put(1) = B_put;
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
V0_call = BSTheory_call(S0,t_0,T,K,r,sigma(1));
V0_put = BSTheory_put(S0,t_0,T,K,r,sigma(1));
V_call(1) = V0_call;
V_put(1) = V0_put;

S(1) = S0;
A_call(1) = normal(d1(S0,t_0,T,K,r,sigma(1))); % delta
P_call(1) = A_call(1)*S(1) + B_call(1);

A_put(1) = normal(d1(S0,t_0,T,K,r,sigma(1)))-1; % delta
P_put(1) = A_put(1)*S(1) + B_put(1);


j = 1;
for k = 2:hedge_freq:discretization_num_t
    hedge_array(j) = k;
    j = j+1;
end;


for i = 2:discretization_num_t
    
    % Stock price simulation
    S(i) = S(i-1)*exp((r-sigma(i)^2/2)*delta_t+sigma(i)*sqrt(delta_t)*randn());
    P_call(i) = S(i)*A_call(i-1)+B_call(i-1)*(1+r*delta_t); % recalculating hedging prtf value using known alpha
    P_put(i) = S(i)*A_put(i-1)+B_put(i-1)*(1+r*delta_t); % recalculating hedging prtf value using known alpha
    
    
    % Recalculating hedging portfolio
    if ismember(i,hedge_array)
    
        A_call(i) = normal(d1(S(i),t(i),T,K,r,sigma(i))); % New alpha recalculation
        A_put(i) = normal(d1(S(i),t(i),T,K,r,sigma(i)))-1; % New alpha recalculation

    else

        A_call(i) = A_call(i-1);
        A_put(i) = A_put(i-1);
     
    end;
    
    % Parameters recalculation
    B_call(i) = P_call(i)-S(i)*A_call(i); % New Cash Deposit value
    B_put(i) = P_put(i)-S(i)*A_put(i); % New Cash Deposit value
    
    % Discounting
    P_call(i) = P_call(i) + (V_call(1) - P_call(1))*exp(r*t(i)); % Hedging Portfolio discounting (for final info)
    P_put(i) = P_put(i) + (V_put(1) - P_put(1))*exp(r*t(i)); % Hedging Portfolio discounting (for final info)
    
    V_call(i) = BSTheory_call(S(i),t(i),T,K,r,sigma(i)); 
    V_put(i) = BSTheory_put(S(i),t(i),T,K,r,sigma(i));
    
    
end;
    
    PLcall = P_call(end)-V_call(end);
    PLput = P_put(end)-V_put(end);

    
return;
    
    
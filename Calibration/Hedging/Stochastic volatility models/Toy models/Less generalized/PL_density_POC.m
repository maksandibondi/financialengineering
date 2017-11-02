%% Variables intervals and discretization
clear; clc;

B0 = 1;

t0 = 0;
T = 5;
r = 0;
K = 10;
S0 = 10;
hedge_freq = 1;
num_of_iter = 50000;

alpha_level = 0.05;

discretization_num_t = 100; % number of discretization

delta_t = (T - t0)/(discretization_num_t+1); % delta y

%% X and T arrays creation

% definition of the t-values on axis
t(1) = t0;

% sigma1(1) = volatility1(); %stochastic vol
sigma1(1) = 0.3; %stochastic vol
sigma(1) = 0.3; % const vol

for q = 2:1:discretization_num_t
    t(q) = t(q-1) + delta_t; 
    % sigma1(q) = volatility1(); %stochastic vol
    % sigma1(q) = volatility2(sigma1(q-1)); %stochastic vol
    sigma1(q) = sigma1(q-1)+sigma1(q-1)*sqrt(delta_t)*randn(); %stochastic vol
    
    sigma(q) = sigma(q-1); % const vol
end;

figure;
plot(t,sigma);
hold on;
plot(t,sigma1);

%% Main  

s = size(hedge_freq,2);

for k = 1:s

    for i = 1:num_of_iter

        % [V_call, V_put, P_call, P_put, A_call, A_put, B_call, B_put, PLcall(i,k), PLput(i,k)] = hedging(B_call0, B_put0, t0, T, discretization_num_t, r, sigma, K, S0, hedge_freq(k));
        [V, P, A, B, PL(i,k)] = hedging_exact_repl(t0, T, discretization_num_t, r, sigma1, K, S0, hedge_freq(k));
        [V_const, P_const, A_const, B_const, PL_const(i,k)] = hedging_exact_repl(t0, T, discretization_num_t, r, sigma, K, S0, hedge_freq(k));
    end;

end;

for k = 1:s
    EV_PL(k) = sum(PL(:,k))/num_of_iter;
    Var_PL(k) = sum(PL(:,k).^2)/num_of_iter - EV_PL(k)^2;
    
    EV_PL_const(k) = sum(PL_const(:,k))/num_of_iter;
    Var_PL_const(k) = sum(PL_const(:,k).^2)/num_of_iter - EV_PL_const(k)^2;
end;

display(EV_PL);
display(Var_PL);
display(EV_PL_const);
display(Var_PL_const);


set_of_PL = -4:0.1:4;
sz = size(set_of_PL,2);
figure;
figure;
% VAR array creation
for k = 1:s
    indic(k) = 0;    
    indic_const(k) = 0;
end; 


for k = 1:s

    for i = 1:sz
    
        counter_density = 0;
        counter_proba = 0;
        counter_density_const = 0;
        counter_proba_const = 0;
        
        for j = 1:num_of_iter
            
            if (i ~= sz)
                
                if (PL(j,k) < set_of_PL(i+1)) && (PL(j,k) >= set_of_PL(i))
                    counter_density = counter_density+1;
                end;
                
                if (PL_const(j,k) < set_of_PL(i+1)) && (PL_const(j,k) >= set_of_PL(i))
                    counter_density_const = counter_density_const+1;
                end;
                
                if (PL(j,k) < set_of_PL(i+1))
                    counter_proba = counter_proba+1;
                end;
                
                if (PL_const(j,k) < set_of_PL(i+1))
                    counter_proba_const = counter_proba_const+1;
                end;              
        
            end;
            
        end;
        
        density(i,k) = (counter_density/num_of_iter)/0.1;
        proba(i,k) = counter_proba/num_of_iter;
        
        density_const(i,k) = (counter_density_const/num_of_iter)/0.1;
        proba_const(i,k) = counter_proba_const/num_of_iter;      
        
        if (proba(i,k)>=alpha_level) && (indic(k)==0)
            VAR(k) = set_of_PL(i);
            indic(k) = indic(k)+1;
        end;
        
        
        if (proba_const(i,k)>=alpha_level) && (indic_const(k)==0)
            VAR_const(k) = set_of_PL(i);
            indic_const(k) = indic_const(k)+1;
        end;
        
    end;

    figure(2);
    hold on;
    plot(set_of_PL(1:end-1), density((1:end-1),k));
    plot(set_of_PL(1:end-1), density_const((1:end-1),k));
    
    
    figure(3);
    hold on;
    plot(set_of_PL(1:end-1), proba((1:end-1),k));
    scatter(VAR(k)*ones(1,20),0.05:0.05:1);
    plot(set_of_PL(1:end-1), proba_const((1:end-1),k));
    scatter(VAR_const(k)*ones(1,20),0.05:0.05:1);
    
    
end;

display(VAR);
display(VAR_const);

figure(2);
legend('dens hedge freq = 1','dens const hedge freq = 1','dens hedge freq = 5','dens const hedge freq = 5');
figure(3);
legend('proba hedge freq = 1','VAR hedge freq = 1','proba const hedge freq = 1','VAR const hedge freq = 1','proba hedge freq = 5','VAR hedge freq = 5','proba const hedge freq = 5','VAR const hedge freq = 5');
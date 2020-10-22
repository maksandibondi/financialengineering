%% Variables intervals and discretization
clear; clc;

B0 = 1;
t0 = 0;
T = 5;
r = 0;
sigma = 0.3;
K = 10;
S0 = 10;
hedge_freq = [1 2 5 10];
num_of_iter = 500;

alpha_level = 0.1;

discretization_num_t = 100; % number of discretization

delta_t = (T - t0)/(discretization_num_t+1); % delta y

%% X and T arrays creation

% definition of the t-values on axis
t(1) = t0;
for q = 2:1:discretization_num_t
    t(q) = t(q-1) + delta_t; 
end;


%% Main  

s = size(hedge_freq,2);

for k = 1:s

    for i = 1:num_of_iter
        [V, P, A, B, PL(i,k)] = hedging_exact_repl(t0, T, discretization_num_t, r, sigma, K, S0, hedge_freq(k));
    end;

end;

for k = 1:s
    EV_PL(k) = sum(PL(:,k))/num_of_iter;
    Var_PL(k) = sum(PL(:,k).^2)/num_of_iter - EV_PL(k)^2;
end;

display(EV_PL);
display(Var_PL);

% plot(hedge_freq, EV_PLcall);
% 
% figure;
% plot(hedge_freq, Var_PLcall);


set_of_PL = -2:0.1:2;
sz = size(set_of_PL,2);
figure;
figure;
% VAR array creation
for k = 1:s
    indic(k) = 0;    
end; 


for k = 1:s

    for i = 1:sz
    
        counter_density = 0;
        counter_proba = 0;
        
        for j = 1:num_of_iter
            
            if (i ~= sz)
                
                if (PL(j,k) < set_of_PL(i+1)) && (PL(j,k) >= set_of_PL(i))
                    counter_density = counter_density+1;
                end;
                
                if (PL(j,k) < set_of_PL(i+1))
                    counter_proba = counter_proba+1;
                end;
        
            end;
            
        end;
        
        density(i,k) = (counter_density/num_of_iter)/0.1;
        proba(i,k) = counter_proba/num_of_iter;
        
        if (proba(i,k)>=alpha_level) && (indic(k)==0)
            VAR(k) = set_of_PL(i);
            indic(k) = indic(k)+1;
        end;
        
    end;

    figure(1);
    hold on;
    plot(set_of_PL(1:end-1), density((1:end-1),k));
    
    figure(2);
    hold on;
    plot(set_of_PL(1:end-1), proba((1:end-1),k));
    scatter(VAR(k)*ones(1,20),0.05:0.05:1);
    
end;

display(VAR);

figure(1);
legend('dens hedge freq = 1','dens hedge freq = 2','dens hedge freq = 5','dens hedge freq = 10');
figure(2);
legend('proba hedge freq = 1','proba hedge freq = 2','proba hedge freq = 5','proba hedge freq = 10');
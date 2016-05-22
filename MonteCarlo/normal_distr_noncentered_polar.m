function[expv1, expv2, variance_exp1, variance_exp2] = normal_distr_noncentered_polar(num_of_iter, num_of_rand, set_of_x)

%% Definition of the loops' size and polars values definition
sz = size(set_of_x, 2);
random_num = rand(num_of_iter,num_of_rand); % creating the matrix of the randoms (2,3 or more uniforms)
q = 1;
while q <= num_of_rand; % transforming the receiving randoms according to the algorithm of the polars
    V(:,q) = random_num(:,q)*2-1;
    q = q + 1;
end; 

%% Main program
i = 1;
while i<=sz
    
    q = 1;
    while q <= num_of_rand
    counter_x_frequency (q) = 0;
    counter_x_density (q) = 0;
    q = q + 1;
    end; 
    
    k = 1;
    while k<=num_of_iter    
    
        q = 1;
        S(k) = 0;
        while q <= num_of_rand
           S(k) = S(k) + V(k,q)^2;
           q = q + 1;
        end;
        
        if S(k) > 1
           q = 1;
           while q <= num_of_rand
                X(k,q) = NaN;
                q = q + 1;
           end;
        else
            q = 1;
            while q <= num_of_rand
                X(k,q) = V(k,q)*(-2*log(S(k))/S(k))^(1/2);
                q = q + 1;
            end;
        end;
           
%% PDF and DDF calculation
        
        if (i ~= sz)
            
            q = 1;
            while q<=num_of_rand
                
            if (X(k,q) < set_of_x(q,i+1)) % if alpha received < x then calculate how many numbers are equal to x therefore know their probability(frequency)
                counter_x_frequency(q) = counter_x_frequency(q) + 1; % adding one more number to current x frequency    
            end;
        
            if (X(k,q) < set_of_x(q,i+1)) && (X(k,q) >= set_of_x(q,i))
                counter_x_density(q) = counter_x_density(q) + 1;
            end;
            
            q = q + 1;
            end;
            
        end;

        k = k+1;
    end;
    
    q = 1;
    while q <= num_of_rand
        probability_of_x(i,q) = counter_x_frequency(q)/num_of_iter; % receive probability from frequency
        density_of_x(i,q) = counter_x_density(q) /num_of_iter; % receive density from counter
        q = q + 1;
    end;
    
    i = i+1;
    
end;

% display (X);
% display(probability_of_x);
% display(density_of_x);

%% Mean and variance calculations
a = X(:,1);
b = X(:,2);

i = 1; m1 = 1; 
while i <= num_of_iter
    if  isnan(a(i))
        i = i + 1;
    else
        X_NEW1(m1) = a(i);
        m1 = m1+1;
        i = i + 1;
    end;
end;
        
i = 1; m2 = 1;    
while i <= num_of_iter
    if isnan(b(i))
        i = i + 1;
    else 
        X_NEW2(m2) = b(i);
        m2 = m2+1;
        i = i + 1;
    end;
end;

% display(X_NEW1);
% display(X_NEW2);
expv1 = mean(X_NEW1);% Calculating the mean of the received Alphas
expv2 = mean(X_NEW2);
variance_exp1 = var(X_NEW1); % Calculating the variance of the received alphas
variance_exp2 = var(X_NEW2);

%% Plotting PDF and DDF

q = 1;
while q<=num_of_rand
scatter(set_of_x(q,:), probability_of_x(:,q), 'r');
hold on;
scatter(set_of_x(q,:), density_of_x(:,q), 'b');
q = q+1;
end;

return 
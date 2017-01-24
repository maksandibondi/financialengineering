function[alpha] = normal_TCL(num_of_iter, mu, sigma)

%% Definition of loops' size

num_of_rand = 10000;

%% Main program
    
    for k = 1:num_of_iter    
    
        Arithmeticmean_ubar = mean(rand(num_of_rand,1)); % Calculating of arithmetic mean
        alpha(k) = (Arithmeticmean_ubar-0.5)*(12*num_of_rand)^(0.5)*sigma + mu; % Looking for alpha non-centered
    
    end;
    
return;
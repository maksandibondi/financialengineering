function[alpha, expv, variance_exp] = discrete_distr (num_of_iter)

%% Definition of loops' size
mass_discrete = [0.1 0.5 1 1.5 2; 1/13 2/13 3/13 4/13 3/13];
sz = size(mass_discrete, 2);
random_num = rand (num_of_iter); 

%% Main program
k = 1;
while k<=num_of_iter    
    
           cumulative_sum = 0;
           i = 1;
           u = random_num(k);
           
           while i<=sz
           
           
           if u <= cumulative_sum + mass_discrete(2,i)
              alpha(k) = mass_discrete(1,i);
              break;
              
           else
              cumulative_sum = cumulative_sum + mass_discrete(2,i);
              i = i+1;
           end;
                
           end;
           
           k = k+1;
end;
        
%% Mean and variance calculations

expv = mean(alpha); % Calculating the mean of the received Alphas
variance_exp = var(alpha); % Calculating the variance of the received alphas

return 
function[variance_exp, expv] = poisson_distr(num_of_iter, lambda, set_of_x)

%% Definition of loops' size
i = 1;
sz = size(set_of_x,2);
random_num = rand(sz,num_of_iter);

%% Main program
while i<=sz
    
    counter_alpha_frequency = 0;
    counter_alpha_density = 0;
    k = 1;
    while k<=num_of_iter
          probability_simulated = random_num(i,k);
          probability_received = exp(-lambda); % initial value for probability
          sum = probability_received; % initial value for sum
          Alpha(i,k) = 0; %initial value for Alpha
          %% 
          while sum<probability_simulated
              Alpha(i,k) = Alpha(i,k) + 1;
              probability_received = probability_received * lambda/Alpha(i,k);
              sum = sum + probability_received;
          end;
          
          
%% PDF and DDF calculation
          if (i ~= sz)
                
            if (Alpha(i,k) < set_of_x(i+1)) % if alpha received < x then calculate how many numbers are equal to x therefore know their probability(frequency)
                counter_alpha_frequency = counter_alpha_frequency + 1; % adding one more number to current x frequency    
            end;
        
            if (Alpha(i,k) < set_of_x(i+1)) && (Alpha(i,k) >= set_of_x(i))
                counter_alpha_density = counter_alpha_density + 1;
            end;
        
            
          end;
          
          k = k+1;
    end;
    
    probability_of_x(i) = counter_alpha_frequency/num_of_iter; % receive probability from frequency
    density_of_x(i) = counter_alpha_density/num_of_iter; % receive density from counter
    
    i = i+1;
    
end
   
%% Mean and variance calculations

% display(Alpha(1:end,1:end));   
% display(probability_of_x);
% display(density_of_x);

expv = mean(Alpha(1,:)); % Calculating the mean of the received Alphas
variance_exp = var(Alpha(1,:)); % Calculating the variance of the received alphas

%% Plotting PDF and DDF
scatter(set_of_x, probability_of_x, 'b');
hold on;
scatter(set_of_x, density_of_x, 'r');

return


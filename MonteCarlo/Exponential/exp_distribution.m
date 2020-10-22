function[expv, variance_exp] = exp_distribution (num_of_iter, theta, set_of_x)

%% Definition of loops' size
u = rand(1,num_of_iter); % Random matrix creation
sz = size(set_of_x,2); 

%% Main program
i = 1; % Walking through discrete values loop
while i<=sz
    
    counter_alpha_frequency = 0; % counters initializing
    counter_alpha_density = 0; % counters initializing
    
    k = 1;
    while k <= num_of_iter
        
        Alpha(i,k) = log(abs(1-u(1,k)))/(-theta); % Formula for exponential distribution

%% PDF and DDF calculation
    if i~=sz % Check the last element in row
        
        if Alpha(i,k) < set_of_x(i+1)
            counter_alpha_frequency = counter_alpha_frequency + 1; % Counting frequency more if element smaller than first x
        end;
        
        if (Alpha(i,k) < set_of_x(i+1)) && (Alpha(i,k) >= set_of_x(i)) % Counting for density if element lays in the interval
            counter_alpha_density = counter_alpha_density + 1;
        end;
        
    end; 
    
    k = k+1;
    
    end;
    
    probability_of_x(i) = counter_alpha_frequency/num_of_iter; % receive probability from frequency
    density_of_x(i) = counter_alpha_density/num_of_iter; % receive density from counter
    
    i = i+1;
    
end;   

%% Mean and variance calculations
expv = mean(Alpha(1,:)); % Calculating the mean of the received Alphas
variance_exp = var(Alpha(1,:)); % Calculating the variance of the received alphas

%% Plotting PDF and DDF
scatter(set_of_x, probability_of_x, 'r');
hold on;
scatter(set_of_x, density_of_x, 'b');


return

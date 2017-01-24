clear; clc;

num_of_iter = 10000;
set_of_x = -3:0.1:3;
sz = size(set_of_x,2);

alpha = normal_polar(num_of_iter);


EV = mean(alpha);
Var = var(alpha);
display(EV);
display(Var);

for i = 1:sz
    
    counter_alpha_frequency = 0;
    counter_alpha_density = 0;
    
    for k = 1:num_of_iter
        
        if (i ~= sz)
                
            if (alpha(k) < set_of_x(i+1)) % if alpha received < x then calculate how many numbers are equal to x therefore know their probability(frequency)
                counter_alpha_frequency = counter_alpha_frequency+1; % adding one more number to current x frequency    
            end;
        
            if (alpha(k) < set_of_x(i+1)) && (alpha(k) >= set_of_x(i))
                counter_alpha_density = counter_alpha_density+1;
            end;
        
        end;
    
    end;
    
    probability_of_x(i) = counter_alpha_frequency/num_of_iter; % receive probability from frequency
    density_of_x(i) = counter_alpha_density/num_of_iter; % receive density from counter
    
    
end;

scatter(set_of_x, probability_of_x, 'r');
hold on;
scatter(set_of_x, density_of_x, 'b');
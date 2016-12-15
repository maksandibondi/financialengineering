clear; clc;

theta = 1;
num_of_iter = 100000;
set_of_x = -5:0.1:5;
sz = size(set_of_x,2);

alpha = normal_distr_reject(theta,num_of_iter);
alpha2 = -normal_distr_reject(theta,num_of_iter);
ALPHA = [alpha alpha2];

EV = mean(ALPHA);
Var = var(ALPHA);
display(EV);
display(Var);

for i = 1:sz
    
    counter_alpha_frequency = 0;
    counter_alpha_density = 0;
    
    for k = 1:2*num_of_iter
        
        if (i ~= sz)
                
            if (ALPHA(k) < set_of_x(i+1)) % if alpha received < x then calculate how many numbers are equal to x therefore know their probability(frequency)
                counter_alpha_frequency = counter_alpha_frequency+1; % adding one more number to current x frequency    
            end;
        
            if (ALPHA(k) < set_of_x(i+1)) && (ALPHA(k) >= set_of_x(i))
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
function [num_of_jumps] = Jumps_t(jump_freq,t)


num_of_jumps = 0;
product = rand();

while product >= exp(-jump_freq*t) 
    
    product = product*rand(); % comes from theory - product of uniform distr
    
    num_of_jumps = num_of_jumps + 1;
    
end;

return;
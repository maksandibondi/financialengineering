function [alpha] = normal_polar(num_of_iter)


%% Generating of U,Y - independent r.v
k = 1;
while k <= num_of_iter
    
        V1 = 2*rand()-1; 
        V2 = 2*rand()-1;
        S = V1^2 + V2^2;
 
        if S < 1 
            alpha(k) = V1*(-2*log(S)/S)^(1/2);
            k = k+1; 
        end;         
        
end;
  
return;
function [alpha] = normal_polar(num_of_iter)


random_num = rand(num_of_iter,2);

for q = 1:2
    V(:,q) = random_num(:,q)*2-1;
end;

%% Generating of U,Y - independent r.v
k = 1;
while k <= num_of_iter
    
        V1 = 2*rand()-1; 
        V2 = 2*rand()-1;
        S = V1^2 + V2^2;
 
        if S < 1 
            alpha(k) = V1*(-2*log(S)/S)^(1/2);
            alpha(k+1) = V2*(-2*log(S)/S)^(1/2);
            k = k+2; 
        end;         
        
end;
  
return;
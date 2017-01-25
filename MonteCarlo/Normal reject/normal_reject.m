function [alpha] = normal_reject(theta,num_of_iter)

%% Generating of U,Y - independent r.v
k = 1;
while k <= num_of_iter
        
        Y(k) = log(abs(1-rand()))/(-theta); % Formula for exponential distribution
        u(k) = rand();
        
        if u(k)<=exp(-((Y(k)-1)^2)/2)
            if rand()<0.5
               alpha(k) = Y(k);
            else
               alpha(k) = -Y(k); 
            end;
            k = k+1;
        end;
             
end;  

return;
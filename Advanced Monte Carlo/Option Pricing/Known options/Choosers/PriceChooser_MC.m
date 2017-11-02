function [chooser] = PriceChooser_MC(S0,K,r,sigma,T_dec,T)

num_of_iter = 500;
num_of_iter2 = 500;


for k = 1:num_of_iter
    S_dec_end(k) = S0*exp((r-(sigma^2)/2)*T_dec+sigma*sqrt(T_dec)*randn());
    
    % Choosing, put or call
    if (S_dec_end(k)>K)
        feature(k) = 1; % call
    else
        feature(k) = 0; % put
    end;
    
    % Simulation of share price and option payoffs between T_dec and T
    optsum = 0;
    if (feature(k) == 1)
        
        for i = 1:num_of_iter2
            S_end(k,i) = S_dec_end(k)*exp((r-(sigma^2)/2)*(T-T_dec)+sigma*sqrt(T-T_dec)*randn());
            optsum = optsum + max(S_end(k,i)-K,0)*exp(-r*(T-T_dec));
        end;
    
    else
        
        for i = 1:num_of_iter2
            S_end(k,i) = S_dec_end(k)*exp((r-(sigma^2)/2)*(T-T_dec)+sigma*sqrt(T-T_dec)*randn());
            optsum = optsum + max(K-S_end(k,i),0)*exp(-r*(T-T_dec));
        end;
        
    end;
    
    % Option price for k-th simulation of S
    price(k) = optsum/num_of_iter2;
        
    

end;

chooser = (sum(price))*exp(-r*T_dec)/num_of_iter;

return;


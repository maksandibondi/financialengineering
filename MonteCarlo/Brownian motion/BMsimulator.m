function[W] = BMsimulator(T,discretization_num_t,method)


delta_t = T/discretization_num_t;
t(1) = 0;
for q = 2:discretization_num_t
    t(q) = t(q-1) + delta_t; 
end;


if (strcmp(method,'TCL'))
    
    for i = 1:discretization_num_t
        
        alpha = normal_TCL(num_of_iter, 0, 1);
        W(i) = sqrt(t(i))*sum(alpha)/sqrt(num_of_iter);
        
    end;
 
elseif (strcmp(method,'Reject'))
    
    for i = 1:discretization_num_t
        
        alpha1 = normal_reject(1,num_of_iter/2);
        alpha2 = -normal_reject(1,num_of_iter/2);
        alpha = [alpha1 alpha2];
        
        W(i) = sqrt(t(i))*sum(alpha)/sqrt(num_of_iter);
        
    end;
    
elseif (strcmp(method,'Polar'))
    
    W(1) = 0;
    
    for i = 2:discretization_num_t
        
        gaussian = normal_polar(i); 
        
        W(i) = delta_t*sum(gaussian);
        
    end;
    
end;
 

return;
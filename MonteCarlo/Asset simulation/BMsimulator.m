function[W] = BMsimulator(T,discretization_num_t,method)


delta_t = T/discretization_num_t;
t(1) = 0;
for q = 2:discretization_num_t
    t(q) = t(q-1) + delta_t; 
end;


if (strcmp(method,'TCL'))
    
    for i = 1:discretization_num_t
        
        W(i) = sqrt(t(i))*sum(normal_TCL(i, 0, 1))/sqrt(i);
        
    end;
 
elseif (strcmp(method,'Reject'))
    
    for i = 1:discretization_num_t
        
        W(i) = sqrt(t(i))*sum(normal_reject(1,i))/sqrt(i);
        
    end;
    
elseif (strcmp(method,'Polar'))
    
    for i = 1:discretization_num_t
        
    alpha = normal_polar(i); 
        
        W(i) = sqrt(t(i))*sum(alpha)/sqrt(i);
        
    end;
    
end;
 

return;
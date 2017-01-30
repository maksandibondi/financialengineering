function[W] = BMsimulator(T,discretization_num_t,method)


delta_t = T/discretization_num_t;
t(1) = 0;
for q = 2:discretization_num_t
    t(q) = t(q-1) + delta_t; 
end;


if (strcmp(method,'TCL'))
    
    W(1) = 0;
    
    gaussian = normal_TCL(discretization_num_t,0,1);
    
    for i = 2:discretization_num_t
        
      W(i) = sqrt(delta_t)*sum(gaussian(1:i));
      % W(i) = sqrt(t(i))*sum(gaussian(1:i))/sqrt(i);
        
    end;
 
elseif (strcmp(method,'Reject'))
    
    W(1) = 0;

    gaussian = normal_reject(1,discretization_num_t,0,1);
    
    for i = 2:discretization_num_t
        
      % W(i) = sqrt(delta_t)*sum(gaussian(1:i));
      W(i) = sqrt(t(i))*sum(gaussian(1:i))/sqrt(i);
        
    end;
    
elseif (strcmp(method,'Polar'))
    
    W(1) = 0;
    
    gaussian = normal_polar(discretization_num_t,0,1);
    
    for i = 2:discretization_num_t
        
        W(i) = sqrt(delta_t)*sum(gaussian(1:i));
        % W(i) = sqrt(t(i))*sum(gaussian(1:i))/sqrt(i);
        
    end;
    
elseif (strcmp(method,'Random walk'))
    
    W(1) = 0;
    r = rand();
    x(1) = (r>0.5)*(-1)+(r<=0.5)*(1);
    for i = 2:discretization_num_t
        
        if rand()>0.5
            x(i) = -1;
        else
            x(i) = 1;
        end;
        
        W(i) = sqrt(t(i))*sum(x(1:i))/sqrt(i);
        % W(i) = sqrt(t(i))*sum(gaussian(1:i))/sqrt(i);
        
    end;
    
end;
 

return;
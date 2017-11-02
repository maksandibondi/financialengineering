function [Y] = discreteY(v1,p1,v2,p2)

    u = rand();
    
    if (u >= p1)
        
        Y = v2;
        
    else 
        
        Y = v1;
        
    end;
    
    

return;


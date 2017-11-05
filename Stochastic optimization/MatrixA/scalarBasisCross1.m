function [integral] = scalarBasisCross1(knots, m, x)
		
integral = 0;
sz = size(x,2);

%alculation of the integral over x values
		
% loop for fm
for i = 1:sz

    if (x(i) <= knots(m - 2)) 
				val(i) = 0;
    elseif (x(i) >= knots(m - 2) && x(i) <= knots(m - 1)) 
				temp = (x(i) - knots(m - 2)) / (knots(m - 1) - knots(m - 2));
				val(i) = (temp*temp*temp) / 6;
    elseif (x(i) >= knots(m - 1) && x(i) <= knots(m)) 
				temp = (x(i) - knots(m - 1)) / (knots(m) - knots(m-1));
				val(i) = (1 / 6)*(1 + 3 * temp + 3 * (temp*temp) - 3 * (temp*temp*temp));
            
    elseif (x(i) >= knots(m) && x(i) <= knots(m+1)) 
				temp = (x(i) - knots(m)) / (knots(m+1) - knots(m));
				val(i) = (1 / 6)*(4 - 6 * (temp*temp) + 3 * (temp*temp*temp));
           
    elseif (x(i) >= knots(m+1) && x(i) <= knots(m+2)) 
				temp = (x(i) - knots(m+1)) / (knots(m+2)- knots(m+1));
				val(i) = (1 / 6)*(1 - temp)^3;
            
    elseif (x(i) >= knots(m+2))
				val(i) = 0;
            
    end;
            
end;

%loop for fm-1
for i = 1:sz

    if (x(i) <= knots(m - 2 - 1)) 
				val(i) = val(i)* 0;
				integral = integral + val(i);
			
            elseif (x(i) >= knots(m - 2 - 1) && x(i) <= knots(m - 1 - 1)) 
				temp = (x(i) - knots(m - 2 - 1)) / (knots(m - 1 - 1) - knots(m - 2 - 1));
				val(i) = val(i)*(temp*temp*temp) / 6;
				integral = integral + val(i);
			
            elseif (x(i) >= knots(m - 1 - 1) && x(i) <= knots(m - 1)) 
				temp = (x(i) - knots(m - 1 - 1)) / (knots(m - 1) - knots(m - 1 - 1) );
				val(i) = val(i)*(1 / 6)*(1 + 3 * temp + 3 * (temp*temp) - 3 * (temp*temp*temp));
				integral = integral + val(i);
			
            elseif (x(i) >= knots(m - 1) && x(i) <= knots(m + 1 - 1)) 
				temp = (x(i) - knots(m - 1)) / (knots(m + 1 - 1) - knots(m - 1));
				val(i) = val(i)*(1 / 6)*(4 - 6 * (temp*temp) + 3 * (temp*temp*temp));
				integral = integral + val(i);
			
            elseif (x(i) >= knots(m + 1 - 1) && x(i) <= knots(m + 2 - 1)) 
				temp = (x(i) - knots(m + 1 - 1)) / (knots(m + 2 - 1) - knots(m + 1 - 1));
				val(i) = val(i)* (1 / 6)*(1 - temp)^3;
				integral = integral + val(i);
		
            elseif (x(i) >= knots(m + 2 - 1)) 
				val(i) = val(i)*0;
				integral = integral + val(i);
    end;
            
end;
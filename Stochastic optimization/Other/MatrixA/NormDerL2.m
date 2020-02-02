function [val] = NormDerL2(knots,m,x)

        val = 0;
		sz = size(x);

		% calculation of the integral of fm^2 over K values
for i = 1:sz 

    if (x(i) < knots(m - 2))
				val = val + 0;
    elseif (x(i) >= knots(m - 2) && x(i) <= knots(m - 1)) 
				temp = (x(i) - knots(m - 2)) / (knots(m - 1) - knots(m - 2));
				val = val + temp*(1/((knots(m - 1) - knots(m - 2))^2));
    elseif (x(i) >= knots(m - 1) && x(i) <= knots(m))
				temp = (x(i) - knots(m - 1)) / (knots(m) - knots(m - 1));
				val = val + (1/(knots(m) - knots(m - 1))^2)-3*temp/((knots(m) - knots(m - 1))^2);
                
    elseif (x(i) >= knots(m) && x(m) <= knots(m + 1)) 
				temp = (x(i) - knots(m)) / (knots(m + 1) - knots(m));
				val = val + (-2/(knots(m+1) - knots(m))^2)+3*temp/((knots(m+1) - knots(m))^2);
                
    elseif (x(i) >= knots(m + 1) && x(i) <= knots(m + 2))
				temp = (x(i) - knots(m + 1)) / (knots(m + 2) - knots(m + 1));
				val = val + (1/(knots(m+2) - knots(m+1))^2)-temp/((knots(m+2) - knots(m+1))^2);
                
    elseif (x(i) >= knots(m + 2)) 
				val = val + 0;
    end;
             
end;
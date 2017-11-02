function [val] = NormL2(knots,m,x)

        val = 0;
		sz = size(x);

		% calculation of the integral of fm^2 over K values
for i = 1:sz 

    if (x(i) < knots(m - 2))
				val = val + 0;
    elseif (x(i) >= knots(m - 2) && x(i) <= knots(m - 1)) 
				temp = (x(i) - knots(m - 2)) / (knots(m - 1) - knots(m - 2));
				val = val + ((temp*temp*temp) / 6)^2;
    elseif (x(i) >= knots(m - 1) && x(i) <= knots(m))
				temp = (x(i) - knots(m - 1)) / (knots(m) - knots(m - 1));
				val = val + ((1 / 6)*(1 + 3 * temp + 3 * (temp*temp) - 3 * (temp*temp*temp)))^2;
    elseif (x(i) >= knots(m) && x(m) <= knots(m + 1)) 
				temp = (x(i) - knots(m)) / (knots(m + 1) - knots(m));
				val = val + ((1 / 6)*(4 - 6 * (temp*temp) + 3 * (temp*temp*temp)))^2;
    elseif (x(i) >= knots(m + 1) && x(i) <= knots(m + 2))
				temp = (x(i) - knots(m + 1)) / (knots(m + 2) - knots(m + 1));
				val = val + ((1 / 6)*((1 - temp))^3)^2;
    elseif (x(i) >= knots(m + 2)) 
				val = val + 0;
    end;
             
end;
function [ out ] = gfunc( y , ymax, y_0)

out = 0.5*log((ymax+y)/(ymax-y))+y_0;

end


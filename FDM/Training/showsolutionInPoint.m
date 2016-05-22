function [solution] = showsolutionInPoint(x,t, u, x_0, delta_x, t_0, delta_t)

x_order = round((x-x_0)/delta_x); % obtaining index by x in u-matrix
t_order = round((t-t_0)/delta_t); % obtaining index by t in u-matrix

if (t_order>0 && x_order>0)
solution = u(t_order, x_order); % receiving element from matrix
elseif (t_order==0 && x_order>0)
solution = u(1, x_order); % receiving element from matrix
elseif (t_order>0 && x_order==0)
solution = u(t_order, 1); % receiving element from matrix       
else
display('mistake');
end;
        

return


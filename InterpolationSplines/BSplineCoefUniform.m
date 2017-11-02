function [val] = BSplineCoefUniform(X,x,m)

if (x <= X(m-2))
  val = 0;
elseif (x>=X(m-2) && x<=X(m-1))
  a = (x-X(m-2))/(X(m-1)-X(m-2));
  val = (a^3)/6;
elseif (x>=X(m-1) && x<=X(m))
  a =(x-X(m-1))/(X(m)-X(m-1));
  val = (1/6)*(1+3*a+3*(a^2)-3*(a^3));
elseif (x>=X(m) && x<=X(m+1))
  a =(x-X(m))/(X(m+1)-X(m));
  val = (1/6)*(4-6*(a^2)+3*(a^3));
elseif (x>=X(m+1) && x<=X(m+2))
  a =(x-X(m+1))/(X(m+2)-X(m+1));
  val = (1/6)*(1-a)^3;
elseif (x>=X(m+2))
  val = 0;
end;

 
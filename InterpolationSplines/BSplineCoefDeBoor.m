function [val] = BSplineCoefDeBoor(X,m,p,u)

if (p == 1)
  if (X(m)<=u && u<X(m+1))
  val = 1;
  else
  val = 0;
  end;
else 
val = ((u-X(m))/(X(m+p)-X(m)))*BSplineCoefDeBoor(X,m,p-1,u) + ((X(m+p+1)-u)/(X(m+p+1)-X(m+1)))*BSplineCoefDeBoor(X,m+1,p-1,u);
end;
 
 
 
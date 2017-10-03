clear; clc;

X = [0:1:10]; % knot vector
x = rand(11,1); % K-vector sz = 101 control points
sz = size(x,1);
p = 4;

% Known form for cubic spline B-spline
res(1) = x(1);
%res(101) = x(end);
K = 0:0.1:10;
iter = 2;
for k = 0.1:0.1:10 % iterate through points
  spl = res(1);
  for m = 3:sz-2
    %cspline = BSplineCoefDeBoor(X,m,p,k);
    cspline = BSplineCoefUniform(X,k,m);
    spl = spl + cspline*x(m);
  endfor;
  res(iter) = spl;
  iter++;
endfor;

figure;
plot(K,res);


% General (non uniform knot) approach using De Boors formula for B-spline
%X = [0 0.8 1.5 3.5 4 5 8 8.5 10 11 13];
X = [0:1:10];
res_(1) = x(1);
%res_(101) = x(end);
K = 0:0.1:10;
iter = 2;
for k = 0.1:0.1:10 % iterate through points
  spl = res(1);
  for m = 1:sz-p-1
    cspline = BSplineCoefDeBoor(X,m,p,k);
    spl = spl + cspline*x(m);
  endfor;
  res_(iter) = spl;
  iter++;
endfor;

size(K)
size(res_)

figure;
plot(K,res_);



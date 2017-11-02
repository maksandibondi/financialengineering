clear; clc;

X = [0:1:10]; % knots sz = 11
x = rand(11,1); % control points (thetas) sz = n-p-1 = 7 values will be used 
sz = size(x,1);
p = 4;

%res(101) = x(end);

K = 0:0.1:10;
%iter = 2;
iter = 1;
for k = 0:0.1:10 % iterate through points where interpolate
  spl = 0;
  for m = 3:sz-2
    if (m==3)
        cspline = 1;
    else
        cspline = BSplineCoefUniform(X,k,m);
    end;
    spl = spl + cspline*x(m);
  end;
  res(iter) = spl;
  iter=iter+1;
end;

figure;
scatter(X,x); hold on;
plot(K,res);


% General (non uniform knot) approach using De Boors formula for B-spline
%X = [0 0.8 1.5 3.5 4 5 8 8.5 10 11 13];
X = [0:1:10];
%res_(1) = x(1);
%res_(101) = x(end);
K = 0:0.1:10;
iter = 1;
for k = 0:0.1:10 % iterate through points
  spl = 0;
  for m = 1:sz-p-1
    cspline = BSplineCoefDeBoor(X,m,p,k);
    spl = spl + cspline*x(m);
  end;
  res_(iter) = spl;
  iter = iter+1;
  end;

size(K)
size(res_)

%figure;
plot(K,res_);



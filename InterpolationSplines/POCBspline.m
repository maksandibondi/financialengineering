clear; clc;

X = [0:1:10]; % knots sz = 11
x = rand(11,1); % control points (thetas) sz = n-p-1 = 7 values will be used 
sz = size(x,1);
p = 4;

% K = 0:0.001:10;
% iter = 1;
% for k = 0:0.001:10 % iterate through points where interpolate
%   spl = 0;
%   for m = 3:sz-2
%     cspline = BSplineCoefUniform(X,k,m);
%     cs(iter,m) = cspline;
%     spl = spl + cspline*x(m);
%   end;
%   res(iter) = spl;
%   iter=iter+1;
% end;
% 
% 
% figure;
% scatter(X,x); hold on;
% plot(K,res); xlabel('knots'); ylabel('sigma(T,x)'); 
% 
% figure;
% for m_ = 3:sz-2
% plot(K,cs(:,m_)); xlabel('knots'); ylabel('Basis(x)'); hold on;
% end;

% General (non uniform knot) approach using De Boors formula for B-spline
X = [0 1.5 3 4.5 6 6.25 6.5 6.75 7 7.5 10 13 14.5 16];
%X = [0:1:10];
%res_(1) = x(1);
%res_(101) = x(end);
K = 3:0.001:13;
iter = 1;
for k = 3:0.001:13 % iterate through points in which interpolate
  spl = 0;
  for m = 1:sz-p-1
    cspline = BSplineCoefDeBoor(X,m,p,k);
    cs(iter,m) = cspline;
    spl = spl + cspline*x(m);
  end;
  res_(iter) = spl;
  iter = iter+1;
end;

size(K)
size(res_)

figure;
plot(K,res_);

figure;
for m_ = 1:sz-p-1
plot(K,cs(:,m_)); hold on;
end;


% Non-uniform knot (program from internet)
X = [0 3 6 6.25 6.5 6.75 7 7.5 10 13 16]; % knots
U = 6.25:0.001:7.5; % points where interpolate
x = rand(1,7); % control points
[C,u] = BSplineCoefDeBoorEx(p,X,x,U);
figure;
plot(U,C);

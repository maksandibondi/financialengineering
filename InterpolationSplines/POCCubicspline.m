clear; clc;

ctrlX = [0:1:5]; % knot vector
ctrlY = [0.5 1 0.5 2 1 0.5]; % K-vector sz = 6 control points
n = size(ctrlX,2);

h = zeros(n-1,1);
for i = 1:n-1
  h(i) = ctrlX(i+1)-ctrlX(i); % size of H is equal to n-1
end;
c1 = 0;
c_last = 0;

c_n = ThomasAlgo(ctrlY, ctrlX, h); % calculating c system using tridiagonal system solver

c = [c1 c_n c_last];

sz = size(c,2);
d_ = zeros(sz-1,1);

%% Calculating vectors b and d from vector c
for i = 1:sz-1
  d_(i) = (c(i+1)-c(i))/(3*h(i));
end;
d = [d_' 0]; % last d = cn/(-3hn) = 0 as last cn = 0

b_ = zeros (sz-1,1);
for i = 1:sz-1
  b_(i) = (ctrlY(i+1)-ctrlY(i))/h(i) - c(i)*h(i) - d(i)*(h(i)^2);
end;
b = [b_' b_(sz-1)];







x = [0:0.1:6];
iter = 1;

for k = 1:n-1 % iterator through knot intervals

for i = 1+10*(k-1):10*k % iterator through polynom points
  S(iter) = ctrlY(k)+b(k+1)*(x(i)-ctrlX(k))+c(k+1)*(x(i)-ctrlX(k))^2+d(k+1)*(x(i)-ctrlX(k))^3; % building a 3-order polynom for each knot interval
  iter = iter+1; 
end;

k = k+1;

end;
size(x)
size(S)

plot(x(1:41),S(1:41));
x = [1 1 1 1 1 1 1; 0 0 1 1 2 2 3; 0 1 0 1 0 1 2];
y = [0 1 0 1 1 0 1];
x = x';
y = y';

m = size(x,1); % number of training examples
order = size(x,2); % number of features

xnorm = ones(m,order);
for k = 1:m % nomalization
    for i = 2:order
    xnorm(k,i) = (x(k,i)-mean(x(:,i)))/(max(x(:,i))-min(x(:,i)));
    end;
end;

display(xnorm);

theta = [1; 1; 1]; % initial guess of theta
learning_rate = 0.01;

gd = [2;2;3];
% while abs(gd) > 0.0001 % for non-linear logarithmic 
% while abs(gd(2)) > 0.0000001 || abs(gd(3))> 0.0000001 % for linear with gradient-vector

while abs(gd(2)) > 0.0000001 || abs(gd(3))> 0.0000001 % for linear with gradient-vector

    gd = gradient (theta,xnorm,y);
    theta = theta-learning_rate*gd';
   
end;

display(theta);
% display(gd);

h = (theta')*xnorm';
figure;
scatter3(xnorm(:,2), xnorm(:,3),y, 'b', 'filled');
xlabel('x1');
ylabel('x2');
zlabel('y');
hold on;
% x1fit = min(xnorm(:,2)):10:max(xnorm(:,2));
% x2fit = min(xnorm(:,3)):10:max(xnorm(:,3));
% [x1fit, x2fit] = meshgrid(x1fit,x2fit);
% mesh(x1fit,x2fit,h);
scatter3(xnorm(:,2), xnorm(:,3),h,'r', 'filled');

newx = [1;3;2];

forecast = (theta')*[newx(1);(newx(2)-mean(x(:,2)))/(max(x(:,2))-min(x(:,2)));(newx(3)-mean(x(:,3)))/(max(x(:,3))-min(x(:,3)))];
display(forecast);

error = sum((y-h').^2)/m;
display(error);


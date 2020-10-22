clear; clc;
X = [-1;1];
precision = 0.001;
% All minimas can be observed going rfom this starting point when x is from
% 1.8 till 2 (the best is obtained for x = 1.95) and y is close to zero. Then we go to minima when function value goes to zero.
syms f(x,y) H G FX FY;

%% Function to optimize
f = (x^2+y-11)^2+(x+y^2-7)^2; % Himmelblow
% f = 10*(x^2-y)^2+(x-1)^2; % Rosenbrock
% f = x^2+y^2+x*y-y+4; first equation
% f = (x^2+y-11)^2+(x+y^2-7)^2; % Himmelblow
% f = (1.5-x+x*y)^2+(2.25-x+x*(y^2))^2+(2.625-x+x*(y^3))^2; % Beale func

%% Symbolic Gradient and Hessian matrix calculations
FX = diff(f,x); % calcul of gradient
FY = diff(f,y); % calcul of gradient
G = [FX;FY];
H = [diff(f,x,2) diff(diff(f,x),y); diff(diff(f,y),x) diff(f,y,2)]; % calcul of hessian

%% Main algorithm
i = 1; eps = 1;
while (eps>= precision & i <=1000)
  
    
    StepValue(:,i) = X;
   
    x = X(1); y = X(2);
    grad = eval(G); % eval of grad
    Hess = eval(H); % eval of hessian
    
    Y = X - Hess\grad;
    eps = abs(Y - X);
    X = X - Hess\grad;
    
    i = i+1;
end;

%% Steps and final optimal function value displaying
display(StepValue);
display(X);
x = X(1); y = X(2);
min = eval(f);
display(min);

%% Griding and ploting the path
[X,Y] = meshgrid(-5:0.01:5, -5:0.01:5);
contour(X,Y,((X.^2)+Y-11).^2+(X+Y.^2-7).^2, [0,0.01,0.03, 0.05,0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 1,2,3, 4, 4.5, 5, 7, 8, 9, 10, 15, 20, 25, 30, 40, 70, 100, 120, 135, 150, 200, 250]); %Himmelblow
% contour(X,Y,(X.^2)+(Y.^2)+X.*Y-Y+4, [0,0.01,0.03, 0.05,0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 1,2,3, 4, 4.5, 5, 7, 8, 9, 10, 15, 20, 25, 30, 40, 70, 100, 120, 135, 150, 200, 250]); % first equation
% contour(X,Y,10*((X.^2)-Y).^2+(X-1).^2, [0,0.01,0.03, 0.05,0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 1,2,3, 4, 4.5, 5, 7, 8, 9, 10, 15, 20, 25, 30, 40, 70, 100, 120, 135, 150, 200, 250]); %rosenbrock
% contour(X,Y,((X.^2)+Y-11).^2+(X+Y.^2-7).^2, [0,0.01,0.03, 0.05,0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 1,2,3, 4, 4.5, 5, 7, 8, 9, 10, 15, 20, 25, 30, 40, 70, 100, 120, 135, 150, 200, 250]); %Himmelblow
% contour(X,Y,(1.5-X+X.*Y).^2+(2.25-X+X.*(Y.^2)).^2+(2.625-X+X.*(Y.^3)).^2, [0,0.01,0.03, 0.05,0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 1,2,5,10]); % Beale
hold on;
plot(StepValue(1,:), StepValue(2,:), '-o');
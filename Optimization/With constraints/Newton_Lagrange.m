clear; clc; 
precision = 0.0001;
X = [10;10;2];
syms f(x,y, lamb1); syms h(x,y,lamb1); syms L; syms HX HY LX LY LLamb;

%% Function to optimize
f = (1/2)*(x^2+y^2-x*y); % problem 2
% f = (x^2+y-11)^2+(x+y^2-7)^2; % Himmelblow
% f = 10*(x^2-y)^2+(x-1)^2; % Rosenbrock
% f = (1/2)*(x^2+y^2); % problem1 
% f = (1/2)*(x^2+y^2-x*y); % problem 2

%% Constraints and lagrangian
h = 4*x+y-20; 
L = f+lamb1*h; % Lagrangian

%% Symbolic Gradient and Hessian matrix calculations for lagrangian and constraints 
LX = diff(L,x); % calcul of gradient of lagrangian symb 
LY = diff(L,y); % calcul of gradient of lagrangian symb
LLamb = diff(L,lamb1); % calcul of gradient of lagrangian symb
HX = diff(h,x); % calcul of gradient of constraint symb
HY = diff(h,y); % calcul of gradient of constraint symb
H = [diff(L,x,2) diff(diff(L,x),y); diff(diff(L,y),x) diff(L,y,2)]; % calcul of hessian of Lagrangian symb

%% Main algorithm
i = 1; eps = 1;
while (eps>= precision & i <=1000)
    
    StepValue(:,i) = X;
    display(StepValue);
    x = X(1); y = X(2); lamb1 = X(3);
    grad = eval([HX;HY]); % eval of grad
    Hess = eval(H); % eval of hessian of Lagrangian
  
    M = [Hess grad; grad' 0];
    V = eval([-LX;-LY;-LLamb]);
 
    Direction = M\V;
    
    Y = X + Direction;
    eps = abs(Y - X);
    X = X + Direction;
    
    i = i+1;
end;

%% Steps and final optimal function value displaying
% display(StepValue);
display(X);
x = X(1); y = X(2);
min = eval(f);
display(min);

%% Griding and ploting the path
[X,Y] = meshgrid(-5:0.01:5, -5:0.01:5);
contour(X,Y,(1/2)*(X.^2+Y.^2-X.*Y), [0,0.01,0.03, 0.05,0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 1,2,3, 4, 4.5, 5, 7, 8, 9, 10, 15, 20, 25, 30, 40, 70, 100, 120, 135, 150, 200, 250]); %Problem2

% contour(X,Y,(1/2)*(X.^2+Y.^2-X.*Y), [0,0.01,0.03, 0.05,0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 1,2,3, 4, 4.5, 5, 7, 8, 9, 10, 15, 20, 25, 30, 40, 70, 100, 120, 135, 150, 200, 250]); %Problem2
% contour(X,Y,10*((X.^2)-Y).^2+(X-1).^2, [0,0.01,0.03, 0.05,0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 1,2,3, 4, 4.5, 5, 7, 8, 9, 10, 15, 20, 25, 30, 40, 70, 100, 120, 135, 150, 200, 250]); %rosenbrock
% contour(X,Y,((X.^2)+Y-11).^2+(X+Y.^2-7).^2, [0,0.01,0.03, 0.05,0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 1,2,3, 4, 4.5, 5, 7, 8, 9, 10, 15, 20, 25, 30, 40, 70, 100, 120, 135, 150, 200, 250]); %Himmelblow
% contour(X,Y,(1/2)*(X.^2+Y.^2), [0,0.01,0.03, 0.05,0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 1,2,3, 4, 4.5, 5, 7, 8, 9, 10, 15, 20, 25, 30, 40, 70, 100, 120, 135, 150, 200, 250]); % problem 1

hold on;
plot(StepValue(1,:), StepValue(2,:), '-o');
a = -5:0.01:5;
plot(a,20-4*a);
    
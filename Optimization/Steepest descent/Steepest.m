clear; clc;
X = [1;-1]; % starting point
precision = 0.001;
syms f(x,y); syms FX FY H gradient;

%% Function to optimize
f = (x^2+y-11)^2+(x+y^2-7)^2; % Himmelblow
% f = 10*(x^2-y)^2+(x-1)^2; % Rosenbrock
% f = x^2+y^2+x*y-y+4; % first equation
% f = (x^2+y-11)^2+(x+y^2-7)^2; % Himmelblow

%% Symbolic gradient and Hessian matrix calculations
FY = diff(f,y); FX = diff(f,x);
gradient = [FX;FY];
H = [diff(FX,x) diff(FX,y); diff(FY,x) diff(FY,y)];

%% Main algorithm
i = 1; % iteration counter
x = X(1); y = X(2);
while sqrt((eval(gradient)')*eval(gradient))>precision
    
    x = X(1); y = X(2);
    grad = eval(-gradient);
    Hess = eval(H);
    
    t = (grad')*grad/((grad')*Hess*grad);
    
    StepValue(:,i) = X;
    X = X + t*grad;
    
    display(StepValue(:,i));
    i = i+1;
end;

%% Steps and final optimal function value displaying 
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
hold on;
plot(StepValue(1,:), StepValue(2,:), '-o');

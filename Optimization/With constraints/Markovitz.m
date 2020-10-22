clear; clc;
precision = 0.0001;
X = [1/3 1/3 1/3 1 1]; % initial weights where 4th and 5th elements are lmbda1 and 2
E = [1 2 3]; % return vector
mu = 1.4; % return wanted
covar = [0.5 0.3 0.05; 0.3 0.3 0.1; 0.05 0.1 0.8]; % covariance matrix historically obtained
syms f(x,y,z); syms h(x,y,z); syms L(x,y,z,lamb1,lamb2); syms HH HH2 LX LY LZ LLamb1 LLamb2; 

%% Function to optimize
f = (1/2)*([x y z]*covar*[x;y;z]); 

%% Constraints and lagrangian
h1 =  x*E(1)+y*E(2)+z*E(3)-mu; % constraint 1
h2 = x+y+z-1; % constraint 2
L = f+lamb1*h1+lamb2*h2; % Lagrangian

%% Symbolic Gradient and Hessian matrix calculations for lagrangian and constraints
LX = diff(L,x); % calcul of gradient of lagrangian symb
LY = diff(L,y); % calcul of gradient of lagrangian symb
LZ = diff(L,z); % calcul of gradient of lagrangian symb
LLamb1 = diff(L,lamb1); % calcul of gradient of lagrangian symb
LLamb2 = diff(L,lamb2); % calcul of gradient of lagrangian symb
HH(1) = diff(h1,x); % calcul of gradient of constraint symb
HH(2) = diff(h1,y); % calcul of gradient of constraint symb
HH(3) = diff(h1,z); % calcul of gradient of constraint symb
HH2(1) = diff(h2,x); % calcul of gradient of constraint symb
HH2(2) = diff(h2,y); % calcul of gradient of constraint symb
HH2(3) = diff(h2,z); % calcul of gradient of constraint symb
H = [diff(LX,x) diff(LX,y) diff(LX,z); diff(LY,x) diff(LY,y) diff(LY,z); diff(LZ,x) diff(LZ,y) diff(LZ,z)]; % calcul of hessian of Lagrangian symb
display(LX); display(LY); display(LZ); display(H);

%% Main algorithm
i = 1; eps = 1;
while (eps >= precision)
    
    StepValue(:,i) = X;
    display(StepValue);
    x = X(1); y = X(2); z = X(3); lamb1 = X(4); lamb2 = X(5);
    grad1 = eval(HH); % eval of grad of constraint 1
    grad2 = eval(HH2); % eval of grad of constraint 2
    Hess = eval(H); % eval of hessian of Lagrangian
  
   % display(grad); display(Hess); display(E);
    
    M = [Hess grad1' grad2'; grad1 0 0; grad2 0 0]; % Matrix M
    display(M);
    V = (eval([-LX -LY -LZ -LLamb1 -LLamb2]))'; % Matrix 
 
    Direction = M\V;
    
    Y = X + Direction';
    eps = abs(Y - X);
    X = X + Direction';
    
    i = i+1;
end;

%% Steps and final optimal function value displaying
% display(StepValue);
display(X);
x = X(1); y = X(2); z = X(3);
min = eval(f);
display(min);

%% Griding and ploting the path
% [X,Y] = meshgrid(-5:0.01:5, -5:0.01:5);
% contour(X,Y,Z,((X.^2)+Y-11).^2+(X+Y.^2-7).^2, [0,0.01,0.03, 0.05,0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 1,2,3, 4, 4.5, 5, 7, 8, 9, 10, 15, 20, 25, 30, 40, 70, 100, 120, 135, 150, 200, 250]); %Markovitz
% % contour(X,Y,10*((X.^2)-Y).^2+(X-1).^2, [0,0.01,0.03, 0.05,0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 1,2,3, 4, 4.5, 5, 7, 8, 9, 10, 15, 20, 25, 30, 40, 70, 100, 120, 135, 150, 200, 250]); %rosenbrock
% % contour(X,Y,((X.^2)+Y-11).^2+(X+Y.^2-7).^2, [0,0.01,0.03, 0.05,0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 1,2,3, 4, 4.5, 5, 7, 8, 9, 10, 15, 20, 25, 30, 40, 70, 100, 120, 135, 150, 200, 250]); %Himmelblow
% 
% hold on;
% plot(StepValue(1,:), StepValue(2,:), '-o');
% a = -5:0.01:5;
% plot(a,1-a);
%     
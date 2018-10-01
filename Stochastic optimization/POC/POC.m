clear; clc;

%% Initial data
K_0 = 0.005; K_l = 200;
T_0 = 0; T_l = 0.5;
discretization_num_K = 200;
discretization_num_T = 30;
S = 100;
r = 0;
% Representation of T vector
T(1) = T_0;
dt = (T_l-T_0)/(discretization_num_T);
for i = 2:discretization_num_T+1
    T(i) = T_0+(i-1)*dt;
end;


%% Simulation of market data
sigma = ones(discretization_num_T+1,discretization_num_K+1)*0.2;
prices = BSPriceMatrixCreator(S,K_0,K_l,T_0,T_l,discretization_num_K, discretization_num_T,r,sigma);

%% Nonuniform discretization grid(the same as in the pricer)
delta_log_K = (log(K_l) - log(K_0)) / (discretization_num_K+3);
y = zeros(1,discretization_num_K+2);

y((discretization_num_K+2)/2) = log(S);
for i = 1:discretization_num_K+2
	y(i) = 0.5*log((log(K_l) + log(K_0)+i*delta_log_K) / (log(K_l) - (log(K_0) + i*delta_log_K))) + y((discretization_num_K+2)/2); 
%     if (y(i)>log(K_l))
%          break;
%     end;
end;

disc_T = 5; disc_K = 10;
ptsToEvalK = exp(y(2:end)) % size = discretization_num_K-2
%ptsToEvalK = exp(y(2:i-1)) % size = discretization_num_K-2
scatter(y(1:end), zeros(1,discretization_num_K+2));
%scatter(y(1:i-1), zeros(1,i-1));

%% Non-uniform2
% delta_log_K = 2*log(K_l) / (discretization_num_K);
% y_0 = log(S);
% 
% for i = 1:discretization_num_K+1
% 	y(i) = 0.5*log((log(K_l) - log(K_l)+i*delta_log_K) / (log(K_l) + (log(K_l) - i*delta_log_K))) + y_0; 
%     if (y(i)>log(K_l))
%         break;
%     end;
% end;
% 
% disc_T = 5; disc_K = 10;
% ptsToEvalK = exp(y(1:i)) % size = discretization_num_K-2
% scatter(y(1:i), zeros(1,i));


%% Creating the sigma surface
%while (dif < 5)

ctrlpts = rand(disc_T+1,disc_K+1);
sigmaaa = SplineLinear2DInterpol(T_0,T_l,K_0,K_l,disc_T, disc_K,ptsToEvalK,T, ctrlpts)
%Spline evaluated at the points at which our non-uniform discretization is
% done. How to give it the size 31*198 as well?


%% Pricing
%[u,K,T] = Pricer(sigma, K_0, K_l, T_0, T_l, discretization_num_K, discretization_num_T, S, r);
[u,K,T] = Pricer(sigmaaa, K_0, K_l, T_0, T_l, discretization_num_K, discretization_num_T, S, r);

fitness = sumOfSqrDif(u(1:end,96:104),prices(1:end,96:104))

%end;

abso = zeros(discretization_num_T+1,discretization_num_K-90);

%% fitting test
for i = 5:discretization_num_T+1
     for j = 1:discretization_num_K-90
         abso(i-4,j) = abs(u(i,j)-prices(i,j))/prices(i,j);
    end;
end;
% 
% %u
% %prices
figure;
surf(K,T,u);
xlabel('K'); ylabel('T'); zlabel('u(K,T)');
% 
size(abso)
figure;
surf(K(1:discretization_num_K-90),T(1:end),abso);
% %abso
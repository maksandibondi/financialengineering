clear; clc;

K_0 = 0.005; K_l = 200;
T_0 = 0; T_l = 0.5;
discretization_num_K = 200;
discretization_num_T = 30;
delta_k = (K_l-K_0)/discretization_num_K;
K(1) = K_0;
for j = 2:discretization_num_K+1
    K(j) = K(j-1) + delta_k;
end;
delta_t = (T_l-T_0)/discretization_num_T;
T(1) = T_0;
for i = 2:discretization_num_T+1
    T(i) = T(i-1) + delta_t;
end;



S = 100;
r = 0;
sigma = ones(discretization_num_T+1,discretization_num_K+1)*0.2;
prices = BSPriceMatrixCreator(S,K,r,sigma,T);



%% Setting initial and final conditions
for j = 2:discretization_num_K+1
    u(1,j) = max(S - (K(1) + (j-1)*delta_k), 0);
end;
for i = 1:discretization_num_T+1
    u(i,1) = S;
end;


[alpha beta gamma] = FillAlphaBetaGamma(sigma,r,K,discretization_num_K+1);


for n = 1:discretization_num_T
    A = TriagonalMatrix(alpha, beta, gamma, n + 1);
    H = A*delta_t+eye(discretization_num_K-1,discretization_num_K-1);
for i = 1:discretization_num_K-1
    B(i) = u(n,i+1)-alpha(n+1,2)*S*delta_t*max(2-i,0);
end;

res = ThomasAlgo(H,B);

for i = 2:discretization_num_K
    u(n+1,i) = res(i-1);
end;

end;



%% fitting test
for i = 5:discretization_num_T+1
    for j = 1:discretization_num_K-90
        abso(i-4,j) = abs(u(i,j)-prices(i,j))/prices(i,j);
    end;
end;

%u
%prices
figure;
surf(K,T,u); 

figure;
surf(K(1:discretization_num_K-90),T(5:end),abso);
xlabel('K'); ylabel('T'); zlabel('(u*(K,T)-market)/market');
%abso
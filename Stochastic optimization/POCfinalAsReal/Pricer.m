function[u, K, T] = Pricer(sigma, K_0, K_l, T_0, T_l, discretization_num_K, discretization_num_T, S, r)

%% Discretization setting
delta_k = (K_l-K_0)/discretization_num_K;
K(1) = K_0;
for j = 2:discretization_num_K+1
    K(j) = K_0 + (j-1)*delta_k;
end;

delta_t = (T_l-T_0)/discretization_num_T;
T(1) = T_0;
for i = 2:discretization_num_T+1
    T(i) = T_0 + (i-1)*delta_t;
end;

%% Setting initial and final conditions
for j = 2:discretization_num_K
    u(1,j) = max(S - (K(1) + (j-1)*delta_k), 0);
end;

for i = 1:discretization_num_T
    u(i,1) = S;
end;

%% Filling the matrix of parameters alpha, beta, gamma
[alpha beta gamma] = FillAlphaBetaGamma(sigma,r,K,discretization_num_K,S); 

%% Pricing
for n = 1:discretization_num_T-1
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



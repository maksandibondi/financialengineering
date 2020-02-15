function[u] = Pricer_dupire(sigma,ptsToEvalK,T, discretization_num_K, discretization_num_T, S, r, discretizationType, nonuniform_method)

%% Discretization setting
delta_k = (ptsToEvalK(end)-ptsToEvalK(1))/discretization_num_K;
delta_t = (T(end)-T(1))/discretization_num_T;

%% Setting initial and final conditions
for j = 2:discretization_num_K %+1
    u(1,j) = max(S - (ptsToEvalK(j)), 0);
end;

for i = 1:discretization_num_T %+1
    u(i,1) = S-ptsToEvalK(1);
end;

%% Filling the matrix of parameters alpha, beta, gamma
[alpha, beta, gamma] = FillAlphaBetaGamma(sigma,r,ptsToEvalK,discretization_num_K,S, discretizationType, nonuniform_method); 

%% Pricing
for n = 1:discretization_num_T-1
    A = TriagonalMatrix(alpha, beta, gamma, n); %+1
    H = A*delta_t+eye(discretization_num_K-2,discretization_num_K-2);
    for i = 1:discretization_num_K-1
            B(i) = u(n,i+1)-alpha(n+1,2)*S*delta_t*max(2-i,0);
    end;

    res = ThomasAlgo(H,B);

    for i = 2:discretization_num_K-1 % no -1
        u(n+1,i) = res(i-1);
    end;

end;



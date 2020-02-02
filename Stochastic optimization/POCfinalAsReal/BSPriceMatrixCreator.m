function [prices] = BSPriceMatrixCreator(S,K_0,K_l,T_0,T_l,discretization_num_K, discretization_num_T,r,sigma, popsize)

%% Discretization setting
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

sz1 = popsize;
sz2 = size(T,2);
sz3 = size(K,2);

%% Filling the matrix

for l = 1:sz1
    for i = 1:sz2
        for j = 1:sz3
        d1 = (1 / (sigma(l,i,j) * sqrt(T(i))))*(log(S / K(j)) + (r + (sigma(l,i,j)^ 2) / 2)*T(i));
        d2 = d1 - sigma(l,i,j) * sqrt(T(i));
        price = normcdf(d1)*S - normcdf(d2)*K(j)*exp(-r*T(i));
        prices(l,i,j) = price;
        end;
    end;
end;
 
return;


end


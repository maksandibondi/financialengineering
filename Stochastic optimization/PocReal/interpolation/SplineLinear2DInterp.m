function [sigma] = SplineLinear2DInterp(T_0,T_l,K_0,K_l, S, disc_T,disc_K,ptsToEvalK,ptsToEvalT, ctrlpts, discretizationType)

%delta_k = (K_l-K_0)/disc_K;
delta_t = (T_l-T_0)/(disc_T-1);
%knots = K_0:delta_k:K_l;
ptsToEvalT = transp(ptsToEvalT);

if (strcmp(discretizationType, 'uniform'))
    %knotsK = getUniformKnots(K_0, K_l, disc_K);
    knotsK = getNonUniformKnots(K_0, K_l, S, disc_K);
else
    knotsK = getNonUniformKnots(0, K_l, S, disc_K); % [0 11 16 18 24 26 30 42 70]
end
%knotsT = getUniformKnots(T_0, T_l, S, disc_T); 
%% B-spline part (K-discretization)
% We evaluate sigma(:,K) in 61 points using spline interpolation and 21 ctrl point by K. have
% to fill zeros with values as well. we obtain sigma(6*62) matrix

if (strcmp(discretizationType, 'uniform'))
    for i = 1:disc_T 
        order = 4;
        sigma_temp(i,:) = Bspline(knotsK,ctrlpts(i,:), order, ptsToEvalK);
        %sigma_temp(i,:) = interp1(knotsK, ctrlpts(i,:), ptsToEvalK, 'linear');
    end;
else
     for i = 1:disc_T
        %sigma(i,:) = interp1(knots,ctrlpts(i,:), ptsToEvalK, 'spline'); 
        order = 4;
        sigma_temp(i,:) = Bspline(knotsK,ctrlpts(i,:), order, ptsToEvalK);
        %sigma_temp(i,:) = interp1(knotsK, ctrlpts(i,:), ptsToEvalK, 'linear', 'extrap');
    end;
end
% figure;
% plot(ptsToEvalK, sigma_temp(i,:));
% hold on;
% plot(knotsK,ctrlpts(i,:));

for j = 1:size(sigma_temp,2)
     knotsT = T_0:delta_t:T_l;
     sigma(:,j) = interp1(knotsT, sigma_temp(:,j), ptsToEvalT);
end;

%figure;
%surf(knots,0:delta_t:0.5,ctrlpts); figure; 
%surf(ptsToEvalK(1:end-20),0:delta_t:0.5,sigma(:,(1:end-20)));

%% Convinient representation ( values far from strike go to 0)
%surf(ptsToEvalK(1:end-50),ptsToEvalT,out(:,(1:end-50)));

end


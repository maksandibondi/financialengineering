function [sigma] = SplineLinear2DInterp(T_0,T_l,K_0,K_l, S, disc_T,disc_K,ptsToEvalK,ptsToEvalT, ctrlpts, discretizationType)

%% Prepare grids
delta_t = (T_l-T_0)/(disc_T-1);
ptsToEvalT = transp(ptsToEvalT);

%% Getting uniform/nonuniform knots
if (strcmp(discretizationType, 'uniform'))
    %knotsK = getUniformKnots(K_0, K_l, disc_K);
    knotsK = getNonUniformKnots(K_0, K_l, S, disc_K);
else
    %knotsK = getUniformKnots(K_0, K_l, disc_K);
    knotsK = getNonUniformKnots(0, K_l, S, disc_K); % [0 11 16 18 24 26 30 42 70]
end

knotsT = T_0:delta_t:T_l;

%% 2D interpolation
% [Xq,Yq] = meshgrid(ptsToEvalK,ptsToEvalT);
% sigma(:,:) = interp2(knotsK, knotsT, ctrlpts, Xq,Yq,  'spline');
% 
% figure;
% plot(ptsToEvalK(1:40), sigma(1,1:40));
% hold on;
% plot(knotsK(1:10), ctrlpts(1,1:10), 'red');
% return;

%% Space discretization
if (strcmp(discretizationType, 'uniform'))
    for i = 1:disc_T 
        order = 4;
        sigma_temp(i,:) = Bspline(knotsK,ctrlpts(i,:), order, ptsToEvalK);
        %sigma_temp(i,:) = interp1(knotsK, ctrlpts(i,:), ptsToEvalK, 'linear');
        %sigma_temp(i,:) = interp1(knotsK, ctrlpts(i,:), ptsToEvalK,'pchip'); % works for uniform knots
    end;
else
     for i = 1:disc_T
        %sigma_temp(i,:) = interp1(knots,ctrlpts(i,:), ptsToEvalK, 'spline'); % works for nonuniform knots
        %sigma_temp(i,:) = interp1(knotsK, ctrlpts(i,:), ptsToEvalK, 'linear', 'extrap');
        order = 4;
        sigma_temp(i,:) = Bspline(knotsK,ctrlpts(i,:), order, ptsToEvalK);
        
    end;
end

%% Time discretization
for j = 1:size(sigma_temp,2)
     sigma(:,j) = interp1(knotsT, sigma_temp(:,j), ptsToEvalT);
end;

end


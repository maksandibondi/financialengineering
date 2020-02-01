function [sigma] = SplineLinear2DInterp(T_0,T_l,K_0,K_l, S, disc_T,disc_K,ptsToEvalK,ptsToEvalT, ctrlpts, discretizationType)

%delta_k = (K_l-K_0)/disc_K;
delta_t = (T_l-T_0)/(disc_T-1);
%knots = K_0:delta_k:K_l;
%order = 4;
ptsToEvalT = transp(ptsToEvalT);

if (strcmp(discretizationType, 'uniform'))
    knotsK = getNonUniformKnots(K_0, K_l, S, disc_K); 
else
    knotsK = getNonUniformKnots(0, K_l, S, disc_K);
end
%knotsT = getUniformKnots(T_0, T_l, S, disc_T); 
%% B-spline part (K-discretization)
% We evaluate sigma(:,K) in 61 points using spline interpolation and 21 ctrl point by K. have
% to fill zeros with values as well. we obtain sigma(6*62) matrix

if (strcmp(discretizationType, 'uniform'))
    for i = 1:disc_T
        %sigma(i,:) = interp1(knots,ctrlpts(i,:), ptsToEvalK, 'spline'); 
        sigma_temp(i,:) = interp1(knotsK, ctrlpts(i,:), ptsToEvalK);
        %sigma(i,:) = Bspline(knots,ctrlpts(i,:), order, ptsToEvalK);

    end;
else
     for i = 1:disc_T
        %sigma(i,:) = interp1(knots,ctrlpts(i,:), ptsToEvalK, 'spline'); 
        sigma_temp(i,:) = interp1(knotsK, ctrlpts(i,:), ptsToEvalK, 'linear', 'extrap');
        %sigma_temp(i,end) =  sigma_temp(i,end-1);
        %sigma(i,:) = Bspline(knots,ctrlpts(i,:), order, ptsToEvalK);

    end;
end
%plot(ptsToEvalK, sigma(i,:));

for j = 1:size(sigma_temp,2)
     knotsT = T_0:delta_t:T_l;
     sigma(:,j) = interp1(knotsT, sigma_temp(:,j), ptsToEvalT);
end;

%figure;
%surf(knots,0:delta_t:0.5,ctrlpts); figure; 
%surf(ptsToEvalK(1:end-20),0:delta_t:0.5,sigma(:,(1:end-20)));

%% Linear discretization by time axis
% We need to evaluate sigma(10,:) in 10 points using linear interpolation in and
% 6 ctrl points by T
% for j=1:size(ptsToEvalT,2)
%     T_inf_idx = 1;
%     T_sup_idx = 1;
%     for m=1:disc_T
%         if (ptsToEvalT(j)>=T_0+(m-1)*delta_t)
%             T_inf_idx = m;
%             continue;
%         else
%             T_sup_idx = m;
%             break;
%         end;
%     end;
%     T_sup = T_0+T_sup_idx*delta_t;
%     T_inf = T_0+T_inf_idx*delta_t;
%     
%     
%     for k = 1:size(ptsToEvalK,2)
%         out(j,k) = (T_sup-ptsToEvalT(j))*sigma(T_inf_idx,k)/(T_sup-T_inf)+(ptsToEvalT(j)-T_inf)*sigma(T_sup_idx,k)/(T_sup-T_inf);
%     end;
%    
% end;


%% Convinient representation ( values far from strike go to 0)
%surf(ptsToEvalK(1:end-50),ptsToEvalT,out(:,(1:end-50)));

end


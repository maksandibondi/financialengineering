function [sigma] = SplineLinear2DInterp(T_0,T_l,K_0,K_l, S, disc_T,disc_K,ptsToEvalK,ptsToEvalT, ctrlpts, discretizationType, interpTypeK)

%% Prepare grids
delta_t = (T_l-T_0)/(disc_T-1);
ptsToEvalT = transp(ptsToEvalT);

%% Getting uniform/nonuniform knots
if (strcmp(discretizationType, 'uniform'))
    knotsK = getUniformKnots(K_0, K_l, disc_K);
else
    knotsK = getNonUniformKnots(0, K_l, S, disc_K); % [0 11 16 18 24 26 30 42 70]
end

knotsT = T_0:delta_t:T_l;

%% Space discretization
if (strcmp(discretizationType, 'uniform'))
    for i = 1:disc_T 
        if (strcmp(interpTypeK, 'bspline'))
            order = 4;
            sigma_temp(i,:) = Bspline(knotsK,ctrlpts(i,:), order, ptsToEvalK); 
            %sigma_temp(i,:) = interp1(knotsK, ctrlpts(i,:), ptsToEvalK, 'pchip');
            %sigma_temp(i,:) = BsplineMama(knotsK,ctrlpts(i,:), ptsToEvalK);
%         if (sigma_temp(i,20) ~= 0)
%             plot([1:46],sigma_temp(i,:));
%         end
%             x = knotsK; v = ctrlpts(i,:); xq=ptsToEvalK;
%             vq2 = interp1(x,v,xq,'spline');
%             hold on
%             plot(x,v,'o',xq,vq2,':.');
%             %xlim([0 2*pi]);
%             title('Spline Interpolation');
%         end

        else
            sigma_temp(i,:) = interp1(knotsK, ctrlpts(i,:), ptsToEvalK, interpTypeK);
        end;
    end;
else
     for i = 1:disc_T
         if (strcmp(interpTypeK, 'bspline'))
            order = 4;
            sigma_temp(i,:) = Bspline(knotsK,ctrlpts(i,:), order, ptsToEvalK);
         else
            sigma_temp(i,:) = interp1(knotsK,ctrlpts(i,:), ptsToEvalK, interpTypeK);
         end;
    end;
end

%% Time discretization
for j = 1:size(sigma_temp,2)
      sigma(:,j) = interp1(knotsT, sigma_temp(:,j), ptsToEvalT, 'spline');
end;


end


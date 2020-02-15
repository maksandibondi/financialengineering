function [ctrl] = getSeedFromImplicitVol(K_0, K_l, T_0, T_l, discretization_num_T, discretization_num_K,disc_T, disc_K,VolImp, knotsK, knotsT)

   ptsToEvalCtrlptsK = getUniformKnots(K_0, K_l, disc_K);
   ptsToEvalCtrlptsT = getUniformKnots(T_0, T_l, disc_T);
   
for i = 1:discretization_num_T 
    order = 4;
    %ctrl_temp(i,:) = Bspline(knotsK,VolImp(i,:), order, ptsToEvalCtrlptsK);
    ctrl_temp(i,:) = interp1(knotsK,VolImp(i,:), ptsToEvalCtrlptsK, 'spline');
end;

for j = 1:size(ctrl_temp,2) % discretization_num_K
      ctrl(:,j) = interp1(knotsT, ctrl_temp(:,j), ptsToEvalCtrlptsT, 'spline');
end;

end
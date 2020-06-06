function [ localVolCalibratedNormalizedScale ] = getInterpolatedLocalVol( localVolCalibrated, T, K, T_normalized, K_normalized, interpMethod )
%GETINTERPOLATEDLOCALVOL Summary of this function goes here
%   Detailed explanation goes here

    disc_T = size(T);

    for i = 1:disc_T 
        localVolCalibratedNormalizedScale_temp(i,:) = interp1(K, localVolCalibrated(i,:), K_normalized, interpMethod, 'extrap');
    end;
    
    for j = 1:size(localVolCalibratedNormalizedScale_temp,2) 
        localVolCalibratedNormalizedScale(:,j) = interp1(T, localVolCalibrated(:,j), T_normalized, interpMethod, 'extrap');
    end;
    
end


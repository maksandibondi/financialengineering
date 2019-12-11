function [ VolImp_filled ] = interpolateImpliedVolVector( VolImp)
%INTERPOLATEIMPLIEDVOL Summary of this function goes here
%   Detailed explanation goes here
 % Interpolate intermediate values
    vectorWithZeros = VolImp;
    x = 1 : length(vectorWithZeros);
    m = vectorWithZeros==0;
    s = interp1(x(~m),vectorWithZeros(~m),x(m));
    %s = spline();
    VolImp(m) = s;
    
    % Interpolate tail values
    vectorWithNans = VolImp;
    x = 1 : length(vectorWithNans);
    % Works only for call option
    if (isnan(vectorWithNans(1)))
        bestGuess = 0.01;
        vectorWithNans(1) = bestGuess;
        VolImp(1) = bestGuess;
    end
    
    if (isnan(vectorWithNans(end)))
        bestGuess = 0.01;
        vectorWithNans(end) = bestGuess;
        VolImp(end) = bestGuess;
    end
  
    m = isnan(vectorWithNans); 
    s = interp1(x(~m),vectorWithNans(~m),x(m));
    VolImp(m) = s;
    
    VolImp_filled = VolImp;


end


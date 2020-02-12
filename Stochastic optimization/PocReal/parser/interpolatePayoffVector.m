function [ Vmarket_filled ] = interpolatePayoffVector( Vmarket, K, S, optionType )
%INTERPOLATEPAYOFFVECTOR Summary of this function goes here
%   Detailed explanation goes here
 
% Interpolate intermediate values
    vectorWithZeros = Vmarket;
    x = 1 : length(vectorWithZeros);
    m = vectorWithZeros==0;
    s = interp1(x(~m),vectorWithZeros(~m),x(m));
    %s = spline();
    Vmarket(m) = s;
    
    % Interpolate tail values
    vectorWithNans = Vmarket;
    x = 1 : length(vectorWithNans);
    % Works only for call option
    if (isnan(vectorWithNans(1)))
        if (optionType == 'C' && K(1) < S)
            bestGuess = S - K(1);
        elseif (optionType == 'C' && K(1) > S)
             bestGuess = 0.01;
        elseif (optionType == 'P' && K(1) < S)
             bestGuess = 0.01;
        elseif (optionType == 'P' && K(1) > S)
             bestGuess = K(1) - S;
        end
        
        vectorWithNans(1) = bestGuess;
        Vmarket(1) = bestGuess;
    end
    
    if (isnan(vectorWithNans(end)))
       if (optionType == 'C' && K(end) < S)
            bestGuess = S - K(end);
        elseif (optionType == 'C' && K(end) > S)
             bestGuess = 0.01;
        elseif (optionType == 'P' && K(end) < S)
             bestGuess = 0.01;
        elseif (optionType == 'P' && K(end) > S)
             bestGuess = K(end) - S;
        end
        vectorWithNans(end) = bestGuess;
        Vmarket(end) = bestGuess;
    end
  
    m = isnan(vectorWithNans); 
    s = interp1(x(~m),vectorWithNans(~m),x(m));
    Vmarket(m) = s;
    
    Vmarket_filled = Vmarket;
end


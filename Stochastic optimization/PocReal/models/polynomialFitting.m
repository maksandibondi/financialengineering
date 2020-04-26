function [ paramReg, rsquared ] = polynomialFitting( x, y, order )
%POLYNOMIALFITTING Summary of this function goes here
%   Detailed explanation goes here

%% Regression
paramReg = polyfit(x,y,order);

%% CALCULATE RESIDUALS AND R2
    y_model = polyval(paramReg,x);
    resid = y - y_model;
    SSresid = sum(resid.^2);
    SStotal = (length(y)-1) * var(y);
    rsquared = 1 - (SSresid/SStotal) * (length(y)-1)/(length(y)-length(paramReg)-1); %% coef of determination

end


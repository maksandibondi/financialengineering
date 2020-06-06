function [ impliedVolCalibrated, impliedVolCalibratedNormalizedScale, S ] = impliedVolCalibratorStockRunner( inputStructure )
%IMPLIEDVOLCALIBRATORSTOCK Summary of this function goes here
%   Detailed explanation goes here
%% Improt folders
addpath('../../parser','-end');
addpath('../../pricing','-end');
addpath('../../pricing/BS','-end');
addpath('../../interpolation', '-end');
addpath('../../utils', '-end');
addpath('../../formatting', '-end');
addpath('../../graphics', '-end');
addpath('..', '-end');
addpath('../Algo', '-end');

inputStruct = inputStructure; %% copy input structure to keep initial setting constant

%% parse market data and format it to prepare for usage here
[K, T, dates, S, Vmarket] = parseStockCSV(inputStruct.marketDataFile, inputStruct.ticker, inputStruct.date, inputStruct.optionType);

impliedVolCalibrated  = calculateImpVolMatrixNewton(T, K, S, inputStruct.r, Vmarket, inputStruct.epsilon_impVol);


if inputStruct.isNormalizedScale
    interpMethod = 'linear';
    impliedVolCalibratedNormalizedScale = getInterpolatedLocalVol(impliedVolCalibrated, T, K, inputStruct.T_normalized, inputStruct.K_normalized, interpMethod);
end;

return;


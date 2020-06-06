function [S, localVolCalibrated, localVolCalibratedNormalizedScale, diffprice, epsilon, weight] = localVolCalibratorStockRunner(inputStructure)
% This module takes data from files for different option historical data
% it builds linear regression on this data

%% Improt folders
addpath('../../parser','-end');
addpath('../../pricing','-end');
addpath('../../interpolation', '-end');
addpath('../../utils', '-end');
addpath('../../formatting', '-end');
addpath('../../graphics', '-end');
addpath('..', '-end');
addpath('../Algo', '-end');

inputStruct = inputStructure; %% copy input structure to keep initial setting constant

%% parse market data and format it to prepare for usage here
[K, T, dates, S, Vmarket] = parseStockCSV(inputStruct.marketDataFile, inputStruct.ticker, inputStruct.date, inputStruct.optionType);
VolImp = ones(size(Vmarket,1), size(Vmarket,2))*0.2;

failure = 1;
iter = 1;
while (failure) %% rerun until success
    if (iter > 1 && iter < 7)
        inputStruct = adaptModelParameters(inputStruct, 'epsilon'); %% increase epsilon
    elseif (iter > 1)
        inputStruct = adaptModelParameters(inputStruct, 'weight'); %% increase weight
    end;
    
    [localVolCalibrated, diffprice, ptsToEvalK, f] = calculateLocalVolGenetic(K, T, S, inputStruct.r, Vmarket, VolImp, inputStruct);
    failure = f; %% true if failed
    iter = iter + 1;
end;

if inputStruct.isNormalizedScale
    interpMethod = 'linear';
    localVolCalibratedNormalizedScale = getInterpolatedLocalVol(localVolCalibrated, T, ptsToEvalK, inputStruct.T_normalized, inputStruct.K_normalized, interpMethod);
end;

epsilon = inputStruct.epsilon;
weight = inputStruct.concentration_weights;

return;
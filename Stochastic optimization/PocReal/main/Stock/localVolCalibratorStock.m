function [S, localVolCalibrated, localVolCalibratedNormalizedScale, diffprice] = localVolCalibratorStock(inputStructure)
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

%% parse market data and format it to prepare for usage here
[K, T, dates, S, Vmarket] = parseStockCSV(inputStructure.marketDataFile, inputStructure.ticker, inputStructure.date, inputStructure.optionType);
VolImp = ones(size(Vmarket,1), size(Vmarket,2))*0.2;
[localVolCalibrated, diffprice, ptsToEvalK] = geneticRunner(K, T, S, inputStructure.r, Vmarket, VolImp, inputStructure);

if inputStructure.isNormalizedScale
    interpMethod = 'linear';
    localVolCalibratedNormalizedScale = getInterpolatedLocalVol(localVolCalibrated, T, ptsToEvalK, inputStructure.T_normalized, inputStructure.K_normalized, interpMethod);
end;

return;
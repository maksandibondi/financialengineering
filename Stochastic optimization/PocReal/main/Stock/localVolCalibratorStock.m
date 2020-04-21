function [localVolCalibrated] = localVolCalibratorStock(inputStructure)
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
%[T, Vmarket, VolImp] = additionalFormattingStock( T,  Vmarket, VolImp, S, K, inputStructure.optionType);
VolImp = ones(size(Vmarket,1), size(Vmarket,2))*0.2;
localVolCalibrated = geneticRunner(K, T, S, inputStructure.r, Vmarket, VolImp, inputStructure);
return;
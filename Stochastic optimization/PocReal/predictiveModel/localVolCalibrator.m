function [localVolCalibrated] = localVolCalibrator(inputStructure)
% This module takes data from files for different option historical data
% it builds linear regression on this data

%% Improt folders
addpath('../parser','-end');
addpath('../pricing','-end');
addpath('../interpolation', '-end');
addpath('../utils', '-end');
addpath('../formatting', '-end');
addpath('../graphics', '-end');

%% parse market data and format it to prepare for usage here
[K, T, S, Vmarket, VolImp] = parseCBOE(inputStructure.marketDataFile, inputStructure.ticker, inputStructure.date, inputStructure.optionType);
[T, Vmarket, VolImp] = additionalFormatting( T,  Vmarket, VolImp, S, K, inputStructure.optionType);
localVolCalibrated = geneticRunner(K, T, S, inputStructure.r, Vmarket, VolImp, inputStructure);
return;
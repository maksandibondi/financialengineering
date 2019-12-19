% This module takes data from files for different option historical data
% it builds linear regression on this data


% Input. Parameters set

% Intermediate res: calculated local vol( K,T ) for every T0 on input 
% Final res: Regression on local vol(K,T) near the money 


% Verifications: does local vol "explain" VIX
% Verifications: does local vol "explain" std deviation
% Verifications: does local vol predicts well new local vol

addpath('../parser','-end');
addpath('../pricing','-end');
addpath('../interpolation', '-end');
addpath('../utils', '-end');


clear;
[K, T, S, Vmarket, VolImp] = parseCBOE('../../../resources/UnderlyingOptionsEODCalcs_2018-11.csv', '^VIX', '2018-11-15', 'C');
geneticRunner(K, T, S, 0, Vmarket, VolImp);
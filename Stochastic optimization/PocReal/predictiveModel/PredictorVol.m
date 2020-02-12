clear;
% Intermediate res: calculated local vol( K,T ) for every T0 on input 
% Final res: Regression on local vol(K,T) near the money 


% Verifications: does local vol "explain" VIX
% Verifications: does local vol "explain" std deviation
% Verifications: does local vol predicts well new local vol

%% Algorithm settings
inputStructure = struct;
inputStructure.T0 = 1000;
inputStructure.t = 0.1;
inputStructure.Nmc = 2^14;
inputStructure.M = 64;
inputStructure.popsize = 20;
inputStructure.disc_T = 6; % reduced dimensionality for T axe
inputStructure.disc_K = 13; % reduced dimensionality for K axe

inputStructure.epsilon = 0.13;
inputStructure.concentration_weights = 0.06;
inputStructure.discretizationType = 'uniform';
inputStructure.interpTypeK = 'bspline';

inputStructure.marketDataFile = '../../resources/UnderlyingOptionsEODCalcs_2018-12.csv';
% '../../resources/UnderlyingOptionsEODCalcs_2018-11.csv'
% '../../resources/UnderlyingOptionsEODCalcs_2018-12.csv'
inputStructure.ticker = '^VIX';
inputStructure.date = '2018-12-24'; %'2018-11-15', '2018-11-30', '2018-12-24'
inputStructure.optionType = 'C';
inputStructure.r = 0;

localVolCalibrator(inputStructure);
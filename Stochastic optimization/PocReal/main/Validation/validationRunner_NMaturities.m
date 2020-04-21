clear;
% Intermediate res: calculated local vol( K,T ) for every T0 on input 
% Final res: Regression on local vol(K,T) near the money 


% Verifications: does local vol "explain" VIX
% Verifications: does local vol "explain" std deviation
% Verifications: does local vol predicts well new local vol

%% Algorithm settings
inputStructure = struct;
inputStructure.epsilon = 0.2; 
inputStructure.discretizationType = 'uniform'; % nonuniform
inputStructure.nonuniform_method = 'user'; % 'log', 'sin'

inputStructure.validationFile = '../../Results/Validation/Validation_2018-11-15.xls';
inputStructure.marketDataFile = '../../../resources/UnderlyingOptionsEODCalcs_2018-11.csv';
% '../../resources/UnderlyingOptionsEODCalcs_2018-11.csv'
% '../../resources/UnderlyingOptionsEODCalcs_2018-12.csv'
inputStructure.ticker = '^VIX';
inputStructure.validationDate = '2018-11-15'; %'2018-11-30', '2018-11-30', '2018-12-24'
inputStructure.optionType = 'C';
inputStructure.r = 0;
inputStructure.vasicek_assumption = 1; % zero

localVolValidator_NMaturities(inputStructure);
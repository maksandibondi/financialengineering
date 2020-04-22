clear;
% Intermediate res: calculated local vol( K,T ) for every T0 on input 
% Final res: Regression on local vol(K,T) near the money 


% Verifications: does local vol "explain" VIX
% Verifications: does local vol "explain" std deviation
% Verifications: does local vol predicts well new local vol

%% Algorithm settings
inputStructure = struct;
inputStructure.T0 = 1000;
inputStructure.t = 0.21; %% Important not to choose too big value so divergence happens, rerun sometimes for different set of random numbes
inputStructure.Nmc = 2^14;
inputStructure.M = 64;
inputStructure.popsize = 20;
inputStructure.disc_T = 6; % reduced dimensionality for T axe
inputStructure.disc_K = 13; % reduced dimensionality for K axe

inputStructure.epsilon = 0.2; % 0.18 for 10-12-2018 ; 1.3 for 2018-12-24  % Changing it allows to explore different solutions
inputStructure.maxIter = 200;
inputStructure.concentration_weights = 0.06;
inputStructure.discretizationType = 'uniform'; % nonuniform
inputStructure.nonuniform_method = 'user'; % 'log', 'sin'
inputStructure.interpTypeK = 'bspline';
%inputStructure.rowByRowMutation = 1;

inputStructure.visualizeResults = 1;
inputStructure.marketDataFile = '../../../resources/UnderlyingOptionsEODCalcs_2019-01.csv';
% '../../resources/UnderlyingOptionsEODCalcs_2018-11.csv'
% '../../resources/UnderlyingOptionsEODCalcs_2018-12.csv'
inputStructure.ticker = '^VIX';
inputStructure.date = '2019-01-15'; %'2018-11-15', '2018-11-30', '2018-12-24'
inputStructure.optionType = 'C';
inputStructure.r = 0;
inputStructure.includeVolAssumption = 'zero'; % zero
inputStructure.vasicek_assumption = 1; % zero
inputStructure.outputdir = fullfile(pwd,'../../Results');

localVolCalibrator(inputStructure);
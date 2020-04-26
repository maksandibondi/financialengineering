clear;
% Intermediate res: calculated local vol( K,T ) for every T0 on input 
% Final res: Regression on local vol(K,T) near the money 


% Verifications: does local vol "explain" VIX
% Verifications: does local vol "explain" std deviation
% Verifications: does local vol predicts well new local vol
addpath('../../models', '-end');


%% Algorithm settings
inputStructure = struct;
inputStructure.T0 = 1000;
inputStructure.t = 0.2; %% Important not to choose too big value so divergence happens, rerun sometimes for different set of random numbes
inputStructure.Nmc = 2^14;
inputStructure.M = 64;
inputStructure.popsize = 20;
inputStructure.disc_T = 6; % reduced dimensionality for T axe
inputStructure.disc_K = 13; % reduced dimensionality for K axe

inputStructure.epsilon = 0.3; % 0.18 for 10-12-2018 ; 1.3 for 2018-12-24  % Changing it allows to explore different solutions
inputStructure.maxIter = 20;
inputStructure.concentration_weights = 1;
inputStructure.discretizationType = 'uniform'; % nonuniform
inputStructure.nonuniform_method = 'user'; % 'log', 'sin'
inputStructure.interpTypeK = 'bspline';
%inputStructure.rowByRowMutation = 1;

inputStructure.marketDataFile = '../../../resources/MSFT.csv';
inputStructure.ticker = 'MSFT';
inputStructure.optionType = 'C';
inputStructure.r = 0;
inputStructure.includeVolAssumption = 'zero'; % zero
inputStructure.vasicek_assumption = 0; % zero
inputStructure.outputdir = fullfile(pwd,'../../Results/Results_MSFT');

inputStructure.visualizeResults = 0;
inputStructure.isNormalizedScale = 1;
inputStructure.T_normalized = [0.05, 0,1, 0.2, 0.5, 1, 1.5, 2];
inputStructure.K_normalized = 100:5:230;

calibrationDates = {'01/03/2020', '01/10/2020', '01/17/2020', '01/24/2020', '01/31/2020', '02/07/2020', '02/14/2020', '02/21/2020', '02/28/2020', '03/06/2020', '03/13/2020', '03/20/2020', '03/27/2020', '04/03/2020'};

%% Calibrate local vol with adaptive param for all dates
for i = 1:size(calibrationDates,2)
    inputStructure.date = cell2mat(calibrationDates(i)); 
    [S, localVolCalibrated, localVolCalibratedNormalizedScale, diffprice, epsilon, weight] = localVolCalibratorStock(inputStructure);
    S0(i) = S;
    eps(i) = epsilon;
    wei(i) = weight;
    %LV(i,:,:) = localVolCalibrated;
    LVN(i,:,:) = localVolCalibratedNormalizedScale;
end;

%% Construct ATM matrix for regression analysis
for k = 1:size(inputStructure.T_normalized, 2)
    for j = 1:size(calibrationDates,2)
        [~,closestIndex] = min(abs(S0(j)-inputStructure.K_normalized'));
        LVATM(k,j) = LVN(j,k,closestIndex); %% get local vol value atm dynamics for all T
    end;
end;



%% Do regression analysis
calibrationDatesNumeric(1) = 0;
for i = 2:size(calibrationDates,2) % get numeric values for calibration dates vector
    calibrationDatesNumeric(i) = daysact(cell2mat(calibrationDates(1)), cell2mat(calibrationDates(i)))/365;
end;

%% (1) Regress ATM Local Volatility on calibration dates for each maturity
for k = 1:size(inputStructure.T_normalized, 2)
    LVATM_ = LVATM(k,:);
    idxnan = isnan(LVATM_); %% find indices with non NAN values
    [prmReg, rsquared] = polynomialFitting(calibrationDatesNumeric(~idxnan),LVATM_(~idxnan),3);
    paramReg(k,:) = prmReg;
    rsq(k) = rsquared;
    LVAMT_model = polyval(paramReg(k,:),calibrationDatesNumeric(~idxnan)); %% values obtained by model

    %% visualize regression results
    figure;
    scatter(calibrationDatesNumeric(~idxnan),LVATM_(~idxnan));
    hold on
    plot(calibrationDatesNumeric(~idxnan),LVAMT_model)
end;
% Conclusion : sometimes negative coef of determination, bad fitting

%% (1) Validate model (1)

%% (2) Regress ATM Local Volatility on implied volatility for each maturity

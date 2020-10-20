% In this file we make a run and study the relation between ATM local vol
% and Realized volatility independent variables in a univariate regression,
% generating reports
clear;
% Intermediate res: calculated local vol( K,T ) for every T0 on input 
% Final res: Regression on local vol(K,T) near the money 


% Verifications: does local vol "explain" VIX
% Verifications: does local vol "explain" std deviation
% Verifications: does local vol predicts well new local vol
import mlreportgen.report.* 
import mlreportgen.dom.* 

addpath('../../models', '-end');
addpath('../../parser', '-end');

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
inputStructure.epsilon_impVol = 0.5;
inputStructure.maxIter = 20;
inputStructure.concentration_weights = 1;
inputStructure.discretizationType = 'uniform'; % nonuniform
inputStructure.nonuniform_method = 'user'; % 'log', 'sin'
inputStructure.interpTypeK = 'bspline';
%inputStructure.rowByRowMutation = 1;

inputStructure.marketDataFile = '../../../resources/MSFT.csv';
inputStructure.spotFile='../../../resources/MSFT_spot.csv';
inputStructure.ticker = 'MSFT';
inputStructure.optionType = 'C';
inputStructure.r = 0;
inputStructure.includeVolAssumption = 'zero'; % zero
inputStructure.vasicek_assumption = 0; % zero
inputStructure.outputdir = fullfile(pwd,'../../Results/Results_MSFT');

inputStructure.visualizeResults = 0;
inputStructure.isNormalizedScale = 1;
inputStructure.T_normalized = [0.05, 0.1, 0.15, 0.2, 0.25, 0.4];
inputStructure.K_normalized = 100:5:230;

calibrationDates = {'01/03/2020', '01/10/2020', '01/24/2020', '01/31/2020', '02/07/2020',  '02/28/2020', '03/06/2020',  '03/20/2020', '03/27/2020', '04/03/2020'};
validationDates = {'01/17/2020','02/14/2020', '02/21/2020', '03/13/2020', '04/09/2020', '04/17/2020'};

%% Prepare numeric calibration and validation dates
calibrationDatesNumeric(1) = 0;
for i = 2:size(calibrationDates,2) % get numeric values for calibration dates vector
    calibrationDatesNumeric(i) = daysact(cell2mat(calibrationDates(1)), cell2mat(calibrationDates(i)))/365;
end;

for i = 1:size(validationDates,2) % get numeric values for calibration dates vector
    validationDatesNumeric(i) = daysact(cell2mat(calibrationDates(1)), cell2mat(validationDates(i)))/365;
end;


%% (3) Calulate realized standard deviation at all maturities w.r.t all calibration dates
% get date format of normalized maturities for each calibration date 
T_normalized_dates = parseSpotCSVDateFromFractionOfYear('../../../resources/MSFT_spot.csv', calibrationDates, inputStructure.T_normalized);
realizedVol = parseSpotCSVRealizedVol('../../../resources/MSFT_spot.csv', calibrationDates, T_normalized_dates);

%% (3) Calibrate Local Volatility with adaptive param for all Calibration dates
for i = 1:size(calibrationDates,2)
    inputStructure.date = cell2mat(calibrationDates(i)); 
    [S, localVolCalibrated, localVolCalibratedNormalizedScale, diffprice, epsilon, weight] = localVolCalibratorStockRunner(inputStructure);
    S0(i) = S;
    eps(i) = epsilon;
    wei(i) = weight;
    %LV(i,:,:) = localVolCalibrated;
    LVN(i,:,:) = localVolCalibratedNormalizedScale;
end;
 % Construct ATM local vol matrix for regression analysis
for k = 1:size(inputStructure.T_normalized, 2)
    for j = 1:size(calibrationDates,2)
        [~,closestIndex] = min(abs(S0(j)-inputStructure.K_normalized'));
        LVATM(k,j) = LVN(j,k,closestIndex); %% get local vol value atm dynamics for all T
    end;
end;

%% (3) Calibrate implied vol with adaptive param for all Calibration dates
 % Calculate implied volatility for all maturities for all calib dates 
for i = 1:size(calibrationDates,2)
    inputStructure.date = cell2mat(calibrationDates(i)); 
    [impliedVolCalibrated, impliedVolCalibratedNormalizedScale, S] = impliedVolCalibratorStockRunner(inputStructure);
    S0(i) = S;
    VolImp(i,:,:) = impliedVolCalibratedNormalizedScale;
end;
 % Construct ATM implied vol matrix for regression analysis
for k = 1:size(inputStructure.T_normalized, 2)
    for j = 1:size(calibrationDates,2)
        [~,closestIndex] = min(abs(S0(j)-inputStructure.K_normalized'));
        VolImpATM(k,j) = VolImp(j,k,closestIndex); %% get imp vol value atm dynamics for all T
    end;
end;

%% (3) Visualize calibration results
for i = 1:size(calibrationDates,2)
    f1(i) = figure; f1(i).Visible = 'off'; f1(i).Tag = 'Local Vol Surface';
    surf(inputStructure.K_normalized, inputStructure.T_normalized, squeeze(LVN(i,:,:)));
    ttl = sprintf('Local volatility surface at t =%s i.e. t=%s', char(calibrationDates(i)), num2str(calibrationDatesNumeric(i)));
    title(ttl);
    f2(i) = figure; f2(i).Visible = 'off'; f2(i).Tag = 'Implied Vol Surface';
    surf(inputStructure.K_normalized, inputStructure.T_normalized, squeeze(VolImp(i,:,:)));
    ttl = sprintf('Implied volatility surface at t =%s i.e. t=%s', char(calibrationDates(i)), num2str(calibrationDatesNumeric(i)));
    title(ttl);    
end;

%% (3) Generate report
report('0.rpt','-oReportCalibRealizedVol.rtf','-frtf');
pause(7);



%% (3) Regress (univariate) ATM Local Volatility on realized T volatility of S for each maturity for all calib dates
for k = 1:size(inputStructure.T_normalized, 2)
    [LVATM_, sortedIndex] = sort(LVATM(k,:));
    RelizedVol_ = realizedVol(:,k)';
    RealizedVol_ = RelizedVol_(sortedIndex);
    idxnan = isnan(LVATM_); %% find indices with non NAN values
    [prmReg, rsquared] = polynomialFitting(LVATM_(~idxnan),RelizedVol_(~idxnan),1);
    [prmReg2, rsquared2] = polynomialFitting(LVATM_(~idxnan),RelizedVol_(~idxnan),3);
    paramReg(k,:) = prmReg;
    paramReg2(k,:) = prmReg2;
    rsq(k) = rsquared;
    rsq2(k) = rsquared2;
    RelizedVol_model = polyval(paramReg(k,:),LVATM_(~idxnan)); %% values obtained by model
    RelizedVol_model2 = polyval(paramReg2(k,:),LVATM_(~idxnan)); %% values obtained by model    

%% visualize regression results
    f(k) = figure; f(k).Visible = 'off'; f(k).Tag = 'ATMonRealized reg';
    s1 = scatter(LVATM_(~idxnan), RelizedVol_(~idxnan));
    hold on
    p1 = plot(LVATM_(~idxnan),RelizedVol_model); ylabel('realized vol'); xlabel('LVATM');
    p2 = plot(LVATM_(~idxnan),RelizedVol_model2); ylabel('realized vol'); xlabel('LVATM');
    ttl = sprintf('Linear regressions of Relized Vol(t) on ATM Local Vol(t) with fixed T=%s',num2str(inputStructure.T_normalized(k)));
    title(ttl);
    legend([s1;p1;p2], 'Realized vol real', 'Realized vol lin regr order1', 'Realized vol lin regr order3');
    legend('boxoff');
end;



%% (3) Validate model on validation points (interm and extrap) with MAE (mean abs error)
% calculate realized vol values on validation points
for i = 1:size(validationDates,2);
    T_normalized_dates = parseSpotCSVDateFromFractionOfYear('../../../resources/MSFT_spot.csv', validationDates, inputStructure.T_normalized);
    realizedVol_vld = parseSpotCSVRealizedVol('../../../resources/MSFT_spot.csv', validationDates, T_normalized_dates);
end;

% calculate local vol real values : Calibrate local vol with adaptive param for all validation dates
for i = 1:size(validationDates,2)
    inputStructure.date = cell2mat(validationDates(i)); 
    [S, localVolCalibrated, localVolCalibratedNormalizedScale, diffprice, epsilon, weight] = localVolCalibratorStockRunner(inputStructure);
    S0_v(i) = S;
    eps_v(i) = epsilon;
    wei_v(i) = weight;
    %LV(i,:,:) = localVolCalibrated;
    LVN_vld(i,:,:) = localVolCalibratedNormalizedScale;
end;
% Construct ATM matrix 
for k = 1:size(inputStructure.T_normalized, 2)
    for j = 1:size(validationDates,2)
        [~,closestIndex] = min(abs(S0_v(j)-inputStructure.K_normalized'));
        LVATM_vld(k,j) = LVN_vld(j,k,closestIndex); %% get local vol value atm dynamics for all T
    end;
end;

% calculate model values
for k = 1:size(inputStructure.T_normalized, 2)
    for i = 1:size(validationDates,2)
        RealizedVol_model_vld(k,i) = polyval(paramReg(k,:), LVATM_vld(k,i));
        RealizedVol_model_vld2(k,i) = polyval(paramReg2(k,:), LVATM_vld(k,i));
    end;
end;
%% visualize validation results localVol/realVol
for k = 1:size(inputStructure.T_normalized, 2)
    [LVATM_vld_, sortedIndex] = sort(LVATM_vld(k,:));
    realizedVol_vld_ = realizedVol_vld(k,:);
    realizedVol_vld_ = realizedVol_vld_(sortedIndex);
    RealizedVol_model_vld_ = RealizedVol_model_vld(k,:);
    RealizedVol_model_vld_ = RealizedVol_model_vld_(sortedIndex);
    RealizedVol_model_vld2_ = RealizedVol_model_vld2(k,:);
    RealizedVol_model_vld2_ = RealizedVol_model_vld2_(sortedIndex);  
    
    f(k) = figure; f(k).Visible = 'off'; f(k).Tag = 'ATMonRealized reg';
    s1_v = scatter(LVATM_vld_, realizedVol_vld_);
    hold on
    p1_v = plot(LVATM_vld_, RealizedVol_model_vld_); ylabel('realized vol'); xlabel('LVATM');
    p2_v = plot(LVATM_vld_, RealizedVol_model_vld2_); ylabel('realized vol'); xlabel('LVATM');
    ttl = sprintf('Validation: Lin regr of Relized Vol(t) on ATM Local Vol(t), fixed T=%s',num2str(inputStructure.T_normalized(k)));
    title(ttl);
    legend([s1_v;p1_v;p2_v], 'Realized vol real', 'Realized vol lin regr order1', 'Realized vol lin regr order3');
    legend('boxoff');
end;
%% visualize validation results valDate/realVol
for l = 1:size(inputStructure.T_normalized, 2)
    fg(l) = figure; fg(l).Visible = 'off'; fg(l).Tag = 'ATMonRealized reg';
    s1_v = scatter(1:size(validationDates,2), realizedVol_vld(l,:));
    hold on
    p1_v = plot(1:size(validationDates,2), RealizedVol_model_vld(l,:)); ylabel('realized vol'); xlabel('validation Dates');
    p2_v = plot(1:size(validationDates,2), RealizedVol_model_vld2(l,:)); ylabel('realized vol'); xlabel('validation Dates');
    ttl = sprintf('Validation: Lin regr of Relized Vol(t) on ATM Local Vol(t), fixed T=%s',num2str(inputStructure.T_normalized(l)));
    title(ttl);
    legend([s1_v;p1_v;p2_v], 'Realized vol real', 'Realized vol lin regr order1', 'Realized vol lin regr order3');
    legend('boxoff');
end;

% Calculate MAE matrix
diff_vld = (abs(realizedVol_vld-RealizedVol_model_vld))./realizedVol_vld;
diff_vld2 = (abs(realizedVol_vld-RealizedVol_model_vld2))./realizedVol_vld;
MSE_vld = sum(sum(diff_vld(:,:)))/(size(diff_vld,1)*size(diff_vld,2));
MSE_vld2 = sum(sum(diff_vld2(:,:)))/(size(diff_vld2,1)*size(diff_vld2,2));


%% (3) Generate report
report('3.rpt','-oReportRealized.rtf','-frtf');
pause(7);

%%  close all; % close all figures
   delete(findall(0));


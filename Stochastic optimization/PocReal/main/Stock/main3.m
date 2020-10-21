% In this file we make a run and study the relation between ATM local vol
% and Realized volatility independent variables in a multivariate regression,
% generating reports
clear;
% Intermediate res: calculated local vol( K,T ) for every T0 on input 
% Final res: Regression on local vol(K,T) near the money 

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
report('0.rpt','-oReportCalibRealizedVolMultivar.rtf','-frtf');
pause(7);



%% (3) Regress (multivariate) ATM Local Volatility on realized T volatility of S for FIXED maturity for all calib dates
RelizedVol_ = realizedVol(:,end); % last value of realized vol
LVATM_ = LVATM(:,:)';
[paramReg, sigma, E,CovB,logL] = mvregress(LVATM_, RelizedVol_); %E = 10*1, sigma = 0,027,

RelizedVol_model = LVATM_*paramReg;
 %% visualize regression results
    f(k) = figure; f(k).Visible = 'off'; f(k).Tag = 'ATMonRealized multivariate reg';
    s1 = scatter([1:size(calibrationDates,2)], RelizedVol_);
    hold on
    p1 = plot([1:size(calibrationDates,2)],RelizedVol_model); ylabel('realized vol'); xlabel('calibration dates');
    ttl = sprintf('Multivariate linear regressions of Relized Vol(t) on ATM Local Vol(t)');
    title(ttl);
    legend([s1;p1], 'Realized vol real', 'Realized vol multivar lin regr order1');
    legend('boxoff');
 %% visualize regression stats
 % plot residuals vs indep vars
 for k = 1:size(inputStructure.T_normalized, 2)
    res1(k) = figure; res1(k).Visible = 'off'; res1(k).Tag = 'ATMonRealized multivariate reg';
    s2 = scatter(LVATM(k,:), E');
    ylabel('residuals'); xlabel('local vol');
    hold on
    ttl = sprintf('Analysis: Multivar regr, residuals on local vol at T=%s',num2str(inputStructure.T_normalized(k)));
    title(ttl);
    legend(s2, 'Residuals');
    legend('boxoff');
 end;
 
 % plot  residuals  vs fitted
 res2 = figure; res2.Visible = 'off'; res2.Tag = 'ATMonRealized multivariate reg';
 s3 = scatter(RelizedVol_model', E');
 ylabel('residuals'); xlabel('fitted values');
 hold on
 ttl = sprintf('Analysis: Multivar regr, residuals on fitted values');
 title(ttl);
 legend(s3, 'Residuals');
 legend('boxoff');
 
 % plot histogram of residuals
 res3 = figure; res3.Visible = 'off'; res3.Tag = 'ATMonRealized multivariate reg';
 s4 = histogram(E');
 ylabel('frequency'); xlabel('residuals');
 hold on
 ttl = sprintf('Analysis: Histogram of residuals');
 title(ttl);
 legend(s4, 'Residuals');
 legend('boxoff');
 
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
RealizedVol_model_vld = LVATM_vld(:,:)'*paramReg;

% Calculate MAE matrix
diff_vld = (abs(realizedVol_vld(:,end)-RealizedVol_model_vld))./realizedVol_vld(:,end);
MSE_vld = sum(sum(diff_vld(:,:)))/(size(diff_vld,1)*size(diff_vld,2));
%% visualize validation results valDate/realVol
    RealizedVol_model_vld_ = RealizedVol_model_vld';
    fg(1) = figure; fg(1).Visible = 'off'; fg(1).Tag = 'ATMonRealized multivariate reg';
    s1_v = scatter(1:size(validationDates,2), realizedVol_vld(1,:));
    hold on
    p1_v = plot(1:size(validationDates,2), RealizedVol_model_vld_(1,:)); ylabel('realized vol'); xlabel('validation Dates');
    ttl = sprintf('Validation: Multivar regr of Relized Vol(t) on ATM Local Vol(t)');
    title(ttl);
    legend([s1_v;p1_v], 'Realized vol real', 'Realized vol multivar lin regr order1');
    legend('boxoff');

%% (3) Generate report
report('4.rpt','-oReportRealizedMultivar.rtf','-frtf');
pause(7);

%%  close all; % close all figures
delete(findall(0));


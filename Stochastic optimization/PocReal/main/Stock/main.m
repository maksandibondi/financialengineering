% In this file we make a run and study the relation between ATM local vol
% and Implied volatility, Calibration Dates independent variables,
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


%% (1, 2) Calibrate Local Volatility with adaptive param for all Calibration dates
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

%% (2) Calibrate implied vol with adaptive param for all Calibration dates
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

%% (1, 2) Visualize calibration results
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

%% (1, 2) Generate report
report('0.rpt','-oReport0.rtf','-frtf');
pause(7);



%% (1) Regress ATM Local Volatility on calibration dates for each maturity 
for k = 1:size(inputStructure.T_normalized, 2)
    LVATM_ = LVATM(k,:);
    idxnan = isnan(LVATM_); %% find indices with non NAN values
    [prmReg, rsquared] = polynomialFitting(calibrationDatesNumeric(~idxnan),LVATM_(~idxnan),1);
    [prmReg2, rsquared2] = polynomialFitting(calibrationDatesNumeric(~idxnan),LVATM_(~idxnan),3);
    paramReg(k,:) = prmReg;
    paramReg2(k,:) = prmReg2; 
    rsq(k) = rsquared;
    rsq2(k) = rsquared2;
    LVAMT_model = polyval(paramReg(k,:),calibrationDatesNumeric(~idxnan)); %% values obtained by model
    LVAMT_model2 = polyval(paramReg2(k,:),calibrationDatesNumeric(~idxnan)); %% values obtained by model
    %% visualize regression results
    f(k) = figure; f(k).Visible = 'off'; f(k).Tag = 'ATMonTime reg';
    s1 = scatter(calibrationDatesNumeric(~idxnan),LVATM_(~idxnan));
    hold on
    p1 = plot(calibrationDatesNumeric(~idxnan),LVAMT_model); xlabel('t'); ylabel('LVATM');
    p2 = plot(calibrationDatesNumeric(~idxnan),LVAMT_model2); xlabel('t'); ylabel('LVATM');
    ttl = sprintf('Linear regressions of ATM Local Vol on calibration date t for fixed T=%s',num2str(inputStructure.T_normalized(k)));
    title(ttl);
    legend([s1;p1;p2], 'LV ATM real', 'LV ATM lin regr order1', 'LV ATM lin regr order3');
end;
% Conclusion : sometimes negative coef of determination, bad fitting

%% (1) Validate model on validation points (interm and extrap) with MAE (mean abs error)
 % calculate model values
for k = 1:size(inputStructure.T_normalized, 2)
    for i = 1:size(validationDates,2)
        LVATM_model_vld(k,i) = polyval(paramReg(k,:), validationDatesNumeric(i));
        LVATM_model_vld2(k,i) = polyval(paramReg2(k,:), validationDatesNumeric(i));
    end;
end;
  % calculate real values : Calibrate local vol with adaptive param for all validation dates
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

 % Calculate MAE matrix
   diff_vld = (abs(LVATM_vld-LVATM_model_vld))./LVATM_vld;
   diff_vld2 = (abs(LVATM_vld-LVATM_model_vld2))./LVATM_vld;
   MSE_vld = sum(sum(diff_vld(:,:)))/(size(diff_vld,1)*size(diff_vld,2));
   MSE_vld2 = sum(sum(diff_vld2(:,:)))/(size(diff_vld2,1)*size(diff_vld2,2));
   
%% (1) Generate report
   report('1.rpt','-oReport1.rtf','-frtf');
   pause(7);
   %close all; % close all figures

   
%% (2) Regress ATM Local Volatility on ATM implied volatility for each maturity for all calib dates
for k = 1:size(inputStructure.T_normalized, 2)
    [LVATM_, sortedIndex] = sort(LVATM(k,:));
    VolImpATM_ = VolImpATM(k,:);
    VolImpATM_ = VolImpATM_(sortedIndex);
    idxnan = isnan(LVATM_); %% find indices with non NAN values
    [prmReg, rsquared] = polynomialFitting(LVATM_(~idxnan),VolImpATM_(~idxnan),1);
    [prmReg2, rsquared2] = polynomialFitting(LVATM_(~idxnan),VolImpATM_(~idxnan),3);
    paramReg(k,:) = prmReg;
    paramReg2(k,:) = prmReg2;
    rsq(k) = rsquared;
    rsq2(k) = rsquared2;
    VolImpATM_model = polyval(paramReg(k,:),LVATM_(~idxnan)); %% values obtained by model
    VolImpATM_model2 = polyval(paramReg2(k,:),LVATM_(~idxnan)); %% values obtained by model    

    %% visualize regression results
    f(k) = figure; f(k).Visible = 'off'; f(k).Tag = 'ATMonImplied reg';
    s1 = scatter(LVATM_(~idxnan),VolImpATM_(~idxnan));
    hold on
    p1 = plot(LVAMT_model,VolImpATM_(~idxnan)); ylabel('impvol ATM'); xlabel('LVATM');
    p2 = plot(LVAMT_model2,VolImpATM_(~idxnan)); ylabel('impvol ATM'); xlabel('LVATM');
    ttl = sprintf('Linear regressions of Implied vol(t) on ATM Local Vol(t) with fixed T=%s',num2str(inputStructure.T_normalized(k)));
    title(ttl);
    legend([s1;p1;p2], 'LV ATM real', 'LV ATM lin regr order1', 'LV ATM lin regr order3');
end;
% Conclusion : sometimes negative coef of determination, bad fitting

%% (2) Validate model on validation points (interm and extrap) with MAE (mean abs error)
% calculate implied vol real values on validation points
for i = 1:size(validationDates,2)
    inputStructure.date = cell2mat(validationDates(i)); 
    [impliedVolCalibrated, impliedVolCalibratedNormalizedScale, S] = impliedVolCalibratorStockRunner(inputStructure);
    S0_v(i) = S;
    VolImp_vld(i,:,:) = impliedVolCalibratedNormalizedScale;
end;
% Construct ATM matrix 
for k = 1:size(inputStructure.T_normalized, 2)
   for j = 1:size(validationDates,2)
        [~,closestIndex] = min(abs(S0(j)-inputStructure.K_normalized'));
        VolImpATM_vld(k,j) = VolImp_vld(j,k,closestIndex); %% get imp vol value atm dynamics for all T
    end;
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
        VolImpATM_model_vld(k,i) = polyval(paramReg(k,:), LVATM_vld(k,i));
        VolImpATM_model_vld2(k,i) = polyval(paramReg2(k,:), LVATM_vld(k,i));
    end;
end;

% Calculate MAE matrix
diff_vld = (abs(VolImpATM_vld-VolImpATM_model_vld))./VolImpATM_vld;
diff_vld2 = (abs(VolImpATM_vld-VolImpATM_model_vld2))./VolImpATM_vld;
MSE_vld = sum(sum(diff_vld(:,:)))/(size(diff_vld,1)*size(diff_vld,2));
MSE_vld2 = sum(sum(diff_vld2(:,:)))/(size(diff_vld2,1)*size(diff_vld2,2));

%% (2) Generate report
   report('2.rpt','-oReport2.rtf','-frtf');
   pause(7);
   
 
   
   
%% close all; % close all figures
   delete(findall(0));









%% (4) Regress ATM Local Volatility on VIX value at t for each maturity for all calib dates


%% (5) Regress ATM Local Volatility on VIX valye at T for each maturity for all calib dates


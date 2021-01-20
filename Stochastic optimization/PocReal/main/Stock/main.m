% In this file we make a run and study the relation between ATM local vol
% and Implied volatility, Calibration Dates independent variables for each maturity,
% generating reports
clear;
% Intermediate res: calculated local vol( K,T ) for every T0 on input 
% Final res: Regression on local vol(K,T) near the money . one regression
% by maturity


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
inputStructure.includeVolAssumption = 'zero'; % implied
inputStructure.vasicek_assumption = 0; % zero
inputStructure.outputdir = fullfile(pwd,'../../Results/Results_MSFT');

inputStructure.visualizeResults = 0;
inputStructure.isNormalizedScale = 1;
inputStructure.T_normalized = [0.05, 0.1, 0.15, 0.2, 0.25, 0.4];
%inputStructure.T_normalized = [0.05, 0.1, 0.15, 0.2, 0.25, 0.4];
inputStructure.K_normalized = 100:5:230;

Stock_normalization_factor = 1000;  % as price is around 100 (for affichage only)

calibrationDates = {'01/03/2020', '01/10/2020', '01/24/2020', '01/31/2020', '02/07/2020',  '02/28/2020', '03/06/2020',  '03/20/2020', '03/27/2020', '04/03/2020'};
validationDates = {'01/17/2020','02/14/2020', '02/21/2020', '03/13/2020', '04/09/2020', '04/17/2020'};

%calibrationDates = {'01/03/2020', '01/10/2020', '01/17/2020', '01/24/2020', '01/31/2020', '02/07/2020', '02/14/2020', '02/21/2020','02/28/2020', '03/06/2020'};
%validationDates = { '03/13/2020', '03/20/2020', '03/27/2020', '04/03/2020', '04/09/2020', '04/17/2020'};

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
    xlabel('Strike'); ylabel('Maturity T'); zlabel('Local vol');
    ttl = sprintf('Local volatility surface at t =%s i.e. t=%s', char(calibrationDates(i)), num2str(calibrationDatesNumeric(i)));
    title(ttl);
    f2(i) = figure; f2(i).Visible = 'off'; f2(i).Tag = 'Implied Vol Surface';
    surf(inputStructure.K_normalized, inputStructure.T_normalized, squeeze(VolImp(i,:,:)));
    xlabel('Strike'); ylabel('Maturity T'); zlabel('Implied vol');
    ttl = sprintf('Implied volatility surface at t =%s i.e. t=%s', char(calibrationDates(i)), num2str(calibrationDatesNumeric(i)));
    title(ttl);    
end;

%% (1, 2) Generate report
report('0.rpt','-oReportCalibLVonDates.rtf','-frtf');
pause(7);


stat = cell(size(inputStructure.T_normalized, 2), 1);
stat2 = cell(size(inputStructure.T_normalized, 2), 1);
anov = cell(size(inputStructure.T_normalized, 2), 1);
anov2 = cell(size(inputStructure.T_normalized, 2), 1);

%% (1) Regress ATM Local Volatility on calibration dates for each maturity 
for k = 1:size(inputStructure.T_normalized, 2)
    LVATM_ = LVATM(k,:);
    idxnan = isnan(LVATM_); %% find indices with non NAN values
    [prmReg, rsquared, resid, mdl] = polynomialFitting(calibrationDatesNumeric(~idxnan),LVATM_(~idxnan),1);
    [prmReg2, rsquared2, resid2, mdl2] = polynomialFitting(calibrationDatesNumeric(~idxnan),LVATM_(~idxnan),3);
    paramReg(k,:) = prmReg;
    paramReg2(k,:) = prmReg2; 
    rsq(k) = rsquared;
    rsq2(k) = rsquared2;
    stat{k} = mdl;
    stat2{k} = mdl2;
    anov{k} = anova(mdl);
    anov2{k} = anova(mdl2);
    [pvaldw{k},DW{k}] = dwtest(mdl,'exact','both');
    [pvaldw2{k},DW2{k}] = dwtest(mdl2,'exact','both');
    LVAMT_model(k,:) = polyval(paramReg(k,:),calibrationDatesNumeric(~idxnan)); %% values obtained by model
    LVAMT_model2(k,:) = polyval(paramReg2(k,:),calibrationDatesNumeric(~idxnan)); %% values obtained by model
    %% visualize regression results
    f1(k) = figure; f1(k).Visible = 'off'; f1(k).Tag = 'ATMonTime reg';
    s1 = scatter(calibrationDatesNumeric(~idxnan),LVATM_(~idxnan));
    hold on
    p1 = plot(calibrationDatesNumeric(~idxnan),LVAMT_model(k,:),'-s'); xlabel('t'); ylabel('LVATM');
    p2 = plot(calibrationDatesNumeric(~idxnan),LVAMT_model2(k,:),'-d'); xlabel('t'); ylabel('LVATM');
    ttl = sprintf('Linear regressions of ATM Local Vol on calibration date t for fixed T=%s',num2str(inputStructure.T_normalized(k)));
    title(ttl);
    legend([s1;p1;p2], 'LV ATM real', 'LV ATM lin regr order1', 'LV ATM lin regr order3');
    legend('boxoff');
    %% visualize regression stats  
 % plot residuals vs indep vars
    f2(k) = figure; f2(k).Visible = 'off'; f2(k).Tag = 'ATMonTime reg';
    s2 = scatter(calibrationDatesNumeric(~idxnan), resid);
    ylabel('residuals'); xlabel('local vol');
    hold on
    ttl = sprintf('Analysis: Multivar regr, residuals on calib dates at T=%s',num2str(inputStructure.T_normalized(k)));
    title(ttl);
    legend(s2, 'Residuals');
    legend('boxoff');
 % plot residuals vs indep vars (3rd order reg)
    f3(k) = figure; f3(k).Visible = 'off'; f3(k).Tag = 'ATMonTime reg';
    s3 = scatter(calibrationDatesNumeric(~idxnan), resid2);
    ylabel('residuals'); xlabel('local vol');
    hold on
    ttl = sprintf('Analysis: Multivar regr, 3rd ordeer residuals on calib dates at T=%s',num2str(inputStructure.T_normalized(k)));
    title(ttl);
    legend(s3, 'Residuals');
    legend('boxoff');
 % plot histogram of residuals
    f4(k) = figure; f4(k).Visible = 'off'; f4(k).Tag = 'ATMonTime reg';
    s4 = histogram(resid,10);
    ylabel('frequency'); xlabel('residuals');
    hold on
    ttl = sprintf('Analysis: Histogram of residuals, T=%s',num2str(inputStructure.T_normalized(k)));
    title(ttl);
    legend(s4, 'Residuals');
    legend('boxoff');
  % plot histogram of residuals 3rd order
    f5(k) = figure; f5(k).Visible = 'off'; f5(k).Tag = 'ATMonTime reg';
    s5 = histogram(resid2,10);
    ylabel('frequency'); xlabel('residuals');
    hold on
    ttl = sprintf('Analysis: Histogram of residuals 3rd order, T=%s',num2str(inputStructure.T_normalized(k)));
    title(ttl);
    legend(s5, 'Residuals');
    legend('boxoff');
  % Normal distribution of residuals
    fig(1) = figure; fig(1).Visible = 'off'; fig(1).Tag = 'ATMonTime reg';
    fig = plotResiduals(mdl,'probability');
    ttl = sprintf('Analysis: Proba plot of residuals, T=%s',num2str(inputStructure.T_normalized(k)));
    title(ttl);
  % Normal distribution of residuals 3rd order
    fig(2) = figure; fig(2).Visible = 'off'; fig(2).Tag = 'ATMonTime reg';
    fig = plotResiduals(mdl2,'probability');
    ttl = sprintf('Analysis: Proba plot of residuals 3rd order, T=%s',num2str(inputStructure.T_normalized(k)));
    title(ttl);    
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
for k = 1:size(inputStructure.T_normalized, 2)
   MSE_vld(k) = sum(diff_vld(k,:))/size(diff_vld,2);
   MSE_vld2(k) = sum(diff_vld2(k,:))/size(diff_vld2,2);
end;
   %% visualize validation results valDate/localVol
for l = 1:size(inputStructure.T_normalized, 2)
    f6(l) = figure; f6(l).Visible = 'off'; f6(l).Tag = 'ATMonTime reg';
    s1_v = scatter(1:size(validationDates,2), LVATM_vld(l,:));
    hold on
    p1_v = plot(1:size(validationDates,2), LVATM_model_vld(l,:), '-s'); ylabel('local vol'); xlabel('validation Dates');
    p2_v = plot(1:size(validationDates,2), LVATM_model_vld2(l,:), '-d'); ylabel('local vol'); xlabel('validation Dates');
    ttl = sprintf('Validation: Lin regr of Local Vol(t) on Dates(t), fixed T=%s',num2str(inputStructure.T_normalized(k)));
    title(ttl);
    legend([s1_v;p1_v;p2_v], 'local vol real', 'local vol lin regr order1', 'local vol lin regr order3');
    legend('boxoff');
end;
%% (1) Generate report
   report('1.rpt','-oReportLVonDates.rtf','-frtf');
   pause(7);
   %close all; % close all figures

   
%% (2) Regress ATM Local Volatility on ATM implied volatility for each maturity for all calib dates
for k = 1:size(inputStructure.T_normalized, 2)
    [LVATM_, sortedIndex] = sort(LVATM(k,:));
    VolImpATM_ = VolImpATM(k,:);
    VolImpATM_ = VolImpATM_(sortedIndex);
    idxnan = isnan(LVATM_); %% find indices with non NAN values
    [prmReg, rsquared, resid, mdl] = polynomialFitting(LVATM_(~idxnan),VolImpATM_(~idxnan),1);
    [prmReg2, rsquared2, resid2, mdl2] = polynomialFitting(LVATM_(~idxnan),VolImpATM_(~idxnan),3);
    paramReg(k,:) = prmReg;
    paramReg2(k,:) = prmReg2;
    rsq(k) = rsquared;
    rsq2(k) = rsquared2;
    stat{k} = mdl;
    stat2{k} = mdl2;
    anov{k} = anova(mdl);
    anov2{k} = anova(mdl2);
    [pvaldw{k},DW{k}] = dwtest(mdl,'exact','both');
    [pvaldw2{k},DW2{k}] = dwtest(mdl2,'exact','both');
    VolImpATM_model(k,:) = polyval(paramReg(k,:),LVATM_(~idxnan)); %% values obtained by model
    VolImpATM_model2(k,:) = polyval(paramReg2(k,:),LVATM_(~idxnan)); %% values obtained by model    

%% visualize regression results
    f7(k) = figure; f7(k).Visible = 'off'; f7(k).Tag = 'ATMonImplied reg';
    s1 = scatter(LVATM_(~idxnan),VolImpATM_(~idxnan));
    hold on
    p1 = plot(LVATM_(~idxnan),VolImpATM_model(k,:),'-s'); ylabel('impvol ATM'); xlabel('LVATM');
    p2 = plot(LVATM_(~idxnan),VolImpATM_model2(k,:),'-d'); ylabel('impvol ATM'); xlabel('LVATM');
    ttl = sprintf('Linear regressions of Implied vol(t) on ATM Local Vol(t) with fixed T=%s',num2str(inputStructure.T_normalized(k)));
    title(ttl);
    legend([s1;p1;p2], 'Implied vol real', 'Implied vol lin regr order1', 'Implied vol lin regr order3');
    legend('boxoff');
     %% visualize regression stats  
 % plot residuals vs indep vars
    f8(k) = figure; f8(k).Visible = 'off'; f8(k).Tag = 'ATMonImplied reg';
    s2 = scatter(LVATM_(~idxnan), resid);
    ylabel('residuals'); xlabel('local vol');
    hold on
    ttl = sprintf('Analysis Imp: Multivar regr, residuals on implied vol at T=%s',num2str(inputStructure.T_normalized(k)));
    title(ttl);
    legend(s2, 'Residuals');
    legend('boxoff');
 % plot residuals vs indep vars (3rd order reg)
    f9(k) = figure; f9(k).Visible = 'off'; f9(k).Tag = 'ATMonImplied reg';
    s3 = scatter(LVATM_(~idxnan), resid2);
    ylabel('residuals'); xlabel('local vol');
    hold on
    ttl = sprintf('Analysis Imp: Multivar regr, 3rd ordeer residuals on implied vol at T=%s',num2str(inputStructure.T_normalized(k)));
    title(ttl);
    legend(s3, 'Residuals');
    legend('boxoff');
 % plot histogram of residuals
    f10(k) = figure; f10(k).Visible = 'off'; f10(k).Tag = 'ATMonImplied reg';
    s4 = histogram(resid,10);
    ylabel('frequency'); xlabel('residuals');
    hold on
    ttl = sprintf('Analysis Imp: Histogram of residuals, T=%s',num2str(inputStructure.T_normalized(k)));
    title(ttl);
    legend(s4, 'Residuals');
    legend('boxoff');
  % plot histogram of residuals 3rd order
    f11(k) = figure; f11(k).Visible = 'off'; f11(k).Tag = 'ATMonImplied reg';
    s5 = histogram(resid2,10);
    ylabel('frequency'); xlabel('residuals');
    hold on
    ttl = sprintf('Analysis Imp: Histogram of residuals 3rd order, T=%s',num2str(inputStructure.T_normalized(k)));
    title(ttl);
    legend(s5, 'Residuals');
    legend('boxoff');
  % Normal distribution of residuals
    figg(1) = figure; figg(1).Visible = 'off'; figg(1).Tag = 'ATMonImplied reg';
    fig = plotResiduals(mdl,'probability');
    ttl = sprintf('Analysis: Proba plot of residuals, T=%s',num2str(inputStructure.T_normalized(k)));
    title(ttl);
  % Normal distribution of residuals 3rd order
    figg(2) = figure; figg(2).Visible = 'off'; figg(2).Tag = 'ATMonImplied reg';
    fig = plotResiduals(mdl2,'probability');
    ttl = sprintf('Analysis: Proba plot of residuals 3rd order, T=%s',num2str(inputStructure.T_normalized(k)));
    title(ttl);    
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
for k = 1:size(inputStructure.T_normalized, 2)
    MSE_vld(k) = sum(diff_vld(k,:))/size(diff_vld,2);
    MSE_vld2(k) = sum(diff_vld2(k,:))/size(diff_vld2,2);
end;
 %% visualize validation results localVol/impliedVol
for l = 1:size(inputStructure.T_normalized, 2)
    [LVATM_vld_, sortedIndex] = sort(LVATM_vld(k,:));
    VolImp_vld_ = VolImp_vld(k,:);
    VolImp_vld_ = VolImp_vld_(sortedIndex);
    
    f12(l) = figure; f12(l).Visible = 'off'; f12(l).Tag = 'ATMonImplied reg';
    s1_v = scatter(LVATM_vld_,VolImp_vld_);
    hold on
    p1_v = plot(LVATM_vld_, VolImpATM_model_vld(l,:),'-s'); ylabel('implied vol'); xlabel('local vol');
    p2_v = plot(LVATM_vld_, VolImpATM_model_vld2(l,:), '-d'); ylabel('implied vol'); xlabel('local vol');
    ttl = sprintf('Validation: Lin regr of  Impied vol(t) on Local Vol(t), fixed T=%s',num2str(inputStructure.T_normalized(k)));
    title(ttl);
    legend([s1_v;p1_v;p2_v], 'local vol real', 'local vol lin regr order1', 'local vol lin regr order3');
    legend('boxoff');
end;
%% (2) Generate report
   report('2.rpt','-oReportImpVolonLV.rtf','-frtf');
   pause(7);
   
   
%% Visualize stock price and aggregated plots of validation/calibration values
[DatesNumeric, sorted_index] = sort([calibrationDatesNumeric, validationDatesNumeric]);
S_unsorted = [S0, S0_v];
S = S_unsorted(sorted_index);
Dates_unsorted = [calibrationDates, validationDates];
Dates = Dates_unsorted(sorted_index);

for k = 1:size(inputStructure.T_normalized, 2)
    % show modeled and real local and implied vols
    localVolUnsorted = [LVATM(k,:), LVATM_vld(k,:)];
    localVolModelUnsorted = [LVAMT_model(k,:), LVATM_model_vld(k,:)];
    localVolModelUnsorted2 = [LVAMT_model2(k,:), LVATM_model_vld2(k,:)];
    impVolUnsorted = [VolImpATM(k,:), VolImp_vld(k,:)];
    impVolModelUnsorted = [VolImpATM_model(k,:), VolImpATM_model_vld(k,:)];
    impVolModelUnsorted2 = [VolImpATM_model2(k,:), VolImpATM_model_vld2(k,:)];
    
    localVOL = localVolUnsorted(sorted_index);
    localVOLModel = localVolModelUnsorted(sorted_index);
    localVOLModel2 = localVolModelUnsorted2(sorted_index);
    impVOL = impVolUnsorted(sorted_index);
    impVOLModel = impVolModelUnsorted(sorted_index);
    impVOLModel2 = impVolModelUnsorted2(sorted_index);

    f13(k) = figure; f13(k).Visible = 'off'; f13(k).Tag = 'AggregatedPlots';
    subplot(2,1,1);
    s1 = scatter(datenum(Dates,'mm/dd/yyyy'), localVOL, 'o');
    hold on
    s2 = plot(datenum(Dates,'mm/dd/yyyy'), localVOLModel, '-s');
    s3 = plot(datenum(Dates,'mm/dd/yyyy'), localVOLModel2, '-d');    
    ylabel('Local Vol'); xlabel('Dates');  set(gca,'xtick',datenum(Dates,'mm/dd/yyyy')); set(gca,'FontSize',4); datetick('x',29,'keepticks');
    ttl = sprintf('Local vol on time models/true with predictions with fixed T=%s', num2str(inputStructure.T_normalized(k)));
    legend([s1;s2;s3], 'Local vol calibrated', 'Local vol modeled', 'Local vol modeled 3rd order');
    title(ttl);
    legend('boxoff');
    
    subplot(2,1,2); 
    p1 = plot(datenum(Dates,'mm/dd/yyyy'), S);
    ylabel('Stcok Spot Price'); xlabel('Dates');  set(gca,'xtick',datenum(Dates,'mm/dd/yyyy')); set(gca,'FontSize',4); datetick('x',29,'keepticks');
    legend('boxoff');
    
    
    f14(k) = figure; f14(k).Visible = 'off'; f14(k).Tag = 'AggregatedPlots';
    subplot(2,1,1);
    s1 = scatter(datenum(Dates,'mm/dd/yyyy'), impVOL);
    hold on
    s2 = plot(datenum(Dates,'mm/dd/yyyy'), impVOLModel,'-s');
    s3 = plot(datenum(Dates,'mm/dd/yyyy'), impVOLModel2, '-d');
    ylabel('Imp Vol'); xlabel('Dates');  set(gca,'xtick',datenum(Dates,'mm/dd/yyyy')); set(gca,'FontSize',4); datetick('x',29,'keepticks');
    ttl = sprintf('Implied vol on LocalVol models/true with predictions with fixed T=%s', num2str(inputStructure.T_normalized(k)));
    legend([s1;s2;s3], 'Imp vol calibrated', 'Imp vol modeled', 'Imp vol modeled 3rd order');
    title(ttl);
    legend('boxoff');
    
    subplot(2,1,2); 
    p1 = plot(datenum(Dates,'mm/dd/yyyy'), S);
    ylabel('Stock spot price'); xlabel('Dates');  set(gca,'xtick',datenum(Dates,'mm/dd/yyyy')); set(gca,'FontSize',4); datetick('x',29,'keepticks');
    legend('boxoff');
end;

report('aggr0.rpt','-oReportAggregatedLVonDatesImp.rtf','-frtf');
pause(7); 
   
   
%% close all; % close all figures
   delete(findall(0));









%% (4) Regress ATM Local Volatility on VIX value at t for each maturity for all calib dates


%% (5) Regress ATM Local Volatility on VIX valye at T for each maturity for all calib dates


function [isLocalVolValidated] = localVolValidator_NMaturities(inputStructure)

%% Improt folders
addpath('../../parser','-end');
addpath('../../pricing','-end');
addpath('../../formatting', '-end');
addpath('../../graphics', '-end');
addpath('../../utils', '-end');


[localVol, K, T, ~] = parseValidationFile(inputStructure.validationFile, inputStructure.validationDate);
[K, T_numeric, S, Vmarket, VolImp] = parseCBOE(inputStructure.marketDataFile, inputStructure.ticker, inputStructure.validationDate, inputStructure.optionType);
[T, Vmarket, VolImp] = additionalFormatting( T,  Vmarket, VolImp, S, K, inputStructure.optionType);

%% Format obtained data
% indexOfT = find(abs(T_numeric-T)<= 0.00001); %% find the index of T in the array T_numeric taking into account float difference
% Vmarket = Vmarket(indexOfT,:);
% localVol = localVol';
% localVol = [zeros(1,size(localVol,2));localVol];
% [T, Vmarket, VolImp] = additionalFormatting( T,  Vmarket, VolImp, S, K, inputStructure.optionType);
% K = K';

%% Get disc
discretization_num_K = size(K,2);
discretization_num_T = size(T,1);
%% Get Vasicek assumption and set r
if (inputStructure.vasicek_assumption)
    r = (20-S)/100;
else
    r = inputStructure.r;
end;

%% Validation
u(1,:,:) = Pricer_dupire(localVol, K, T, discretization_num_K, discretization_num_T, S, r, inputStructure.discretizationType, inputStructure.nonuniform_method);
[fitness, ~, ~] = sumOfSqrDif_(u(1,:,:), Vmarket(:,:), S, K, inputStructure.epsilon, 1); % cost funtion for n-th member of population
 

if fitness < 0.3
     isLocalVolValidated = 1;
 else
     isLocalVolValidated = 0;
 end;
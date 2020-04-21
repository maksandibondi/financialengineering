function[K, T_numeric_sorted, T_sorted, S, Vmarket] = parseStockCSV(filename, symbol, quoteDate, optionType)

% Input filename, symbol (VIX), quoteDate (every date on which calculate
% vol), optionType = C

% Output SELECT strike, expiration (fraction of year), close(K,T), impliedVol(K,T) FROM table WHERE quoteDate.
%

addpath('../../utils', '-end');

%% Parse csv file and format data for certain date

allData = readtable(filename); 

rows = strcmp(allData.symbol, symbol) & strcmp(allData.call_put, optionType) & strcmp(allData.date, quoteDate);
vars = {'expiration','strike','ask', 'bid', 'adjustedClose'};
workTable = allData(rows, vars); % Filtered table of values for given quoted date and instrument, option type

%% Select expiration
T_non_filtered = allData(rows, 'expiration');
T = cell2mat(table2array(unique(T_non_filtered)));
numOfT = size(T,1);
T_numeric = daysact(quoteDate, T)/365;
T_numeric_sorted = sort(T_numeric);
for idx = 1:size(T_numeric,1)
    %% used to know which non-numerical value corresponds to each of sorted values
    original_index (idx) = find(abs(T_numeric - T_numeric_sorted(idx))<0.00001); 
    T_sorted(idx, :) = T(original_index (idx), :);
end;
display(numOfT);

%% Select strike values and do uniform grid 
allK = workTable(:, {'strike'});
max_el = max(table2array(allK));
min_el = min(table2array(allK));
step = getStep(min_el, max_el);
K = transp(min_el : step : max_el) ;
numOfK = size(K,1);
display(numOfK);

%% Select underlying value (same for all maturities as spot)
S = workTable(1, {'adjustedClose'}).('adjustedClose');
display(S);

%% Select close prices and implied volatilities
Vmarket = zeros(numOfT, numOfK);
%VolImp = zeros(numOfT, numOfK);
for i = 1 : numOfT
    % Select strike/close pair table for needed maturity
    rows = strcmp(workTable.expiration, T_sorted(i,:));
    nonuniformK = workTable(rows, {'strike'});
    ask_price = workTable(rows, {'ask'});
    bid_price = workTable(rows, {'bid'});
    %vols_implied = workTable(rows, {'implied_volatility_1545'});
    for j = 1 : numOfK
        [lia, Index] = ismember(K(j,1), table2array( nonuniformK)); % any element in nonuniforK, get 'close'
        if (Index ~= 0)
           Vmarket(i,j) =  bid_price(Index,1).bid + (ask_price(Index,1).ask - bid_price(Index,1).bid)/2; %select a price average between bid/ask
           %VolImp(i,j) = vols_implied(Index,1).implied_volatility_1545;
        end
    end
    
   Vmarket(i,:) = interpolatePayoffVector(Vmarket(i, :), K(:,1), S, optionType);
   %VolImp(i,:) = interpolateImpliedVolVector(VolImp(i,:));
end


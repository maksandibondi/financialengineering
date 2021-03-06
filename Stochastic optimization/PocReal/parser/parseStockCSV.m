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
    nonuniformK = table2array(workTable(rows, {'strike'}));
    ask_price = table2array(workTable(rows, {'ask'}));
    bid_price = table2array(workTable(rows, {'bid'}));
    %vols_implied = workTable(rows, {'implied_volatility_1545'});
    
    price = (bid_price) + (ask_price - bid_price)/2;  
    %Clean zero values from price vector
    for l = 2:size(price)-1
        if (price(l) == 0 || isnan(price(l)))
            price(l) = (price(l+1)-price(l-1))/2;
        end;
    end;
    
    %% Provide K(0) and K(end) values to get the same scale whatever maturity we use and interpolate values in betweem
    if (optionType == 'C')
        if abs(nonuniformK(1)-K(1,1))>0.0001
            price = cat(1,S-K(1,1),price); % add initial value S-K
            nonuniformK = cat(1,K(1,1),nonuniformK); %add initial K = Kmin << S
        end;
        if abs(nonuniformK(end)-K(end,1))>0.0001
            price = cat(1,price,0.0001); % add final value 0
            nonuniformK = cat(1,nonuniformK,K(end,1)); %add final K = Kmax >> S
        end;
    elseif (optionType == 'P')
        if abs(nonuniformK(1)-K(1,1))>0.0001
            price = cat(1,0.0001,price);  % add initial value 0
            nonuniformK = cat(1,K(1,1),nonuniformK); %add initial K = Kmin << S
        end;
        if abs(nonuniformK(end)-K(end,1))>0.0001
            price = cat(1,price,K(end,1)-S); % add final value K-S
            nonuniformK = cat(1,nonuniformK,K(end,1)); %add final K = Kmax >> S
        end; 
    end;
    
    Vmarket(i,:) = interp1(nonuniformK, price, K(:,1), 'linear', 'extrap');
end


function[K, T_numeric, S, Vmarket, VolImp] = parseCBOE(filename, symbol, quoteDate, optionType)

% Input filename, symbol (VIX), quoteDate (every date on which calculate
% vol), optionType = C

% Output SELECT strike, expiration (fraction of year), close(K,T), impliedVol(K,T) FROM table WHERE quoteDate.
%

%% Parse cboe file and format data for certain date

allData = readtable(filename); 

rows = strcmp(allData.underlying_symbol, symbol) & strcmp(allData.option_type, optionType) & strcmp(allData.quote_date, quoteDate);
vars = {'expiration','strike','option_type', 'close', 'active_underlying_price_1545', 'implied_volatility_1545'};
workTable = allData(rows, vars); % Filtered table of values for given quoted date and instrument, option type

%% Select expiration
T_non_filtered = allData(rows, 'expiration');
T = cell2mat(table2array(unique(T_non_filtered)));
numOfT = size(T,1);
T_numeric = daysact(quoteDate, T)/365;
display(numOfT);

%% Select strike values and do uniform grid 
allK = workTable(:, {'strike'});
max_el = max(table2array(allK));
min_el = min(table2array(allK));
K = transp(min_el : 1 : max_el) ;
numOfK = size(K,1);
display(numOfK);

%% Select underlying value (same for all maturities as spot)
S = workTable(1, {'active_underlying_price_1545'}).active_underlying_price_1545;
display(S);

%% Select close prices and implied volatilities
Vmarket = zeros(numOfT, numOfK);
VolImp = zeros(numOfT, numOfK);
for i = 1 : numOfT
    % Select strike/close pair table for needed maturity
    rows = strcmp(workTable.expiration, T(i,:));
    nonuniformK = workTable(rows, {'strike'});
    prices = workTable(rows, {'close'});
    vols_implied = workTable(rows, {'implied_volatility_1545'});
    for j = 1 : numOfK
        [lia, Index] = ismember(K(j,1), table2array( nonuniformK)); % any element in nonuniforK, get 'close'
        if (Index ~= 0)
           Vmarket(i,j) = prices(Index,1).close;
           VolImp(i,j) = vols_implied(Index,1).implied_volatility_1545;
        end
    end
    
   Vmarket(i,:) = interpolatePayoffVector(Vmarket(i, :), K(:,1), S, optionType);
   VolImp(i,:) = interpolateImpliedVolVector(VolImp(i,:));
end


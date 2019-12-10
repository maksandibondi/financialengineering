%function[K, T, V(K,T), sigmaImp] = parseCBOE(filename, symbol, quoteDate, optionType)

% Input filename, symbol (VIX), quoteDate (every date on which calculate
% vol), optionType = C

% Output SELECT strike, expiration (fraction of year), close(K,T), impliedVol(K,T) FROM table WHERE quoteDate.
%


% Parse cboe file and format data for certain date

allData = readtable('../../resources/UnderlyingOptionsEODCalcs_2018-11.csv'); 

rows = strcmp(allData.underlying_symbol,'^VIX') & strcmp(allData.option_type,'C') & strcmp(allData.quote_date, '2018-11-15');
vars = {'expiration','strike','option_type', 'close', 'active_underlying_price_1545', 'implied_volatility_1545'};
workTable = allData(rows, vars); % Filtered table of values for given quoted date and instrument, option type

% Select expiration
T_non_filtered = allData(rows, 'expiration');
T = unique(T_non_filtered);
numOfT = size(T,1);
display(numOfT);

% Select strike values and do uniform grid 
%rows = strcmp(workTable.expiration, T(1,1).expiration);
allK = workTable(:, {'strike'});
max_el = max(table2array(allK));
min_el = min(table2array(allK));

K = min_el : 1 : max_el ;
numOfK = size(K,2);
display(numOfK);

Vmarket = zeros(numOfT, numOfK);
% Select close prices
for i = 1 : numOfT
    for j = 1 : numOfK
        %display(j+(i-1)*size(K,1));
        Vmarket(i,j) = interpolateMarketPrice(i,j,
        %Vmarket(i,j) = workTable(j + (i-1)*size(K,1), :).close; % function that generates interpolated or not value for given T,K
    end
end


% Choose All different expirations to T
% Choose all different strikes to K
% Choose all different Close price to matrix V(K,T). Interpolate where no price available
% Choose implied volatility to matrix Imp(K,T) to use as seed . Interpolate where no vol  available
% Choose underlying asset price as S0 scalar
% Test function with real input quote_date


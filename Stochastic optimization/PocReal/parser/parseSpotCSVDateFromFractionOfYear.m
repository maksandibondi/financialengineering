function [T_dates] = parseSpotCSVDateFromFractionOfYear(filename, startDates, endDatesFraction)
addpath('../../utils', '-end');

allData = readtable(filename); 
formatOut = 'yyyy-mm-dd';

for i = 1 : size(startDates,2)
    dates(i,:) = datestr(startDates(i),formatOut);
end;

for i = 1 : size(dates,1)
    for j = 1 : size(endDatesFraction, 2)   
        endDate = addtodate(datenum(dates(i,:),formatOut), 365*endDatesFraction(j), 'day');
        mat = datestr(endDate, formatOut);
        rows = strcmp(allData.Date, mat);
        %if (~ismember(1,rows))
        [~,ind1] = min(abs(datenum(allData.Date,formatOut)-endDate));
        %vars = {'Date'};
        %foundRows = allData(rows, vars); % Filtered table of values for given quoted date and instrument, option type
        T_dates(i,j,:) = allData.Date(ind1, :);
    end;
end;



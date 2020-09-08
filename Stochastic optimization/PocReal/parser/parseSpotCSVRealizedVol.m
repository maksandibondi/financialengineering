function [vol] = parseSpotCSVRealizedVol(filename, dates, maturities)
addpath('../../utils', '-end');

allData = readtable(filename); 
formatOut = 'yyyy-mm-dd';
for i = 1:size(dates,2)
    datesnew(i,:) = datestr(dates(i),formatOut);
end;
dates = datesnew;

for i = 1:size(maturities,1)
    for j = 1:size(maturities,2)
        maturitiesnew(i,j,:) = datestr(maturities(i,j,:),formatOut);
    end;
end;
maturities = maturitiesnew;

%size(maturities,1) == size(dates,2)
for i = 1 : size(maturities,1)
    for j = 1 : size(maturities,2)    
        dat = dates(i,:);
        rows = strcmp(allData.Date, dat);
        vars = {'Date','Close'};
        foundRows = allData(rows, vars); % Filtered table of values for given quoted date and instrument, option type
        closeAtDate = foundRows.Close;
       
        mat = maturities(i,j,:);
        rows = strcmp(allData.Date, mat);
        vars = {'Date','Close'};
        foundRows = allData(rows, vars); % Filtered table of values for given quoted date and instrument, option type
        closeAtMat = foundRows.Close;
        
        vol(i,j) = abs(closeAtMat - closeAtDate)/closeAtDate;
    end;
end;

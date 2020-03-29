function [data, K, T_numeric, T] = parseValidationFile(validationFile, validationDate)
    allData = readtable(validationFile);
    rows = strcmp(allData.validationDate, validationDate);
    
    K = table2array(allData(rows, 'K'));
    data = table2array(allData(rows, 'LocalVol'));
    T_numeric = table2array(unique(allData(rows, 'T')));
    T = table2array(unique(allData(rows, 'validatedMaturity')));
end


function [ T, Vmarket, VolImp ] = additionalFormatting( T,  Vmarket, VolImp, S, K, type )
%ADDITIONALFORMATTING Summary of this function goes here
%   Detailed explanation goes here

T = [zeros(1,1); T];

if (strcmp(type, 'C'))
    for j = 1:size(K)
        newRow(j) = max(S-K(j), 0);
        newRowVol(j) = 0;
    end;
    Vmarket = [newRow; Vmarket];
    VolImp = [newRowVol; VolImp];
end;



end


function [ triagonal ] = TriagonalMatrix( alpha, beta, gamma, timeIndex )

    sz = size(alpha,2)-2;

        
    triagonal(1, 1) = beta(timeIndex, 2);
	triagonal(1, 2) = gamma(timeIndex, 2); 
	triagonal(sz, sz) = beta(timeIndex, sz+1); %% 59 , 59
	triagonal(sz, sz - 1) = alpha(timeIndex, sz+1); %% 59 , 58

		for i = 2:sz-1
			triagonal(i, i - 1) = alpha(timeIndex, i + 1); 
			triagonal(i, i) = beta(timeIndex, i + 1);
			triagonal(i, i + 1) = gamma(timeIndex, i + 1);
        end;

		return;

end


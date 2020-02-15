function [dif, difT, dpret] = sumOfSqrDif_(u,prices, S, ptsToEvalK, epsilon, concentration_weights)
sz1 = size(u,2);
sz2 = size(u,3);
dpret = 0;
sumall = 0;
for i = 2:sz1
    sum(i) = 0;
    for j = 1:sz2
        x = (ptsToEvalK(j)-S)/ptsToEvalK(end);
        lambda = concentration_weights; %defines how concentrated weights are around S
        weight(i,j) = (1/(lambda*sqrt(2*pi)))*exp(-0.5*(x/lambda)^2);
        sum(i) = sum(i) + weight(i,j);
    end;
    sumall = sumall + sum(i);
end;

% weight(:,:) = 0;
% weight(2,27) = 5;
% weight(3,27) = 4;
% weight(4,27) = 3;
% weight(5,27) = 2;
% sumall = 14;

weight(:,:) = 0;
sumall = 10;
for ii = 2:11
    weight(ii, 27) = 1;
end;

dif = 0;
for i = 2:sz1
    difT(i) = 0;
    for j = 1:sz2
        %dif = dif + abs((prices(i,j)-u(1,i,j)))*weight(i,j)/sum;
        %dif = dif + (prices(i,j)-u(1,i,j))^2*weight(i,j)/sum(i);
        
        %% Good measure would be sum of abs differences 
        %dif = dif + (((prices(i,j)-u(1,i,j)))/prices(i,j))^2*weight(i,j)/sumall;
        dp(i,j) = abs((prices(i,j)-u(1,i,j))/prices(i,j));
        dif = dif + dp(i,j)*weight(i,j)/sumall;
        difT(i) = difT(i) + dp(i,j)*weight(i,j)/sum(i);
    end;
end;

 if (dif < 0.2) %% return matrix of abs differences only on last iter
     difftest = dp;
     mtest = squeeze(u(1,:,:));
 end
 if (dif < epsilon) %% return matrix of abs differences only on last iter
     dpret = dp;
 end

end


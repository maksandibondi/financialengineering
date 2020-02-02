function [dif] = sumOfSqrDif_(u,prices, S, K)
sz1 = size(u,2);
sz2 = size(u,3);

sumall = 0;
for i = 1:sz1
    sum(i) = 0;
    for j = 1:sz2
        x = (K(j)-S)/K(end);
        lambda = 0.05; %defines how concentrated weights are around S
        weight(i,j) = (1/(lambda*sqrt(2*pi)))*exp(-0.5*(x/lambda)^2);
        sum(i) = sum(i) + weight(i,j);
    end;
    sumall = sumall + sum(i);
end;

dif = 0;
for i = 1:sz1
    for j = 1:sz2
        %dif = dif + abs((prices(i,j)-u(1,i,j)))*weight(i,j)/sum;
        %dif = dif + (prices(i,j)-u(1,i,j))^2*weight(i,j)/sum(i);
        
        %% Good measure would be sum of abs differences 
        dif = dif + sqrt((((prices(i,j)-u(1,i,j)))/prices(i,j))^2)*weight(i,j)/sumall;
        
    end;
end;
 %dif = sqrt(dif);

end


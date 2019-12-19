function [dif] = sumOfSqrDif_(u,prices)
sz1 = size(u,2);
sz2 = size(u,3);

% sum = 0;
for i = 1:sz1
    sum(i) = 0;
    for j = 1:sz2
        x = (j-sz2/2)/sz2;
        weight(i,j) = (1/(0.1*sqrt(2*pi)))*exp(-0.5*(x/0.1)^2);
        % sum = sum + weight(i,j);
        sum(i) = sum(i) + weight(i,j);
    end;
end;

dif = 0;
for i = 1:sz1
    for j = 1:sz2
        dif = dif + (prices(i,j)-u(1,i,j))^2*weight(i,j)/sum(i);
        %dif = dif + abs((prices(1,i,j)-u(1,i,j)))*weight(i,j)/sum;
    end;
end;

end


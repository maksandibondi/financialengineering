function [dif] = sumOfSqrDif(u,prices)
sz1 = size(u,1);
sz2 = size(u,2);

dif = 0;
for i = 1:sz1
    for j = 1:sz2
        dif = dif + (prices(i,j)-u(i,j))^2;
    end;
end;

end


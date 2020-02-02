clear; clc;
% Weight function such that sum weights = 1; 
sz1 = 5;
sz2 = 20;
sum = 0;
for i = 1:sz1
    for j = 1:sz2
        x = (j-sz2/2)/(sz2);
        weight(i,j) = (1/(0.1*sqrt(2*pi)))*exp(-0.5*(x/0.1)^2);
        sum = sum + weight(i,j);
    end;
end;

scatter(1:1:20,weight(1,:)/sum);
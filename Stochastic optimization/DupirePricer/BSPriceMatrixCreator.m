function [prices] = BSPriceMatrixCreator(S,K,r,sigma,T)
sz1 = size(T,2);
sz2 = size(K,2);

for i = 1:sz1
    for j = 1:sz2
       d1 = (1 / (sigma(i,j) * sqrt(T(i))))*(log(S / K(j)) + (r + (sigma(i,j)^ 2) / 2)*T(i));
       d2 = d1 - sigma(i,j) * sqrt(T(i));
       price = normcdf(d1)*S - normcdf(d2)*K(j)*exp(-r*T(i));
       prices(i,j) = price;
    end;
end;
 
return;


end


S0 = 5430.3; 
t = 0;
T = 4/12;
r = 0.05;
K = [5125 5225 5325 5425 5525 5625 5725 5825];
M = [475 405 340 280.5 226 179.5 139 105];
precision = 0.0001;
sz = size(K,2);

for j = 1:sz

eps = 1;
sigma = sqrt(2*abs((log(S0/K(j))+r*T)/T));


while (abs(eps)>= precision)

eps = BSTheory(sigma,S0,r,T,K(j))-M(j);
sigma = sigma - eps/Vega(sigma,S0,r,T,K(j));
    
end;

sigma_imp(j) = sigma;
display(sigma_imp(j));

end;

plot(K,sigma_imp);



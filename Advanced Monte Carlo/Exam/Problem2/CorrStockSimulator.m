function [path1,path2] = CorrStockSimulator(T,N,S10, S20, r,sigma1, sigma2, ro)

delta_t = T/N;
path1(1) = S10;
path2(1) = S20;


W1 = BMsimulator(T,N,0);
W2 = BMsimulator(T,N,0);

bmpath1 = W1;
bmpath2 = ro*W1+sqrt(1-ro^2)*W2;

for i = 2:N
    
path1(i) = path1(i-1) + path1(i-1)*(r*delta_t+sigma1*(bmpath1(i)-bmpath1(i-1)));
path2(i) = path2(i-1) + path2(i-1)*(r*delta_t+sigma2*(bmpath2(i)-bmpath2(i-1)));

end;
 
return;
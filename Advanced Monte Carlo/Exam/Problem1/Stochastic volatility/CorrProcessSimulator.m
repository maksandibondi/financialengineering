function [path1,path2] = CorrProcessSimulator(T,N,S0,r,sigma0,ro)

delta_t = T/N;
path1(1) = sigma0;
path2(1) = S0;


W1 = BMsimulator(T,N,0);
W2 = BMsimulator(T,N,0);

bmpath1 = W1;
bmpath2 = ro*W1+sqrt(1-ro^2)*W2;

for i = 2:N
  
% volatility path    
path1(i) = path1(i-1) + (1/2)*sqrt(abs(path1(i-1)))*(bmpath2(i)-bmpath2(i-1));
path2(i) = path2(i-1) + path2(i-1)*((r*delta_t+path1(i-1)*(bmpath1(i)-bmpath1(i-1))));

end;
 
return;
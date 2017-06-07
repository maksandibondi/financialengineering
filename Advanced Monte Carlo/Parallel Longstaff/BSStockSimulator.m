function[S] = BSStockSimulator(S0, mu, sigma, T, delta_t)

 
 discretization_num_t = T/delta_t; 
 
 
 t(1) = 0; 
 for i = 2:discretization_num_t 
     t(i) = t(i-1)+delta_t; 
 end; 
 
 
 W = BMsimulator(T,discretization_num_t,'Reject'); 
 S(1) = S0; 
 for i = 2:discretization_num_t 
     S(i) = S(i-1) + S(i-1)*((mu-sigma^2/2)*delta_t+sigma*(W(i)-W(i-1))); 
 end; 
 
return;

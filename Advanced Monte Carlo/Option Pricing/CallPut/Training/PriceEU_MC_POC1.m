%% Main
for i = 1:41
    S0(i) = 0.5*(i-1);
    [priceC(i),priceP(i)] = PriceEU_MC(S0(i),10,0.01,0.3,0.5);
end;

%% Initial condition
plot(S0,priceC);
hold on;
plot(S0,priceP);

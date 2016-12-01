[priceC1,priceP1,S1,t1] = Price_EU_creator(0,20,0,0.5,0.1,0.3,10);

[priceC2,priceP2,S2,t2] = Price_EU_creator(0,20,0,0.25,0.1,0.3,10);

[priceC3,priceP3,S3,t3] = Price_EU_creator(0,20,0,0,0.1,0.3,10);

plot(S1,priceC1(end,:));
hold on;
plot(S2,priceC2(end,:));
plot(S3,priceC3(end,:));
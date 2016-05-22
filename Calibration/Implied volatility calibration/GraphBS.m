function[V]=GraphBS(sigma,r,T,K)
    for j=1:20
        asset_price(j)=j;
        V(j)=blackscholestheory(sigma,asset_price(j),r,T,K);
    end
    plot(asset_price,V);
end
function [out] = payoff_butterfly(S,K)


if (S>=K) && (S<=2*K)
    out = S-K;
elseif (S>=2*K) && (S<=3*K)
    out = 3*K-S;
else
    out = 0;
end;

return;
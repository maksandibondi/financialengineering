function[out] = payoff(S,K)

if S<=K
    out = S-K;
elseif S>K && S<=2*K
    out = 0;
elseif S>2*K
    out = S-2*K;
end;

return;
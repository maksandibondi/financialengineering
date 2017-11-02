function [out] = payoff_straddle(S,K)

if S<=K
    out = K-S;
else
    out = S-K;
end;

return;
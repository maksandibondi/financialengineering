function [payoff] = payoff_EU(S,K,type)
if (strcmp(type,'Put'))
payoff = max(K-S,0);
elseif (strcmp(type,'Call'))
payoff = max(S-K,0);
end;

return; 
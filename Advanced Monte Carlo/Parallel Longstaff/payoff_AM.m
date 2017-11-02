function [payoff] = payoff_AM(S,K,type)
if (strcmp(type,'Put'))
payoff = K-S;
else 
payoff = S-K;
end;

return; 
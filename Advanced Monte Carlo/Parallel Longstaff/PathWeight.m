function [weight] = PathWeight(K,S,type)
if (strcmp(type,'Put'))
  if (K-S)>0
    weight = 1;
  else 
   weight = 0;
  end;
  
else 
if (S-K)>0
    weight = 1;
  else 
   weight = 0;
  end;
  
end;

return; 
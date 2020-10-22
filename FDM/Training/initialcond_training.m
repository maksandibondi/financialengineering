function[out] = initialcond_training(x)


out = x*sin(pi*x); % training 2
% out = sin(2*pi*x); % training 1


% if (x<=0)
%     
%     out = 0;
%     
% elseif (0<=x && x<=1/2)
%     
%     out = 2*x;
%     
% elseif (1/2<=x && x<=1)
%     
%     out = -2*x+2;
%     
% else
%     out = 0;
% end;
    
return
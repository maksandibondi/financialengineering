function [out] =  initialcond_bs(x,strike)
out =  max(strike-x,0); % for call neumann and dirichlet
return

% max(x-strike,0); % for call neumann and dirichlet
% max(strike-x,0); % for put neumann
%                ; % for butterfly
% 
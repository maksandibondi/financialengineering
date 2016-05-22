function [out] =  boundaryfinal_bs(t,strike, rate, t_final, x_final)
out = 0; % call dirichlet;
return

% x_final - strike*exp(-rate*(t_final-t))% call dirichlet;
% 1; % call neumann;
% 0; % put neumann
%                ; % for butterfly

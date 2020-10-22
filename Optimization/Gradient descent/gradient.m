function[out] = gradient(theta, x, y) % costfunction for multivariate linear regression

h = ((theta')*x')'; % hypothesis linear

% h = ones(size(x,1))/(ones(size(theta,2))'+exp(-(theta')*x')); % sigmoid function hypothesis for nonlinear logistic
% h = ((theta')*x')'; % hypothesis linear
% h = ((theta')*(x.^2)')'; % hypothesis non-linear squared
% h = ((theta')*(x.^3)')'; % hypothesis non-linear cubic
% h = ((theta')*(x.^4)')'; % hypothesis non-linear care

m = size(x,1); % size of learning set

    out = (1/m)*(h-y)'*x; % derivative of cost function 1/m(h-y)^2 linear
    
    % out = (1/m)*(h-y)'*x; % derivative of cost function 1/m(h-y)^2 linear
    % out = -(1/m)*(y*log(h)+(1-y)*log(1-h); % for nonlinear logistic
    % out = (1/m)*(h-y)'*(x.^2); % for non-linera hypothesis squared
    % out = (1/m)*(h-y)'*(x.^3); % for non-linera hypothesis cubic
    % out = (1/m)*(h-y)'*(x.^4); % for non-linera hypothesis care
    
    
return

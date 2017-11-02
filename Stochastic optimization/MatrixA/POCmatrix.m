clear; clc;
T = [0:0.125:0.5];
logK = [-5.29:((5.29+5.29)/53):5.29]; % axis [-log(k); log(k)]
knots = [-5.29:((5.29+5.29)/5):5.29]; % interpolating knots (each form M-4 knots has corresponding control point theta)
discretization_num_K = 54;


% Find x - array of non-uniform distanced values of logK
sz = size(logK,2)
delta_log_K = (2*(logK(end))) / (discretization_num_K-1)
x(1) = -(logK(end));

for i = 1:sz
	x(i) = 0.5*log(((logK(end)) - (logK(end))+i*delta_log_K) / ((logK(end)) - (-(logK(end)) + i*delta_log_K))) + x(1);
end;


x;
scatter(x,zeros(54,1));
hold on;
scatter(knots, zeros(6,1));

A = FillMatrixA(T,logK,knots)
%B = chol(A,'lower')

% now we have to take control points theta_new = theta + B*normal() and get new
% values

% now new logK presented as cubic spline using theta_new (POC.m)
function[out] = velocity(x, x_0, A, sigma, Q)

% out = 0;
out = A*exp(-((x-x_0)^2)/(2*sigma^2))*(((x-x_0)^2)/(sigma^2)*cos(Q*(x-x_0))*cos(Q/2)+2*sin(Q*(x-x_0))*sin(Q/2));


return


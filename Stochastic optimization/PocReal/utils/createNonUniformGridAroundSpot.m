function [grid, y] = createNonUniformGridAroundSpot(Kmax, Kmin, discretization_num_K, S, method, finness)
    if nargin < 6
    finness = 10;
    end
    y = 0;
    
if (strcmp(method,'log'))
    %%  non-uniform
    ymax = log(Kmin);
    y_0 = log(S);
    h = 2*ymax/(discretization_num_K+1);
    for i = 1:discretization_num_K+1
         y(i) = gfunc(-ymax+i*h, ymax, y_0); 
    end;
    grid = exp(y(1:end-1));
    return;
end;

if (strcmp(method, 'user'))
    %%  non-uniform
    grid = getNonUniformKnots( Kmin, Kmax, S, discretization_num_K);
    y = log(grid);
    return;
end;
    
if (strcmp(method,'log2'))
    %%  non-uniform
    ymax = log(Kmax);
    ymin = -ymax;
    center = log(S); % point where we want concentrated values
    h_uniform = (ymax-ymin) / (discretization_num_K);

    for i = 1 : discretization_num_K
        y(i) = nonUniTransform(ymin+i*h_uniform, ymax, center);
        %h(i-1) = y(i)-y(i-1);
    end;
    grid = exp(y);
    return;
end;

if (strcmp(method,'sin'))
    %%  non-uniform
    alpha = finness;
    c1 = (sinh((Kmin - S)/alpha))^(-1);
    c2 = (sinh((Kmax - S)/alpha))^(-1);
    for i = 1:discretization_num_K
        grid(i) = S + alpha*sinh(c2*i/discretization_num_K + c1*(1-i/discretization_num_K));
    end;
    grid(1) = Kmin;
    grid(discretization_num_K) = Kmax;
    return;
end;

if (strcmp(method,'gauss'))
    %%  non-uniform
    sigma = finness;
    mean = (S-Kmin)/(Kmax-Kmin);
    Xnonsorted = normrnd(mean, sigma, [discretization_num_K*10,1]);
    X = sort(Xnonsorted);
    grid(1) = Kmin;
    grid(discretization_num_K) = Kmax;
    centralElements = X((end-discretization_num_K)/2+1:(end+discretization_num_K)/2);
    
    for i = 2:discretization_num_K-1
        grid(i) = grid(1) + centralElements(i)*(Kmax-Kmin);
    end;
    return;
end;

if (strcmp(method,'quadratic'))
    %%  non-uniform
    point = (S-Kmin)/(Kmax-Kmin);
    numOfPointsAtLeft = ceil(discretization_num_K*point);
    numOfPointsAtRight = discretization_num_K - numOfPointsAtLeft;
   
    dh = (Kmax-Kmin)/(2*discretization_num_K);
    h = -Kmax:dh:Kmax;
    
    a = 0.5;
    y = a*(h.^2)+S;
    expected = Kmax-Kmin;
    while (abs(expected-sum(y(numOfPointsAtRight:numOfPointsAtRight+discretization_num_K)))>0.05)
           a_old = a;
           z_old = sum(y(numOfPointsAtRight:numOfPointsAtRight+discretization_num_K));
           a = a - (expected-sum(y(numOfPointsAtRight:numOfPointsAtRight+discretization_num_K)));
           y = a*(h.^2)+S;
    end
    
    
    
    
    for i = 2:discretization_num_K-1
       % grid(i) = grid(1) + centralElements(i)*(Kmax-Kmin);
    end;
    return;
end;

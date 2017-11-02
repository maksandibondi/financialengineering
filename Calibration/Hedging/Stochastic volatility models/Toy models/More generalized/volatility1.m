function [volat] = volatility1()

    u = rand();
    
    if u>=0.99
        volat = 0.5;
    else
        volat = 0.3;
    end;

return;


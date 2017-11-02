function [volat] = volatility2(volatility_known)

    if volatility_known ~= 0.3
        volat_new = 0.3;
    else
        volat_new = 0.5;
    end;
    
    u = rand();
    
    if u>=0.95
        volat = volat_new;
    else
        volat = volatility_known;
    end;

return;
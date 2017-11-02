function [R] = Ornstein_ITOSimulator(R0,a,b,sigma,T,delta_t)

discretization_num_t = T/delta_t;

t(1) = 0;
for i = 2:discretization_num_t
    t(i) = t(i-1)+delta_t;
end;

W = BMsimulator(T,discretization_num_t,'Polar');
R(1) = R0;
for i = 2:discretization_num_t
    
    for s = 2:i
        integral(s) = exp(a*(t(s)-t(i)))*(W(s)-W(s-1));
    end;
    R(i) = R(1)*exp(-a*t(i))+b*(1-exp(-a*t(i)))+sigma*sum(integral);
end;

return;
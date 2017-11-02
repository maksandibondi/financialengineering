function[Amp] = Amplitude_discrete(jump_freq,t,v1,p1,v2,p2)

Amp = 1;

q = Jumps_t(jump_freq,t);

for i = 1:q
    
    Amp = Amp*(1+discreteY(v1,p1,v2,p2));
    
end;

return
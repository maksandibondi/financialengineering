function[Amp] = Amplitude(m_y,sigma_j,jump_freq,t)

Amp = 1;

q = Jumps_t(jump_freq,t);

for i = 1:q
    
    Amp = Amp*(1+m_y)*exp(-sigma_j^2/2+sigma_j*randn());
    
end;

return
function [path] = BMsimulator(T,N,initial_point)


delta_t = T/N;
t_vector(1) = 0;

for k = 2:N
    t_vector(k) = t_vector(k-1)+delta_t;
end;

path(1) = initial_point;
for n = 2:N
    path(n) = path(n-1)+sqrt(delta_t)*randn();
end;


return;
function [A] = FillMatrixA(T,x,knots)


Nknots = size(knots,2)
M = Nknots-4
N = size(T,2)
sz = N*M;
A = zeros(sz,sz);

% filling terms with fm norm
for i = 1:N-1
	for m = 3:Nknots-2
        f = NormL2(knots, m, x); % squared integral (sum of fm for all K)
		A((i-1)*M + m - 2, (i-1)*M + m - 2) = A((i-1)*M + m - 2, (i-1)*M + m - 2)+ f*T(i)/((T(i + 1) - T(i))^2);
		A((i)*M + m - 2, (i)*M + m - 2) = A((i)*M + m - 2, (i)*M + m - 2) + f*T(i)/ ((T(i + 1) - T(i))^2);
		A((i-1)*M + m - 2, (i)*M + m - 2) = A((i-1)*M + m - 2, (i)*M + m - 2) -2 * f*T(i) / ((T(i + 1) - T(i))^2);
    end;
end;


% filling terms with second derivative of fm norm
for i = 1:N-1
	for m = 3:Nknots-2
        f = NormDerL2(knots, m, x); % squared integral (sum of fm for all K)
		A((i-1)*M + m - 2, (i-1)*M + m - 2) = A((i-1)*M + m - 2, (i-1)*M + m - 2)+ f*((T(i)^3)/3-T(i+1)*T(i)^2+T(i+1)^2*T(i))/((T(i + 1) - T(i))^2);
		A((i)*M + m - 2, (i)*M + m - 2) = A((i)*M + m - 2, (i)*M + m - 2) + f*((T(i)^3)/3-T(i+1)*T(i)^2)/ ((T(i + 1) - T(i))^2);
		A((i-1)*M + m - 2, (i)*M + m - 2) = A((i-1)*M + m - 2, (i)*M + m - 2) +2 * f*((T(i)^3)/3)/ ((T(i + 1) - T(i))^2);
    end;
end;


% filling terms with fm , fm-1 scalar
for  i = 1:N-1
    for m = 4:Nknots-3 

				f = scalarBasisCross1(knots, m, x);
				A((i-1)*M + m - 2, (i-1)*M + m - 2 - 1) = A((i-1)*M + m - 2, (i-1)*M + m - 2 - 1)+ f*T(i) / (T(i + 1) - T(i))^2;
				A((i)*M + m - 2, (i)*M + m - 2 - 1) = A((i)*M + m - 2, (i)*M + m - 2 - 1) + f*T(i) /(T(i + 1) - T(i))^2;
                
                A((i-1)*M + m - 2, (i)*M + m - 2 - 1) = A((i-1)*M + m - 2, (i)*M + m - 2 - 1) + f*T(i) / (T(i + 1) - T(i))^2;
                A((i-1)*M + m - 2 - 1, (i)*M + m - 2) = A((i-1)*M + m - 2 - 1, (i)*M + m - 2) + f*T(i) / (T(i + 1) - T(i))^2;
                A((i-1)*M + m - 2, (i)*M + m - 2 - 1)
                A((i-1)*M + m - 2 - 1, (i)*M + m - 2)
    end;     
end;


% filling terms with fm , fm+1 scalar
for  i = 1:N-1

    for m = 4:Nknots-3 
				
                f = scalarBasisCross2(knots, m, x);
				A((i-1)*M + m - 2, (i-1)*M + m - 2 + 1) = A((i-1)*M + m - 2, (i-1)*M + m - 2 + 1) + f*T(i) / (T(i + 1) - T(i))^2;
				A((i)*M + m - 2, (i)*M + m - 2 + 1) = A((i)*M + m - 2, (i)*M + m - 2 + 1) + f*T(i) /(T(i + 1) - T(i))^2;
                
                A((i-1)*M + m - 2 + 1, (i)*M + m - 2) = A((i-1)*M + m - 2 + 1, (i)*M + m - 2) + f*T(i) / (T(i + 1) - T(i))^2
				A((i-1)*M + m - 2, (i)*M + m - 2 + 1) = A((i-1)*M + m - 2, (i)*M + m - 2 + 1) + f*T(i) / (T(i + 1) - T(i))^2
    end;
            
end;



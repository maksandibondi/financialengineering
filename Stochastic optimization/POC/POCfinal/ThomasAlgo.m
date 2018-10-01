function [res] = ThomasAlgo(A, B) % ctrl y - values from y1 to yn+1
n = size(A,1);
  
c_(1) = A(1, 2) / A(1, 1);
d_(1) = B(1) / A(1, 1);

		% forward walk: looking for transformed A matrix members
		for i = 2:n-1 
			c_(i) = A(i, i + 1) / (A(i, i) - c_(i - 1) * A(i, i - 1));
			d_(i) = (B(i) - d_(i - 1) * A(i, i - 1)) / (A(i, i) - c_(i - 1) * A(i, i - 1));
        end;
    
    
	d_(n) = (B(n) - d_(n - 1) * A(n, n - 1 )) / (A(n, n) - c_(n-1) * A(n, n - 1));
		
		% backward walk: finding the solutions
		res(n) = d_(n);

		for i = n:-1:2 
			res(i - 1) = d_(i - 1) - res(i) * c_(i - 1); % c is res
        end;
    
    
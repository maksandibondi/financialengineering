function [c] = ThomasAlgo(ctrlY, ctrlX, h) % ctrl y - values from y1 to yn+1

n = size(ctrlX,2); 


B = zeros(n-2,1);
for i = 1:n-2
  B(i) = (ctrlY(i+2)-ctrlY(i+1))/h(i+1) -  (ctrlY(i+1)-ctrlY(i))/h(i);
end;


A = zeros(n-2, n-2);
A(1,1) = 2*(h(1)+h(2));
A(1,2) = h(2);
A(n-2,n-3) = 2*(h(n-2)+h(n-1));
A(n-2,n-2) = h(n-1);

for i = 2:n-3
    A(i,i-1) = h(i);
    A(i,i) = 2*(h(i)+h(i+1));
    A(i,i+1) = h(i+1);
end;
  
  A
  
c_(1) = A(1, 2) / A(1, 1);
d_(1) = B(1) / A(1, 1);

		% forward walk: looking for transformed A matrix members
		for i = 2:n-3 
			c_(i) = A(i, i + 1) / (A(i, i) - c_(i - 1) * A(i, i - 1));
			d_(i) = (B(i) - d_(i - 1) * A(i, i - 1)) / (A(i, i) - c_(i - 1) * A(i, i - 1));
        end;
    
    
	d_(n-2) = (B(n-2) - d_(n - 3) * A(n - 2, n - 3)) / (A(n-2, n-2) - c_(n-3) * A(n-2, n - 3));
		
		% backward walk: finding the solutions
		c(n-2) = d_(n-2);

		for i = n-2:-1:2 
			c(i - 1) = d_(i - 1) - c(i) * c_(i - 1); % c is res
        end;
    
    
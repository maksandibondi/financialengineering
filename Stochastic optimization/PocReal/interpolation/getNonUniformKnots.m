function [ knots ] = getNonUniformKnots( X_0, X_last, concentrationPoint, numOfPoints )
concentration_force = 0.4;

% gets interpolation knots concentrated around S
knots(1) = X_0;
knots(numOfPoints) = X_last;
S = concentrationPoint;

if (floor(numOfPoints/2) == (numOfPoints/2))
    isPair = 1;
else
    isPair = 0;
    knots(floor(numOfPoints/2) + 1) = S;
end


for j = 2:floor(numOfPoints/2)
    knots(j) = S - (1/(concentration_force*(j^2))) * (S - X_0);
end;

for j = (floor(numOfPoints/2) + 1 + 1 * (1 - isPair)) : numOfPoints-1
    knots(j) = S + (1/(concentration_force*((numOfPoints+1-j)^2))) * (X_last - S);
end;

end


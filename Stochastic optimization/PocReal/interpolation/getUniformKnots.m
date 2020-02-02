function [ knots ] = getUniformKnots( X_0, X_last, numOfPoints )
% gets interpolation knots concentrated around S
knots(1) = X_0;
knots(numOfPoints) = X_last;

delta_X = (X_last-X_0)/numOfPoints;


for j = 2:numOfPoints-1
    knots(j) = knots(j-1)+delta_X;
end;

end


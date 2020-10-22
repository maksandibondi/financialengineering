function[res] = Bspline(knots,ctrlpts,order,ptsToEval)

sz = size(knots,2);
p = order;
sz2 = size(ptsToEval,2);

iter = 1;
for i = 1:sz2 % iterate through points where interpolate
  spl = 0;
  for m = 3:sz-2
    cspline = BSplineCoefUniform(knots,ptsToEval(i),m);
    spl = spl + cspline*ctrlpts(m);
  end;
  res(iter) = spl;
  iter=iter+1;
end;



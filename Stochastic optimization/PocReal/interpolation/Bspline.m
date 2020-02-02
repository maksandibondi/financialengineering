function[res] = Bspline(knots,ctrlpts,order,ptsToEval)

sz = size(knots,2);
p = order;
sz2 = size(ptsToEval,2);


%% Uniform knots
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


% %% Non-uniform knots
% iter = 1;
% for i = 1:sz2 % iterate through points where interpolate (ptsToEval)
%   spl = 0;
%   for m = 3:sz-order-1 % iterate through knots
%     cspline = BSplineCoefDeBoor(knots,m,order,ptsToEval(i));
%     cs(iter,m) = cspline;
%     spl = spl + cspline*ctrlpts(m);
%   end;
%   res(iter) = spl;
%   iter=iter+1;
% end;
% 
% 
% i = 1;
% for j = 1:sz2
%     if (ptsToEval(j)<knots(order))
%         continue;
%     elseif (ptsToEval(j)<knots(end-order-1))
%         newPtsToEval(i) = ptsToEval(j);
%         i = i+1;
%     else
%         break;
%     end;
% end;
% newPtsToEval
% knots
% ctrlpts
% res = BSplineCoefDeBoorEx(p,knots,ctrlpts(3:end-2),newPtsToEval);

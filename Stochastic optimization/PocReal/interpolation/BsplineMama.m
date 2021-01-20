function[res] = BsplineMama(knots,ctrlpts, ptsToEvalK)

for k = 1:size(ptsToEvalK,2)
    
    for i = 2:size(knots,2)
        if (ptsToEvalK(k)==knots(i-1))
            res(k) = ctrlpts(i-1);
            break;
        end;

        if (ptsToEvalK(k)==knots(i))
            res(k) = ctrlpts(i);
            break;
        end;

        if (ptsToEvalK(k) > knots(i-1) && ptsToEvalK(k) < knots(i))
            a = ctrlpts(i-1);
            b = diff(ctrlpts(i-1:i))./diff(knots(i-1:i));
            %c = diff(b)./diff(x(2:length(x)));
            %d = diff(c)./diff(x(3:length(x)));

            res(k) = a + b*(ptsToEvalK(k)-knots(i-1));% + c + d;
            break;
        else    
            continue;
        end;
        res(k) = 0;
    end;
end;


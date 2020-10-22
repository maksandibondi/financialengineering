function [i] = ChoosenIdx( rand, distmass )

    cumul = distmass(1);
    for i = 2:size(distmass,2)
        if (rand<cumul)
            return;
        else 
            cumul = cumul + distmass(i);
        end;
    end;

end


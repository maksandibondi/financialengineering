function [alpha, beta, gamma] = FillAlphaBetaGamma(sigma, r, K, discretization_num_K, S, discretizationType,ptsToEvalK)

d1 = size(sigma,1); % discretization_num_T
d2 = size(sigma,2); % discretization_num_K

    if (strcmp(discretizationType,'uniform'))
        % uniform grid on -ymax:ymax
        y(1) = log(K(1));
        ymax = log(K(end));
        h = (ymax-y(1))/discretization_num_K;

        for i = 2:d2
            y(i) = y(i-1)+i*h;
        end;

        % find alpha, beta, gamma
        I = ones(d1,d2);

        alpha = -((sigma) .^ 2) / (2 * h ^ 2) - (1 / (2 * h)) * (r + ((sigma) .^ 2) / 2);

        beta = ((sigma) .^ 2) / (h ^ 2);

        gamma = -((sigma) .^ 2) / (2 * h ^ 2) + (1 / (2 * h))*(r + ((sigma) .^ 2) / 2);

        return;
    
    else

            delta_log_K = (log(K(end)) - log(K(1))) / (discretization_num_K); %+2);
           
            
            y = zeros(1,d2);
            y(1) = log(S);   
            for i = 2:d2+1
                y(i) = 0.5*log((log(K(end)) + log(K(1))+i*delta_log_K) / (log(K(end)) - (log(K(1)) + i*delta_log_K))) + y(1); 
            end;
            %m = log(K(1))+20*delta_log_K
            %log(K(end))

            % matrix h_i
            h = zeros(d1,d2);
            for n = 1:d1
                for j = 1:d2-1
                    h(n, j) = y(j+1) - y(j);
                end;
            end;

            % matrix h_(i-1)
            h_ = zeros(d1,d2);
            for n = 1:d1
                for j = 2:d2
                    h_(n, j) = h(n,j-1);
                end;
            end;

            % matrix h_i + h_(i-1)
            %real(h)
            %real(h_)
            h_h = h + h_;
            %real(h_h) 


        % find alpha, beta, gamma
        I = ones(d1,d2);

        alpha = ((sigma) .^ 2).*(I * (-1) ./ (h_h .* h_) - I ./ (I * 4 .* h_)) - I*r ./ (I*2 .* h_);

        beta = ((sigma) .^ 2) .* (I ./ (h_h.*h) + I ./ (h_h.*h_) + I*0.25.*(I ./ h_ - I ./ h)) + I*0.5*r.*(I ./ h_ - I ./ h);

        gamma = ((sigma) .^ 2).*(I * (-1) ./ (h_h .* h) + I ./ (I * 4 .* h)) + I*r ./ (I * 2 .* h);

        return;

    end

end





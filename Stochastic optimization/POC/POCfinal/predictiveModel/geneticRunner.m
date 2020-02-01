function [] = geneticRunner(K, T, S, r, Vmarket, VolImp, discretizationType)
% Troubles:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2) The same non-uniform discetization to be used for market data
% simulation (?), PDE resolution (function invoked) !!!!!!!!!
%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3) spline creation (knots non-uniform) using Working De Boor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initial data

K_0 = K(1,1); K_l = K(end,1);
T_0 = T(1,1); T_l = T(end,1);
discretization_num_K = size(K,1);
discretization_num_T = size(T,1);

T0 = 1000;
t = 0.07;
Nmc = 2^14;
M = 64;
popsize = 20;


if (strcmp(discretizationType, 'uniform'))
    ptsToEvalK = transp(K(:,1));
else 
    %% Nonuniform discretization grid(the same as in the pricer)
    ymax = log(K_l);
    ymin = -ymax;
    y_0 = log(S);
    h = 2*ymax/(discretization_num_K+1);
    for i = 1:discretization_num_K+1
         y(i) = gfunc(-ymax+i*h, ymax, y_0); 
    end;
    ptsToEvalK = exp(y(1:end-1));
    %% Non uniform discretiztion
    % creates non-uniform grid including points K_0, K_l. other values will be
    % interpolated later
    %ptsToEvalK = createNonUniformGridAroundSpot(K_l, K_0, discretization_num_K, S, 'log', 6);
    %ptsToEvalK = getNonUniformKnots(K_0, K_l, S, discretization_num_K);
end

%% Getting non uniform market prices
if (~strcmp(discretizationType, 'uniform'))
    knotsK = transp(K(:,1));%getNonUniformKnots(0, K_l*10, S, discretization_num_K);
    for i = 1:discretization_num_T
        Vmarket_temp(i,:) = interp1(knotsK, Vmarket(i,:), ptsToEvalK, 'linear' ,'extrap');
    end;
    Vmarket = Vmarket_temp;
end


iter = 0;

%% Algo
for k = 1:Nmc

    for n = 1:popsize
        disc_T = 6; disc_K = 10;
        if (k == 1)
        ctrlpts = rand(disc_T,disc_K);
        %Spline evaluated at the points at which our non-uniform discretization is
        % done. How to give it the size 31*198 as well?
        else
        ctrlpts(:,:) = ctr(n,:,:);   
        end;
        % We have to get interpolated sigma by 6*21 in 10*61 points to eval
        sigmaaa = SplineLinear2DInterp(T_0,T_l,K_0,K_l,S,disc_T, disc_K, ptsToEvalK, T, ctrlpts, discretizationType);

        %% Pricing
        u(n,:,:) = Pricer_dupire(sigmaaa, K, T, discretization_num_K, discretization_num_T, S, r, discretizationType, ptsToEvalK);

        %fitness(n) = sumOfSqrDif(u(n,10:20,96:104),prices(n,10:20,96:104)); % cost
        fitness(n) = sumOfSqrDif_(u(n,:,:), Vmarket(:,:), S, K); % cost funtion for n-th member of population

        sig(n,:,:) = sigmaaa; 
        ctr(n,:,:) = ctrlpts;
    end;

    %% Genetic part
        [minfitness, index_best] = min(fitness);
        if (minfitness>7)
            T_ = T0*(1-k*(mod(k,M)==0)/Nmc)^4;

            % Selection
            for n = 1:popsize  
                if (rand() > exp(-fitness(n)/T_)) 
                     for l = 1:popsize
                           % Proba of chosing another candidate. sum of these proba == 1
                           p(l) = exp(-fitness(l)/T_)/sum(exp(-fitness/T_)); 
                     end;
                     idx = ChoosenIdx(rand(),p);
                     ctr(n,:,:) = ctr(idx,:,:);
                end;
            end;

            % Mutation ( applied independently to all ctrl points (6*21))
            for n = 1:popsize
            sz2 = size(ctr,2);
            sz3 = size(ctr,3);
            for m = 1:sz2
                for g = 1:sz3
                    init = ctr(n,m,g);
                    ctr(n,m,g) = ctr(n,m,g) + t*(2*rand()-1);
                    if (ctr(n,m,g)<0 || ctr(n,m,g)>1) % check the constraints
                        ctr(n,m,g) = init;
                    end;
                end;
            end;
            end;

        else 
            break;
        end;

        iter = iter+1
        fit = min(fitness)

end;


x1 = 1:1:6;
x2 = 1:1:21;
Z = squeeze(sig(index_best,:,:));
figure;
surf(K,T,Z);
hold on;
xlabel('K'); ylabel('T'); zlabel('sigma(K,T)');
figure;
hold on;

plot(T,sig(index_best,:,9),'g');
plot(T,sig(index_best,:,10), 'r');
plot(T,sig(index_best,:,11), 'b');
plot(T,sig(index_best,:,12), 'o');



hold on;
figure;
surf(K,T,squeeze(u(n,:,:)));
hold on;
%surf(K,T,ones(size(T,2),size(K,2))*0.6, 'FaceColor', [1 0 1]);
xlabel('K'); ylabel('T'); zlabel('u(K,T)');
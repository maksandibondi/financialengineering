function [localVolCalibrated] = geneticRunner(K, T, S, r, Vmarket, VolImp, inputStructure)

%% Initial parameters
K_0 = K(1,1); K_l = K(end,1);
T_0 = T(1,1); T_l = T(end,1);
discretization_num_K = size(K,1);
discretization_num_T = size(T,1); 

disc_T = inputStructure.disc_T;
disc_K = inputStructure.disc_K;
T0 = inputStructure.T0;
t = inputStructure.t;
Nmc = inputStructure.Nmc;
M = inputStructure.M;
popsize = inputStructure.popsize;
epsilon = inputStructure.epsilon;
concentration_weights = inputStructure.concentration_weights;
discretizationType = inputStructure.discretizationType;
interpTypeK = inputStructure.interpTypeK;
nonuniform_method = inputStructure.nonuniform_method;
includeVolAssumption1 = inputStructure.includeVolAssumption1;
includeVolAssumption2 = inputStructure.includeVolAssumption2;
includeVolAssumptionImplicit = inputStructure.includeVolAssumptionImplicit;
%rowByRowMutation = inputStructure.rowByRowMutation;

rootdir = fullfile(pwd, '..');
outputFile = strcat(rootdir, '/Results/', inputStructure.date, '_', discretizationType, '_', interpTypeK, '_', 'concentration_weights_', strrep(num2str(concentration_weights),'.',''), 'epsilon_', strrep(num2str(epsilon), '.', ''), '.xls');
outputFig = strcat(rootdir, '/Results/figures/', inputStructure.date, '_', discretizationType, '_', interpTypeK, '_', 'concentration_weights_', strrep(num2str(concentration_weights),'.',''), 'epsilon_', strrep(num2str(epsilon), '.', ''), '_placeholder', '.png');

%% Creating uniform/non-uniform grid
if (strcmp(discretizationType, 'uniform'))
    %% Uniform disc grid
    ptsToEvalK = transp(K(:,1));
else 
    %% Nonuniform discretization grid(the same as in the pricer)
    [ptsToEvalK, ~] = createNonUniformGridAroundSpot(K_l, K_0, discretization_num_K, S, nonuniform_method);
end

%% Getting non uniform market prices by interpolating uniform if needed
if (strcmp(discretizationType, 'nonuniform'))
    knotsK = transp(K(:,1));
    for i = 1:discretization_num_T
        Vmarket_temp(i,:) = interp1(knotsK, Vmarket(i,:), ptsToEvalK, 'linear' ,'extrap');
        %Vmarket_temp(i,:) = interp1(knotsK, Vmarket(i,:), ptsToEvalK, 'spline');
    end;
    Vmarket = Vmarket_temp;
end


iter = 0;

%% Algo
for k = 1:Nmc

    for n = 1:popsize
        if (k == 1)
            %ctrlpts = squeeze(VolImp(1:disc_T, 1:disc_K)); %% if we want a good seed as implied vol
            %% Seed is very significant but can be fought with mutation coefficient high
            rng('shuffle');
            ctrlpts = rand(disc_T,disc_K);
            if includeVolAssumption1 %% Assumtion that local vol is bigger for bigger maturity
                for i = 1 : disc_T
                    ctrlpts(i,:) = ctrlpts(i,:)*(i/disc_T);
                end;
            end;
            if includeVolAssumption2 %% Assumtion that local vol is bigger for lower maturity
                for i = 1 : disc_T
                    ctrlpts(i,:) = ctrlpts(i,:)*disc_T/i;
                end;
            end;
            if includeVolAssumptionImplicit
                ctrlpts = getSeedFromImplicitVol(K_0, K_l, T_0, T_l, discretization_num_T, discretization_num_K,disc_T, disc_K,VolImp, ptsToEvalK, T);
            end;
            
        else
            %% Ctrlpts are the last ctr fenerated in mutation phase
            ctrlpts(:,:) = ctr(n,:,:);   
        end;
        
        % We have to get interpolated sigma by interpolating points to eval
        % from disc_T*disc_K knots
        sigmaaa = SplineLinear2DInterp(T_0,T_l,K_0,K_l,S,disc_T, disc_K, ptsToEvalK, T, ctrlpts, discretizationType, interpTypeK);
%%      Show vol splined
%         if  n==1
%                  plot(ptsToEvalK, sigmaaa(1,:));
%                  %plot(ptsToEvalK, sigmaaa(5,:));
%                  hold on;
%         end;
        
        %% Pricing
        u(n,:,:) = Pricer_dupire(sigmaaa, ptsToEvalK, T, discretization_num_K, discretization_num_T, S, r, discretizationType, nonuniform_method);
        [fitness(n), difT, ~] = sumOfSqrDif_(u(n,:,:), Vmarket(:,:), S, ptsToEvalK, epsilon, concentration_weights); % cost funtion for n-th member of population

        sig(n,:,:) = sigmaaa; 
        ctr(n,:,:) = ctrlpts;
        
%%         if rowByRowMutation
%             [~, worstindexT(n)] = max(difT);
%         end;
    end;

    [minfitness, index_best] = min(fitness);
    %% Show vol best
%          if minfitness > 0.16
%             plot(ptsToEvalK, squeeze(sig(index_best,1, :)), 'r');
%             hold on;
%          elseif minfitness < 0.155 && minfitness > 0.1355
%             plot(ptsToEvalK, squeeze(sig(index_best,1, :)), 'b');
%             hold on;
%           elseif minfitness < 0.13 && minfitness > 0.11
%             plot(ptsToEvalK, squeeze(sig(index_best,1, :)), 'k');
%             hold on; 
%            elseif minfitness < 0.11 
%             plot(ptsToEvalK, squeeze(sig(index_best,1, :)), 'g');
%             hold on; 
%          end;
         
         %% Genetic part
         if (minfitness>epsilon)
            T_ = T0*(1-k*(mod(k,M)==0)/Nmc)^4;
           
         %% Selection
            rng('shuffle');
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
            
         %% Mutation ( applied independently to all ctrl points (6*21))
            for n = 1:popsize
                sz2 = size(ctr,2);
                sz3 = size(ctr,3);
%%                 if rowByRowMutation
%                     [~,m] = min(abs([1:size(ctr,2)] - worstindexT(n)));
%                     for g = 1:sz3
%                         init = ctr(n,m,g);
%                         ctr(n,m,g) = ctr(n,m,g) + t*(2*rand()-1);
%                         if (ctr(n,m,g)<0 || ctr(n,m,g)>1) % check the constraints
%                            ctr(n,m,g) = init;
%                         end;
%                     end; 
%%                 else 
                    for m = 1:sz2
                        for g = 1:sz3
                            init = ctr(n,m,g);
                            ctr(n,m,g) = ctr(n,m,g) + t*(2*rand()-1);
                            if (ctr(n,m,g)<0 || ctr(n,m,g)>1) % check the constraints
                               %ctr(n,m,g) = init;
                            end;
                        end;
                    end;
%                 end;
            end;

         else 
             break;
         end; 
         
         iter = iter+1;
         display(minfitness);
end;
localVolCalibrated = squeeze(sig(index_best,:,:));


%%  Visual part

%% Write axes to file
xlswrite(outputFile, ptsToEvalK, 1, 'B1');
xlswrite (outputFile, T, 1, 'A2');
xlswrite(outputFile, ptsToEvalK, 1, strcat('B', num2str(discretization_num_T+2,2)));
xlswrite (outputFile, T, 1, strcat('A', num2str(discretization_num_T+3,2)));

%% Write diff in prices into file
[~, ~, diffprice] = sumOfSqrDif_(u(index_best,:,:), Vmarket(:,:), S, ptsToEvalK, epsilon, concentration_weights);
xlswrite(outputFile, diffprice, 1, 'B2');

%% Write local volatility into file
xlswrite(outputFile, localVolCalibrated, 1, strcat('B', num2str(discretization_num_T+3,2)));

%% find mediane
mediane = size(find(ptsToEvalK < S),2); %% Find index of last in(out) the money element

%% Draw sigma best
figure;
surf(ptsToEvalK,T,localVolCalibrated);
hold on;
title('3D plot of local volatility');
xlabel('K'); ylabel('T'); zlabel('localVolCalibrated(K,T)');
saveas(gcf, strrep(outputFig,'placeholder','2'));
insertPictureToExcel( outputFile, strrep(outputFig,'placeholder','2'), 2 );

%% Draw sigma best near median
rangeOfInterest = mediane-3:mediane+3;
Z = localVolCalibrated(:,rangeOfInterest);
figure;
surf(ptsToEvalK(rangeOfInterest),T,Z);
hold on;
title('3D plot of local volatility reduced to near the money');
xlabel('K'); ylabel('T'); zlabel('localVolCalibrated(K,T)');
saveas(gcf, strrep(outputFig,'placeholder','3'));
insertPictureToExcel( outputFile, strrep(outputFig,'placeholder','3'), 3 );

%% BS surface
figure;
surf(ptsToEvalK,T,squeeze(u(index_best,:,:)));
hold on;
title('3D plot of BS option prices obtained with local volatility');
xlabel('K'); ylabel('T'); zlabel('u(K,T)');
saveas(gcf, strrep(outputFig,'placeholder','4'));
insertPictureToExcel( outputFile, strrep(outputFig,'placeholder','4'), 4 );

%% Draw diff of local vol obtained with implied vol
Z = abs(localVolCalibrated - VolImp);
figure;
surf(ptsToEvalK,T,Z);
hold on;
title('3D plot of diffbetween local Vol Calibrated and implied vol');
xlabel('K'); ylabel('T'); zlabel('diff(K,T)');
saveas(gcf, strrep(outputFig,'placeholder','5'));
insertPictureToExcel( outputFile, strrep(outputFig,'placeholder','5'), 5 );

%% Draw diff of local vol obtained with implied vol plot
if (strcmp(discretizationType,'nonuniform'))
    %% nonuniform
    figure;
    hold on;
    plot(T, sig(index_best,:,30), 'r');
    plot(T, sig(index_best,:,31), 'b');
    xlabel('T'); ylabel('sigma'); 
    figure;
    hold on;
    plot(T, Vmarket(:,31), 'r');
    plot(T, u(index_best,:,31), 'b');
    xlabel('T'); ylabel('u');
    figure;
    surf(ptsToEvalK(29:33),T,(squeeze(u(index_best,:,29:33)) - Vmarket(:,29:33))./Vmarket(:,29:33));
    hold on;
    xlabel('K'); ylabel('T'); zlabel('(u(K,T) - Vmarket(K,T))/Vmarket');
else
    %% uniform   
    figure;
    hold on;
    plot(T, localVolCalibrated(:,mediane -3 ), 'c');
    plot(T, localVolCalibrated(:,mediane -2 ), 'r');
    plot(T, localVolCalibrated(:,mediane -1 ), 'y');
    plot(T, localVolCalibrated(:,mediane), 'b');
    plot(T, localVolCalibrated(:,mediane + 1), 'k');
    plot(T, localVolCalibrated(:,mediane + 2), 'g');
    plot(T, localVolCalibrated(:,mediane +3 ), 'm');
    xlabel('T'); ylabel('localVolCalibrated');
    title('2D plot of sigma values for different K');
    legend (strcat('K=',num2str(ptsToEvalK(mediane)-3,2)),strcat('K=',num2str(ptsToEvalK(mediane)-2,2)),strcat('K=',num2str(ptsToEvalK(mediane)-1,2)), strcat('K=',num2str(ptsToEvalK(mediane),2)), strcat('K=',num2str(ptsToEvalK(mediane+1),2)), strcat('K=',num2str(ptsToEvalK(mediane+2),2)), strcat('K=',num2str(ptsToEvalK(mediane+3),2)));
    saveas(gcf, strrep(outputFig,'placeholder','6'));
    insertPictureToExcel( outputFile, strrep(outputFig,'placeholder','6'), 6 );
    
    figure;
    hold on;
    plot(T, Vmarket(:,mediane), 'r');
    plot(T, u(index_best,:,mediane), 'b');
    xlabel('T'); ylabel('u'); 
    title(strcat('2D plot comparing market price with obtained price at K=', num2str(ptsToEvalK(mediane))));
    legend (strcat('Vmarket at K=',num2str(ptsToEvalK(mediane),2)),strcat('u at K=',num2str(ptsToEvalK(mediane),2)));
    saveas(gcf, strrep(outputFig,'placeholder','7'));
    insertPictureToExcel( outputFile, strrep(outputFig,'placeholder','7'), 7 );
    
    figure;
    rangeOfInterest = mediane-3:mediane+3;
    surf(ptsToEvalK(rangeOfInterest), T, diffprice(:,rangeOfInterest));
    hold on;
    xlabel('K'); ylabel('T'); zlabel('(u(K,T) - Vmarket(K,T))/Vmarket');
    title('3D plot showing absolute diff between market price with obtained price near the money');
    saveas(gcf, strrep(outputFig,'placeholder','8'));
    insertPictureToExcel( outputFile, strrep(outputFig,'placeholder','8'), 8);
end;



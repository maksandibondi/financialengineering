function [localVolCalibrated, diffprice, ptsToEvalK, failure] = calculateLocalVolGenetic(K, T, S, r, Vmarket, VolImp, inputStructure)

%% Initial parameters
failure = 0;
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
maxIter = inputStructure.maxIter;
epsilon = inputStructure.epsilon;
concentration_weights = inputStructure.concentration_weights;
discretizationType = inputStructure.discretizationType;
interpTypeK = inputStructure.interpTypeK;
nonuniform_method = inputStructure.nonuniform_method;
includeVolAssumption = inputStructure.includeVolAssumption;
%rowByRowMutation = inputStructure.rowByRowMutation;
vasicek_assumption = inputStructure.vasicek_assumption;
outputdir = inputStructure.outputdir;
visualizeResults = inputStructure.visualizeResults;

outputFile = strcat(outputdir, '/', inputStructure.ticker, '_', strrep(inputStructure.date,'/','-'), '_', discretizationType, '_', interpTypeK, '_', 'concwgh_', strrep(num2str(concentration_weights),'.',''), 'eps_', strrep(num2str(epsilon), '.', ''), '_mut_',strrep(num2str(t),'.',''),  '.xls');
outputFig = strcat(outputdir, '/figures/', inputStructure.ticker, '_', strrep(inputStructure.date,'/','-'), '_', discretizationType, '_', interpTypeK, '_', 'concwgh_', strrep(num2str(concentration_weights),'.',''), 'eps_', strrep(num2str(epsilon), '.', ''), '_mut_',strrep(num2str(t),'.',''), '_pholder', '.png');

%% Set vol implied to zero if not passed as parameter
if (VolImp == 0)
    VolImp = zeros(size(Vmarket,1), size(Vmarket,2));
end;

%% Setting up pseudo vasicek assumtion
if (vasicek_assumption)
    r = (20-S)/100;
end;



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

%% Getting non uniform implied vol by interpolating uniform if needed
if (strcmp(discretizationType, 'nonuniform'))
    knotsK = transp(K(:,1));
    for i = 1:discretization_num_T
        VolImp_temp(i,:) = interp1(knotsK, VolImp(i,:), ptsToEvalK, 'linear' ,'extrap');
        %Vmarket_temp(i,:) = interp1(knotsK, Vmarket(i,:), ptsToEvalK, 'spline');
    end;
    VolImp = VolImp_temp;
end



%% Algo
iter = 0;
isFound = 0;
for k = 1:Nmc

    for n = 1:popsize
        if (k == 1)
            %% Seed is very significant but can be fought with mutation coefficient high
            ctrlpts = generateCtrlpts(K_0, K_l, T_0, T_l,discretization_num_T, discretization_num_K,disc_T, disc_K, VolImp, ptsToEvalK, T, includeVolAssumption,discretizationType, S);
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
                            if (ctr(n,m,g)<0) %%|| ctr(n,m,g)>1) % check the constraints
                               ctr(n,m,g) = init;
                            end;
                        end;
                    end;
%                 end;
            end;

         else 
             isFound = 1;
             break;
         end; 
        
         
         iter = iter+1;
         if (iter > maxIter)
             failure = 1;
             break;
         end;
         minfitnessList(iter) = minfitness;
         display(minfitness);
end;
localVolCalibrated = squeeze(sig(index_best,:,:));
[~, ~, diffprice] = sumOfSqrDif_(u(index_best,:,:), Vmarket(:,:), S, ptsToEvalK, epsilon, concentration_weights);

%%  Visual part
if visualizeResults

    if isFound

        %% Write axes to file
        xlswrite(outputFile, ptsToEvalK, 1, 'B1');
        xlswrite (outputFile, T, 1, 'A2');
        %% Write diff in prices into file
        xlswrite(outputFile, diffprice, 1, 'B2');

        rownum = discretization_num_T+3;
        xlswrite(outputFile, ptsToEvalK, 1, strcat('B', num2str(rownum,2)));
        rownum = rownum+1;
        xlswrite (outputFile, T, 1, strcat('A', num2str(rownum,2)));
        %% Write local volatility into file
        xlswrite(outputFile, localVolCalibrated, 1, strcat('B', num2str(rownum,2)));

        rownum = rownum + discretization_num_T+3;
        xlswrite(outputFile, ptsToEvalK, 1, strcat('B', num2str(rownum,2)));
        rownum = rownum + 1;
        xlswrite (outputFile, T, 1, strcat('A', num2str(rownum,2)));
        %% Write implied volatility into file
        xlswrite(outputFile, VolImp, 1, strcat('B', num2str(rownum,2)));

        %% find mediane
        mediane = size(find(ptsToEvalK < S),2); %% Find index of last in(out) the money element

        %% Draw sigma best
        figure;
        surf(ptsToEvalK,T,localVolCalibrated);
        hold on;
        title('3D plot of local volatility');
        xlabel('K'); ylabel('T'); zlabel('localVolCalibrated(K,T)');
        saveas(gcf, strrep(outputFig,'pholder','2'));
        insertPictureToExcel( outputFile, strrep(outputFig,'pholder','2'), 2 );

        %% Draw sigma best near median
        rangeOfInterest = mediane-3:mediane+3;
        Z = localVolCalibrated(:,rangeOfInterest);
        figure;
        surf(ptsToEvalK(rangeOfInterest),T,Z);
        hold on;
        title('3D plot of local volatility reduced to near the money');
        xlabel('K'); ylabel('T'); zlabel('localVolCalibrated(K,T)');
        saveas(gcf, strrep(outputFig,'pholder','3'));
        insertPictureToExcel( outputFile, strrep(outputFig,'pholder','3'), 3 );

        %% BS surface
        figure;
        surf(ptsToEvalK,T,squeeze(u(index_best,:,:)));
        hold on;
        title('3D plot of BS option prices obtained with local volatility');
        xlabel('K'); ylabel('T'); zlabel('u(K,T)');
        saveas(gcf, strrep(outputFig,'pholder','4'));
        insertPictureToExcel( outputFile, strrep(outputFig,'pholder','4'), 4 );

        %% Draw diff of local vol obtained with implied vol
        Z = abs(localVolCalibrated - VolImp);
        figure;
        surf(ptsToEvalK,T,Z);
        hold on;
        title('3D plot of diffbetween local Vol Calibrated and implied vol');
        xlabel('K'); ylabel('T'); zlabel('diff(K,T)');
        saveas(gcf, strrep(outputFig,'pholder','5'));
        insertPictureToExcel( outputFile, strrep(outputFig,'pholder','5'), 5 );

        %% Draw diff of local vol obtained with implied vol plot 
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
        saveas(gcf, strrep(outputFig,'pholder','6'));
        insertPictureToExcel( outputFile, strrep(outputFig,'pholder','6'), 6 );

        figure;
        hold on;
        plot(T, Vmarket(:,mediane), 'r');
        plot(T, u(index_best,:,mediane), 'b');
        xlabel('T'); ylabel('u'); 
        title(strcat('2D plot comparing market price with obtained price at K=', num2str(ptsToEvalK(mediane))));
        legend (strcat('Vmarket at K=',num2str(ptsToEvalK(mediane),2)),strcat('u at K=',num2str(ptsToEvalK(mediane),2)));
        saveas(gcf, strrep(outputFig,'pholder','7'));
        insertPictureToExcel( outputFile, strrep(outputFig,'pholder','7'), 7 );

        figure;
        rangeOfInterest = mediane-3:mediane+3;
        surf(ptsToEvalK(rangeOfInterest), T, diffprice(:,rangeOfInterest));
        hold on;
        xlabel('K'); ylabel('T'); zlabel('(u(K,T) - Vmarket(K,T))/Vmarket');
        title('3D plot showing absolute diff between market price with obtained price near the money');
        saveas(gcf, strrep(outputFig,'pholder','8'));
        insertPictureToExcel( outputFile, strrep(outputFig,'pholder','8'), 8);

        %% Draw convergence
        figure;
        hold on;
        plot(1:iter, minfitnessList);
        title('2D plot of  convergence of genetic algorithm');
        xlabel('iter'); ylabel('fitness');
        saveas(gcf, strrep(outputFig,'pholder','9'));
        insertPictureToExcel( outputFile, strrep(outputFig,'pholder','9'), 9);
    else 
        %% Draw convergence
        figure;
        hold on;
        plot(1:iter, minfitnessList);
        title('2D plot of  convergence of genetic algorithm');
        xlabel('iter'); ylabel('fitness');
        saveas(gcf, strrep(outputFig,'pholder','1'));
        xlswrite(outputFile, 0, 1, 'B1');
        insertPictureToExcel( outputFile, strrep(outputFig,'pholder','1'), 1);
    end;

end;


clear; clc;
% Troubles:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2) The same non-uniform discetization to be used for market data
% simulation (?), PDE resolution (function invoked) !!!!!!!!!
%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3) spline creation (knots non-uniform) using Working De Boor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initial data
K_0 = 10; K_l = 70;
T_0 = 0.0164; T_l = 0.4192;
discretization_num_K = 61;
discretization_num_T = 10;
S = 20.46;
r = 0;
% Representation of T vector
T(1) = T_0;
dt = (T_l-T_0)/(discretization_num_T-1);
for i = 2:discretization_num_T
    T(i) = T(i-1) + dt;
end;


T0 = 1000;
t = 0.05;
Nmc = 2^14;
M = 64;
popsize = 20;


%% Simulation of market data
for (i = 1:popsize)
    prices(i,:,:) = cell2mat(struct2cell(load('Vmarket.mat'))); 
end;

%% Nonuniform discretization grid(the same as in the pricer)
delta_log_K = (log(K_l) - log(K_0)) / (discretization_num_K+3); % discr .. K + 3
y = zeros(1,discretization_num_K+1);
y_0 = log(S);
for i = 2:discretization_num_K
	y(i) = 0.5*log((log(K_l) + log(K_0)+i*delta_log_K) / (log(K_l) - (log(K_0) + i*delta_log_K))) + y_0; 
end;
disc_T = 5; disc_K = 20;
ptsToEvalK = exp(y(1:end-1)) % size = discretization_num_K+1
iter = 0;

%% Algo
for k = 1:Nmc

for n = 1:popsize
    
if (k == 1)
%ctrlpts = rand(disc_T+1,disc_K+1-4);
ctrlpts = rand(disc_T+1,disc_K+1);
sigmaaa = SplineLinear2DInterpol(T_0,T_l,K_0,K_l,disc_T, disc_K,ptsToEvalK,T, ctrlpts);
%sigmaaa = ones(discretization_num_T+1,discretization_num_K+1)*0.1;
%Spline evaluated at the points at which our non-uniform discretization is
% done. How to give it the size 31*198 as well?
else
ctrlpts(:,:) = ctr(n,:,:);
    sigmaaa = SplineLinear2DInterpol(T_0,T_l,K_0,K_l,disc_T, disc_K,ptsToEvalK,T, ctrlpts);    
end;

%% Pricing
[u(n,:,:),K,T] = Pricer(sigmaaa, K_0, K_l, T_0, T_l, discretization_num_K, discretization_num_T, S, r);



%fitness(n) = sumOfSqrDif(u(n,10:20,96:104),prices(n,10:20,96:104)); % cost
fitness(n) = sumOfSqrDif(u(n,:,:), prices(n,:,:)); % cost

sig(n,:,:) = sigmaaa;
ctr(n,:,:) = ctrlpts;
end;

%% Genetic part
    [minfitness, index_best] = min(fitness);
    if (minfitness>2)
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
    
        % Mutation ( applied independently to all ctrl points)
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
surf(K,T,ones(size(T,2),size(K,2))*0.2, 'FaceColor', [1 0 1]);
xlabel('K'); ylabel('T'); zlabel('sigma(K,T)');

figure;

surf(K(1:25),T, squeeze(u(index_best,:,1:25) - prices(index_best,:,1:25)./u(index_best,:,1:25)));
hold on;
xlabel('K'); ylabel('T'); zlabel('u(K,T) - Vmarket(K,T)');

figure;
surf(K,T,squeeze(u(index_best,:,:)));
hold on;
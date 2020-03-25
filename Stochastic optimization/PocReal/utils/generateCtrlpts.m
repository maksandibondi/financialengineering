function [ctrlpts] = generateCtrlpts(K_0, K_l, T_0, T_l,discretization_num_T, discretization_num_K,disc_T, disc_K, VolImp, ptsToEvalK, T,  includeVolAssumption,discretizationType, S)

rng('shuffle');
            ctrlpts = rand(disc_T,disc_K);
            
            if (strcmp(includeVolAssumption, 'zero'))
                ctrlpts = ctrlpts*0; % zero
                return;
            end;
            
            if (strcmp(includeVolAssumption, 'implicit'))
                ctrlpts = 0.5*getSeedFromImplicitVol(K_0, K_l, T_0, T_l, discretization_num_T, discretization_num_K,disc_T, disc_K,VolImp, ptsToEvalK, T, discretizationType, S);
            end;
            return;
end
function [ VolImp ] = calculateImpVolMatrixNewton( T, K, S, r, Vmarket, eps )
%CALCULATEIMPVOLNEWTON Summary of this function goes here
%   Detailed explanation goes here

for i = 1:size(T,1)
    for j = 1:size(K,1)
        %% Arbitrage condition
        if (max(S-K(j)*exp(-r*T(i)),0) < Vmarket(i,j) < S) && (T(i)~=0) % for call
        % if max(K*exp(-r*T)-S0,0) < M < S0 % for put       

            sigma = sqrt(2*abs((log(S/K(j))+r*T(i))/T(i)));

            %% Newton algorithm
            it=0;
            while abs(BSTheory(sigma,S,r,T(i),K(j))-Vmarket(i,j)) > eps
                it = it + 1;
                sigma = sigma - (BSTheory(sigma,S,r,T(i),K(j)) - Vmarket(i,j))/Vega(sigma,S,r,T(i),K(j));
%                 if it>100 && it<110
%                     X = ['i is ',num2str(i)]
%                     Y = ['j is ',num2str(j)]
%                     Z = ['sigma is ',num2str(sigma)]
%                     BSTH = ['bs theory is ' , num2str(BSTheory(sigma,S,r,T(i),K(j)))]
%                     vm = ['vmarket is ', num2str(Vmarket(i,j))]
%                     vegaa = ['vega is ', num2str(Vega(sigma,S,r,T(i),K(j)))]
%                 end;
            end

            VolImp(i,j) = sigma;

        else

            VolImp(i,j) = 0;   

        end;
    end; 
end;

return;


function [FKRateMax] = FKR_Ht_Genetic(Nrounds,AlgVars,GeneticVars)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This code calculates the lower bound for the key rate in the finite
%%% block size regime.
%%%
%%% INPUT: Nrounds   - No. of rounds employed by Alice and Bob
%%%        Ddelta    - No. of modules
%%%        Maxf,Minf - bounds for the min-tradeoff function 
%%%        low_bound - min-tradeoff function (NOT the DW key rate!)
%%%        Es        - Smoothing parameter
%%%        EphS      - Smoothing parameter of the physical protocol
%%%        EcPE      - Completeness tolerance for PE
%%%        Etom      - Tolerance for tomography
%%%        EphNA     - Nonabortion tolerance for the physical protocol
%%%        Eec       - Error correction tolerance
%%%
%%% OUTPUT: FRate - Finite key rate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Data preparation


Ppoint   = AlgVars.probdist;        % Primal point
Dpoint   = AlgVars.dualdist;        % Dual point
R_DW     = AlgVars.R_inf;           % Key rate
L        = GeneticVars.L;           % Distance
Ppe      = AlgVars.Ppe;             % PE probability
MaxMinf  = AlgVars.MaxMinf;         % Spread of dualvars
SDPbound = R_DW + AlgVars.deltaEC;  % SDP bound
Ddelta   = GeneticVars.Ddelta;      % No. of subsectors


% Quick check
if L<7
    if MaxMinf > 360 || R_DW <= 0
        FKRateMax = 0;
        return
    end
else
    if MaxMinf > 150 || R_DW <= 0
        FKRateMax = 0;
        return
    end
end

%% Security parameters


Epsilons = GeneticVars.Epsilons;
Es       = Epsilons(1);
EcPE     = Epsilons(2);
EnA      = Epsilons(3);
Eec      = Epsilons(4);
Etom     = Epsilons(5);

% EphS is limited by Lemma 2 of the paper
EphS    = Es*1.00001 + sqrt(2*Etom/EnA);

% Validity check
if Es>1-sqrt(2*Etom/EnA)
    error('Smoothing Es is too big!')
elseif 2*Etom>EnA
    error('Etom is bigger than EnA')
end

% Function Gamma(x)
if Es^2 < eps
    % Taylor approximation due to floating point precision
    Ga = @(x) -log2(0.5*x.^2);
else
    Ga = @(x) -log2(1-sqrt(1-x.^2));
end

    

%% Optimization of the finite key rate


% Dimensions of output O, key register and PE+Tomo registers
dOut= 17*(4*(Ddelta+1) + 1)*5*5;
dZ  = 4 + 1;
dXZ = 17*(4*(Ddelta+1) + 1);

    
% Grid search for the optimal scaling of 
% the key rounds and the EAT parameter a

% Rough bounds for a-1
denom = Nrounds*(R_DW-sqrt(3*log2(2/Eec)/Nrounds)*log(dZ+3))-...
    2*log2(1/EphS)-2*Ga(Es/4)-3*Ga(Es/16);
minimum = (log2(1/(EnA-Etom))+Ga(Es/16))/denom;

jj=0;
FKRateMax = 0; 
FKRate    = -1;
stalling  = 0;

for scale=minimum:minimum/10:1

    a = 1+scale;
    jj = jj + 1;
    if mod(jj,500)==0 
        fprintf('jj=%d \n',jj); 
        if FKRate<0
            fprintf('So far: FKRate = %.8e \n',FKRate);
        end
    end

    if mod(jj,500)==0
        return
    end
    
    % p^key = 1 - Nrounds^(-b)
    % p^PE  = Nrounds^(-b) Ppe
    % p^tom = Nrounds^(-b) (1-Ppe)

    % Here we optimize the scaling b wrt the value of a
    A = (a-1)*log(2)*MaxMinf^2 /2;
    B = (a-1)*log(2)*log2(2*dOut^2 +1);
    C = SDPbound + Ppe*log2(5) + log2(dXZ);
    
    clear bta
    syms  bta
    bta = vpasolve(Nrounds^(2*bta)* (A+ ...
        (B*MaxMinf^2)/(2*sqrt(2+(Nrounds^bta)*MaxMinf^2)))-C==0);
    b = double(bta);
    [dTolPE,dTolTom] = DeltaTol(Nrounds,b,Ppe,Ppoint,Dpoint,EcPE,Etom);


    % Quick check that the value of b is useful
    if R_DW < C*Nrounds^(-b) + A*Nrounds^b +...
            B*sqrt(2+ Nrounds^b *MaxMinf^2) +...
            (2+log2(2*dOut^2+1)^2)*(a-1)*log(2)/2 +dTolPE+dTolTom
        continue
    end
    
    % Calculation of the coeffs. for the finite key rate
    One = A*Nrounds^b + B*sqrt(2 + Nrounds^b *MaxMinf^2) + ...
        (2+log2(2*dOut^2 + 1)^2)*(a-1)*log(2)/2;
    
    K_exp = (a-1)*(2*log2(dOut)+MaxMinf);
    K_num = log(2^(2*log2(dOut) + MaxMinf) + exp(2))^3 * 2^(K_exp);
    K_den = 6*log(2)*(2-a)^3;
    Two   = (K_num*(a-1)^2)/K_den;
    
    Three = C*Nrounds^(-b);
    
    Four = Nrounds^(-1/2) *(sqrt(log(32/((EnA-Etom)*Es^2))/2)*log2(5) +...
        log2(3+dZ)*sqrt(3*log2(2/Eec)) +...
        sqrt(log(512/((EnA-Etom)*Es^2))/2)*log2(dXZ));
    
    Five = Nrounds^(-1) *((Ga(Es/16) + a*log(1/(EnA-Etom)))/(a-1) + ...
        2*log(1/EphS) + 2*Ga(Es/4) + 3*Ga(Es/16));
    
    % Finite key rate
    FKRate = R_DW - One - Two - Three - Four - Five -dTolPE -dTolTom;
    
    % Check the value of the key rate        
    if FKRate<FKRateMax
        % Convergence criterion - find the peak, and once it is
        % verified that it is the peak, break the simulation
        
        if FKRateMax>0
            stalling=stalling+1;
            if stalling>10
                fprintf('  Optimality reached: %.2e \n',FKRate);
                return
            end
        else
            stalling=0;
        end
        continue
    else
        if mod(jj,150)==0
            fprintf('Success at %.6e, %.6f, %d, %e \n',a-1,b,jj,FKRate);
        end
        % Write the new optimal key rate 
        FKRateMax = FKRate;
    end
end
    
  

end
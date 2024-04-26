%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code performs the optimization for the asymptotic secret key
% via a Frank-Wolfe algorithm. This is also a previous step towards
% distilling the key in the finite regime.
%
% REQUIREMENTS
% - YALMIP, with solvers MOSEK and SDPT3
% - CVX
% - Quantinf, the MATLAB package by Toby Qubitt
%
% IN - Set of initial variables, given by a basis generated with
%      OpGenerator_Ht.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Load set of matrices
load('Nc12D09d09.mat');

BasisVars = [];
BasisVars.Nc       = Nc;
BasisVars.Delta    = Delta;
BasisVars.Ddelta   = Ddelta;
BasisVars.GammaRaw = GammaRaw;
BasisVars.GammaBar = GammaBar;
BasisVars.alphaBar = alphaBar;
BasisVars.LambdaA  = LambdaA;


    
% Optimal amplitudes for each distance (only for the Shannon limit)
Opt_Ampli = [1.07,1.05,1.03,1.01,0.99, 0.96,0.93,0.9,0.87,0.85, ...
    0.83,0.83,0.83,0.83,0.83, 0.82,0.82,0.82,0.82,0.82, ...
    0.81,0.81,0.8,0.8,0.79, 0.79,0.78,0.78,0.77,0.77, ...
    0.77,0.77,0.77,0.76,0.76, 0.75,0.75,0.75,0.74,0.74, ...
    0.74,0.74,0.74,0.74,0.74, 0.73,0.73,0.73,0.73,0.73, ...
    0.73,0.73,0.73,0.73,0.73, 0.72,0.72,0.72,0.72,0.72, ...
    0.72,0.72,0.72,0.72,0.72, 0.71,0.71,0.71,0.71,0.71,...
    0.71]; % For Length>70 the optimal amplitude is constant


% First run
L0     = 0;
Length = L0; 
if Length<70
    amplitude = Opt_Ampli(Length+1);
else
    amplitude = 0.71;
end

% Excess noise
xi = 0.01;

% Preparation of variables
InitialVars = InitialCode_Ht(Length,amplitude,xi,BasisVars);

% Optimization algorithm
[AlgVars] = Alg1_Ht(BasisVars,BasisAB,InitialVars,xi,Length);


% Extraction of results
R_inf1  = AlgVars.R_inf1;
R_inf   = AlgVars.R_inf;
MaxMinf = AlgVars.MaxMinf;
deltaEC = AlgVars.deltaEC;
epsilon = AlgVars.epsilon;
T       = AlgVars.T;
Ppe     = AlgVars.Ppe;

probdist = AlgVars.probdist;
dualdist = AlgVars.dualdist;


% Results with heading
FILE1 = fopen(['Rate_xi' ...
    num2str(xi*100) 'Nc' num2str(Nc) 'D0' ...
    num2str(Delta*10) 'd0' num2str(delta*10) '.csv'],'a');
fprintf(FILE1,'Nc, xi, Delta, delta \n');
fprintf(FILE1,'%d, %.2f, %.1f, %.1f \n\n',Nc,xi,Delta,delta);
fprintf(FILE1,'L, Amp, R_1, R_CVX, Maxf - Minf, deltaEC, eps, T, Ppe \n');
fprintf(FILE1,'%d, %.2f, %.10f, %.10f, %.10f, %.10f, %d, %d, %.15f \n',...
    Length, amplitude, R_inf1, R_inf, MaxMinf, deltaEC, epsilon, T, Ppe);
fclose(FILE1);

% Primal variables
FILE2 = fopen(['VarsP_xi' ...
    num2str(xi*100) 'Nc' num2str(Nc) 'D0' ...
    num2str(Delta*10) 'd0' num2str(delta*10) '.csv'],'a');
fprintf(FILE2,'L, Amp, Ppe, pVars,\n');
fprintf(FILE2,'%d, %.2f, %.15f, ',Length,amplitude,Ppe); 
fprintf(FILE2,'%.16e, ',probdist(:));
fprintf(FILE2,'\n');
fclose(FILE2);

% Dual variables
FILE3 = fopen(['VarsD_xi' ...
    num2str(xi*100) 'Nc' num2str(Nc) 'D0' ...
    num2str(Delta*10) 'd0' num2str(delta*10) '.csv'],'a');
fprintf(FILE3,'L, Amp, Ppe, dVars,\n');
fprintf(FILE3,'%d, %.2f, %.15f, ',Length,amplitude,Ppe); 
fprintf(FILE3,'%.16e, ',dualdist(:));
fprintf(FILE3,'\n');
fclose(FILE3);


%% General loop

for i = L0+1:100
    % Take a new amp and run again
    if i<71
	    Length    = Length+1;
        amplitude = Opt_Ampli(Length+1);
    else
        Length    = Length+5;
        amplitude = 0.71;
    end

    % Stop the loop when the distance is too large
    if Length>=230
        fprintf('Process finished \n');
        break
    end
    
    % Preparation of variables
    InitialVars = InitialCode_Ht(Length,amplitude,xi,BasisVars);
    
    % First algorithm
    [AlgVars] = Alg1_Ht(BasisVars,BasisAB,InitialVars,xi,Length);

    % Extraction of results
    R_inf1  = AlgVars.R_inf1;
    R_inf   = AlgVars.R_inf;
    MaxMinf = AlgVars.MaxMinf;
    deltaEC = AlgVars.deltaEC;
    epsilon = AlgVars.epsilon;
    T       = AlgVars.T;
    Ppe     = AlgVars.Ppe;
    
    probdist = AlgVars.probdist;
    dualdist = AlgVars.dualdist;
     
    % Save the results in the .csv file
    FILE1 = fopen(['Rate_xi' ...
    num2str(xi*100) 'Nc' num2str(Nc) 'D0' ...
    num2str(Delta*10) 'd0' num2str(delta*10) '.csv'],'a');
    fprintf(FILE1,'%d, %.2f, %.10f, %.10f, %.10f, %.10f, %d, %d, %.15f \n',...
        Length, amplitude, R_inf1, R_inf, MaxMinf, deltaEC, epsilon, T, Ppe);
    fclose(FILE1);

    % Primal variables
    FILE2 = fopen(['VarsP_xi' ...
    num2str(xi*100) 'Nc' num2str(Nc) 'D0' ...
    num2str(Delta*10) 'd0' num2str(delta*10) '.csv'],'a');
    fprintf(FILE2,'%d, %.2f, %.15f, ',Length,amplitude,Ppe); 
    fprintf(FILE2,'%.16e, ',probdist(:));
    fprintf(FILE2,'\n');
    fclose(FILE2);
    
    % Dual variables
    FILE3 = fopen(['VarsD_xi' ...
    num2str(xi*100) 'Nc' num2str(Nc) 'D0' ...
    num2str(Delta*10) 'd0' num2str(delta*10) '.csv'],'a');
    fprintf(FILE3,'%d, %.2f, %.15f, ',Length,amplitude,Ppe); 
    fprintf(FILE3,'%.16e, ',dualdist(:));
    fprintf(FILE3,'\n');
    fclose(FILE3);

end    
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This code performs the genetic optimization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('Nc12D09d09.mat');


BasisVars = [];
BasisVars.Nc       = Nc;
BasisVars.Delta    = Delta;
BasisVars.Ddelta   = Ddelta;
BasisVars.GammaRaw = GammaRaw;
BasisVars.GammaBar = GammaBar;
BasisVars.alphaBar = alphaBar;
BasisVars.LambdaA  = LambdaA;


% Excess noise and no. of rounds
xi = 0.01;
nr = '5e11';
Nr = str2double(nr);

% epsilon-parameters
Es   = 1e-10;
EcPE = 1e-10;
EnA  = 1e-4;
Eec  = 1e-10;
Etom = 1e-10;

% Create file and header
NAME = ['FKRate_n' nr 'xi' num2str(xi*100) '.csv' ];
FILE = fopen(NAME,'a');

fprintf(FILE, 'L, FKRate, Zeta \n');
fclose(FILE);

% Optimal amplitudes for distances [0,70]
Opt_Ampli = [1.07,1.05,1.03,1.01,0.99, 0.96,0.93,0.9,0.87,0.85, ...
    0.83,0.83,0.83,0.83,0.83, 0.82,0.82,0.82,0.82,0.82, ...
    0.81,0.81,0.8,0.8,0.79, 0.79,0.78,0.78,0.77,0.77, ...
    0.77,0.77,0.77,0.76,0.76, 0.75,0.75,0.75,0.74,0.74, ...
    0.74,0.74,0.74,0.74,0.74, 0.73,0.73,0.73,0.73,0.73, ...
    0.73,0.73,0.73,0.73,0.73, 0.72,0.72,0.72,0.72,0.72, ...
    0.72,0.72,0.72,0.72,0.72, 0.71,0.71,0.71,0.71,0.71,...
    0.71];

for L=1:2:40

    % Load the initial variables
    amplitude   = Opt_Ampli(L+1);
    InitialVars = InitialCode_Ht(L,amplitude,xi,BasisVars);

    GeneticVars          = [];
    GeneticVars.L        = L;
    GeneticVars.Epsilons = [Es,EcPE,EnA,Eec,Etom];
    GeneticVars.GammaRaw = GammaRaw;
    GeneticVars.Ddelta   = Ddelta;

    % Call the algorithm 
    GeneticOpt(xi,InitialVars,GeneticVars,Nr,NAME)
end

fprintf('Genetic subroutine finished \n');
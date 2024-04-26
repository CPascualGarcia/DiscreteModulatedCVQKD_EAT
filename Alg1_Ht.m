%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code performs the minimization of the Devetak-Winter rate
% according to Winnick et al. and the FR provided by
% Hu et al. The minimization is based on an SDP solved with the Frank-Wolfe
% algorithm.
% IN
%  BasisAB - Matrix basis for Alice&Bob's Hilbert space
%  GammaBar- Orthonormalized version of Bob's operators
%  gammaBar- orthonormalized version of Bob's statistical distribution
%  gammaRaw- primal variables (Bob's stat. distribution for PE)
%  Gr - Facial reduced version of the coherent key map
%  Zr - Facial reduced version of the dephasing map for the key register
%
% OUT
%  skipAlg2 - This indicates early convergence for the whole minimization
%  rho      - State shared by Alice and Bob, suboptimal according to the SDP
%  epsilon  - Numerical error of the SDP constraints
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [AlgVars] = Alg1_Ht(BasisVars,BasisAB,InitialVars,xi,Length)


%% Load the variables


GammaBar = BasisVars.GammaBar;
GammaRaw = BasisVars.GammaRaw;
gammaRaw = InitialVars.gammaRaw;
Gr       = InitialVars.Gr;
Zr       = InitialVars.Zr;
gammaBar = InitialVars.gammaBar;


%% Generation of the seed

% We start by generating the seed - an initial state to proceed with the
% Frank-Wolfe algorithm. For this we use Yalmip to calculate a density
% matrix that satisfies the constraints gammaBar within a certain margin of
% tolerance

% In the optimization problem, we generate a seed with a perturbation for 
% the constraint. If \rho is not >= 0, the algorithm relaxes the margin 
% and tries again

[Dim,~,S] = size(GammaBar);
rhozero   = zeros(Dim);
for i = 1:S
    rhozero = rhozero + GammaBar(:,:,i)*gammaBar(i);
end

% Variables
omega = sdpvar(1,Dim^2 - S,'full','real');
lambd = sdpvar(1,1);
Omega = BasisAB(:,:,S+1:end);
DeltaRho = squeeze(sum(Omega.*repmat(reshape(omega,1,1,Dim^2 - S),[Dim,Dim,1]),3));

% Optimization
Optimizable = [ rhozero + DeltaRho>= lambd*eye(Dim)];

% Remember to set the solver: sdpsettings('solver','sdpt3|mosek')
diagnostics = optimize(Optimizable,-lambd,sdpsettings('solver','sdpt3'));
DelRho      = value(DeltaRho);
rho         = rhozero + DelRho;

% Hermiticity
if ~ishermitian(rho); rho = (rho + rho')/2; end
% Slight correction to avoid neg. eigenv.
rho = rho + eye(Dim)*1e-17;

[~,c2] = chol(rho); % Check PSDness
H = 0; 
while c2 ~= 0
    ei = eig(rho);
    rho = rho - eye(Dim)*ei(1);
    [~,c2] = chol(rho);

    if H > 50
        break
    end
    H = H +1;
end
% Again to avoid small neg. eig.
rho = rho + eye(Dim)*1e-17;

numError = zeros(1,S);
for i = 1:S
    numError(1,i) = abs(trace(rho*GammaRaw(:,:,i)) - gammaRaw(i));
end
epsilon = max(numError);
fprintf('Max error is given by eps = %d \n', epsilon);


fprintf('Initial density matrix defined \n');


%% Optimization loop

fprintf('Starting FW algorithm \n');


SOLVTAG = 1;
SOLVER  = 'sdpt3';
ROUNDS  = 300;
for i = 1:ROUNDS
    fprintf('Round %d ---------------\n',i);
    % Gradient for the Relative Entropy
    GradD = Grad_Ht(rho,Gr,Zr);
    
    % Variables of the optimization
    omega = sdpvar(1,Dim^2 - S,'full','real');
    Omega = BasisAB(:,:,(S+1):end);
    DeltaRho = squeeze(sum(Omega.*repmat(reshape(omega,1,1,Dim^2 - S),[Dim,Dim,1]),3));
    
    % Optimization
    % Notice how we have a constraint to ensure PSDness
    constraint  = DeltaRho + rho;
    Optimizable = [constraint >=0];
    diagnostics = optimize(Optimizable,real(trace(transpose(DeltaRho)*GradD)),sdpsettings('solver',SOLVER,'verbose',0));
    DelRho      = value(DeltaRho);
    fprintf('Value of the trace: %d\n',real(trace(transpose(DelRho)*GradD)));
    fprintf('Error in the minimization: i*%d\n',imag(trace(transpose(DelRho)*GradD)));
    fprintf(' ----------- \n Perfoming the line search \n');
    
    % Minimization f(rho_i + lambda*DeltaRho) wrt lambda in [0,1]
    [lambda_m,NewD] = MinSearch(rho,DelRho,Gr,Zr);
    
    fprintf('Minimum at lambda = %d \n',lambda_m);
    OldD = RelEnt_Ht(rho,Gr,Zr);   
    fprintf('New ObjVal: %.14f, gap ObjVal: %d \n',NewD,NewD-OldD);
    
    % If the algorithm reaches a stagnation point, the solver is changed
    if abs(OldD-NewD)<1e-14 && SOLVTAG == 0
        SOLVER  = 'sdpt3';
        SOLVTAG = 1;
	elseif abs(OldD-NewD)<1e-14 && SOLVTAG == 1
        SOLVER  = 'mosek';
        SOLVTAG = 0;
    end
    
    % If the coeff. of the minimization is too small, we stop the loop
    if lambda_m <= 1e-10; break; end
    
    % Apply the min. and start over
    rho = rho + lambda_m*DelRho;
    
    % Check whether the lower bound is satisfactory
    if mod(i,5) == 0
        [AlgVars] = Alg2_Ht(GammaRaw,rho,xi,InitialVars);
        if abs((AlgVars.R_inf1 - AlgVars.R_inf)/AlgVars.R_inf1) < 0.02
            fprintf('Optimization finished \n');
            break
        end
    end
end

% In the end we obtain a density matrix rho, which will be used for the 
% finite-key analysis

%% Save the suboptimal state
AlgVars.rho = rho;

writematrix(rho,['Finite_Analysis/Matrices/BestRho_L' num2str(Length) 'xi' ...
    num2str(100*xi) '.csv'])

fprintf('-----------------------\nExecution completed\n');

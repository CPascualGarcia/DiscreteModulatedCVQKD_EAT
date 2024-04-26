%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code performs the lower-bounding of the Devetak-Winter rate
% according to Winnick et al. as well as the basic characterization of the
% min-tradeoff function (i.e. the dual variables of the SDP, plus the
% extremals of the function).
% IN
%  BasisAB - Matrix basis for Alice&Bob's Hilbert space
%  GammaBar- Orthonormalized version of Bob's operators
%  gammaBar- orthonormalized version of Bob's statistical distribution
%  gammaRaw- primal variables (Bob's stat. distribution for PE)
%  Gr      - Facial reduced version of the coherent key map
%  Zr      - Facial reduced version of the dephasing map for the key register
%
% OUT
%  rho      - State shared by Alice and Bob, suboptimal according to the SDP
%  epsilon  - Numerical error of the SDP constraints
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [AlgVars] = Alg2_Ht(GammaRaw,rho,xi,InitialVars)


%% Load variables

amp      = InitialVars.amp;
pA       = InitialVars.pA;
eta      = InitialVars.eta;
gammaRaw = InitialVars.gammaRaw;
Gr       = InitialVars.Gr;
Zr       = InitialVars.Zr;

%% Constants

% Gradient and Relative Entropy
GradD = Grad_Ht(rho,Gr,Zr);
D     = RelEnt_Ht(rho,Gr,Zr);

% Constant G0
G0 = real(D - trace(transpose(rho)*GradD));

fprintf('Grad(D) and D defined \n');


%% Preparation of the maximization in beta

% Calculation of parameter epsilon, which denotes the 
% numerical imprecision in rho and the constraints

% Modification of the set of constraints. We remove the GGM tomography to
% place the ICPOVM (which is l.d. on the GGMs and the PE operators)
[Dim,~,S] = size(GammaRaw);

% ICPOVM for the tomography
[ICPOVM,gammaIC] = IC_POVM(amp,pA);

numError = zeros(S+4,1);
for i = 1:S-12
    numError(i) = abs(trace(rho*GammaRaw(:,:,i))-gammaRaw(i));
end
for j = 1:16
    numError(i+j) = abs(trace(rho*kron(ICPOVM(:,:,j),eye(Dim/4)))-gammaIC(j));
end
epsilon = max(numError);
epsilon = epsilon*1.001; %This ensures the state is inside our state space
fprintf('Max error is given by epsilon = %d \n', epsilon);

fprintf('Elements of the maximization in beta defined \n');


%% Maximization with CVX

probdist   = [gammaRaw(1:end-12)';gammaIC(:)];
sum_Gammas = 0;
cvx_precision best
cvx_solver sdpt3
cvx_begin
    variable nu(S+4)
    variable mu(S+4)
    for i=1:S-12
        sum_Gammas = sum_Gammas +nu(i)*transpose(GammaRaw(:,:,i));
    end
    for i=1:16
        sum_Gammas = sum_Gammas +nu(S-12+i)*transpose(kron(ICPOVM(:,:,i),eye(Dim/4)));
    end
    maximize (dot(probdist,nu)-epsilon*dot(ones(1,S+4),mu))
    subject to % Set of constraints
        GradD - sum_Gammas == hermitian_semidefinite(Dim)
        for i=1:S+4
            abs(nu(i)) <= mu(i)
        end
cvx_end

% Improvement of the maximization
% This must be done by hand bc. CVX is not precise
% enough for this maximization
for i=1:S+4
	if abs(nu(i))*1.001 <= mu(i)
    	mu(i) = abs(nu(i))*1.001;
	end
end

% Max-Min terms for the finite key
maximum = dot(probdist,nu) - epsilon*dot(ones(1,S+4),mu);
[MaxMinf,Ppe,~] = MaxMin(nu);

% Checkings
T = -18;
[~,c2] = chol(GradD - sum_Gammas);

% If the maximization was not achieved properly, we 
% repeat it with a perturbation T for the >=0
% constraint (usual source of errors)
while c2 ~= 0
    
    T = T + 1;
    sum_Gammas = 0;
    
    cvx_begin
        cvx_precision high
        variable nu(S+4)
        variable mu(S+4)
        for i=1:S-12
            sum_Gammas = sum_Gammas + nu(i)*transpose(GammaRaw(:,:,i));
        end
        for i=1:16
            sum_Gammas = sum_Gammas + nu(i)*transpose(kron(ICPOVM(:,:,i),eye(Dim/4)));
        end
        maximize (dot(probdist,nu) - epsilon*dot(ones(1,S+4),mu))
        subject to % Set of constraints
            GradD - sum_Gammas - 10^(T)*eye(Dim) == hermitian_semidefinite(Dim)
            for i=1:S+4
                abs(nu(i)) <= mu(i)
            end
    cvx_end
    
    % Improvement of the maximization
    for i=1:S+4
        if abs(nu(i))*1.001 <= mu(i)
            mu(i) = abs(nu(i))*1.001;
        end
    end
    maximum = (dot(probdist,nu) - epsilon*dot(ones(1,S+4),mu));
    [MaxMinf,Ppe,~] = MaxMin(nu);

    if T == -9
	% A perturbation was not enough to achieve the 
	% maximization - we have to start over or abort
        fprintf('MAXIMIZATION FAILED------------------- \n')
        maximum = 0;
        break
    end
    [~,c2] = chol(GradD - sum_Gammas);
end

% When no perturbation at all was used, we set it to -30
if T == -18; T = -30; end

if strcmp(cvx_status,'Inaccurate/Solved') == 1 || strcmp(cvx_status,'Solved') == 1
    % Lower bound
    low_bound = G0 + maximum;
    low_bound_Alg1 = real(D);
else
    % Failed lower bound
    MaxMinf = 99; Ppe = 0.5; Ptom = 0.5;
    low_bound = 0;
    low_bound_Alg1 = real(D);
end

fprintf('Bound on the min. tradeoff func = %d \n',low_bound);
fprintf('With Alg1: %d \n', low_bound_Alg1);
fprintf('-----------------------\nExecution completed\n');


%% Final bound for the min. tradeoff function

% Note that we consider EC at the Shannon limit. Later on, the inefficiency
% can be taken into account by multiplying the deltaEC by a factor f >= 0.

% Conditional probability p(z|x)
Pzx = zeros(1,16);
j = 0;
for ze=1:4
    for x=1:4
        j = j+1;
        Pro = @(g,t) g.*exp(-(abs(g.*exp(1i*t)-sqrt(eta).*amp(x)).^2)./(1 + eta*xi/2))./(pi*(1 + eta*xi/2));
        Pzx(j)=integral2(Pro,0,inf,pi*(2*ze-3)/4,pi*(2*ze-1)/4);       
    end
end

% Renormalization (actually only necessary when there's postselection)
p_passzx = sum(Pzx)/4;
Pzx      = Pzx/p_passzx;

% Mutual information of EC cost
MInfo = 0;
Pz    = zeros(1,4);
for ze=1:4
    for j=1:4
        Pz(ze) = Pz(ze) + Pzx(1,4*(ze-1) + j)*pA;
    end

    for x = 1:4
        MInfo = MInfo + pA*Pzx(4*(ze-1) + x)*log2(Pzx(4*(ze-1) + x)/Pz(ze));
    end
end

% Information leakage per signal
deltaEC = VNent(diag(Pz)) - MInfo;

% Asymptotic key rate
R_inf = low_bound - p_passzx*deltaEC;
fprintf('Reliable key rate: \n  R_inf >= %d \n------------\n',R_inf);
R_inf1 = real(D) - p_passzx*deltaEC;
fprintf('With algorithm 1: R_inf >= %d \n',R_inf1);


%% Results
AlgVars = [];
AlgVars.epsilon  = epsilon;
AlgVars.probdist = probdist;
AlgVars.dualdist = nu;
AlgVars.T        = T;
AlgVars.MaxMinf  = MaxMinf;
AlgVars.Ppe      = Ppe;
AlgVars.deltaEC  = deltaEC;
AlgVars.R_inf1   = R_inf1;
AlgVars.R_inf    = R_inf;


end
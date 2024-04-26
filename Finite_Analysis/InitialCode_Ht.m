%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code prepares all the parameters to obtain the curves
% for our paper with a heterodyne setup AND including rho_A. 
% IN
%  Length    - distance between Alice and Bob
%  amplitude - amplitude of Alice's coherent states
%  xi        - Excess noise at the opt. fibers
%  Ddelta    - coefficient Delta/delta of Bob's modules
%  GammaRaw  - Bob's operators for PE
%  GammaBar  - Orthonormalized set of GammaRaw
%  alphaBar  - Norm of the Orthogonal operators GammaRaw
%  LambdaA   - Matrix basis for Alice's tomography of her marginal
% OUT
%  amp      - Set of states employed by Alice 
%  pA       - vector of probabilities for the coherent states
%  gammaRaw - Primal variables of the SDP (Bob's measurement results)
%  R        - Set of region operators for Bob's key map
%  G        - Map of the coherent key generation
%  Z        - Pinching map for the key generation
%  Gr       - Facial reduction of G
%  Zr       - Idem for Z
%  gammaBar - Orthonormalized set of primal variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [InitialVars] = InitialCode_Ht(Length,amplitude,xi,BasisVars)


%% Load variables of the basis

Nc       = BasisVars.Nc;
Delta    = BasisVars.Delta;
Ddelta   = BasisVars.Ddelta;
GammaRaw = BasisVars.GammaRaw;
GammaBar = BasisVars.GammaBar;
alphaBar = BasisVars.alphaBar;
LambdaA  = BasisVars.LambdaA;

%% Parameters of the simulations

% Amplitude of the coherent state
amp = [amplitude, - amplitude, 1i*amplitude, - 1i*amplitude];

% probabilities for each state
pA = 0.25; % Note that all states have the same probability

% Transmittance 
eta = 10^(-0.02*Length);

% Definition of constraint operators and expectation values
gammaRaw = zeros(1,16*(Ddelta+1) + 12);

fprintf('Parameters defined \n');


%% Region operators for the key map

% Here we calculate the region operators

R = zeros(Nc+1,Nc+1,4);
for n = 0:Nc
    for m = n:Nc
        integ = @(g) (g.^(m+n+1)).*exp(- g.^2)/(pi*sqrt(factorial(n)*factorial(m)));
        radial = integral(integ,0,inf);
        ang = @(t) exp(1i.*t.*(n-m));
        
        R(n+1,m+1,1) = radial*integral(ang,-pi/4,pi/4);
        
        R(n+1,m+1,2) = radial*integral(ang,pi/4,3*pi/4);
        
        R(n+1,m+1,3) = radial*integral(ang,3*pi/4,5*pi/4);
        
        R(n+1,m+1,4) = radial*integral(ang,5*pi/4,7*pi/4);
    end
end
% Completition of the operators
for k = 1:4; R(:,:,k) = R(:,:,k) + triu(R(:,:,k),1)'; end

    
fprintf('Region operators with postselection defined \n');


%% Postprocessing map

% Square roots of the region operators
sqR = zeros(Nc+1,Nc+1,4);

for i = 1:4
    sq = sqrtm(R(:,:,i));
    
    % Here we make sure the operator is hermitian
    if ~ishermitian(sq); sq = (sq + sq')/2; end 
    sqR(:,:,i) = sq;
end


%% Postprocessing maps G and Z

% Maps G and Z
G = cell(1);
Z = cell(1,4);

% Definition of register z, which keeps the bits of the key, and the
%superoperators of G and Z
[Dim,~,~,] = size(GammaRaw);
K = zeros(Dim,Dim,4);
for i = 1:4
    K(:,:,i) = kron(eye(4),sqR(:,:,i));
    
    B = zeros(4,4);
    B(i,i) = 1;
    Z{1,i} = kron(kron(B,eye(4)),eye(Nc+1));
end

% Definition of the whole postprocessing map
G{1} = [K(:,:,1);K(:,:,2);K(:,:,3);K(:,:,4)];

% Reduced operators
[Gr,Zr] = FR(G,Z);

fprintf('Maps Gr and Zr defined \n');


%% Coefficients gammaRaw & gammaBar

% The expectation values of the constraints are given by the
% coefficients calld gammaRaw, while the gammaBars will define the same 
% expectation values after the Gram-Schmidt process

delta = Delta/Ddelta;
for r = 1:(Ddelta+1)
    for i = 0:3
        for j = 1:4
            % Definition of the conditional probability distribution for gammaRaw
            Pro = @(g,t) g.*exp(-(abs(g.*exp(1i*t)-sqrt(eta).*amp(j)).^2)./(1 + eta*xi/2))./(pi*(1 + eta*xi/2));
        
            if r == (Ddelta+1)
                gammaRaw(1,16*(r-1) + 4*i + j) = pA*integral2(Pro,Delta,inf,pi*(2*i-1)/4,pi*(2*i + 1)/4);
            else
                gammaRaw(1,16*(r-1) + 4*i + j) = pA*integral2(Pro,delta*(r-1),delta*r,pi*(2*i-1)/4,pi*(2*i + 1)/4);
            end
        end
    end
end

% Here the same for the constraint \rho_A
for x= 1:12
    MatrixA = LambdaA(:,:,x+4);
    ancilla = 0;
    for i = 1:4
        for j = 1:4
            ovrlp = exp(1i*imag(amp(i)*conj(amp(j))) - (abs(amp(i)-amp(j))^2)/2);
            ancilla = ancilla + MatrixA(j,i)*pA*ovrlp;
        end
    end
    gammaRaw(1,16*(Ddelta+1) + x) = ancilla;
end

% Construction of the gammaBars with a chain
gammaBar = zeros(1,16*(Ddelta+1) + 12);
for i = 1:(16*(Ddelta+1) + 12)   
    pj = 0;
    
    % Loop of the construction
    for j = 1:i-1
        pj = pj + gammaBar(1,j)*real(trace(GammaBar(:,:,j)*GammaRaw(:,:,i)));
    end
    
    gammaBar(1,i) = (gammaRaw(1,i) - pj)/alphaBar(1,i);
end

fprintf('Coefficients gammaBar defined \n');

%% Results
InitialVars = [];
InitialVars.amp      = amp;
InitialVars.pA       = pA;
InitialVars.eta      = eta;
InitialVars.gammaRaw = gammaRaw;
InitialVars.Gr       = Gr;
InitialVars.Zr       = Zr;
InitialVars.gammaBar = gammaBar;

fprintf('-----------------------\nInitial variables defined\n');
end

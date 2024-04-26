%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code prepares all the constraints to obtain the curves 
% for a heterodyne setup.
%
% IN
%  Nc - Value of the cutoff
%  Delta,delta - Choice for the modulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Parameters

% Cutoff and size of rho
Nc = 12;
Dim = 4*Nc + 4;

% Modules and postselection for the constraints
Delta = 0.9;
delta = 0.9;
if mod(Delta,delta) ~= 0
   error('2*Delta and delta should be multiples'); 
end
Ddelta = round(Delta/delta);

% Detection parameter
cut = 1e2;

% Basis of constraint operators
GammaRaw = zeros(Dim,Dim, (Ddelta +1)*16 + 12);


%% Definition of the constraint operators


for r = 1:(Ddelta+1)
    R = zeros(Nc+1,Nc+1,4);
    
    for n=0:Nc
        for m = n:Nc
            integ = @(g) (g.^(m+n+1)).*exp(- g.^2)/(pi*sqrt(factorial(n)*factorial(m)));
            if r == (Ddelta+1)
                radial = integral(integ,Delta,inf);
            else
                radial = integral(integ,delta*(r-1),delta*r);
            end        
            ang = @(t) exp(1i.*t.*(n-m));
        
            for k = 1:4
                angular = integral(ang, pi*(2*k-3)/4, pi*(2*k-1)/4);
                R(n+1,m+1,k) = radial*angular;
            end
        end
    end
    % Completition of the operators
    for k=1:4; R(:,:,k) = R(:,:,k) + triu(R(:,:,k),1)'; end
    
    for i = 0:3
        for j = 1:4
            A = zeros(4,4);
            A(j,j) = 1;
            GammaRaw(:,:,16*(r-1) + 4*i + j) = kron(A,R(:,:,i+1));
        end
    end
end

% Removal of unnecessary variables
clear r R n m integ radial ang k angular i j A;

fprintf('Region operators for the constraints defined \n');


%% GellMann Bases

% We need to define a matrix bases for subspaces A and B
% For this we use the GellMann matrices normalized to their dimensions

% GGM basis in A
LambdaA = zeros(4,4,16); % We need 15 matrices (plus the identity)

% Identity and diagonal matrices
LambdaA(:,:,1) = eye(4); 
LambdaA(:,:,2) = sqrt(2)*[1,0,0,0;0,-1,0,0;0,0,0,0;0,0,0,0];
LambdaA(:,:,3) = sqrt(2/3)*[1,0,0,0;0,1,0,0;0,0,-2,0;0,0,0,0];
LambdaA(:,:,4) = sqrt(1/3)*[1,0,0,0;0,1,0,0;0,0,1,0;0,0,0,-3];

% Symmetric matrices
LambdaA(:,:,5) = sqrt(2)*[0,1,0,0;1,0,0,0;0,0,0,0;0,0,0,0];
LambdaA(:,:,6) = sqrt(2)*[0,0,1,0;0,0,0,0;1,0,0,0;0,0,0,0];
LambdaA(:,:,7) = sqrt(2)*[0,0,0,1;0,0,0,0;0,0,0,0;1,0,0,0];
LambdaA(:,:,8) = sqrt(2)*[0,0,0,0;0,0,1,0;0,1,0,0;0,0,0,0];
LambdaA(:,:,9) = sqrt(2)*[0,0,0,0;0,0,0,1;0,0,0,0;0,1,0,0];
LambdaA(:,:,10) = sqrt(2)*[0,0,0,0;0,0,0,0;0,0,0,1;0,0,1,0];

% Antisymmetric matrices
LambdaA(:,:,11) = sqrt(-2)*[0,-1,0,0;1,0,0,0;0,0,0,0;0,0,0,0];
LambdaA(:,:,12) = sqrt(-2)*[0,0,-1,0;0,0,0,0;1,0,0,0;0,0,0,0];
LambdaA(:,:,13) = sqrt(-2)*[0,0,0,-1;0,0,0,0;0,0,0,0;1,0,0,0];
LambdaA(:,:,14) = sqrt(-2)*[0,0,0,0;0,0,-1,0;0,1,0,0;0,0,0,0];
LambdaA(:,:,15) = sqrt(-2)*[0,0,0,0;0,0,0,-1;0,0,0,0;0,1,0,0];
LambdaA(:,:,16) = sqrt(-2)*[0,0,0,0;0,0,0,0;0,0,0,-1;0,0,1,0];

fprintf('Basis LambdaA given by GGMs defined \n')

% GGM basis in B
LambdaB = zeros(Nc+1,Nc+1,(Nc+1)*(Nc+1));

% Diagonal GGMs
LambdaB(:,:,1) = eye(Nc + 1); % Identity matrix in B
for r = 1:Nc
    ancilla = zeros(Nc+1,Nc + 1);
    ancilla(r+1,r+1) = -r;
    for j = 1:r
        ancilla(j,j) = 1;
    end
    LambdaB(:,:,1 + r) = sqrt((Nc + 1)/(r*(r+1)))*ancilla;
end

% Symmetric GGMs
h = 0; % Auxiliary index to order the matrices
for i = 1:Nc
    for j = i+1:Nc+1
        h = h+1;
        ancilla = zeros(Nc+1,Nc + 1);
        ancilla(i,j) = 1;
        ancilla(j,i) = 1;
        LambdaB(:,:,Nc+1 + h) = ancilla*sqrt((Nc+1)/2);
    end
end

% Antisymmetric GGMs
h = 0;
for i = 1:Nc
    for j = i+1:Nc+1
        h = h + 1;
        ancilla = zeros(Nc+1,Nc + 1);
        ancilla(i,j) = -1i;
        ancilla(j,i) = 1i;
        LambdaB(:,:,Nc+1 + (Nc*(Nc+1)/2) + h) = ancilla*sqrt((Nc+1)/2);
    end
end

% Completion of the initial set of operators
for j = 1:12
    GammaRaw(:,:,16*(Ddelta+1) + j) = kron(LambdaA(:,:,j+4),eye(Nc+1));
end

% Removal of unnecessary variables
clear r ancilla j h i;

fprintf('Basis LambdaB for the complementary operators defined \n');

%% Gram-Schmidt process for the operators

% Here we use a Gram-Schmidt process to build an orthonormal set of
% matrices from the initial set of interval operators

% Complete basis of ALL matrices
BasisAB = zeros(Dim,Dim,Dim^2);

% GS process for the basis of constraint operators
GammaBar = zeros(Dim,Dim,16*(Ddelta+1) + 12);
alphaBar = zeros(1,16*(Ddelta+1) + 12); % 'norm' of the GammaBars

% Here we take the number of GammaRaws and the first GammaBar
[~,~,S] = size(GammaRaw);
alphaBar(1,1) = sqrt(trace(GammaRaw(:,:,1)*GammaRaw(:,:,1)'));
GammaBar(:,:,1) = GammaRaw(:,:,1)/alphaBar(1,1);
BasisAB(:,:,1) = GammaBar(:,:,1);

% General loop
for k=2:S
    GammaBar(:,:,k) = GammaRaw(:,:,k);
    for j = 1:(k-1)
        GammaBar(:,:,k) = GammaBar(:,:,k) - GammaBar(:,:,j)*real(trace(GammaBar(:,:,k)*GammaBar(:,:,j)')/trace(GammaBar(:,:,j)*GammaBar(:,:,j)'));
    end

    % Checkup of linear independence
    check1 = fix(GammaBar(:,:,k)*10^6)/(10^6);
	check2 = all(check1 == 0);
    if sum(check2) == Dim
        fprintf('  Lin. dep. matrix at k=%d \n',k);
    end

	% By eq. (30) of Winnick et al. the GammaBars need to be normalized to 1
    alphaBar(1,k) = sqrt(trace(GammaBar(:,:,k)*GammaBar(:,:,k)'));
    GammaBar(:,:,k) = GammaBar(:,:,k)/alphaBar(1,k); 
    % And finally we add the new operator to the basis
    BasisAB(:,:,k) = GammaBar(:,:,k);
end

% Removal of unnecessary variables
clear k j check1 check2;

fprintf('Operators GammaBar defined and orthonormal\n');

%% Completion of the basis

% For this task we need two indices - one that runs through the indices of
%the GGMs and another one over the operators GammaBar. The reason why
%we need two indices is because some GGMs will be linearly dep. on the
%previously defined matrices, so these must be discarded.

% Notice also that by the definition rho = sum \bar{gamma} \bar{Gamma}+...
%we are actually performing a sort of Bloch decomposition which means that 
%either the coefficients \bar{gamma} bear the normalization factor
%or the GammaBar's are normalized to the unit. The latter is the approach 
%taken here.

h = 0; % With this index we skip the linearly dependent GGMs
J = zeros(1,(Nc+1)^2); % This simplifies the detection of l.d. matrices
for i = 1:16 % Index for GGMs in A
    for j = 1:((Nc+1)*(Nc+1)) % Index for GGMs in B
        
        if mod(j,100)==0
            fprintf('j = %d\n',j);
        end
        
        % Shortcut of the GS process
        if isdiag(LambdaA(:,:,i)) == 0 && j > 1
            % Normalization to the unit
            projection = kron(LambdaA(:,:,i),LambdaB(:,:,j))/sqrt(Dim);

            % Here we make sure the matrix is hermitian
            if ~ishermitian(projection)
                projection = (projection + projection')/2;
            end
            
            %The matrix is now properly defined, and we add it to the basis
            BasisAB(:,:,(Nc^2 +2*Nc +1)*(i-1) + j + S - h) = projection;  
            continue    
        end
        
        % Sortcut for l.d. matrices
        if J(j) > 0
            fprintf('  Lin. dep. matrix at i=%d , j=%d, h=%d \n',i,j,h+1);
            h = h + 1;
            continue
        end 
        
        % Sum of projections
        pj = kron(LambdaA(:,:,i),LambdaB(:,:,j));
        for k=1:((Nc^2 +2*Nc +1)*(i-1) + j + (S-1) - h)
            pj = pj - real(trace(pj*BasisAB(:,:,k)')/trace(BasisAB(:,:,k)*BasisAB(:,:,k)'))*BasisAB(:,:,k);
        end      
        
        % Here we check if the matrix is l.d. which requires a truncation
        check1 = fix(pj*cut)/(cut);
        check2 = all(check1 == 0);
        if sum(check2) == Dim
            fprintf('  Lin. dep. matrix at i=%d , j=%d, h=%d \n',i,j,h+1);
            h = h + 1;
            J(j) = 1;
        else
            % As the matrix is not l.dep. we use it for the basis, which
            % requires the matrix to be normalized to the unit
            pj = pj/sqrt(trace(pj*pj'));
            
            % Here below we make sure the matrix is hermitian
            if ~ishermitian(pj); pj = (pj+pj')/2; end
            
            %The matrix is now properly defined, we add it to the basis
            BasisAB(:,:,(Nc^2 +2*Nc +1)*(i-1) + j + S - h) = pj;
        end
    end
end

fprintf('No of GammaBars: %d \n',S);
fprintf('Discarded matrices: %d \n', h);
if h ~= S
    error('Not all l.d. matrices were detected!');
end

% Removal of unnecessary variables
clear Dim cut S h i j projection LambdaB t pj k check1 check2;

fprintf('Complete matrix basis properly defined \n');          
fprintf('-----------------------\nExecution completed\n');

 function GradD = Grad_Ht(rho,Gr,Zr)

fprintf('Calculation of Grad(D) \n');

%% Action of channels

% Action of channel G and hermiticity
Gmap = applychan(Gr,rho,'kraus');
if ~ishermitian(Gmap); Gmap = (Gmap + Gmap')/2; end

% Check that Gmap is >=0
[~,check] = chol(Gmap);
if check ~=0
    error('Gmap is not > 0');
    %fprintf('Gmap is not >0. Numerical help used \n');
    %Gmap = Gmap*(1-1e-17)+eye(length(Gmap))*1e-17;
end

% Action of channel Z and hermiticity
Zmap = applychan(Zr,rho,'kraus');
if ~ishermitian(Zmap); Zmap = (Zmap + Zmap')/2; end

% Check that Zmap is >=0
[~, check] = chol(Zmap);
if check ~= 0
	error('Zmap is not > 0');
    %fprintf('Zmap is not > 0. Numerical help used \n');
    %Zmap = Zmap*(1-1e-15)+eye(length(Zmap))*1e-15;
end

%% Matrix Gradient for RelEnt

% Definition of adjoint maps
Gradj = cell(1);
Gradj{1} = Gr{1}';
Zradj = cell(1,4);
for i=1:4
    Zradj{i} = Zr{i}';
end

% Logarithms and correction of residuals
% Notice the ./log(2) to calculate the logm in base 2 instead of e
logG = logm(Gmap)./log(2);
if ~ishermitian(logG); logG = (logG + logG')/2; end
logZ = logm(Zmap)./log(2);
if ~ishermitian(logZ); logZ = (logZ + logZ')/2; end

% Portion corresponding to G
grad1 = applychan(Gradj,logG,'kraus') + ...
    applychan(Gradj,eye(size(Gradj{1},2)),'kraus');

% Portion corresponding to Z
grad2 = applychan(Zradj,logZ,'kraus') + ...
    applychan(Zradj,eye(size(Zradj{1},2)),'kraus');

% Gradient and correction
GradDT = grad1 - grad2;
% Here we take the transpose, which is the actual gradient, and 
% ensure the output is hermitian
GradD = transpose(GradDT);
if ~ishermitian(GradD); GradD = (GradD+GradD')/2; end

end

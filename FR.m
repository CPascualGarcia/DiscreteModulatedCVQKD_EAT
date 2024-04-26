function [Gr,Zr] = FR(G,Z)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function performs facial reduction
% for the maps G and Z
% Careful because the new map Z is NOT
% anymore a projector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Reduction of the map G

G0 = applychan(G,eye(size(G{1},2)),'kraus');
if ~ishermitian(G0); G0 = (G0 + G0')/2; end
[V,D] = eig(G0);

% Application of the reduction
ranG = rank(D);
if ranG < length(G0)
    V = V(:,end-ranG+1:end); % FR applied
    fprintf('Rank of G reduced from %d to %d \n',length(G0),ranG);
else
    V = eye(length(G0)); % No FR possible
end
Gr = cell(1);
Gr{1} = V'*G{1};


%% Reduction of the map Z

G1 = applychan(Gr,eye(size(G{1},2)),'kraus');
if ~ishermitian(G1); G1 = (G1 + G1')/2; end
Z0 = applychan(Z,V*G1*V','kraus');
if ~ishermitian(Z0); Z0 = (Z0 + Z0')/2; end
[U,T] = eig(Z0);

% Application of the reduction
ranZ = rank(T);
U = U(:,end-ranZ+1:end);
Zr = cell(1,4);
% We perform a rotation even in the event of no FR because this makes the
% distribution of eigenvalues more uniform - as no eigenvalues tend to
% zero, there is a smaller risk of having neg. eig. due to numerical
% imprecisions (and remind that the VNentropy is invariant under unitaries)
if ranZ < length(Z0)
    for i = 1:length(Z)
        Zr{i} = U'*Z{i}*G{1};
        fprintf('Rank of Z reduced from %d to %d \n',length(Z0),ranZ);
    end
else
    for i = 1:length(Z)
        Zr{i} = U'*Z{i}*G{1}; % No FR possible
    end
end


end
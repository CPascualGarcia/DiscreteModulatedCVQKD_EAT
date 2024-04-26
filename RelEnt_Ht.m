function D = RelEnt_Ht(rho,Gr,Zr)

% If the state is not >=0 we discard it for the
% minimization by setting the RelEnt. to the max.
[~,c] = chol(rho);
if c ~= 0
    D = 2;
    return
end

% Action of channel G, and hermicity
Gmap = applychan(Gr,rho,'kraus');
if ~ishermitian(Gmap); Gmap = (Gmap+Gmap')/2; end

% And an extra check to verify that the final G is >= 0
[~,check] = chol(Gmap);
if check ~= 0
    D = 2; % This ensures that this RelEnt will be neglected in the minimization
    return
end

% Action of channel Z
Zmap = applychan(Zr,rho,'kraus');

% Here we force Z's output to be hermitian
if ~ishermitian(Zmap); Zmap = (Zmap + Zmap')/2; end
[~, check] = chol(Zmap);
if check ~= 0
    D = 2;
    return
end

% Relative entropy
D = VNent(Zmap) - VNent(Gmap);
end



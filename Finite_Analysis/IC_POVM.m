
function [ICPOVM,gammaIC] = IC_POVM(amp,pA)

% 2D vectors that define the ICPOVM for qubits
maat = zeros(2,4);
maat(:,1) = [1;0];
maat(:,2) = [1;sqrt(2)]/sqrt(3);
maat(:,3) = [1;sqrt(2)*exp(1i*2*pi/3)]/sqrt(3);
maat(:,4) = [1;sqrt(2)*exp(1i*4*pi/3)]/sqrt(3);


% Here we use their outer products to create the ICPOVM vectors in 4D
phi = zeros(4,16);
jj = 0;
for i=1:4
    for j = 1:4
        jj = jj + 1;
        phi(:,jj) = kron(maat(:,i),maat(:,j));
    end
end

% ICPOVM vectors in 4D
ICPOVM = zeros(4,4,16);
summa = 0;
for i = 1:16
    ICPOVM(:,:,i) = 0.25*phi(:,i)*phi(:,i)';
    if ~ishermitian(ICPOVM(:,:,i)); ICPOVM(:,:,i)=(ICPOVM(:,:,i)+ICPOVM(:,:,i)')/2;end
    % Force PSDness if necessary
    while min(eig(ICPOVM(:,:,i)))<0
        ICPOVM(:,:,i)=ICPOVM(:,:,i)+1e-17*eye(4);
    end
    summa = summa + ICPOVM(:,:,i);
end


[~,~,S] = size(ICPOVM);
alpha = zeros(S,1);
alpha(1) = sqrt(trace(ICPOVM(:,:,1)*ICPOVM(:,:,1)'));
Bar(:,:,1) = ICPOVM(:,:,1)/alpha(1);
ICPOVMBar(:,:,1) = Bar(:,:,1);

% General loop
for k=2:S
    Bar(:,:,k) = ICPOVM(:,:,k);
    for j = 1:(k-1)
        Bar(:,:,k) = Bar(:,:,k) - Bar(:,:,j)*real(...
            trace(Bar(:,:,k)*Bar(:,:,j)')/trace(Bar(:,:,j)*Bar(:,:,j)'));
    end

    % Checkup of linear independence
    check1 = fix(Bar(:,:,k)*1e12)/1e12;
	check2 = all(check1 == 0);
    if sum(check2) == 4; fprintf('  Lin. dep. matrix at k=%d \n',k); end

	% By eq. (30) of Winnick et al. the GammaBars need to be normalized to 1
    alpha(k) = sqrt(trace(Bar(:,:,k)*Bar(:,:,k)'));
    Bar(:,:,k) = Bar(:,:,k)/alpha(k); 
    % And finally we add the new operator to the basis
    ICPOVMBar(:,:,k) = Bar(:,:,k);
end

% for k=1:S
%     if ~ishermitian(ICPOVMBar(:,:,k));
%         display(k);
%     end
% end


% Expectation values
gammaIC=zeros(S,1);
for x= 1:S
    MatrixA = ICPOVM(:,:,x);
    ancilla = 0;
    for i = 1:4
        for j = 1:4
            ovrlp = exp(1i*imag(amp(i)*conj(amp(j)))-(abs(amp(i)-amp(j))^2)/2);
            ancilla = ancilla+MatrixA(j,i)*pA*ovrlp;
        end
    end
    gammaIC(x)=real(ancilla);
end


end
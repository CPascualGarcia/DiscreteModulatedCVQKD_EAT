%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This code calculates the tolerance margins for the segments of the
%%% min-tradeoff function related to PE and tomography
%%%
%%% IN
%%%  Nrounds - No. of rounds employed for the key distillation
%%%  b       - Scaling of test rounds
%%%  Ppe     - Probability of having a test round employed for PE
%%%  Ppoint  - Prob. dist. that defines the primal point of the SDP
%%%  Dpoint  - Dual point for the lower bound of the SDP
%%%  Etom    - Upper bound for the deviations in Alice's marginal
%%%  EcPE    - Upper bound in the physical protocol for observing
%%%            unacceptable deviations in the PE section of the
%%%            min-tradeoff function for an honest implementation
%%% OUT
%%%  dTolPE  - Margin of deviations for the PE section of the min-tradeoff
%%%            function
%%%  dTolTom - Margin of deviations for the tomography section of the 
%%%            min-tradeoff function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [dTolPE,dTolTom] = DeltaTol(Nrounds,b,Ppe,Ppoint,Dpoint,EcPE,Etom)

%% Notation and preliminaries

n = Nrounds;
TOffVars = [Dpoint(1:end-16)/Ppe;Dpoint(end-15:end)/(1-Ppe)];
Maxf = max(TOffVars);
% Minf = min(TOffVars);
% error('Sanity check Maxminf=Maxf-Minf');

PiPE  = [Ppoint(1:end-16)*n^(-b);1-sum(Ppoint(1:end-16))*n^(-b)];
% MuPE  = Dpoint(1:end-16);
hPE   = [-Maxf+TOffVars(1:end-16);0]*(n^b -1);
DPE   = abs(max(hPE)-min(hPE));

PiTom = [Ppoint(end-15:end)*n^(-b);1-sum(Ppoint(end-15:end))*n^(-b)];
% MuTom = Dpoint(end-15:end);
hTom  = [-Maxf+TOffVars(end-15:end);0]*(n^b -1);
DTom  = abs(max(hTom)-min(hTom));


%% Margin for PE


Bpe = 0;
% g = zeros(length(MuTom),1);
% c = zeros(length(MuTom),1);
for i=1:length(PiTom)-1
    gi = PiPE(i)*(1-sum(PiPE(1:i)))/(1-sum(PiPE(1:i-1)));
    ci = hPE(i) - PiPE(i+1:end)'*hPE(i+1:end)/(1-sum(PiPE(1:i)));
%     g(i) = gi; c(i) = ci;
    Bpe = Bpe+gi*ci^2;
end

dTolPE = 2*sqrt(log2(n/EcPE)*Bpe/n) + 3*DPE*log2(2/EcPE)/n;


%% Margin for the tomography


Btom = 0;
% g = zeros(length(MuTom),1);
% c = zeros(length(MuTom),1);
for i=1:length(PiTom)-1
    gi = PiTom(i)*(1-sum(PiTom(1:i)))/(1-sum(PiTom(1:i-1)));
    ci = hTom(i) - PiTom(i+1:end)'*hTom(i+1:end)/(1-sum(PiTom(1:i)));
%     g(i) = gi; c(i) = ci;
    Btom = Btom+gi*ci^2;
end

dTolTom = 2*sqrt(log2(2*n/Etom)*Btom/n) + 3*DTom*log2(4/Etom)/n;



%% Improvement based on simplified bounds

hPE   = [TOffVars(1:end-16);0]*(n^b -1);
hTom  = [TOffVars(end-15:end);0]*(n^b -1);

Bpe = 0;
for i=1:length(PiTom)-1
    gi = PiPE(i)*(1-sum(PiPE(1:i)))/(1-sum(PiPE(1:i-1)));
    ci = hPE(i) - PiPE(i+1:end)'*hPE(i+1:end)/(1-sum(PiPE(1:i)));
    Bpe = Bpe+gi*ci^2;
end
dTolPE2 = 2*sqrt(log2(n/EcPE)*Bpe/n) + 3*DPE*log2(2/EcPE)/n;

Btom = 0;
for i=1:length(PiTom)-1
    gi = PiTom(i)*(1-sum(PiTom(1:i)))/(1-sum(PiTom(1:i-1)));
    ci = hTom(i) - PiTom(i+1:end)'*hTom(i+1:end)/(1-sum(PiTom(1:i)));
    Btom = Btom+gi*ci^2;
end
dTolTom2 = 2*sqrt(log2(2*n/Etom)*Btom/n) + 3*DTom*log2(4/Etom)/n;

if (dTolPE+dTolTom)>(dTolPE2+dTolTom2)
    dTolPE  = dTolPE2;
    dTolTom = dTolTom2;
end


end








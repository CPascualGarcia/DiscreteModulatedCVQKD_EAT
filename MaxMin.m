function [MaxMinf,Ppe,Ptom] = MaxMin(y)

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This code calculates the max. and min. bounds for the min-tradeoff
%%% function according to the dual maximization.
%%%
%%% IN: y - variables of the dual maximization
%%%        
%%% OUT: MaxMinf   - Difference between the bounds: Maxf - Minf
%%%      Ppe, Ptom - Renormalized probabilities for tomography and PE
%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Split the set of variables into PE variables and tomography variables
yPE  = y(1:end-16);
ytom = y(end-15:end);

% Extremal values of each set
Max1 = max([max(yPE),max(ytom)]); % max(yPE) usually
Max2 = min([max(yPE),max(ytom)]); % max(ytom)
Min1 = min([min(yPE),min(ytom)]); % min(yPE) usually
Min2 = max([min(yPE),min(ytom)]); % min(ytom)

% Check the indices of the extremals and calculate the optimal renormalized
% probability distributions
IndMax=find(y==Max1); IndMin=find(y==Min1);

% Both extremals belong to PE variables
if IndMax <= length(yPE) && IndMin <= length(yPE)
    cvx_precision best
    cvx_begin
        variable x
        minimize ((Max1-Min1)*inv_pos(x))
        subject to
            x <= 1
            x >= 0
            (1-x)*Max1 - Max2*x >= 0
            (1-x)*Min1 - Min2*x <= 0
    cvx_end
    Ppe = x; Ptom = 1-x;
    MaxMinf = (Max1-Min1)/Ppe;

% Both extremals belong to tomography variables
elseif IndMax > length(yPE) && IndMin > length(yPE)
    cvx_precision best
    cvx_begin
        variable x
        minimize ((Max1-Min1)*inv_pos(x))
        subject to
            x <= 1
            x >= 0
            (1-x)*Max1 - Max2*x >= 0
            (1-x)*Min1 - Min2*x <= 0
    cvx_end
    Ptom = x; Ppe = 1-x;
    MaxMinf = (Max1-Min1)/Ptom;
    
% The extremals belong to different sets
else
    cvx_precision best
    cvx_begin
        variable x
        minimize (Max1*inv_pos(x) - Min1*inv_pos(1-x))
        subject to
            x <= 1
            x >= 0
            (1-x)*Max1 - Max2*x >= 0 % x <= Max1/(Max1+Max2)
            x*Min1 - (1-x)*Min2 <= 0
	cvx_end
    
    if IndMax <= length(yPE) % The max. is related to PE rounds
        Ppe = x; Ptom = 1-x;
        MaxMinf = Max1/Ppe - Min1/Ptom;
    else % Max. related to tomography rounds
        Ptom = x; Ppe = 1-x;
        MaxMinf = Max1/Ptom - Min1/Ppe;
    end
end

end


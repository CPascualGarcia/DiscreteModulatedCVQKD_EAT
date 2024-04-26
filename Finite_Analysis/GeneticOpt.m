
function [] = GeneticOpt(xi,InitialVars,GeneticVars,Nr,NAME)

    % No. of vectors kappa called per iteration
    N_ELEMs = 100;

    % Load the variables of interest
    L        =    GeneticVars.L;
    Epsilons = GeneticVars.Epsilons;
    GammaRaw = GeneticVars.GammaRaw;

    % Load the matrix
    NAME1 = ['Matrices\BestRho_L' num2str(L) 'xi' num2str(xi*100) '.csv'];
    rho   = readmatrix(NAME1);

    % Evaluate the rates without perturbation
    AlgVars = Alg2_Ht_Genetic(GammaRaw,rho,xi,InitialVars,zeros(1,3));
    MaxMin0 = AlgVars.MaxMinf;
    
    [R0] = FKR_Ht_Genetic(Nr,AlgVars,GeneticVars);

    % Create the perturbations for the dual
    kappa = rand(N_ELEMs,3)/1e3;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% GENETIC SUBROUTINE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    N_ITER    = 5;
    ListRates = zeros(N_ITER,1);
    ListKappas = zeros(N_ITER,3);
    for k=1:N_ITER

        if Nr>6e12 % No genetic optimization necessary
            ListRates(1)   = R0;
            ListKappas(1,:) = zeros(1,3);
            break
        end

        fprintf('ROUND %d_________________________________\n',k);
        Results = [R0,MaxMin0];

        for z=1:N_ELEMs
            if mod(z,5)==0
                fprintf('--------------------------------- z: %d\n',z);
            end

            % Evaluate the key rates for the round
            AlgVars = Alg2_Ht_Genetic(GammaRaw,rho,xi,InitialVars,kappa(z,:));
            MaxMin  = AlgVars.MaxMinf;

            % Fitness function
            [R_fin] = FKR_Ht_Genetic(Nr,AlgVars,GeneticVars);

            % Save the results
            Results = [Results;R_fin,MaxMin];

        end

        % Include the original bound
        kappa = [zeros(1,3);kappa];

        %%% SELECTION
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Save the best outcomes
        ListRates(k)   = max(Results(:,1));
        AUX            = Results(:,1)==max(Results(:,1));
        fprintf('Best result of the round: %d \n',ListRates(k));
        if sum(AUX)>1
            tmp      = find(AUX,1,"first");
            AUX      = logical(zeros(length(AUX),1));
            AUX(tmp) = true; 
        end
        ListKappas(k,:) = kappa(AUX,:);

        % Eliminate the worst-performing values
        Ind1 = Results(:,1)>=quantile(Results(:,1),0.1);
        Ind2 = Results(:,2)<=quantile(Results(:,2),0.95);
        Ind  = Ind2(Ind1);


        % Filter the resulting vectors
        R1     = Results(Ind1,:);
        Kappa1 = kappa(Ind1,:);
        R2     = R1(Ind,:);
        Kappa2 = Kappa1(Ind,:);

        % Create the new population
        newGen = zeros(N_ELEMs,3);

        % ELISTISM - directy choose the 10% best performing
        % vectors and include them in the new population
        tmp  = R2(:,1)       >=quantile(Results(:,1),0.9);
        newGen(1:nnz(tmp),:) = Kappa2(tmp,:);

        for JJ=1:N_ELEMs-nnz(tmp)
            Pick1 = randi([1 length(Kappa2)]);
            Pick2 = randi([1 length(Kappa2)]);
            
            parent1 = Kappa2(Pick1,:);
            score1  = R2(Pick1,1);
    
            parent2 = Kappa2(Pick2,:);
            score2  = R2(Pick2,1);
    
            offsp = Mate(parent1,score1,parent2,score2);
            newGen(nnz(tmp)+JJ,:) = offsp;
        end

        % Take the new vectors and start over
        kappa = newGen;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% END OF THE GENETIC ALGORITHM
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Extract the final key rate
    FinalKey   = max(ListRates);
    FinalKappa = ListKappas(ListRates==max(ListRates),:);
    if size(FinalKappa,1)>1
        FinalKappa = FinalKappa(1,:);
    end

    FILE2 = fopen(NAME,'a');
    fprintf(FILE2,'%d, %.6e, ',L, FinalKey);
    fprintf(FILE2,'%.16e, ',FinalKappa(1:end-1));
    fprintf(FILE2,'%.16e \n',FinalKappa(length(FinalKappa)));
    fclose(FILE2);
end

function offsp = Mate(parent1,score1,parent2,score2)

offsp = zeros(1,length(parent1));

for k =1:length(parent1)

    mutation = rand(1);
    A = 1e-8;
    B = 1e-1;

    if mutation>0.9
        offsp(k) = A + (B-A)*rand(1);
    elseif score1>0 || score2>0
        if binornd(1,score1/(score1+score2))<=0.5
            offsp(k) = parent1(k);
        else
            offsp(k) = parent2(k);
        end
    else
        if mutation<0.45
            offsp(k) = parent1(k);
        else
            offsp(k) = parent2(k);
        end
    end
end
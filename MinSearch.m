function [lambda_m,after] = MinSearch(rho,DelRho,Gr,Zr)

% Preliminaries
fun = @(t) RelEnt_Ht(rho+t*DelRho,Gr,Zr);
before = RelEnt_Ht(rho,Gr,Zr);
upb = 1;
p = 0;

% Minimization
opts = optimset('TolX',1e-13);
while p == 0
    lambda_m = fminbnd(fun,0,upb,opts);
    after = RelEnt_Ht(rho + lambda_m*DelRho,Gr,Zr);
    if after > before
        upb = upb*0.75;
    else
        break
    end
end
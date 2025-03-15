function F = surf_modelfun(x, Akobs, Aksyn, indextable, Thetaij, v, R, casemodel)

Nk = size(indextable, 1); % number of data
Ns = max(indextable(:, 1)); % number of sensor
Nb = max(indextable(:, 2)); % number of sensor

f = cell(Nk, 1);

for k = 1:Nk
    i = indextable(k, 1);
    j = indextable(k, 2);
    theta = Thetaij(i, j); % incidence angle in degree
    Si = x(i);

    if casemodel == "directionality"
        TR = x(Ns+1);
%         b = x(Ns+2);
%         Ak = Aksyn(k)*Si*exp(-a*abs(theta)^b);
        Ak = Aksyn(k)*Si*k_surfeffect(TR, theta, v, R);
        
    elseif casemodel == "mixed"
        Tj = x(Ns+j);
        TR = x(Ns+Nb+1);
%         b = x(Ns+Nb+2);
%         Ak = Aksyn(k)*Si*Tj*exp(-a*abs(theta)^b);
        Ak = Aksyn(k)*Si*Tj*k_surfeffect(TR, theta, v, R);
        
    else
        error("case model is unknown.");
    end

    f{k} = Ak - Akobs(k);
end

F = horzcat(f{:});

end
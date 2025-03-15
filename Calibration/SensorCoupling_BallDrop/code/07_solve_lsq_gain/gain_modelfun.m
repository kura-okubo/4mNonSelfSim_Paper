function F = gain_modelfun(x, Akobs, Aksyn, indextable)

Nk = size(indextable, 1); % number of data
Ns = max(indextable(:, 1)); % number of sensor

f = cell(Nk, 1);

for k = 1:Nk
    i = indextable(k, 1);
    j = indextable(k, 2);
    Si = x(i);
    Tj = x(Ns+j);
    Ak = Aksyn(k)*Si*Tj;
    f{k} = Ak - Akobs(k);
%     f{k} = sqrt(w(k)) * (Ak - Akobs(k)); %with weight 
end

F = horzcat(f{:});

end

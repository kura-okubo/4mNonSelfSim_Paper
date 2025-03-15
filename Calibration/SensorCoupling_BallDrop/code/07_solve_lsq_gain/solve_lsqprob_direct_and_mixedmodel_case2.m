% compute non-linear least square method for sensor calibration
% extended to directionality and mixed
clear all;
set(0,'DefaultTextFontsize',16, ...
    'DefaultTextFontname','Arial', ...
    'DefaultTextFontWeight','normal', ...
    'DefaultTextFontname','Arial', ...
    'DefaultAxesFontsize',16, ...
    'DefaultAxesFontname','Arial', ...
    'DefaultLineLineWidth', 1.5)

set(0,'defaulttextinterpreter','none')
figdir="./figure";
if ~exist(figdir) mkdir(figdir); end

%% Parameters
mingain = 0.0; % lower bound of gain
maxgain = 10.0; % upper bound of gain
minbdcoef = -10.0; % lower bound of gain
maxbdcoef = 10.0; % upper bound of gain
minTR = 0.0; % lower bound of AEsurf coef. TR
maxTR = 10.0; % upper bound of AEsurf coef. TR
% minb =-10.0; % lower bound of directionality coef. b
% maxb =10.0; % upper bound of directionality coef. b
initgain = 0.6; % initial gain
initbdcoef = 1.0; % initial ball-drop constant
initTR = 5e-6; %[s] initial constant of rise timeof P pulse
% initb = 1.0;

v = 6200; %[m/s] p wave velocity of rock medium
R = 6.35e-3; %8.0e-3%6.5e-3; %[m] Radius of AE sensor surface

%% read PAZ file to read the gain
% PAZ = load("../../data/AE_resp_dataandcoef_fronttop.mat");
% preamp = 40; % [dB] amplification of preamp included in the LDV calibration
% PAZ_gain = PAZ.k/PAZ.u_normfact/(10^(preamp/20)); % [V/(m/s)] gain obtained from LDV calibration

%%
datarootdir = "../../data/DATA_surfaceeffect_Aij_Case2";

Aijobs  = importdata(fullfile(datarootdir, 'Aijobs.csv'));
Aijsyn  = importdata(fullfile(datarootdir, 'Aijsyn.csv'));
Thetaij = importdata(fullfile(datarootdir, 'Thetaij.csv'));
Alphaij = importdata(fullfile(datarootdir, 'Alphaij.csv'));
rij     = importdata(fullfile(datarootdir, 'rij.csv'));
lij     = importdata(fullfile(datarootdir, 'lij.csv'));
lrij     = importdata(fullfile(datarootdir, 'loverrij.csv'));
Caseij  = readtable(fullfile(datarootdir, 'Caseij.csv'), 'FileType', 'delimitedtext', 'Delimiter',...
    ',', 'ReadVariableNames', false);

%% Check consistency of matrix
maxdist = 400; %mm
assert(max(max(lij)) < maxdist);
maxtheta = max(max(Thetaij));

%% Vectorize Matrix
Nd = nnz(lij); % number of ball drop data used to inversion
[Ns, Nb] = size(Caseij); % number of sensors and balldrop event

Akobs = zeros(Nd, 1); % vectorized amplitude
Aksyn = zeros(Nd, 1);
Slin  = zeros(Ns, 1); % gain inverted from linear lsq
Tlin  = ones(Nb, 1); % source coefficient is fixed to one for the case of linear lsq
Snonlin = zeros(Ns, 1); % gain inverted from non-linear lsq
Tnonlin = zeros(Nb, 1); % source coefficient for the case of non-linear lsq
indextable = zeros(Nd,2); % index table for k and (i, j)

kcount = 0;
for i = 1:Ns
    for j = 1:Nb
        if lij(i, j) ~= 0
            kcount = kcount + 1;
            Akobs(kcount) = Aijobs(i, j);
            Aksyn(kcount) = Aijsyn(i, j);
            indextable(kcount, :) = [i, j];
        end
    end
end

%%
% test gain_modelfun for mixed model
x0dir = zeros(Ns+1, 1);
x0dir(1:Ns) = initgain;
x0dir(Ns+1) = initTR;
% x0dir(Ns+2) = initb;
F = surf_modelfun(x0dir, Akobs, Aksyn, indextable, Thetaij, v, R, "directionality");
% F = dir_modelfun(x0dir, Akobs, Aksyn, indextable, Thetaij, "directionality");
norm(F)
%%
x0mix = zeros(Ns+Nb+1, 1);
x0mix(1:Ns) = initgain;
x0mix(Ns+1:Ns+Nb) = initbdcoef;
x0mix(Ns+Nb+1) = initTR;
% x0mix(Ns+Nb+2) = initb;

% F = dir_modelfun(x0mix, Akobs, Aksyn, indextable, Thetaij, "mixed");
F = surf_modelfun(x0mix, Akobs, Aksyn, indextable, Thetaij, v, R, "mixed");
norm(F)
%% Case 3. model with directionality
lbdir = zeros(Ns+1, 1);
ubdir = zeros(Ns+1, 1);
lbdir(1:Ns) = mingain;
ubdir(1:Ns) = maxgain;
lbdir(Ns+1) = minTR;
ubdir(Ns+1) = maxTR;

% fdir= @(x) dir_modelfun(x, Akobs, Aksyn, indextable, Thetaij, "directionality");
fdir= @(x) surf_modelfun(x, Akobs, Aksyn, indextable, Thetaij, v, R, "directionality");

options = optimoptions('lsqnonlin','OptimalityTolerance', 1e-12, 'StepTolerance', 1e-12, ...
    'FunctionTolerance', 1e-12, 'Display', 'iter', 'Algorithm', 'trust-region-reflective'); %'levenberg-marquardt');

[xdir,resnormdir,residualdir,exitflagdir,outputdir] = lsqnonlin(fdir,x0dir,lbdir,ubdir, options);

%% Case 4. mixed model

lbmix = zeros(Ns+Nb+1, 1);
ubmix = zeros(Ns+Nb+1, 1);
lbmix(1:Ns) = mingain;
ubmix(1:Ns) = maxgain;
lbmix(Ns+1:Ns+Nb) = minbdcoef;
ubmix(Ns+1:Ns+Nb) = maxbdcoef;
lbmix(Ns+Nb+1) = minTR;
ubmix(Ns+Nb+1) = maxTR;

% fmix = @(x) dir_modelfun(x, Akobs, Aksyn, indextable, Thetaij, "mixed");
fmix = @(x) surf_modelfun(x, Akobs, Aksyn, indextable, Thetaij, v, R, "mixed");

options = optimoptions('lsqnonlin','OptimalityTolerance', 1e-12, 'StepTolerance', 1e-12, ...
    'FunctionTolerance', 1e-12, 'Display', 'iter', 'Algorithm', 'trust-region-reflective'); %'levenberg-marquardt');

[xmix,resnormmix,residualmix,exitflagmix,outputmix] = lsqnonlin(fmix,x0mix,lbmix,ubmix, options);

%% Plot result
fig = figure(1);
fig.Units = 'normalized';
fig.Position = [0 1 0.4 0.7];
clf(fig,'reset'); cla(fig,'reset'); hold on;clf

th = linspace(1, 85,101);
% ydir = exp(-xdir(Ns+1).*th .^xdir(Ns+2));
% ymix = exp(-xmix(Ns+Nb+1).*th .^xmix(Ns+Nb+2));
ydir = zeros(length(th), 1);
ymix = zeros(length(th), 1);
for i = 1:length(th)
    ydir(i) = k_surfeffect(xdir(Ns+1), th(i), v, R);
    ymix(i) = k_surfeffect(xmix(Ns+Nb+1), th(i), v, R);
end

subplot(5, 1, 1); hold on;
plot(xdir(1:Ns), 'k+-');
plot(xmix(1:Ns), 'ro-');
xlim([1, Ns]);
ylim([0, 2.0]);
title("Si")

box on;

subplot(5, 1, 2); hold on;
plot(xmix(Ns+1:Ns+Nb), 'ko-');
xlim([1, Nb]);
ylim([0, 2.0]);
title("Tj")
box on;

subplot(5, 1, 3); hold on;
plot(th, ydir, 'k+-');
plot(th, ymix, 'ro-');
title("beta aperture effect")
% xlim([0, 90]);
% ylim([0.0, 1.0]);
box on;

subplot(5, 1, 4); hold on;
plot(1, xdir(Ns+1), 'ks');
plot(2, xmix(Ns+Nb+1), 'rv');
xlim([0.0, 3]);
title(sprintf("TR aperture effect: %4.2f us average",...
    (xmix(Ns+Nb+1) + xdir(Ns+1))*1e6/2))
box on;


subplot(5, 1, 5); hold on;
plot(residualdir, 'k-');
plot(residualmix, 'r-');
% ylim([0.0, 1.5]);
title("Redidual")
box on;

fprintf("dir residual, mix residual = %f, %f\n", resnormdir, resnormmix);

%% Post-analysis
k_highres = find(abs(residualdir)>2*std(residualdir));
ij_highres = indextable(k_highres, :);

%% Compute AIC

AICdir = Nd*log(resnormdir/Nd) + 2*(Ns+1); % updated 2024.1.24; removed Nd*(log(2*pi) + 1)
AICmix = Nd*log(resnormmix/Nd) + 2*(Ns+Nb+1);

fprintf("AIC directionality model: %f, AIC mix model: %f\n", AICdir, AICmix);

%% Compute AICc
kdir = Ns+1;
kmix = Ns+Nb+1;
AICcdir = AICdir + (2*kdir*(kdir+1))/(Nd - kdir - 1);
AICcmix = AICmix + (2*kmix*(kmix+1))/(Nd - kmix - 1);

fprintf("AICc directionality model: %f, AICc mix model: %f\n", AICcdir, AICcmix);

%% Save results into csv
Akest_dir = zeros(Nd, 1);
Akest_mix= zeros(Nd, 1);
Theta_k = zeros(Nd, 1);
Alpha_k = zeros(Nd, 1);
r_k = zeros(Nd, 1);
l_k = zeros(Nd, 1);
lr_k = zeros(Nd, 1);
Caseid_k = cell(Nd, 1);

for k = 1:Nd
    i = indextable(k, 1);
    j = indextable(k, 2);
    Sidir = xdir(i);
    Simix = xmix(i);
    Tjmix = xmix(Ns+j);
    theta = Thetaij(i, j);
%     alphadir = exp(-xdir(Ns+1) .* abs(theta) .^ xdir(Ns+2));
%     alphamix = exp(-xmix(Ns+Nb+1) .* abs(theta) .^ xmix(Ns+Nb+2));
    Akest_dir(k) = Akobs(k)/Sidir/k_surfeffect(xdir(Ns+1), theta, v, R);
    Akest_mix(k) = Akobs(k)/Simix/Tjmix/k_surfeffect(xmix(Ns+Nb+1), theta, v, R);
    Theta_k(k) =  Thetaij(i, j);
    Alpha_k(k) =  Alphaij(i, j);
    r_k(k) =  rij(i, j);
    l_k(k) =  lij(i, j);
    lr_k(k) =  lrij(i, j);
    Caseid_k(k) =  table2array(Caseij(i,j));
end

datadir = "../../data/DATA_surfaceeffect_Aij_Case2";
writematrix(indextable, fullfile(datadir, "directionality_indextable_mix.csv"));
writematrix([Akest_dir Akest_mix], fullfile(datadir, "directionality_Akest_mix.csv"));
writematrix([xdir(1:Ns) xmix(1:Ns)], fullfile(datadir, "directionality_Si_mix.csv"));
writematrix([xmix(Ns+1:Ns+Nb)], fullfile(datadir, "directionality_Tj_mix.csv"));
writematrix([xdir(Ns+1), xmix(Ns+Nb+1)], fullfile(datadir, "directionality_TR_mix.csv"));
writematrix([r_k l_k lr_k Theta_k Alpha_k], fullfile(datadir, "directionality_distandtheta_mix.csv"));
writematrix([residualdir' residualmix'], fullfile(datadir, "directionality_residual_mix.csv"));
writematrix([AICdir, AICmix, AICcdir, AICcmix], fullfile(datadir, "directionality_AIC_mix.csv"));
writematrix([xdir(Ns+1), xmix(Ns+Nb+1)], fullfile(datadir, "directionality_TR_mix.csv"));
writecell(Caseid_k, fullfile(datadir, "directionality_caseid.csv"));

%% Another plot
fig = figure(2);clf;
fig.Units = 'normalized';
fig.Position = [0.5 1 0.4 0.6];
clf(fig,'reset'); cla(fig,'reset'); clf;

subplot(1, 1, 1); hold on;
plot(abs(Alpha_k), Aksyn.*r_k, "+");
plot(abs(Alpha_k), Akest_dir.*r_k, "v");
plot(abs(Alpha_k), Akest_mix.*r_k, "o");
box on;
daspect([5e3, 1, 1])

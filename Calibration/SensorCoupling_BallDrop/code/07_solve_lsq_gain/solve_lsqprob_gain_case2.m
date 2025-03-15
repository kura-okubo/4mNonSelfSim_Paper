% compute non-linear least square method for sensor calibration
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
mingaincoef = 0.0; % lower bound of gain
maxgaincoef = 10.0; % upper bound of gain
minbdcoef = -10.0; % lower bound of Tj
maxbdcoef = 10.0; % upper bound of Tj
initgain = 0.3; % initial coefficient of gain
initbdcoef = 1.0; % initial ball-drop constant

%% read PAZ file to read the gain
% PAZ = load("../../data/AE_resp_dataandcoef_fronttop.mat");
% preamp = 40; % [dB] amplification of preamp included in the LDV calibration
% PAZ_gain = PAZ.k/PAZ.u_normfact/(10^(preamp/20)); % [V/(m/s)] gain obtained from LDV calibration

%%
% datarootdir = "../../data/DATA_surfaceeffect_Aij";
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
% test gain_modelfun
x0 = zeros(Ns+Nb, 1);
x0(1:Ns) = initgain;
x0(Ns+1:Ns+Nb) = initbdcoef;

F = gain_modelfun(x0, Akobs, Aksyn, indextable);

norm(F)

%% Case 1. Linear least square problem
C = zeros(Nd, Ns); % Ndata x Ngain, where Tj = 1.0

for k = 1:Nd
    i = indextable(k, 1);
    j = indextable(k, 2);
    assert(Aksyn(k) == Aijsyn(i, j));
    C(k, i) = Aksyn(k);
end

d = Akobs;
A = []; b = []; Aeq = []; beq = [];
lblin = zeros(Ns, 1);
ublin = zeros(Ns, 1);
lblin(1:Ns) = mingaincoef;
ublin(1:Ns) = maxgaincoef;
[xlin,resnormlin,residuallin,exitflaglin,outputlin,lambdalin] = lsqlin(C,d,A,b,Aeq,beq,lblin,ublin);

%% Case 2. Non-linear least square problem
lbnonlin = zeros(Ns+Nb, 1);
lbnonlin(1:Ns) = mingaincoef;
ubnonlin(1:Ns) = maxgaincoef;
lbnonlin(Ns+1:Ns+Nb) = minbdcoef;
ubnonlin(Ns+1:Ns+Nb) = maxbdcoef;

%% Case 2. Non-linear least square problem
options = optimoptions('lsqnonlin','OptimalityTolerance', 1e-12, 'StepTolerance', 1e-12, ...
    'FunctionTolerance', 1e-12, 'Display', 'iter', 'Algorithm', 'trust-region-reflective'); %'levenberg-marquardt');

fnonlin = @(x) gain_modelfun(x, Akobs, Aksyn, indextable);

[xnonlin,resnormnonlin,residualnonlin,exitflagnonlin,outputnonlin] = lsqnonlin(fnonlin,x0,lbnonlin,ubnonlin, options);

%% Plot result
fig = figure(1);
fig.Units = 'normalized';
fig.Position = [0 1 0.4 0.7];
clf(fig,'reset'); cla(fig,'reset'); hold on;clf

subplot(3, 1, 1); hold on;
plot(xlin(1:Ns), 'k+-');
plot(xnonlin(1:Ns), 'ro-');
xlim([1, Ns]);
ylim([0, 2.0]);
box on;

subplot(3, 1, 2); hold on;
plot(ones(Nb), 'k+-');
plot(xnonlin(Ns+1:Ns+Nb), 'ro-');
xlim([1, Nb]);
ylim([0.0, 2.0]);
box on;

subplot(3, 1, 3); hold on;
plot(residuallin, 'k-');
plot(residualnonlin, 'r-');
% ylim([0.0, 1.5]);
box on;

fprintf("lin residual, nonlin-residual = %f, %f\n", resnormlin, resnormnonlin);

%% Post-analysis
k_highres = find(abs(residuallin)>2*std(residuallin));
ij_highres = indextable(k_highres, :);

%% Compute AIC

AIClin = Nd*log(resnormlin/Nd) + 2*(Ns);  % updated 2024.1.24; removed Nd*(log(2*pi) + 1)
AICnonlin = Nd*log(resnormnonlin/Nd) + 2*(Ns+Nb);

fprintf("AIC linear: %f, AIC nonlinear: %f\n", AIClin, AICnonlin);

%% Compute AICc
klin = Ns;
knonlin = Ns+Nb;
AICclin = AIClin + (2*klin*(klin+1))/(Nd - klin - 1);
AICcnonlin = AICnonlin + (2*knonlin*(knonlin+1))/(Nd - knonlin - 1);

fprintf("AICc directionality model: %f, AICc mix model: %f\n", AICclin, AICcnonlin);

%% Save results into csv
Akest_lin = zeros(Nd, 1);
Akest_nonlin = zeros(Nd, 1);
Theta_k = zeros(Nd, 1);
Alpha_k = zeros(Nd, 1);
r_k = zeros(Nd, 1);
l_k = zeros(Nd, 1);
lr_k = zeros(Nd, 1);
Caseid_k = cell(Nd, 1);

for k = 1:Nd
    i = indextable(k, 1);
    j = indextable(k, 2);
    Silin = xlin(i);
    Sinonlin = xnonlin(i);
    Tjnonlin = xnonlin(Ns+j);
    Akest_lin(k) = Akobs(k)/Silin/1.0;
    Akest_nonlin(k) = Akobs(k)/Sinonlin/Tjnonlin;
    Theta_k(k) =  Thetaij(i, j);
    Alpha_k(k) =  Alphaij(i, j);
    r_k(k) =  rij(i, j);
    l_k(k) =  lij(i, j);
    lr_k(k) =  lrij(i, j);
    Caseid_k(k) =  table2array(Caseij(i,j));
end

datadir = "../../data/DATA_surfaceeffect_Aij_Case2";
writematrix(indextable, fullfile(datadir, "directionality_indextable.csv"));
writematrix([Akest_lin Akest_nonlin], fullfile(datadir, "directionality_Akest.csv"));
writematrix([xlin(1:Ns) xnonlin(1:Ns)], fullfile(datadir, "directionality_Si.csv"));
writematrix([xnonlin(Ns+1:Ns+Nb)], fullfile(datadir, "directionality_Tj.csv"));
writematrix([r_k l_k lr_k Theta_k Alpha_k], fullfile(datadir, "directionality_distandtheta.csv"));
writematrix([residuallin residualnonlin'], fullfile(datadir, "directionality_residual.csv"));
writematrix([AIClin, AICnonlin, AICclin, AICcnonlin], fullfile(datadir, "directionality_AIC.csv"));
writecell(Caseid_k, fullfile(datadir, "directionality_caseid.csv"));

%% Another plot
fig = figure(2);
fig.Units = 'normalized';
fig.Position = [0.5 1 0.4 0.6];
clf(fig,'reset'); cla(fig,'reset'); clf;

subplot(1, 1, 1); hold on;
plot(abs(Alpha_k), Aksyn.*r_k, "+");
plot(abs(Alpha_k), Akest_lin.*r_k, "v");
plot(abs(Alpha_k), Akest_nonlin.*r_k, "o");
xvec = linspace(0, 90, 51);
plot(xvec, 0.036*cosd(xvec), "b--");
box on;
daspect([4e3, 1, 1])
% xlim([50, 85])
% ylim([0, 0.07])

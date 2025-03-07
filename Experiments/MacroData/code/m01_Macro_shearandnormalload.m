% Script to plot the evolution of slip during the stick-slip events
clear all;
set(0,'DefaultTextFontsize',16, ...
    'DefaultTextFontname','Arial', ...
    'DefaultTextFontWeight','normal', ...
    'DefaultTextFontname','Arial', ...
    'DefaultAxesFontsize',16, ... 
    'DefaultAxesFontname','Arial', ...
    'DefaultLineLineWidth', 1.0)

figdir="../figure";
if ~exist(figdir) mkdir(figdir); end
addpath("../../../utils/matlabcode_biax_v03");

%% Load cases

M = struct();

for expr_id = [86, 87, 93, 94] %, 99] % 93-94 same condition with 86-87 with GP
    expr_id
    % pname = sprintf('/Volumes/NIEDSSD01/FB03data/fb03-%03d/biax/', expr_id);
    pname = sprintf('/Volumes/Okuboetal2025_masterHDD/FB03data/fb03-%03d/biax/', expr_id);
    runID = sprintf('fb03-%03d', expr_id);

    pname, runID
    [tmat,DX,NP,SS,FCM]=fmbiaxmacroS3(pname,runID);

    M.(sprintf("FB03_%03d", expr_id)).tmat=tmat;
    M.(sprintf("FB03_%03d", expr_id)).DX=DX;
    M.(sprintf("FB03_%03d", expr_id)).NP=NP;
    M.(sprintf("FB03_%03d", expr_id)).SS=SS;
    M.(sprintf("FB03_%03d", expr_id)).FCM=FCM;

end

%
save("../data/MacroData_raw.mat", "M"); 
%%
% load("../data/MacroData_raw.mat", "M");

%% Plot Main stick-slip experiments
fig = figure(1);
fig.Units = 'point';
fig.Position = [0 500 800 450];
clf(fig,'reset'); cla(fig,'reset'); hold on; box on; 

tshift = 20.2;
plot(M.FB03_087.tmat- tshift, M.FB03_087.SS, "k-", "DisplayName", "FB03-087");
plot(M.FB03_094.tmat, M.FB03_094.SS, "DisplayName", "FB03-094");

xlim([0, 250]);
ylim([0.5, 0.7]);
xlabel("Time [s]");
ylabel("Shear stress [MPa]");

legend('Location', 'northeast');
exportgraphics(gcf, sprintf("../figure/Macro_shear_stress.png"), "Resolution",300);


%%
fig = figure(2);
fig.Units = 'point';
fig.Position = [0 500 800 450];
clf(fig,'reset'); cla(fig,'reset'); hold on; box on; 

tshift = 20.2;
plot(M.FB03_087.tmat- tshift, M.FB03_087.FCM, "k-", "DisplayName", "FB03-087");
plot(M.FB03_094.tmat, M.FB03_094.FCM, "DisplayName", "FB03-094");

xlim([0, 250]);
ylim([0.3, 0.45]);
xlabel("Time [s]");
ylabel("Friction Coefficient");

legend('Location', 'northeast');
exportgraphics(gcf, sprintf("../figure/Macro_frictioncoef.png"), "Resolution",300);


%%
fig = figure(3);
fig.Units = 'point';
fig.Position = [0 500 800 450];
clf(fig,'reset'); cla(fig,'reset'); hold on; box on; 

tshift = 20.2;
h1 = plot(M.FB03_087.tmat- tshift, M.FB03_087.NP, "k-");
h2 = plot(M.FB03_094.tmat, M.FB03_094.NP, "r-", "DisplayName", "FB03-094");

% Plot mean value
h3 = plot(M.FB03_087.tmat- tshift, mean(M.FB03_087.NP, 2), "k--", "LineWidth", 3);
h4 = plot(M.FB03_094.tmat, mean(M.FB03_094.NP, 2), "r--", "LineWidth", 2);

xlim([0, 250]);
ylim([0.5, 4]);
xlabel("Time [s]");
ylabel("Normal stress [MPa]");

legend([h1(1), h2(1), h3(1), h4(1)], 'FB03-087', 'FB03-094', 'mean', 'mean', 'Location', 'northeast');
exportgraphics(gcf, sprintf("../figure/Macro_normalload.png"), "Resolution",300);


%% Plot preparation stage
fig = figure(4);
fig.Units = 'point';
fig.Position = [0 500 800 450];
clf(fig,'reset'); cla(fig,'reset'); hold on; box on; 

tshift = 80;
plot(M.FB03_086.tmat - tshift, M.FB03_086.SS, "k-", "DisplayName", "FB03-086");
plot(M.FB03_093.tmat, M.FB03_093.SS, "DisplayName", "FB03-093");

xlim([345, 420]);
% ylim([0.5, 0.7]);
xlabel("Time [s]");
ylabel("Shear stress [MPa]");

legend('Location', 'northwest');
exportgraphics(gcf, sprintf("../figure/Pre_Macro_shear_stress.png"), "Resolution",300);


%%
fig = figure(5);
fig.Units = 'point';
fig.Position = [0 500 800 450];
clf(fig,'reset'); cla(fig,'reset'); hold on; box on; 

tshift = 80;
plot(M.FB03_086.tmat - tshift, M.FB03_086.FCM, "k-", "DisplayName", "FB03-086");
plot(M.FB03_093.tmat, M.FB03_093.FCM, "DisplayName", "FB03-093");

xlim([345, 420]);
% ylim([0.5, 0.7]);
xlabel("Time [s]");
ylabel("Friction Coefficient");

legend('Location', 'northwest');
exportgraphics(gcf, sprintf("../figure/Pre_Macro_frictioncoef.png"), "Resolution",300);

%%
fig = figure(6);
fig.Units = 'point';
fig.Position = [0 500 800 450];
clf(fig,'reset'); cla(fig,'reset'); hold on; box on; 

tshift = 80;
h1 = plot(M.FB03_086.tmat- tshift, M.FB03_086.NP, "k-");
h2 = plot(M.FB03_093.tmat, M.FB03_093.NP, "r-", "DisplayName", "FB03-094");

% Plot mean value
h3 = plot(M.FB03_086.tmat- tshift, mean(M.FB03_086.NP, 2), "k--", "LineWidth", 3);
h4 = plot(M.FB03_093.tmat, mean(M.FB03_093.NP, 2), "r--", "LineWidth", 2);

xlim([0, 420]);
% ylim([0.5, 4]);
xlabel("Time [s]");
ylabel("Normal stress [MPa]");

legend([h1(1), h2(1), h3(1), h4(1)], 'FB03-086', 'FB03-093', 'mean', 'mean', 'Location', 'northwest');
exportgraphics(gcf, sprintf("../figure/Pre_Macro_normalload.png"), "Resolution",300);




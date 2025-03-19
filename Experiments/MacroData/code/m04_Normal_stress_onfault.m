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

%% Load decimated data

D = struct();

for expr_id = [86, 87, 93, 94]
    expr_id
    pname = sprintf('/Volumes/Okuboetal2025_masterHDD/FB03data/fb03-%03d/biax/', expr_id);
    runID = sprintf('fb03-%03d', expr_id);
    
    pname, runID
    [tmat,Snmat,Taumat3,Taumat2,Spmat]=fmbiaxlocalSTS3(pname,runID);
    
    D.(sprintf("FB03_%03d", expr_id)).tmat=tmat;
    D.(sprintf("FB03_%03d", expr_id)).Snmat=Snmat;
    D.(sprintf("FB03_%03d", expr_id)).Taumat3=Taumat3;
    D.(sprintf("FB03_%03d", expr_id)).Taumat2=Taumat2;
    D.(sprintf("FB03_%03d", expr_id)).Spmat=Spmat;

end

%%
save("../data/DecimatedData_raw_paper.mat", "D");
%%
load("../data/DecimatedData_raw_paper.mat", "D");

%% Plot normal stress
% Coordinate of strain gauges
% No correction with initial rock specimen location in this plot for simplicity.
SG3_x = [170:240:4000, 230:240:4000];
SG2_x = [290:240:4000, 110:240:3900];

%% Remove the offset
offsetrange = 10;

for expr_id = [86, 87, 93, 94]

    D.(sprintf("FB03_%03d", expr_id)).Snmat_offsetremoved = D.(sprintf("FB03_%03d", expr_id)).Snmat - mean(D.(sprintf("FB03_%03d", expr_id)).Snmat(1:offsetrange, :)); 
    D.(sprintf("FB03_%03d", expr_id)).Spmat_offsetremoved = D.(sprintf("FB03_%03d", expr_id)).Spmat - mean(D.(sprintf("FB03_%03d", expr_id)).Spmat(1:offsetrange, :)); 
    D.(sprintf("FB03_%03d", expr_id)).Taumat2_offsetremoved = D.(sprintf("FB03_%03d", expr_id)).Taumat2 - mean(D.(sprintf("FB03_%03d", expr_id)).Taumat2(1:offsetrange, :)); 
    D.(sprintf("FB03_%03d", expr_id)).Taumat3_offsetremoved = D.(sprintf("FB03_%03d", expr_id)).Taumat3 - mean(D.(sprintf("FB03_%03d", expr_id)).Taumat3(1:offsetrange, :)); 

end

%% Preparation
fig = figure(1);
fig.Units = 'point';
fig.Position = [0 500 1000 600];
clf(fig,'reset'); cla(fig,'reset'); hold on; box on; 

ylimit = [-1, 6];
sn_apply = 2; 
bar_width = 0.5;

hold on; box on; grid on;

% nt = 490000;

bar(SG3_x(1:16), D.FB03_086.Snmat_offsetremoved(end, 1:16), bar_width);
bar(SG3_x(17:end), D.FB03_086.Snmat_offsetremoved(end, 17:end), bar_width);

% plot([0, 4100], [sn_apply, sn_apply], "k--");
% plot flat jack
% plot(NP_x*1e3, D.NP(end, :), "ro-", "LineWidth",2);

title(sprintf("FB03-086 %4.1f[s]", D.FB03_086.tmat(end)));
xlabel("Location [mm]");
ylabel("Normal stress [MPa]");
xlim([0, 4100]);
ylim(ylimit);
legend({"South 1-16", "North 17-32"}, "Location", "northwest");

set(gcf, 'Color', 'w');
exportgraphics(gcf, sprintf("../figure/normalstress_FB03_086.png"), "Resolution", 150);

%%
fig = figure(2);
fig.Units = 'point';
fig.Position = [0 500 1000 600];
clf(fig,'reset'); cla(fig,'reset'); hold on; box on; 

ylimit = [-1, 6];
sn_apply = 2; 
bar_width = 0.5;

hold on; box on; grid on;

% nt = 490000;

bar(SG3_x(1:16), D.FB03_093.Snmat_offsetremoved(end, 1:16), bar_width);
bar(SG3_x(17:end), D.FB03_093.Snmat_offsetremoved(end, 17:end), bar_width);

% plot([0, 4100], [sn_apply, sn_apply], "k--");
% plot flat jack
% plot(NP_x*1e3, D.NP(end, :), "ro-", "LineWidth",2);

title(sprintf("FB03-093 %4.1f[s]", D.FB03_093.tmat(end)));
xlabel("Location [mm]");
ylabel("Normal stress [MPa]");
xlim([0, 4100]);
ylim(ylimit);
legend({"South 1-16", "North 17-32"}, "Location", "northwest");

set(gcf, 'Color', 'w');
exportgraphics(gcf, sprintf("../figure/normalstress_FB03_093.png"), "Resolution", 150);


%% Main stick-slip
fig = figure(3);
fig.Units = 'point';
fig.Position = [0 500 1000 600];
clf(fig,'reset'); cla(fig,'reset'); hold on; box on; 

ylimit = [-1, 6];
sn_apply = 2; 
bar_width = 0.5;

hold on; box on; grid on;

nt = 97658;

bar(SG3_x(1:16), D.FB03_087.Snmat(nt, 1:16), bar_width);
bar(SG3_x(17:end), D.FB03_087.Snmat(nt, 17:end), bar_width);

% plot([0, 4100], [sn_apply, sn_apply], "k--");
% plot flat jack
% plot(NP_x*1e3, D.NP(end, :), "ro-", "LineWidth",2);

title(sprintf("FB03-087 %4.1f[s]", D.FB03_087.tmat(nt)));
xlabel("Location [mm]");
ylabel("Normal stress [MPa]");
xlim([0, 4100]);
ylim(ylimit);
legend({"South 1-16", "North 17-32"}, "Location", "northwest");

set(gcf, 'Color', 'w');
exportgraphics(gcf, sprintf("../figure/normalstress_FB03_087.png"), "Resolution", 150);

%%
fig = figure(4);
fig.Units = 'point';
fig.Position = [0 500 1000 600];
clf(fig,'reset'); cla(fig,'reset'); hold on; box on; 

ylimit = [-1, 6];
sn_apply = 2; 
bar_width = 0.5;

hold on; box on; grid on;

nt = 97658;

bar(SG3_x(1:16), D.FB03_094.Snmat(nt, 1:16), bar_width);
bar(SG3_x(17:end), D.FB03_094.Snmat(nt, 17:end), bar_width);

% plot([0, 4100], [sn_apply, sn_apply], "k--");
% plot flat jack
% plot(NP_x*1e3, D.NP(end, :), "ro-", "LineWidth",2);

title(sprintf("FB03-094 %4.1f[s]", D.FB03_094.tmat(nt)));
xlabel("Location [mm]");
ylabel("Normal stress [MPa]");
xlim([0, 4100]);
ylim(ylimit);
legend({"South 1-16", "North 17-32"}, "Location", "northwest");

set(gcf, 'Color', 'w');
exportgraphics(gcf, sprintf("../figure/normalstress_FB03_094.png"), "Resolution", 150);



% 
% 
% 
% %% Plot Main stick-slip experiments
% fig = figure(1);
% fig.Units = 'point';
% fig.Position = [0 500 800 450];
% clf(fig,'reset'); cla(fig,'reset'); hold on; box on; 
% 
% tshift = 20.2;
% plot(M.FB03_087.tmat- tshift, M.FB03_087.SS, "k-", "DisplayName", "FB03-087");
% plot(M.FB03_094.tmat, M.FB03_094.SS, "DisplayName", "FB03-094");
% 
% xlim([0, 250]);
% ylim([0.5, 0.7]);
% xlabel("Time [s]");
% ylabel("Shear stress [MPa]");
% 
% legend('Location', 'northeast');
% exportgraphics(gcf, sprintf("../figure/Macro_shear_stress.png"), "Resolution",300);
% 
% 
% %%
% fig = figure(2);
% fig.Units = 'point';
% fig.Position = [0 500 800 450];
% clf(fig,'reset'); cla(fig,'reset'); hold on; box on; 
% 
% tshift = 20.2;
% plot(M.FB03_087.tmat- tshift, M.FB03_087.FCM, "k-", "DisplayName", "FB03-087");
% plot(M.FB03_094.tmat, M.FB03_094.FCM, "DisplayName", "FB03-094");
% 
% xlim([0, 250]);
% ylim([0.3, 0.45]);
% xlabel("Time [s]");
% ylabel("Friction Coefficient");
% 
% legend('Location', 'northeast');
% exportgraphics(gcf, sprintf("../figure/Macro_frictioncoef.png"), "Resolution",300);
% 
% 
% %%
% fig = figure(3);
% fig.Units = 'point';
% fig.Position = [0 500 800 450];
% clf(fig,'reset'); cla(fig,'reset'); hold on; box on; 
% 
% tshift = 20.2;
% h1 = plot(M.FB03_087.tmat- tshift, M.FB03_087.NP, "k-");
% h2 = plot(M.FB03_094.tmat, M.FB03_094.NP, "r-", "DisplayName", "FB03-094");
% 
% % Plot mean value
% h3 = plot(M.FB03_087.tmat- tshift, mean(M.FB03_087.NP, 2), "k--", "LineWidth", 3);
% h4 = plot(M.FB03_094.tmat, mean(M.FB03_094.NP, 2), "r--", "LineWidth", 2);
% 
% xlim([0, 250]);
% ylim([0.5, 4]);
% xlabel("Time [s]");
% ylabel("Normal stress [MPa]");
% 
% legend([h1(1), h2(1), h3(1), h4(1)], 'FB03-087', 'FB03-094', 'mean', 'mean', 'Location', 'northeast');
% exportgraphics(gcf, sprintf("../figure/Macro_normalload.png"), "Resolution",300);
% 
% 
% %% Plot preparation stage
% fig = figure(4);
% fig.Units = 'point';
% fig.Position = [0 500 800 450];
% clf(fig,'reset'); cla(fig,'reset'); hold on; box on; 
% 
% tshift = 80;
% plot(M.FB03_086.tmat - tshift, M.FB03_086.SS, "k-", "DisplayName", "FB03-086");
% plot(M.FB03_093.tmat, M.FB03_093.SS, "DisplayName", "FB03-093");
% 
% xlim([345, 420]);
% % ylim([0.5, 0.7]);
% xlabel("Time [s]");
% ylabel("Shear stress [MPa]");
% 
% legend('Location', 'northwest');
% exportgraphics(gcf, sprintf("../figure/Pre_Macro_shear_stress.png"), "Resolution",300);
% 
% 
% %%
% fig = figure(5);
% fig.Units = 'point';
% fig.Position = [0 500 800 450];
% clf(fig,'reset'); cla(fig,'reset'); hold on; box on; 
% 
% tshift = 80;
% plot(M.FB03_086.tmat - tshift, M.FB03_086.FCM, "k-", "DisplayName", "FB03-086");
% plot(M.FB03_093.tmat, M.FB03_093.FCM, "DisplayName", "FB03-093");
% 
% xlim([345, 420]);
% % ylim([0.5, 0.7]);
% xlabel("Time [s]");
% ylabel("Friction Coefficient");
% 
% legend('Location', 'northwest');
% exportgraphics(gcf, sprintf("../figure/Pre_Macro_frictioncoef.png"), "Resolution",300);
% 
% %%
% fig = figure(6);
% fig.Units = 'point';
% fig.Position = [0 500 800 450];
% clf(fig,'reset'); cla(fig,'reset'); hold on; box on; 
% 
% tshift = 80;
% h1 = plot(M.FB03_086.tmat- tshift, M.FB03_086.NP, "k-");
% h2 = plot(M.FB03_093.tmat, M.FB03_093.NP, "r-", "DisplayName", "FB03-094");
% 
% % Plot mean value
% h3 = plot(M.FB03_086.tmat- tshift, mean(M.FB03_086.NP, 2), "k--", "LineWidth", 3);
% h4 = plot(M.FB03_093.tmat, mean(M.FB03_093.NP, 2), "r--", "LineWidth", 2);
% 
% xlim([0, 420]);
% % ylim([0.5, 4]);
% xlabel("Time [s]");
% ylabel("Normal stress [MPa]");
% 
% legend([h1(1), h2(1), h3(1), h4(1)], 'FB03-086', 'FB03-093', 'mean', 'mean', 'Location', 'northwest');
% exportgraphics(gcf, sprintf("../figure/Pre_Macro_normalload.png"), "Resolution",300);
% 
% 
% 

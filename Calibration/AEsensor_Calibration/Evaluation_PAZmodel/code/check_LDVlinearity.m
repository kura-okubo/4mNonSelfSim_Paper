% Script to check the linearity of LDV test for sensor calibration
%2022/09/17 Kurama Okubo

clear all;
set(0,'DefaultTextFontsize',14, ...
    'DefaultTextFontname','Arial', ...
    'DefaultTextFontWeight','normal', ...
    'DefaultTextFontname','Arial', ...
    'DefaultAxesFontsize',14, ... 
    'DefaultAxesFontname','Arial', ...
    'DefaultLineLineWidth', 1.5)
set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'})

addpath("../../../../utils/matlabcode_biax_v03/");
%% read linear data
% input_ohm = 50; % input inpedence

V = [10, 25, 50, 75, 100, 125, 150, 175, 200]; %[V]
fs = 1e7;

D = importdata('../../AEsensor_Calibration_DataConversion/data/LDV_linearprog.csv');

%%
tmat = D.data(:, 1);
datmat = D.data(:, 2:end);

% trim data
xlimit = [0, 0.001];
trimind = find(tmat >= xlimit(2), 1);
tmat_trim = tmat(1:trimind);
dat_trim = datmat(1:trimind, :);

%% Plot scatter
fig = figure(1); clf;
fig.Units = 'point';
fig.Position = [0 800 700 700];
box on; axis equal; hold on;

sc = viridis(9);

figdir_linear = "../figure/linear_all/";
if ~isfolder(figdir_linear); mkdir(figdir_linear); end

dat_baseind = find(V==100); % base is V = 100V;
dat_base = datmat(1:trimind, dat_baseind);

for i = 1:9
    clf; 
    hold on; box on; grid on;
    S = scatter(dat_base*1e3, dat_trim(:, i)*1e3, 10, sc(i, :), ...
        "MarkerEdgeColor", "None", "MarkerFaceColor", sc(i, :),...
        "MarkerFaceAlpha", 0.8, "Marker", "o");
    xlim([-0.5, 0.5]);
    ylim([-0.5, 0.5]);
    xlabel("Amplitude with source 100V");
    ylabel(sprintf("Amplitude with source %dV", V(i)));
    exportgraphics(gcf, figdir_linear+sprintf("linearity_%dV.png",V(i)) , "Resolution", 150);

end




%% Plot associated with 100V vs 200V
fig = figure(2); clf;
fig.Units = 'point';
fig.Position = [0 800 700 700];
box on; axis equal; hold on;

subplot(3, 1, 3); hold on;

tr_100 = datmat(:, dat_baseind);
tr_200 = datmat(:, 9);
tr_100 = tr_100 - mean(tr_100(1:100));
tr_200 = tr_200 - mean(tr_200(1:100));


plot(tmat*1e3, tr_100*1e3, "k-", "LineWidth",1.5, "DisplayName","V100");
plot(tmat*1e3, tr_200*1e3/1.7911, "r-", "LineWidth",1.5, "DisplayName","V200 / 1.79");

xlim([0.01, 0.04]);
ylim([-1.2, 1.2]);
xlabel("Time [ms]");
ylabel("Amplitude [V]");

box on;
legend("Location","northeast");

subplot(3, 1, [1, 2]); hold on;
axis equal; box on; grid on;

% compute linear regression
X = [ones(length(tr_100),1) tr_100];

b = X \ tr_200;
yCalc2 = X*b;
Rsq2 =  1 - sum((tr_200 - yCalc2).^2)/sum((tr_200 - mean(tr_200)).^2);

S = scatter(tr_100*1e3, tr_200*1e3, 10, sc(2, :), ...
    "MarkerEdgeColor", "None", "MarkerFaceColor", sc(2, :),...
    "MarkerFaceAlpha", 0.8, "Marker", "o");

plot(tr_100*1e3, yCalc2*1e3, "r-");
text(0.1, -0.8, sprintf("k=%4.4f, R^2=%4.2f", b(2), Rsq2));

xlim([-1, 1]);
ylim([-1, 1]);
xticks(-1:0.5:1);
yticks(-1:0.5:1);
xlabel("Amplitude with V100 [V]");
ylabel("Amplitude with V200 [V]");

set(gcf, 'Color', 'w');
exportgraphics(gcf, "../figure/linear_all/linearity_V100vsV200.png", "Resolution", 150);

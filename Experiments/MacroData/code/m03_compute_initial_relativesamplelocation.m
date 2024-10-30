% Script to plot macroscopic data
clear all;
set(0,'DefaultTextFontsize',16, ...
    'DefaultTextFontname','Arial', ...
    'DefaultTextFontWeight','normal', ...
    'DefaultTextFontname','Arial', ...
    'DefaultAxesFontsize',16, ... 
    'DefaultAxesFontname','Arial', ...
    'DefaultLineLineWidth', 1.0)

set(0,'defaulttextinterpreter','none')
figdir_root="../figure";
if ~exist(figdir_root) mkdir(figdir_root); end

addpath("../../../utils/matlabcode_biax_v03");

%% read data

expr_id = 87;

pname = sprintf('/Volumes/4mGouge_WorkHDD/FB03data/fb03-%03d/biax/', expr_id);
runID = sprintf('fb03-%03d',expr_id);

[tmat,DX,NP,SS,FCM]=fmbiaxmacroS3_initrelativeloc(pname,runID);

%%
fig = figure(1); clf; box on; hold on;
fig.Units = 'point';
fig.Position = [0 800 800 450];


plot(tmat, DX(:, 1), "-", "DisplayName","DX west");
plot(tmat, DX(:, 2), "-", "DisplayName","DX east");

xlabel("Time [s]");
ylabel("Absolute value of laser disp sensor [mm]");

u0 = -3.97; %[mm] absolute laser disp when the east side of bottom rock is aligned to the top
ue = mean(DX(1:10, 2));
xtop_shift = 100-(ue-u0);

title(sprintf("%s U0 east: %4.2fmm u0east:%4.2fmm \nxtop_shift:%4.2fmm Init cylinder length: %4.2fmm", runID, ue, u0, xtop_shift, 100-xtop_shift));
legend("Location","best");

figname = sprintf("../figure/%s_initial_relativesampleloc.png", runID);

exportgraphics(fig, figname);

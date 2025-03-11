% Script to plot AIC
% 2022/09/09 Kurama Okubo
% update: for master script
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

%%
figdir_AIC = "../figure/AIC/";
if ~isfolder(figdir_AIC); mkdir(figdir_AIC); end

%% Read input (LDV) and output (AEsensor) data

casename_list = ["frontcenter", "fronttop", "sidecenter"];

ic = 2;

casename = casename_list(ic);

A = load(sprintf("../data/AE_resp_dataandcoef_%s.mat", casename));

%% Plot AIC
fig = figure(1); clf;
fig.Units = 'point';
fig.Position = [0 800 700 600];
% clf(fig,'reset'); cla(fig,'reset'); hold on;
box on; grid on;

x = A.all_AIC(:, 1); %na
y = A.all_AIC(:, 2); %nb
% z = A.all_AIC(:, 3); %AIC
z = normalize(A.all_AIC(:, 3), 'range'); % normalized AIC

N = size(A.AIC_mask, 1); % Number Of Points Desired
xv = linspace(min(x), max(x), N);
yv = linspace(min(y), max(y), N);
[X,Y] = ndgrid(xv, yv);
Z = griddata(x, y, z, X, Y);

%%
% create mask patch
mask = ones(N+2, N+2); % extend mask to make patch
mask(2:N+1, 2:N+1) = A.AIC_mask;
[mX,mY] = ndgrid(0:N+1, 0:N+1);

%%


hold on;

% contourf(X, Y, Z, 100,'edgecolor','none');
% plot AIC
h=imagesc(xv, yv, Z');
colormap("viridis");
caxis([0.0, 0.1]);
% caxis([0.031, 0.032]);

hc = colorbar;
hc.Label.String = 'normalized AIC';


% % plot mask
for i = 1:N
    for j = 1:N
        if A.AIC_mask(i, j) == 1           
            plot(X(i, j), Y(i, j), "x", "MarkerSize", 12, "Color", "w");
        end
    end
end


% set(h, 'EdgeColor', 'none');
set(gca, 'YDir', 'normal');
axis('equal')

% plot edgelines
hold on;
edge_x = repmat((0:N)+0.5, N+1,1);
edge_y = repmat((0:N)+0.5, N+1,1).';
plot(edge_x ,edge_y, 'k', "LineWidth", 0.1) % vertical lines
plot(edge_x.', edge_y.', 'k', "LineWidth", 0.1) % horizontal lines

% plot best AIC
hold on;
% plot(A.Na_best, A.Nb_best, "o", "MarkerSize",12,"MarkerFaceColor","w","MarkerEdgeColor","k");
plot(A.Na_best, A.Nb_best, "p", "MarkerSize",20,"MarkerFaceColor","y","MarkerEdgeColor","k");

% % plot the combination used for the analysis: Na=20, Nb=6;
% plot(20, 6, "p", "MarkerSize",20,"MarkerFaceColor","y","MarkerEdgeColor","k");


xlabel("m");
ylabel("n");
zlabel("AIC");

xlim([0.5, 30.5]);
ylim([0.5, 30.5]);


set(gcf, 'Color', 'w');
figname = sprintf("../figure/AIC/FigS_AIC_all_%s.eps", casename);
exportgraphics(fig, figname);

figname = sprintf("../figure/AIC/FigS_AIC_all_%s.png", casename);
exportgraphics(fig, figname, "Resolution", 80);

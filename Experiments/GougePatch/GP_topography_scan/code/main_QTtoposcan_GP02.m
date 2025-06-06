% Script to compute topology
% 2023/02/23 Kurama Okubo
% 2023/12/21 update for rectangular host rock
clear all;
set(0,'DefaultTextFontsize',10, ...
    'DefaultTextFontname','Arial', ...
    'DefaultTextFontWeight','normal', ...
    'DefaultTextFontname','Arial', ...
    'DefaultAxesFontsize',10, ... 
    'DefaultAxesFontname','Arial', ...
    'DefaultLineLineWidth', 1.0)
set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'})

addpath("../../../../utils/matlabcode_biax_v03");

%%

param.l_space = 27.0;
param.ds = 0.05;
param.fs = 20;
param.fs_TEAC = 100; % frequency of TEAC recording
param.vpps = 500;

param.R_sample = 40.0;
param.center_x1 = 13.35;
param.center_x2 = 13.65;

param.Vout_thresh = 0.1;

param.low_clipval = -2.0; %[um] lower valuer of clip and interpolation for missing data point
param.gsize = 20;
param.sig = 10; 

%% Read data

ifLoadfromOriginal = false;

if ifLoadfromOriginal
    
    pname = '../data/';
    fname = 'GPTOPO-02_001';
    
    [tmat, Zmat, Trigmat]=TEACtopo(pname,fname);
    
    figure(1); clf;
    subplot(2, 1, 1);
    plot(tmat, Zmat, "+-");
    
    subplot(2, 1, 2);
    plot(tmat, Trigmat);
    
    %% Trim the data
    % We have an artifact at the end of the data due to the stop signal from TEAC.
    % Thus, trim the signal
    et_ind = find(tmat>1.8e4, 1);
    tmat = tmat(1:et_ind);
    Zmat = Zmat(1:et_ind);
    Trigmat = Trigmat(1:et_ind);
    
    %% Downsampling the data
    % TEAC data is sampled high-frequency enough to avoid the aliasing along
    % scan trace.
    param.rdecim = ceil(param.fs_TEAC/param.fs);
    Zmat_downsampled = decimate(Zmat,param.rdecim, 8);
    tmat_downsampled = downsample(tmat,param.rdecim);
    Trigmat_downsampled = decimate(Trigmat,param.rdecim);
    
    figure(9); clf; hold on;
    plot(tmat, Zmat, "+-");
    plot(tmat_downsampled, Zmat_downsampled, "o-");
    xlim([3720, 3750]);
    
    %%
    figure(10); clf; hold on;
    plot(tmat_downsampled(1:end-1), diff(Trigmat_downsampled), "o-");
    % xlim([3720, 3750]);
    
    %% Plot raw data
    fig=figure(11); clf; hold on;
    fig.Units = 'point';
    fig.Position = [0 500 800 600];
    clf(fig,'reset'); cla(fig,'reset'); 
    hold on; box on; axis equal;
    
    xlimit=[0, 1.7786e4];
    
    subplot(211); hold on; box on;
    plot(tmat_downsampled, Zmat_downsampled, "ko-");
    title("Laser Z raw data");
    xlim(xlimit);
    xlabel("Time [s]");
    ylabel("u [μm]");
    
    subplot(212); hold on; box on;
    plot(tmat_downsampled(1:end-1), diff(Trigmat_downsampled), "bo-");
    title("Trigger signal raw data");
    xlim(xlimit);
    xlabel("Time [s]");
    ylabel("Volt [V]");
    
    exportgraphics(gcf, sprintf('../figure/Ztopo_rawdata_%s.png', fname), "Resolution", 150);
    
    
    %% Compute topography from the Zmat time series
    [Ztopo, x1vec, x2vec, I] = compute_QTtopo_rect(Zmat_downsampled, tmat_downsampled, Trigmat_downsampled, param);

    %% Save data
    save(sprintf("../data/Ztopo_data_%s.mat", fname), "x1vec", "x2vec", "Ztopo", "I");
else
    % loading the preprocessed 2D topography data
    fname = 'GPTOPO-02_001';
    load(sprintf("../data/Ztopo_data_%s.mat", fname));
end

%% Plot result
fig = figure(2);
fig.Units = 'point';
fig.Position = [0 500 300 300];
clf(fig,'reset'); cla(fig,'reset'); 
hold on; box on; axis equal;
colormap("viridis");

xshift = 0.21; % adjust center x

[C, h] = contourf(x1vec-xshift, x2vec, Ztopo, 151, 'Tag', "topo");
% pcolor(x1vec, x2vec, Ztopo, 'Tag', "topo");
set(h,'LineColor','none');
% set(h, 'EdgeColor', 'none');

cb = colorbar;
cb.Label.String = "Height [μm]";

%https://jp.mathworks.com/matlabcentral/answers/1449709-adjusting-width-of-horizontal-colorbar#answer_783834
ax1 = gca();
ax1Pos = ax1.Position;
pos = cb.Position;
pos(3) = 0.8*pos(3);
cb.Position = pos;
% The above automatically change the ax1.Position
% We restore the origninal ax1.Position
ax1.Position = ax1Pos;


% shading interp
shading flat
clim([0, 60]); %[0, 80]);

xlim([-8, 8]);
ylim([-8, 8]);

xticks([-8, -4, 0, 4, 8 ]);
yticks([-8, -4, 0, 4, 8 ]);

xlabel("x [mm]");
ylabel("z [mm]");

% plot scale
plot([-4, 4], [-5, -5], "w-")

% convert contourf as raster

exportgraphics(gca, '../figure/Ztopo_processed.png', "Resolution", 70);

% Save the contourf with raster
% We modified vecrast.m for the GP topography plot.
% Theo Michelis (2025). vecrast (https://www.mathworks.com/matlabcentral/fileexchange/154667-vecrast), MATLAB Central File Exchange. Retrieved March 18, 2025.
% Copyright (c) 2023, Theodoros Michelis
% All rights reserved.


% To save the mixed raster and vector file, change the following section in
% vecrast:

% 1. 
% for i = 1:length(axesHandleVector)
%     axesPosition{i} = get(axesHandleVector(i), 'position'); % Fix: get axes size
%     toRemove = [...
%         findobj(axesHandleVector(i).Children, 'Tag','topo');... % update
%         ];
%     set(toRemove, 'visible', 'off');
% end

% 2. 
% for i = 1:length(axesHandleRaster)
%     axesPosition{i} = get(axesHandleRaster(i), 'position'); % Fix: get axes size
%     contents = findall(axesHandleRaster(i));
%     toKeep = [...
%         findobj(axesHandleRaster(i).Children, 'Tag','topo');... % update
%         ];
%     toRemove = setxor(contents, toKeep);
%     set(toRemove, 'visible', 'off');
% end

%3.
% if isempty(findall(rasterFigure, 'type', 'surface'))
% if isempty(findall(rasterFigure, 'type', 'contour'))

% Remove comment out after downloading and modifying vecrast
vecrast_topo(fig, '../figure/Ztopo_processed', 300, 'bottom', 'eps');

% We swhitched to the one-way scanning with trigger signals.
%%% Note: There is a limitation in the horizontal stripes in the figure
%%% which could be the unbalance between the RL and LR scanning.
%%% We need to either one-way scan or use 2D laser displacement sensor.

%% Plot surface
fig = figure(3);
fig.Units = 'point';
fig.Position = [0 500 1200 800];
clf(fig,'reset'); cla(fig,'reset'); 
hold on; box on; axis equal;
colormap("viridis"); grid on;

h = surf(x1vec, x2vec, Ztopo);
set(h, 'EdgeColor', 'none');

cb = colorbar;
cb.Label.String = "Height [um] ";
shading interp
% shading flat
% shading faceted 
clim([0, 80]);

xlim([-15, 15]);
ylim([-15, 15]);
zlim([-5, 100]);

% view(30,30);
view(-30,40);

xlabel("x1 [mm]");
ylabel("x2 [mm]");
zlabel("Height [μm]");

daspect([1 1 5e1]);

exportgraphics(gca, '../figure/Ztopo_surface.png', "Resolution", 70);

% % % Note: There is a limitation in the horizontal stripes in the figure
% % % which could be the unbalance between the RL and LR scanning.
% % % We need to either one-way scan or use 2D laser displacement sensor.

%% Plot surface for schematic
fig = figure(12);
fig.Units = 'point';
fig.Position = [0 500 1200 800];
clf(fig,'reset'); cla(fig,'reset'); 
hold on; box on; axis equal;

colormap("gray"); grid on;

h = surf(x1vec, x2vec, Ztopo);
set(h, 'EdgeColor', 'none');

% cb = colorbar;
% cb.Label.String = "Height [um] ";
shading interp
shading flat
% shading faceted 
clim([-20, 110]);

xlim([-15, 15]);
ylim([-15, 15]);
zlim([-5, 100]);

% view(30,30);
view(-40,25);

xlabel("x1 [mm]");
ylabel("x2 [mm]");
zlabel("Height [μm]");

daspect([1 1 15e1]);
axis off;

exportgraphics(gca, '../figure/Ztopo_surface_schematic.eps');

%% Debug 
fig = figure(9);
fig.Units = 'point';
fig.Position = [0 500 1200 800];
clf(fig,'reset'); cla(fig,'reset'); 
hold on; box on; axis equal;
colormap("viridis");
h = surf(x1vec, x2vec, I.Zmat_topo);
set(h, 'EdgeColor', 'none');

% caxis([0, 10]);
title("Zmat topo before plane removal")
xlim([-15, 15]);
ylim([-15, 15]);
% zlim([-160e-3, -100e-3]);
cb = colorbar;
cb.Label.String = "Height [um] ";

%% Plot crosssection

xbounds = [-5, 5];
plotinds = find((x2vec > xbounds(1)) & (x2vec < xbounds(2)));

plotind_centerind = find(x2vec >= 0, 1);
%%
fig = figure(4);
fig.Units = 'point';
fig.Position = [0 500 1200 400];
clf(fig,'reset'); cla(fig,'reset'); 
hold on; box on; grid on;
colormap("viridis");

% for ii = plotinds(1:plotstep:end)
%     lh = plot(x1vec, Ztopo(ii, :), "-");
%     lh.Color = [lh.Color 0.5];
% end

lh1 = plot(x1vec, Ztopo(plotind_centerind, :), "k-o");
% lh2 = plot(x1vec, Ztopo(:, plotind_centerind), "b-");

xlim([-12.5, 12.5]);
% ylim([-3, 15]);
ylim([-3, 200]);

xlabel("x1 [mm]");
ylabel("Topography [μm]");

exportgraphics(gca, '../figure/Ztopo_crosssection.png', "Resolution", 70);

%% Compute mean in 1D Line

Ztopo_line = Ztopo(plotind_centerind, :);
Ztopo_line(Ztopo_line<10) = NaN;
mean_topo = mean(Ztopo_line, "omitnan");
max_topo = max(Ztopo_line);
fprintf("Line mean topo: %4.2fµm, max topo: %4.2fµm\n", mean_topo, max_topo);

%% Compute mean and max height in 2D
Ztopo_meanmax = Ztopo;

Nx1 = size(Ztopo_meanmax, 1);
Nx2 = size(Ztopo_meanmax, 2);

x1_center = 0; %[mm] shift the center
x2_center = 0; %[mm] shift the center

patchradius = 4.2;

for i = 1:Nx1
    for j = 1:Nx2
        x1 = x1vec(i);
        x2 = x2vec(j);
        R = norm([(x1-x1_center), (x2-x2_center)]);
        if R > patchradius
            Ztopo_meanmax(i,j) = NaN;
        end
    end
end

fig = figure(15); clf;
surf(x1vec, x2vec, Ztopo_meanmax);
shading interp
colorbar;

mean_topo_2D = mean(Ztopo_meanmax, "all", "omitnan");
medoan_topo_2D = median(Ztopo_meanmax, "all", "omitnan");
max_topo_2D = max(max(Ztopo_meanmax));
fprintf("2D mean topo: %4.2fµm, median: %4.2fµm, max topo: %4.2fµm\n", mean_topo_2D,medoan_topo_2D,max_topo_2D);
fo = fopen('../data/mean_topo.txt','w');
fprintf(fo, "2D mean topo: %4.2fµm, median: %4.2fµm, max topo: %4.2fµm\n", mean_topo_2D,medoan_topo_2D,max_topo_2D);
fclose(fo);

%% PSD of topography

window = 2^4;
noverlap = round(window/2);
nfft = 2^7;

plotstep = 1;
pinds_stepped = plotinds(1:plotstep:end);
pxx_all = zeros(length(pinds_stepped), nfft/2+1);
itrace = 0;

for ii = pinds_stepped
    tr = Ztopo(ii, :);
    tr = tr(~isnan(tr));
    [pxx,f] = pwelch(tr,window,noverlap,nfft, 1/param.ds);
    itrace = itrace+1;
    pxx_all(itrace, :) = pxx;
end

%%

fig = figure(5);
fig.Units = 'point';
fig.Position = [0 500 600 500];
clf(fig,'reset'); cla(fig,'reset'); 
hold on; box on; grid off;

lh = plot(f, pxx_all, "k-");

for l = 1:length(lh)
    lh(l).Color = [lh(l).Color 0.2];
end

% stacked trace
plot(f, mean(pxx_all, 1), "r-", "LineWidth", 3);

%     tr = Ztopo(plotind_centerind, :);

% plot self-similar psd
% Dunham et al. 2011b
Pm = @(k, alpha) (2*pi)^3 * alpha^2 ./ abs(k);

% plot grid
alpha1 = [10^(-6), 10^(-2), 10^(-1)];
for i = 1:3
    pxxm = Pm(f(2:end)*1e3, alpha1(i));
    plot(f(2:end), 5e14* pxxm, "--");
end

xlim([0.1, 10]);
ylim([1e-4, 1e4]);

set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');

xlabel("k_{x1} [1/mm]");
ylabel("PSD of topography [mm^3]");

exportgraphics(gca, '../figure/Ztopo_PSD.png', "Resolution", 70);

%%


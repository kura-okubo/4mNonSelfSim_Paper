% Script to read the output of prescale mobile and plot the data and
% convert to mat file.
% Updated for GPPS-02
% 2023/08/23 Kurama Okubo

% 2025.03.18 update for master plot
% Added gray background with 2 MPa

clear all; close all;
set(0,'DefaultTextFontsize',10, ...
    'DefaultTextFontname','Arial', ...
    'DefaultTextFontWeight','normal', ...
    'DefaultTextFontname','Arial', ...
    'DefaultAxesFontsize',10, ... 
    'DefaultAxesFontname','Arial', ...
    'DefaultLineLineWidth', 1.0)

set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'})
addpath("../../../../utils/matlabcode_biax_v03");

%% Parameters
% if ~isfolder("../figure/GPPS-02-pressure");  mkdir("../figure/GPPS-02-pressure"); end

finame = "20230822202542/GPPS-02-LW-20230822201834_pressData.tsv"; % transfer the data from iPad.

Dout = 25.06;  % diameter of rock to rock
patchradius = 4; % [mm] size of gauge patch.

A = importdata("../data/"+finame);
P = A.data;
metadata = split(A.textdata(2));

%%
% A3 size tentative

Npx_x = str2double(cell2mat(metadata(2)));
Npx_y = str2double(cell2mat(metadata(3)));
ds = str2double(cell2mat(metadata(4)));
ds = 0.355*0.66 %66; % NOTE: Need to check the resolution; this is emperically adjusted such that the diameter is 8mm.
xvec = (0:Npx_x-1) * ds;
yvec = (0:Npx_y-1) * ds;

%% Trim gauge patch

figure(9);
[~, h] = contourf(xvec, yvec, P, 101);
set(h,'LineColor','none')
set(gca, 'YDir','reverse')

%% Read the center coordinate from file
C = readmatrix("../data/GPPS-02_LW_centercoords.csv");

% select from 5, 6, 8, 9, 11
gp_id = 6; % We select case 3 of 2.0 MPa for master case

coords = C(C(:, 1)==gp_id, 2:3);
xc = coords(1);
yc = coords(2);

%%

windowsize=30; %mm

xvec_trim_ind = find(xvec>xc-windowsize/2 & xvec<xc+windowsize/2);
yvec_trim_ind = find(yvec>yc-windowsize/2 & yvec<yc+windowsize/2);

P_trim = P(yvec_trim_ind, xvec_trim_ind);

%%

xvec_trim = (0:length(xvec_trim_ind)-1) * ds;
yvec_trim = (0:length(yvec_trim_ind)-1) * ds;

% shift half so that the middle is 0
xvec_trim = xvec_trim - xvec_trim(end)/2;
yvec_trim = yvec_trim - yvec_trim(end)/2;


%% Update: set NaN for outside of base rock cylinder
for i = 1:length(xvec_trim)
    for j = 1: length(yvec_trim)
        if norm([xvec_trim(i), yvec_trim(j)], 2) > Dout/2
            P_trim(j, i) = NaN;
        end
    end
end


%%
fig = figure(1);
fig.Units = 'point';
fig.Position = [0 500 300 300];
clf(fig,'reset'); cla(fig,'reset'); hold on; box on;

c = colormap("magma");
c(1:1, :) = 0.8*ones(1, 3); % color on the base rock sample

set(gca,'Color','white')

R = Dout/2;
theta = linspace(0, 2*pi, 360);

xouter = R * cos(theta);
youter = R * sin(theta);

xgc = patchradius * cos(theta);
ygc = patchradius * sin(theta);

colormap(c);
[~, h] = contourf(xvec_trim, yvec_trim, P_trim, 101, 'Tag', 'pressure');
% h = pcolor(xvec_trim, yvec_trim, P_trim, 'FaceColor', 'interp', 'Tag', "pressure", EdgeColor="none");
set(h,'LineColor','none')

plot(xouter, youter, "k--", "LineWidth", 1);
plot(xgc, ygc, "k:");

axis equal;
% caxis([0, 10]);
clim([2.0, 6.5]);
cb = colorbar;
cb.Label.String = "Normal pressure [MPa]";

daspect([1 1 1]);
ax = gca();
ax.Layer = 'top';

xlabel("x [mm]");
ylabel("z [mm]");

% xlim([-windowsize/2, windowsize/2]);
% ylim([-windowsize/2, windowsize/2]);
xlim([-13, 13]);
ylim([-13, 13]);


set(gca, 'YDir','reverse')

% Compute mean and max
max_P = max(max((P_trim)));
P_trim_pad = P_trim;
P_trim_pad(P_trim_pad<2.0) = NaN; % ignore the small values
mean_P = mean(P_trim_pad, "all","omitnan");


% Save the contourf with raster
% We modified vecrast.m for the GP topography plot.
% Theo Michelis (2025). vecrast (https://www.mathworks.com/matlabcentral/fileexchange/154667-vecrast), MATLAB Central File Exchange. Retrieved March 18, 2025.
% Copyright (c) 2023, Theodoros Michelis
% All rights reserved.

% title(sprintf("GPPS-02-%03d mean:%2.1fMPa max:%2.1fMPa", gp_id, mean_P, max_P));
title(sprintf("\\sigma_n^{macro} = 2.0MPa, mean:%2.1fMPa, max:%2.1fMPa", mean_P, max_P));
% exportgraphics(gcf, sprintf("../figure/gaugepatch_fromPrescalemobile_GPPS-02-%03d.png", gp_id), "Resolution",70);
% exportgraphics(gcf, sprintf("../figure/gaugepatch_fromPrescalemobile_GPPS-02-%03d.eps", gp_id), 'ContentType', 'vector');

% remove comment out after modifying the vecrast; see GP_topography_scan/code/main_QTtoposcan_GP02.m
vecrast_pressure(fig, sprintf('../figure/gaugepatch_fromPrescalemobile_GPPS-02-%03d', gp_id), 300, 'bottom', 'eps');



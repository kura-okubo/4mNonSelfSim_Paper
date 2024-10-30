% Script to compute topology
% 2023/02/23 Kurama Okubo
% 2023/12/21 update for rectangular host rock & add noise level of Laser
% displacement sensor

clear all;
set(0,'DefaultTextFontsize',16, ...
    'DefaultTextFontname','Arial', ...
    'DefaultTextFontWeight','normal', ...
    'DefaultTextFontname','Arial', ...
    'DefaultAxesFontsize',16, ... 
    'DefaultAxesFontname','Arial', ...
    'DefaultLineLineWidth', 1.0)

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
pname = '../data/';
fname = 'GPTOPO-02_001';

%%
load(sprintf("../data/Ztopo_data_%s.mat", fname));

%% Select the traces
% we select the traces on and off the patch

tr_len = 8.0; %[mm]
tr_npts = round(tr_len/param.ds);

tr_x_init = find(-4<x1vec, 1);
tr_on_y_ind = find((-3.0<x2vec)&(x2vec<3.0));
tr_off_y_ind = find(((-7.5<x2vec)&(x2vec<-5.5) | (5.5<x2vec)&(x2vec<7.5)));

%%
fig = figure(1);
fig.Units = 'point';
fig.Position = [0 500 600 600];
clf(fig,'reset'); cla(fig,'reset'); 
hold on; box on; axis equal;
colormap("viridis");

[C, h] = contourf(x1vec, x2vec, Ztopo, 151);
set(h,'LineColor','none');
% set(h, 'EdgeColor', 'none');

cb = colorbar;
cb.Label.String = "Height from host rock [μm] ";
cb.Label.FontSize = 16;
% shading interp
shading flat
caxis([0, 60]);

xlim([-10, 10]);
ylim([-10, 10]);

xlabel("x [mm]");
ylabel("y [mm]");

% % trace off the patch
% for i = tr_off_y_ind(1:7:end)
%     plot(x1vec([tr_x_init, tr_x_init+tr_npts]), x2vec(i) * ones(2), "y-");
% end
% 
% % trace on the patch
% for i = tr_on_y_ind(1:7:end)
%     plot(x1vec([tr_x_init, tr_x_init+tr_npts]), x2vec(i) * ones(2), "r-");
% end

% plot range of PSD measurement
tr_on_init = x2vec(tr_on_y_ind(1));
tr_on_end = x2vec(tr_on_y_ind(end));
x1_init = x1vec(tr_x_init);
x1_end = x1vec(tr_x_init+tr_npts);
plot([x1_init, x1_init+0.2], [tr_on_init, tr_on_init], "r-");
plot([x1_init, x1_init+0.2], [tr_on_end, tr_on_end], "r-");
plot([x1_end, x1_end+0.2], [tr_on_init, tr_on_init], "r-");
plot([x1_end, x1_end+0.2], [tr_on_end, tr_on_end], "r-");


tr_off_init1 = x2vec(tr_off_y_ind(1));
tr_off_init2 = x2vec(tr_off_y_ind(ceil(length(tr_off_y_ind)/2)));
tr_off_end1 = x2vec(tr_off_y_ind(floor(length(tr_off_y_ind)/2)));
tr_off_end2 = x2vec(tr_off_y_ind(end));
plot([x1_init, x1_init+0.2], [tr_off_init1, tr_off_init1], "g-");
plot([x1_init, x1_init+0.2], [tr_off_init2, tr_off_init2], "g-");
plot([x1_init, x1_init+0.2], [tr_off_end1, tr_off_end1], "g-");
plot([x1_init, x1_init+0.2], [tr_off_end2, tr_off_end2], "g-");

plot([x1_end, x1_end+0.2], [tr_off_init1, tr_off_init1], "g-");
plot([x1_end, x1_end+0.2], [tr_off_init2, tr_off_init2], "g-");
plot([x1_end, x1_end+0.2], [tr_off_end1, tr_off_end1], "g-");
plot([x1_end, x1_end+0.2], [tr_off_end2, tr_off_end2], "g-");

axis off;

exportgraphics(gca, '../figure/Ztopo_processed_wpsdtraces.png', "Resolution", 600);

% Plot axis for the master plot

%%%%%%%%%%%%%%%%%%%
fig = figure(1);
fig.Units = 'point';
fig.Position = [0 500 600 600];
clf(fig,'reset'); cla(fig,'reset'); 
hold on; box on; axis equal;
colormap("viridis");

% [C, h] = contourf(x1vec, x2vec, Ztopo, 151);
% set(h,'LineColor','none');
% set(h, 'EdgeColor', 'none');

cb = colorbar;
cb.Label.String = "Height from host rock [μm] ";
cb.Label.FontSize = 16;
% shading interp
shading flat
caxis([0, 60]);

xlim([-10, 10]);
ylim([-10, 10]);

xlabel("x [mm]");
ylabel("y [mm]");

% % trace off the patch
% for i = tr_off_y_ind(1:7:end)
%     plot(x1vec([tr_x_init, tr_x_init+tr_npts]), x2vec(i) * ones(2), "y-");
% end
% 
% % trace on the patch
% for i = tr_on_y_ind(1:7:end)
%     plot(x1vec([tr_x_init, tr_x_init+tr_npts]), x2vec(i) * ones(2), "r-");
% end

% plot range of PSD measurement
% tr_on_init = x2vec(tr_on_y_ind(1));
% tr_on_end = x2vec(tr_on_y_ind(end));
% x1_init = x1vec(tr_x_init);
% plot([x1_init, x1_init+0.2], [tr_on_init, tr_on_init], "r-");
% plot([x1_init, x1_init+0.2], [tr_on_end, tr_on_end], "r-");
% 
% 
% tr_off_init1 = x2vec(tr_off_y_ind(1));
% tr_off_init2 = x2vec(tr_off_y_ind(ceil(length(tr_off_y_ind)/2)));
% tr_off_end1 = x2vec(tr_off_y_ind(floor(length(tr_off_y_ind)/2)));
% tr_off_end2 = x2vec(tr_off_y_ind(end));
% plot([x1_init, x1_init+0.2], [tr_off_init1, tr_off_init1], "g-");
% plot([x1_init, x1_init+0.2], [tr_off_init2, tr_off_init2], "g-");
% plot([x1_init, x1_init+0.2], [tr_off_end1, tr_off_end1], "g-");
% plot([x1_init, x1_init+0.2], [tr_off_end2, tr_off_end2], "g-");

axis on;

exportgraphics(gca, '../figure/Ztopo_processed_wpsdtraces_axis.eps', "Resolution", 150);

%% PSD of topography

window = 2^5;
noverlap = round(window/2);
nfft = 2^7;

pxx_all_on = zeros(length(tr_on_y_ind), nfft/2+1);
pxx_all_off = zeros(length(tr_off_y_ind), nfft/2+1);

%%
itrace_on = 0;
for ion = tr_on_y_ind
    tr_on_all = Ztopo(ion, :) * 1e-3;  %[um] -> [mm]
    tr_on = tr_on_all(tr_x_init:tr_x_init+tr_npts);% bounds within the target length
    [pxx_on,f_on] = pwelch(tr_on,window,noverlap,nfft, 1/param.ds);
    itrace_on = itrace_on+1;
    pxx_all_on(itrace_on, :) = pxx_on;
end

itrace_off = 0;
% figure(9);clf;box on; hold on; % debug
for ioff = tr_off_y_ind
    tr_off_all = Ztopo(ioff, :) * 1e-3; %[um] -> [mm]
    tr_off = tr_off_all(tr_x_init:tr_x_init+tr_npts);
    [pxx_off,f_off] = pwelch(tr_off,window,noverlap,nfft, 1/param.ds);
    itrace_off = itrace_off+1;
    pxx_all_off(itrace_off, :) = pxx_off;
    % if max(tr_off)>0.03
    %     ioff
    %     plot(tr_off);
    % end
end

assert(all(f_on==f_off));
f = f_on;


%% Compute noise level of Laser displacement sensor
% GPTOPO-02_002 records the noise of LDS by just hold the laser on the host
% rock for 1 minute. Read the data and compute its spectrum as the
% background noise.

fname_noise = 'GPTOPO-02_002';
[tmat_noise, Zmat_noise, ~]=TEACtopo(pname,fname_noise);


noise_inds = find(tmat_noise<60);
tmat_noise = tmat_noise(noise_inds);
Zmat_noise = Zmat_noise(noise_inds) - mean(Zmat_noise(noise_inds));

figure(9); clf; box on; hold on;
subplot(2, 1, 1);
plot(tmat_noise, Zmat_noise, "k-");
xlabel("Time [s]");
ylabel("Height [μm]");
title("Noise level of LDT")
exportgraphics(gca, '../figure/LDT_noiselevel.png', "Resolution", 150);

%%
subplot(2, 1, 2);
plot(x1vec, I.Zmat_gaussfiltered(:, tr_off_y_ind(5)), "b-");
xlabel("Location [mm]");
ylabel("Height [μm]");
title("Topography off patch")
%%
% To compute the noise level, we assume the LDT moves with param.ds on the
% perfectly smoothed surface such that only the noise of LDT appears on the data.
tr_noise_all = Zmat_noise * 1e-3;  %[um] -> [mm]
[pxx_noise,f_noise] = pwelch(tr_noise_all,window,noverlap,nfft, 1/param.ds);


%% plot PSD
fig = figure(3);
fig.Units = 'point';
fig.Position = [0 500 600 500];
clf(fig,'reset'); cla(fig,'reset'); 
hold on; box on; grid on;

lon = plot(f, pxx_all_on, "r-");
for l = 1:length(lon)
    lon(l).Color = [lon(l).Color 0.3];
end

loff = plot(f, pxx_all_off, "-", "Color", [0.1,0.1,0.1]);
for l = 1:length(loff)
    loff(l).Color = [loff(l).Color 0.1];
end

% stacked trace
pxx_stack_on = mean(pxx_all_on, 1);
pxx_stack_off = mean(pxx_all_off, 1);
lc_stack = "k"; % [0.4660 0.6740 0.1880];
plot(f, pxx_stack_on, "-", "LineWidth", 2.5, "Color", lc_stack);
plot(f, pxx_stack_off, "-", "LineWidth", 2.5, "Color", lc_stack);

xline(1/param.ds/2,':','Nyquist wavenumber', "FontSize", 14);

text(1.5e-1, 2e-3, "On gauge patch", "FontSize", 16);
text(1.5e-1, 7e-5, "On host rock", "FontSize", 16);

hleglines = [lon(1) loff(1)];
% legend(hleglines, 'On patch', 'Off patch');

%---plot noise level---%
plot(f_noise, pxx_noise, "--", "LineWidth", 1.5, "Color", "k");
text(1.5e-1, 1.2e-8, "Noise level of LDT", "FontSize", 16);

%---plot PSD of 1D von Karman type ACF---%
% Ref: Sato et al. eq. (2.13b)
vonKarman1D = @(m, epsi, a, kappa) (2 * sqrt(pi) * gamma(kappa+0.5) * epsi^2 * a)./(gamma(kappa)*(1+a^2*m.^2).^(kappa+0.5));
gaussian1D = @(m, epsi, a) epsi^2*pi^(1/2)*a*exp(-m.^2*a^2/4); % 1D Fourier transform of gaussian 

% We apply the normalization of residual

f_fit_max_on = length(f); % use all the datapoint
f_fit_max_off = find(f>1, 1); % use up to noise level

obj_fun_on = @(x) norm((vonKarman1D(f(1:f_fit_max_on), x(1), x(2), x(3))-pxx_stack_on(1:f_fit_max_on)));
obj_fun_off = @(x) norm((vonKarman1D(f(1:f_fit_max_off), x(1), x(2), x(3))-pxx_stack_off(1:f_fit_max_off)));
obj_fun_off_gauss = @(x) norm((gaussian1D(f(1:f_fit_max_off), x(1), x(2))-pxx_stack_off(1:f_fit_max_off)));
% 
x0_on=[0.03, 1.0, 0.4];
x0_off=[0.0025, 2.0, 0.95];
x0_off_gauss=[6e-3, 3.5];

options = optimset('Display','iter','TolFun', 1e-8,'TolX', 1e-4);

x_on  = fminsearch(obj_fun_on, x0_on, options);
x_off_von = fminsearch(obj_fun_off, x0_off, options);
x_off_gauss = fminsearch(obj_fun_off_gauss, x0_off_gauss, options);

lc_vk = [0.4660 0.6740 0.1880];
lc_g = [0.9290 0.6940 0.1250];

plot(f, vonKarman1D(f, x_on(1), x_on(2), x_on(3)), "--", "LineWidth", 3, "Color", "b");
% plot(f, vonKarman1D(f, x_off_von(1), x_off_von(2), x_off_von(3)), ":", "LineWidth", 3, "Color","b");
% plot(f, vonKarman1D(f, 0.0025, 1.0, 1.0), ":", "LineWidth", 3, "Color","b");
% plot(f, gaussian1D(f, x_off_gauss(1), x_off_gauss(2)), ":", "LineWidth", 3, "Color", lc_g);

%---------------------------------%
% decorate plot
% xline(f(f_fit_max),'--','Fit max range', "FontSize", 14);


xlim([0.1, 15]);
ylim([1*10^(-9.), 1e-2]);

set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');

xlabel("k_{x} [1/mm]");
ylabel("PSD of topography [mm^3]");
 
fo = fopen("../data/vonKalman_fit.txt", "w");
fprintf(fo, "on gauge von Kalman: epsi, a[mm], kappa=%f, %f, %f\n", x_on(1), x_on(2), x_on(3));
fprintf(fo, "off gauge von Kalman: epsi, a[mm], kappa=%f, %f, %f\n", x_off_von(1), x_off_von(2), x_off_von(3));
fprintf(fo, "off gauge Gaussian: epsi, a[mm]=%f, %f\n", x_off_gauss(1), x_off_gauss(2));
fclose(fo);
% 
exportgraphics(gca, '../figure/Ztopo_PSD_comparison.png', "Resolution", 300);
exportgraphics(gca, '../figure/Ztopo_PSD_comparison.eps');



% plot event data: stick-slip event 29
% render image and vector file for master figure
% update: 2025.02.12 merge vector and raster figure
% update: 2025.04.25 updated ylabels.

clear all; close all;

set(0,'DefaultTextFontsize',10, ...
    'DefaultTextFontname','Helvetica', ...
    'DefaultTextFontWeight','normal', ...
    'DefaultTextFontname','Helvetica', ...
    'DefaultAxesFontsize',10, ...
    'DefaultAxesFontname','Helvetica', ...
    'DefaultLineLineWidth', 1.0)
set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'})

figdir="../figure";
if ~exist(figdir) mkdir(figdir); end
addpath("../../utils/matlabcode_biax_v03");

%% Load the compiled data

ifSavefig = true;
ifPlotAE = true; % true if plotting AE
ifPlotStrain = true;

% ifPlotVector = true;
% ifPlotRaster = true;

% Range to compute rupture velocity
x_bound_typecheck = [1.25, 2.0];
xrange_east = [1.0, 3.0];
xrange_west = [0.5, 2.0];

expr_id = 87;
runID = sprintf('FB03_%03d', expr_id);
T = readtable(sprintf("../../Experiments/DetectEvent/data/p02_eventtype_fb03-%03d.csv", expr_id));

event_type = zeros(length(size(T, 1)), 3); %1. event type %2. rupture vel %3. event start time from strain
event_id = 29;
event_datdir="../data/eventdata_master";

load(event_datdir + sprintf("/eventdata_FB03_%03d_event%02d.mat", expr_id, event_id));

% apply low-frequency band
lowpass_winlen = 1e-3;
SS_HFS_filtered = movmean_causal(SSmacro, round(lowpass_winlen*Fs_strain));

%% Plot gouge event from visually picked datasheet
A = readtable("../../Experiments/DetectEvent/data/p06_visual_pick_gougeevents_merged.csv", 'NumHeaderLines',5);
Aevents = A((A.Var2==event_id) & (string(A.Var5)=={'Ordinary'}), :);

%% dump the plot parameter to the file
fo_param =  fopen(sprintf('../data/plotparam_eventid%03d.txt', event_id), 'w');
fprintf(fo_param, "lowpass Macroscopic shear: %.1f kHz\n", lowpass_winlen*Fs_strain/1e3);

%% Shift the start time
% Trim the time before the event nucleation
Tshift = 22e-3;
Tstart = Tstart + Tshift;

tmat_AE_event = tmat_AE_event - Tshift;
tmat_slip_event = tmat_slip_event - Tshift;
tmat_strain_event = tmat_strain_event - Tshift;
tmat_SS = tmat_macro - Tshift;
% Plot parameter
xlimit=[0, 30];

% if ~ifPlotStrain
%     lc_AE = "k";  %sns_colorpalette(1, 1);
% else
%     lc_AE = "k"; %[0.2,0.2,0.2,0.2];
% end

lc_AE = [0, 0, 0];
lw_strain = 0.9;
lc_strain = [200/255, 0, 0]; %[180/255, 28/255, 57/255]; %[200/255, 0, 0]; %sns_colorpalette(1, 4);
mc_strainpeak = [103, 171, 209]./255; %sns_colorpalette(1, 1); %"r";
lc_rupture = 'k';
mc_gougepatch = 'w'; sns_colorpalette(1, 9);
%--------------------------------%
%% 1. Plot accumulation of slip
%--------------------------------%

fig = figure(1); clf; hold on; box on;
fig.Units = 'point';
fig.Position = [0 800 700 400];

Tlength = event_winlen;

%--------------------------------------------%
%% 1. Plot AE and strain with rupture velocity
%--------------------------------------------%

subplot(12, 1, 1:10); hold on; box on;

% load event type
assert(Tstart-Tshift == T.event_starttime(event_id));

%--Plot AE---%
%%
ifBandpass = true;
bandpass_order = 3;
fmin = 0.1e6; %0.6e5;
fmax = 1.0e6; % 0.6e6
ampnorm_AE = 25; %10; %5e-5;

fprintf(fo_param, "Bandpass AE sensor: %.2f - %.2f MHz\n", fmin/1e6, fmax/1e6);
fprintf(fo_param, "Bandpass AE butterworth order: %d\n", bandpass_order);

% apply bandpass filter
if ifBandpass
    f_nyq = Fs_AE/2;
    [b,a] = butter(bandpass_order, [fmin fmax]/f_nyq, "bandpass");
    %     freqz(b, a, [], fs_read);
    %     ax=gca();
    %     set(ax, 'XScale', 'log');
    AEdatmat = filtfilt(b,a,AEdatmat);
end

Nsensors_AE = size(AEdatmat, 2);
[~, AE_sortedinds] = sort(AEsensor_x);

for i = 1:Nsensors_AE
    ii = AE_sortedinds(i);
    y_shift = AEsensor_x(ii)/1e3;

    if ifPlotAE % & ifPlotVector
        p_AE = plot(tmat_AE_event*1e3, AEdatmat(:, ii)/ampnorm_AE + y_shift, "-", "Color", lc_AE, "LineWidth", 0.5, "Tag", "AEwaveforms");
    end
end

% Plot scale of volt
% if ifPlotRaster
scale_x = xlimit(2) * 0.12; %22.5;
scale_y = 0.08;
scale_len = 5.0; %1; %[V]
plot([scale_x, scale_x], [scale_y-(scale_len/ampnorm_AE)/2, scale_y+(scale_len/ampnorm_AE)/2], "-", "Color", lc_AE, "LineWidth", 1);
text(scale_x+0.2, scale_y, sprintf("%.1fV", scale_len));
% end

%--Plot strain--%
% NOTE: we recompute the polyfit to superimpose the linear regression of rupture
% velocity to the plot.
% Then, we renew the rupture velocity of nucleation.

%%
NSGB = size(Taumat2, 2);
NSGT = size(Taumat3, 2);

%% Remove offset of strain
offsetind = 100;
Taumat2_removeoffset = Taumat2 - mean(Taumat2(1:offsetind, :));
Taumat3_removeoffset = Taumat3 - mean(Taumat3(1:offsetind, :));

%% Apply moving window average
SG_smooth_winlen = 100;
SG_smooth_winlen_t = SG_smooth_winlen/Fs_strain;
fpass = 1/SG_smooth_winlen_t;

fprintf(fo_param, "Lowpass strain gauge : %.2f MHz\n", fpass/1e6);

% Taumat2_removeoffset_smoothed = movmean(Taumat2_removeoffset, SG_smooth_winlen, 1);
% Taumat3_removeoffset_smoothed = movmean(Taumat3_removeoffset, SG_smooth_winlen, 1);

% Apply low-pass filter
[b, a] = butter(bandpass_order, fpass/(Fs_strain/2), "low");
fprintf("Strain low-pass filter with %4.4f[ms]\n" + ...
    "low-pass: %4.4f[kHz]\n", SG_smooth_winlen/Fs_strain*1e3, 1/SG_smooth_winlen_t/1e3);

% lowcutfreq_strain = 5.0; %[Hz]
% fprintf("Strain moving window smoothing with %4.4f[ms]\n" + ...
%     "bandpass: %4.4f-%4.4f[kHz]\n", SG_smooth_winlen/Fs_strain*1e3, lowcutfreq_strain*1e-3, 1/SG_smooth_winlen_t/1e3);

% Apply band-pass filter
% [b, a] = butter(2, [lowcutfreq_strain fpass]/(Fs_strain/2), "bandpass");

% Replace NaN to avoid the filter error
SGB_nanchan = isnan(Taumat2_removeoffset(1, :));
SGT_nanchan = isnan(Taumat3_removeoffset(1, :));
Taumat2_removeoffset(:, SGB_nanchan) = 0;
Taumat3_removeoffset(:, SGT_nanchan) = 0;

Taumat2_removeoffset_smoothed = filtfilt(b, a, Taumat2_removeoffset .* tukeywin(size(Taumat2_removeoffset, 1),0.05));
Taumat3_removeoffset_smoothed = filtfilt(b, a, Taumat3_removeoffset .* tukeywin(size(Taumat3_removeoffset, 1),0.05));

% refill NaN
Taumat2_removeoffset_smoothed(:, SGB_nanchan) = NaN;
Taumat3_removeoffset_smoothed(:, SGT_nanchan) = NaN;

%% Compute rupture pattern and rupture velocity
% 1. compute timing of maximum shear strain at trace
% 2. qualify if the S/N of max strain is larger than threshold
% 3. linear regression to evaluate the rupture initiation and its velocity

SNratio_threshold = 100;

[R, pvel, ruptype, event_type] = rupvel_linearfit(event_id, Tstart, SGB_x, SGT_x,...
    Taumat2_removeoffset_smoothed, Taumat3_removeoffset_smoothed,...
    tmat_strain_event, x_bound_typecheck, xrange_east, xrange_west, SNratio_threshold);


%% For fb03-087 event 29: additional rupture fitting

xrange_east_2 = [0.25, 0.75];

[~, pvel_2, ruptype_2, event_type_2] = rupvel_linearfit(event_id, Tstart, SGB_x, SGT_x,...
    Taumat2_removeoffset_smoothed, Taumat3_removeoffset_smoothed,...
    tmat_strain_event, x_bound_typecheck, xrange_east_2, xrange_west, SNratio_threshold);


%% Check if the estimated rupture velocity is close to p02
% NOTE: As we trimmed the strain traces and the moving window average
% results in a slightly different smoothing causing the shift of imax, the
% slope of polyfit is not identical to the p02.
% However, it should be close each other.
assert(T.eventtype(event_id) == ruptype)
if ~isnan(pvel(1))
    assert(abs(T.nuc_rupturevel(event_id) - pvel(1))/abs(T.nuc_rupturevel(event_id)) < 0.3); % assert if the difference is less than 30 %.
end
%--------------------------------------------%

%% Update 2024.06.08 Compute color contour of normalized shear stress
% To improve the visuallization, we compute the contour of shear strain
% release the memory
clear Taumat2 Taumat2_removeoffset;
clear Taumat3 Taumat3_removeoffset;

Z_strain = zeros(NSGB+NSGT, length(tmat_strain_event));

%%
[sensorloc_vec, I] = sort([SGB_x, SGT_x]);
% Normalize the amplitude by its absolute maximum value
for i = 1:length(I)
    ind = I(i);
    if ind >= 33
        % sensor of SGT
        strain_maxamp = max(abs(Taumat3_removeoffset_smoothed(:, ind-32)));
        Z_strain(i, :) = Taumat3_removeoffset_smoothed(:, ind-32)/strain_maxamp;
    else
        % sensor of SGB
        strain_maxamp = max(abs(Taumat2_removeoffset_smoothed(:, ind)));
        Z_strain(i, :) = Taumat2_removeoffset_smoothed(:, ind)/strain_maxamp;
    end
end

%% remove the nan sensor
nansensor = isnan(Z_strain(:, 1));

dx_sensor = mean(diff(sensorloc_vec));
nanlocs = sensorloc_vec(nansensor);

%%
Z_strain(nansensor, :) = [];
sensorloc_vec(nansensor) = [];


%% Make interpolate object
F = griddedInterpolant({sensorloc_vec, tmat_strain_event}, Z_strain, 'linear');

%%
plotstep = 10;

fprintf(fo_param, "Strain contour plot decimate: %d\n", plotstep);

xcq = sensorloc_vec;
ycq = tmat_strain_event(1:plotstep:end);
vcq = F({xcq, ycq});

%% debug plot
% fig = figure(2); clf; hold on; box on;
% fig.Units = 'point';
% fig.Position = [0 800 800 600];
% pcolor(yq*1e3, xq, vq);
%
% shading interp
% % mask nan sensors
% for i = 1:length(nanlocs)
%     rectangle('Position',[0 nanlocs(i)-dx_sensor/2 xlimit(2) dx_sensor],...
%         'FaceColor',[ 0.6, 0.6 ,0.6 ], 'EdgeColor','None') %[x y w h]
% end
% xlim(xlimit);
% colorbar;

%%
ampnorm = 0.4; %5e-5;

% if ifPlotRaster

% plot rupture velocity
% for i = 1:size(R, 1)
%     plot(R(i, 2)*1e3, R(i, 1), "s", "MarkerSize", 12, "MarkerFaceColor", mc_strainpeak, "MarkerEdgeColor",'k');
%     %     plot(R(i, 2)*1e3, R(i, 1) + R(i, 3)/ampnorm, "bo");
% end

% plot peak strain
% scatter(R(:, 2)*1e3, R(:, 1), 80, "MarkerFaceColor", mc_strainpeak, ...
% "MarkerEdgeColor", 'k', "Marker","square");

for ii = 1:NSGB
    y_shift = SGB_x(ii)/1e3;
    % p_SGB = plot(tmat_strain_event*1e3, Taumat2_removeoffset_smoothed(:, ii)/ampnorm + y_shift, "-", "Color", lc_strain, "LineWidth", lw_strain);
end

% plot SGT
for ii = 1:NSGT
    y_shift = SGT_x(ii)/1e3;
    % plot(tmat_strain_event*1e3, Taumat3_removeoffset_smoothed(:, ii)/ampnorm + y_shift, "-", "Color", lc_strain, "LineWidth", lw_strain);
end

rupture_propagation_slopeshift = 0.3;

if ruptype <= 2 % skip ruptyre == 3
    %     xq = linspace(min(R(:, 2)), max(R(:, 2)), 101);
    xq = linspace(12*1e-3, 21.1*1e-3, 101);
    yq = polyval(pvel, xq);
    plot(xq*1e3, yq-rupture_propagation_slopeshift, "-", "Color", lc_rupture, "LineWidth", 1);
end

% text(15, 1.958, sprintf("%4.1fm/s", abs(pvel(1))), "FontSize",17);
fprintf(fo_param, "first rupture velocity: %.1f m/s\n", abs(pvel(1)));

% Plot second rupture velocity
if ruptype <= 2 % skip ruptyre == 3
    xq_2 = linspace(21.4*1e-3, 22.1*1e-3, 101);
    yq_2 = polyval(pvel_2, xq_2);
    plot(xq_2*1e3, yq_2-rupture_propagation_slopeshift-0.3, "-", "Color", lc_rupture, "LineWidth", 1);
end
% text(20.6, 0.1, sprintf("%4.1fm/s", abs(pvel_2(1))), "FontSize",17);
fprintf(fo_param, "second rupture velocity: %.1f m/s\n", abs(pvel_2(1)));

% end

% xlabel("Time [ms]");
% ylabel("Distance from western edge of fault [m]");
ylabel("Easting [m]");
xlim(xlimit);
ylim([-0.1, 4.1]);

if ifPlotStrain
    % Plot scale of shear stress
    scale_x = xlimit(2) * 0.025; % 0.835; %23.5;
    scale_y = 0.1;
    scale_len = 0.1; %[MPa]
    % plot([scale_x, scale_x], [scale_y-(scale_len/ampnorm)/2, scale_y+(scale_len/ampnorm)/2], "-", "Color", lc_strain, "LineWidth", 2);
    % text(scale_x+0.1, scale_y, "0.1MPa");
end

%% Plot contour

% if ifPlotVector
pcolor(ycq*1e3, xcq/1e3, vcq, "Tag", "stress_contour");
% end

colormap(parula)
shading interp
% mask nan sensors
for i = 1:length(nanlocs)
    rectangle('Position',[0 (nanlocs(i)-dx_sensor/2)/1e3 xlimit(2)*1e3 dx_sensor/1e3],...
        'FaceColor',[ 0.6, 0.6 ,0.6 ], 'EdgeColor','None',  "Tag", "rec_nanstress") %[x y w h]
end

hax = gca;
hax.Children = circshift(hax.Children, -2);

% if ifPlotRaster
%     axis on;
%     xlabel("");
%     ylabel("");
%     xticks([]);
%     yticks([]);
%     set(gcf, 'color', 'none'); set(gca, 'color', 'none');
% end

% if ifPlotRaster
%% colorbar
cb = colorbar(hax);
cb.Limits=[-1, 1];
cb.Ticks = [-1, 0, 1];
cb.Label.String = "Normalized shear stress change";
pos = cb.Position;
pos(1) = 0.915;
pos(3) = 0.015;
pos(4) = 0.2;
cb.Position = pos;
%----------------------------------%
%----------------------------------%
%----------------------------------%

%% Plot gouge patch location

gougepatch_loc = 750 + (0:6)*500;
gouge_rec_x = 1.5;

% Plot fault rectangle for visuallization
faultwidth = 0.7;
faultlength = 4.1;
rec_fault = rectangle('Position',[gouge_rec_x-faultwidth/2 0 faultwidth faultlength],...
    'FaceColor', [0.6, 0.6, 0.6, 0.6], 'EdgeColor','k', "LineWidth", 0.5);

p_gouge = scatter(ones(7, 1)*gouge_rec_x, gougepatch_loc/1e3, 25, "Marker", "o", "MarkerEdgeColor","k",...
    "MarkerFaceColor", mc_gougepatch, "LineWidth",0.75);
% Child = get(p_gouge, 'Children');
% set(Child, 'Clipping', 'on');

%% Plot velocity reference
vp = 6.2;
vs = 3.6;
vscale_len = 2.2;
v_dt = 5; %[ms]
v_t0= 6.0; %[ms]
v_x0= 3.2; %[mm]

dx_vp = vp*v_dt;
dx_vs = vs*v_dt;
dx_vr = 0.92*vs*v_dt;

th_vp = atan2(dx_vp, v_dt);
th_vs = atan2(dx_vs, v_dt);
% th_vr = atan2(dx_vr, v_dt);

plot([v_t0, v_t0+vscale_len/tan(th_vp)], [v_x0, v_x0-vscale_len], "k-", "LineWidth", 0.75)
plot([v_t0, v_t0+vscale_len/tan(th_vs)], [v_x0, v_x0-vscale_len], "k-", "LineWidth", 0.75)
% plot([v_t0, v_t0+vscale_len/tan(th_vr)], [v_x0, -(v_x0+vscale_len)], "k-")

text(5.8, 0.92, "{\it c_p}", "FontSize",10);
text(6.5, 0.92, "{\it c_s}", "FontSize",10);


%% Annotate the gouge-mediated event
% annotation('arrow',[0.6923, 0.6755],[0.6135, 0.549]);
% annotation('arrow',[0.700581395348837 0.679263565891472],...
% [0.612094395280236 0.564896755162241],'HeadStyle','cback2');

% annotation('arrow',[0.456395348837209 0.469961240310078],...
%     [0.83 0.79]); %, 'color', 'b');
%
%
% annotation('arrow',[0.577519379844961 0.593023255813953],...
%     [0.765961651917404 0.72]); %, 'color', 'r');
%
% annotation('arrow',[0.652727272727273 0.661909090909091],...
%     [0.627686327077748 0.56],'LineWidth',1,'HeadWidth',8,...
%     'HeadLength',8);
%

annotation('arrow',[0.451076320939334 0.465037508153946],...
    [0.842787682333873 0.801729567029404],'LineWidth',0.75,'HeadWidth',4,...
    'HeadLength',4);

annotation('arrow',[0.521526418786691 0.535487606001303],...
    [0.760129659643434 0.719071544338966],'LineWidth',0.75,'HeadWidth',4,...
    'HeadLength',4);

annotation('arrow',[0.6445 0.655],...
    [0.599761013218786 0.558702897914318],'LineWidth',0.75,'HeadWidth',4,...
    'HeadLength',4);

Aevent = Aevents(2, :);
gougeevent_ID = "ID:"+string(Aevent.Var7);
text(18.5, 2.105, gougeevent_ID, "FontWeight", "bold");


%% Annotate event waveforms for zoomed panel
Torigin = 130.73454890; %gougeevetn id 61
rec_margin = 0.1;
Tevent = (Torigin - Tstart) * 1e3;
sensor_xmin = 1.2; % synchronize with the zoomed figure
sensor_xmax = 2.2;

rectangle('Position',[Tevent-rec_margin sensor_xmin 0.12+2*rec_margin sensor_xmax-sensor_xmin],...
    'FaceColor',"None", 'EdgeColor','k') %[x y w h]
%
% pause

%%

% for master figure
if ifBandpass
    % titlestr = sprintf("FB03-%03d event %02d: T%4.4f-%4.4f[s] Bandpass filtered: %.0f-%.0fkHz",...
    %     expr_id, event_id, Tstart, Tstart+xlimit(2)*1e-3, fmin/1e3, fmax/1e3);
    titlestr = sprintf("Stick-slip event %02d: T%4.4f-%4.4f [s]", event_id, Tstart, Tstart+xlimit(2)*1e-3);
else
    titlestr = " "; %sprintf("FB03-%03d event %02d: T%4.4f-%4.4f[s]",...
    %expr_id, event_id, Tstart, Tstart+xlimit(2)*1e-3);
end

title(titlestr, "FontWeight", "normal");

% if ifPlotStrain
%     legend([p_AE, p_SGB], ["AE waveform", "Shear stress change"], "Location", "northwest");
% else
%     legend([p_AE], ["AE waveform"], "Location", "northwest");
% end

% gca();
% pos = get(gca, 'Position');
% pos(1) = pos(1) + 0.03;
% set(gca, 'Position', pos)

set(gca,'XTickLabels',[])

subplot(12, 1, 11:12); hold on; box on;

plot(tmat_SS*1e3, SS_HFS_filtered, "k-", "LineWidth", 1);
xlabel("Time [ms]");
ylabel({"Macroscopic"; "shear stress [MPa]"});
xlim(xlimit);
ylim([0.638, 0.662]);
yticks([0.64, 0.65, 0.66]);

%% plot vline for the onset of main shock
% xline(23.525, "k--", "LineWidth",1.0);

init_SS_k = round(Tshift*Fs_strain)+1;
init_SS = SS_HFS_filtered(init_SS_k);

grad_SS = gradient(SS_HFS_filtered(init_SS_k:end));

if ismember(event_id, [1])
    onset_slope = 6e-6;
elseif ismember(event_id, [35, 41, 48])
    onset_slope = 5e-6;
else
    onset_slope = 3e-6;
end

onset_ind = init_SS_k+find(abs(grad_SS)>onset_slope, 1);
xline(tmat_SS(onset_ind)*1e3, "k--", "LineWidth",1.0);

%% Plot timing of mainshock
xline(hax, tmat_SS(onset_ind)*1e3, "k--", "LineWidth",1.0);

%% Plot annotation of stress drop
SS_init_ind = find((xlimit(1) <= tmat_SS), 1);
taup_macro = mean(SSmacro(SS_init_ind:SS_init_ind+100));
plot([27.1, 27.9], [taup_macro, taup_macro], "k-", "LineWidth", 0.75);
annotation('arrow',[0.84 0.84], [0.209 0.139],'LineWidth', 0.75,'HeadWidth',4,...
    'HeadLength',4);
text(27.65, 0.652, "\Delta\tau^{macro}", "Interpreter", "tex")
%% Add figure annotation
% annotation('textbox', [0.065, 0.95, 0, 0], 'string', '(a)', "FontWeight","bold", "FontSize", 17);
% annotation('textbox', [0.065, 0.26, 0, 0], 'string', '(b)', "FontWeight","bold", "FontSize", 17);

%% Add east and west annotation
annotation('textbox', [0.178, 0.928, 0, 0], 'string', 'E', "FontWeight","normal", "FontSize", 9.5);
annotation('textbox', [0.178, 0.295, 0, 0], 'string', 'W', "FontWeight","normal", "FontSize", 9.5);

% end
%% Save the figure separating the vector and raster images
% To decrease the size of figure, we save the contour and AE waveforms as
% raster.

% The following scripts are modified version originally from vecrast:
% Theo Michelis (2025). vecrast (https://www.mathworks.com/matlabcentral/fileexchange/154667-vecrast), MATLAB Central File Exchange. Retrieved February 12, 2025.
% Copyright (c) 2023, Theodoros Michelis
% All rights reserved.

figname_all = sprintf(figdir+'/plot_eventexample_master_fb03-%03d_event%02d_%d_wshearloading_filter%d_forFig2_master', expr_id, event_id, ifPlotStrain, ifBandpass);

figureHandle = fig;
stack = 'bottom';
% Set figure units to points
set(figureHandle, 'units', 'points');

% get(findall(figureHandle, 'type', 'axes'), 'position') % debug

% Ensure figure size and paper size are the same
figurePosition = get(figureHandle, 'Position');
set(figureHandle, 'PaperUnits', 'points', 'PaperSize', [figurePosition(3) figurePosition(4)])
set(figureHandle, 'PaperPositionMode', 'manual', 'PaperPosition', [0 0 figurePosition(3) figurePosition(4)]);

% get(findall(figureHandle, 'type', 'axes'), 'position') % debug

% Create a copy of the figure and remove smoothness in raster figure
rasterFigure = copyobj(figureHandle, groot);
vectorFigure = copyobj(figureHandle, groot);

% Fix vector image axis limits based on the original figure
% (this step is necessary if these limits have not been defined)
axesHandle = findall(vectorFigure, 'type', 'axes');
for i = 1:length(axesHandle)
    xlim(axesHandle(i), 'manual');
    ylim(axesHandle(i), 'manual');
    zlim(axesHandle(i), 'manual');
end

% Create axis in vector figure to fill with raster image
rasterAxis = axes(vectorFigure, 'color', 'none', 'box', 'off', 'units', 'points');
set(rasterAxis, 'position', [0 0 figurePosition(3) figurePosition(4)]);
uistack(rasterAxis, stack);


% Ensure fontsizes are the same in all figures
figText = findall(figureHandle, 'type', 'text');
rastText = findall(rasterFigure, 'type', 'text');
vecText = findall(vectorFigure, 'type', 'text');
for i=1:length(figText)
    set(rastText(i), 'FontSize', get(figText(i), 'FontSize'));
    set(vecText(i), 'FontSize', get(figText(i), 'FontSize'));
end

%%
% Raster Figure ----------------------------------------------------------
% Select what to remove from raster figure
axesHandleRaster = findall(rasterFigure, 'type', 'axes');
set(axesHandleRaster, 'color', 'none');

%%
% Here we select the stress contour and the AE waveforms to be plotted in raster

% Correct the size of contour axis considering the colorbar
correct1 = get(axesHandleRaster(1), 'position');
correct2 = get(axesHandleRaster(2), 'position');
correct2(3) = correct1(3); % set same size as the axis 1.
set(axesHandleRaster(2), 'position', correct2);

for i = 1:length(axesHandleRaster)
    axesPosition{i} = get(axesHandleRaster(i), 'position'); % Fix: get axes size
    contents = findall(axesHandleRaster(i));
    if i==2
        toKeep=[...
            findobj(axesHandleRaster(i).Children, 'Tag','AEwaveforms'); ...
            findobj(axesHandleRaster(i).Children, 'Tag','stress_contour'); ...
            findobj(axesHandleRaster(i).Children, 'Tag','rec_nanstress') ...            
            ];
    else
        toKeep = [];
    end
    toRemove = setxor(contents, toKeep);
    set(toRemove, 'visible', 'off');
end

%%
% Switch transparency of figure if a surface object is found in the figure.
% (This fix solves the raster print black edge lines of pcolor)
if isempty(findall(rasterFigure, 'type', 'surface'))
    set(rasterFigure, 'GraphicsSmoothing', 'off', 'color', 'none');
else
    set(rasterFigure, 'GraphicsSmoothing', 'off', 'color', 'w');
end

% Restore original axes size and font
% axesPosition
for i = 1:length(axesHandleRaster)
	set(axesHandleRaster(i), 'position', axesPosition{i});
end

% Remove all annotations, colorbars and legends from raster figure
annotations = findall(rasterFigure, 'Tag', 'scribeOverlay');
colorbarHandle = findall(rasterFigure, 'type', 'colorbar');
legendHandle = findall(rasterFigure, 'type', 'legend');
set([annotations; colorbarHandle; legendHandle], 'visible', 'off');

%%
% Print rasterFigure on temporary .png
% ('-loose' ensures that the bounding box of the figure is not cropped)
resolution = 300;
print(rasterFigure, [figname_all+'Temp.png'], '-dpng', ['-r' num2str(resolution) ], '-loose', '-image');
% close(rasterFigure);

%%

% Vector Figure -----------------------------------------------------------
% Select what to keep in vector figure
axesHandleVector = findall(vectorFigure, 'type', 'axes');
set(axesHandleVector, 'color', 'none');
for i = 1:length(axesHandleVector)
    axesPosition{i} = get(axesHandleVector(i), 'position'); % Fix: get axes size
    toRemove = [...
            findobj(axesHandleVector(i).Children, 'Tag','AEwaveforms'); ...
            findobj(axesHandleVector(i).Children, 'Tag','stress_contour'); ...
            findobj(axesHandleVector(i).Children, 'Tag','rec_nanstress') ...            
            ];
    set(toRemove, 'visible', 'off');
end

%%
% Restore original axes size
for i = 1:length(axesHandleVector) 
	set(axesHandleVector(i), 'position', axesPosition{i});
end

% Ensure Caxis limits and fonts match the original figure -----------------
figureAxesHandle = findall(figureHandle, 'type', 'axes');
for i = 1:length(figureAxesHandle)
    climOriginal = figureAxesHandle(i).CLim;
    axesHandle(i).CLim = climOriginal;

    zlimOriginal = figureAxesHandle(i).ZLim;
    axesHandle(i).ZLim = zlimOriginal;

    axesHandleVector(i).YAxis.FontSize = figureAxesHandle(i).YAxis.FontSize;
    axesHandleVector(i).XAxis.FontSize = figureAxesHandle(i).YAxis.FontSize;
    axesHandleVector(i).YLabel.FontSize = figureAxesHandle(i).YLabel.FontSize;
    axesHandleVector(i).XLabel.FontSize = figureAxesHandle(i).XLabel.FontSize;
end

% Insert Raster image into the vector figure
[A, ~, alpha] = imread([figname_all+'Temp.png']);

if isempty(alpha)==1
    imagesc(rasterAxis, A);
else
    imagesc(rasterAxis, A, 'alphaData', alpha);
end
axis(rasterAxis, 'off');

% Bring all annotations on top
% NOTE: THIS FOR SOME REASON CHANGES THE FONTS AGAIN
% annotations = findall(vectorFigure, 'Tag', 'scribeOverlay');
% for i = 1:length(annotations)
%     uistack(annotations(i), 'top');
% end

% Ensure figure has finished drawing
% drawnow;

% Finalise ----------------------------------------------------------------
% Remove raster image from directory
% delete([figname_png+'Temp.png']); % COMMENT THIS IF YOU WANT TO KEEP PNG FILE

%%
% Print and close the combined vector-raster figure
if ifSavefig
    print(vectorFigure, figname_all+".pdf", '-dpdf', '-loose', '-vector', ['-r' num2str(resolution) ]);
    % print(vectorFigure, figname_all+".eps", '-deps', '-loose', '-vector', ['-r' num2str(resolution) ]);
    print(vectorFigure, figname_all+".png", '-dpng', '-loose', '-vector', ['-r' num2str(300) ]);
end

% close(vectorFigure);
% plot event data: stick-slip all events
% 2025.3.19 update for master plot

clear all;
set(0,'DefaultTextFontsize',14, ...
    'DefaultTextFontname','Arial', ...
    'DefaultTextFontWeight','normal', ...
    'DefaultTextFontname','Arial', ...
    'DefaultAxesFontsize',14, ...
    'DefaultAxesFontname','Arial', ...
    'DefaultLineLineWidth', 1.0)

set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'})

figdir_jpg="../figure/jpg";
if ~exist(figdir_jpg) mkdir(figdir_jpg); end

addpath("../../../utils/matlabcode_biax_v03");

%% Load the compiled data

ifSavefig = true;
ifPlotAE = true; % true if plotting AE
ifPlotStrain = true;

%% Plot gouge event from visually picked datasheet
A = readtable("../../../Experiments/DetectEvent/data/p06_visual_pick_gougeevents_merged.csv", 'NumHeaderLines',5);

%% Range to compute rupture velocity
x_bound_typecheck = [1.25, 2.0];
xrange_east = [2.0, 3.0];
xrange_west = [0.5, 2.0];

SNratio_threshold = 300; % threshold of peak strain detection

expr_id = 87;
runID = sprintf('FB03_%03d', expr_id);
T = readtable(sprintf("../../../Experiments/DetectEvent/data/p02_eventtype_fb03-%03d.csv", expr_id));

for event_id = 1:size(T, 1) % loop for all the event

% event_id = 3;

fprintf("Start processing stick-slip ID %02d\n", event_id);

event_datdir = "/Volumes/Okuboetal2025_masterHDD/4mBIAX_eventdata_master/p03_eventdata_FB03_087";

load(event_datdir + sprintf("/eventdata_FB03_%03d_event%02d.mat", expr_id, event_id));

% apply low-frequency band
lowpass_winlen = 1e-3;
SS_HFS_filtered = movmean_causal(SSmacro, round(lowpass_winlen*Fs_strain));


%% Shift the start time
% Trim the time before the event nucleation
Tshift = 0;
Tstart = Tstart + Tshift;

tmat_AE_event = tmat_AE_event - Tshift;
tmat_slip_event = tmat_slip_event - Tshift;
tmat_strain_event = tmat_strain_event - Tshift;
tmat_SS = tmat_macro - Tshift;
% Plot parameter

% if ~ifPlotStrain
%     lc_AE = "k";  %sns_colorpalette(1, 1);
% else
%     lc_AE = "k"; %[0.2,0.2,0.2,0.2];
% end

lc_AE = "k";
lw_strain = 0.9;
lc_strain = [200/255, 0, 0]; %[180/255, 28/255, 57/255]; %[200/255, 0, 0]; %sns_colorpalette(1, 4);
mc_strainpeak = [103, 171, 209]./255; %sns_colorpalette(1, 1); %"r";
lc_rupture = 'k';
mc_gougepatch = 'w'; sns_colorpalette(1, 9);

gougeid_fontcolor = "b";

fig = figure(1); clf; hold on; box on;
fig.Units = 'point';
fig.Position = [0 800 1200 900];

%--------------------------------------------%
%% 1. Plot AE and strain with rupture velocity
%--------------------------------------------%
subplot(12, 1, 1:10); hold on; box on;

% load event type
assert(Tstart-Tshift == T.event_starttime(event_id));

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
SG_smooth_winlen = 1000;
SG_smooth_winlen_t = SG_smooth_winlen/Fs_strain;
fpass = 1/SG_smooth_winlen_t;

% Apply low-pass filter
[b, a] = butter(4, fpass/(Fs_strain/2), "low");
fprintf("Strain moving window smoothing with %4.4f[ms]\n" + ...
    "low-pass: %4.4f[kHz]\n", SG_smooth_winlen/Fs_strain*1e3, 1/SG_smooth_winlen_t/1e3);

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

% UPDATE 2024.03.14 use the movmean just for identifying the rupture type
% same with the p02_evalevent_strain.m
%% apply moving window average
smooth_winlen = 100;
Taumat2_removeoffset_smoothed_ruptype = movmean(Taumat2_removeoffset, smooth_winlen, 1);
Taumat3_removeoffset_smoothed_ruptype = movmean(Taumat3_removeoffset, smooth_winlen, 1);

[R, pvel, ruptype, event_type] = rupvel_linearfit(event_id, Tstart, SGB_x, SGT_x,...
    Taumat2_removeoffset_smoothed_ruptype, Taumat3_removeoffset_smoothed_ruptype,...
    tmat_strain_event, x_bound_typecheck, xrange_east, xrange_west, SNratio_threshold);

%% Check if the estimated rupture velocity is close to p02
% NOTE: As we trimmed the strain traces and the moving window average
% results in a slightly different smoothing causing the shift of imax, the
% slope of polyfit is not identical to the p02.
% However, it should be close each other.

assert(T.eventtype(event_id) == ruptype)
if ~isnan(pvel(1))
    assert(abs(T.nuc_rupturevel(event_id) - pvel(1))/abs(T.nuc_rupturevel(event_id)) < 0.5); % assert if the difference is less than 50 %.
end
%--------------------------------------------%


%% Evaluate the time shift
% For the supplementary plots of all stick-slip events, we shift the start 
% time to the earliest peak strain timing.

plot_pretrig_t = 20e-3; % plot from the plot_pretrig_t before the onset of preslip
onset_t_shift = (min(R(:, 2))-plot_pretrig_t)*1e3;

if onset_t_shift<0
    onset_t_shift = 5.0; % for some cases, the pretrigger should be short.
end

%--Plot AE---%
%%
ifBandpass = true;
bandpass_order = 3;
fmin = 0.1e6;
fmax = 1.0e6;
ampnorm_AE = 15; %25; %5e-5;

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

    if ifPlotAE
        p_AE = plot(tmat_AE_event*1e3-onset_t_shift, AEdatmat(:, ii)/ampnorm_AE + y_shift, "-", "Color", lc_AE, "LineWidth", 0.7);
    end
end

%% Plot strain
ampnorm = 0.4; %5e-5;

if ifPlotStrain

    % plot rupture velocity
    % for i = 1:size(R, 1)
    %     plot(R(i, 2)*1e3, R(i, 1), "s", "MarkerSize", 12, "MarkerFaceColor", mc_strainpeak, "MarkerEdgeColor",'k');
    %     %     plot(R(i, 2)*1e3, R(i, 1) + R(i, 3)/ampnorm, "bo");
    % end

    % R_plot = R((0.5<R(:, 1)) & (R(:, 1)<3.5), :); % bound the plotting peak strain
    % plot peak strain location
    R_plot = R;
    scatter(R_plot(:, 2)*1e3-onset_t_shift, R_plot(:, 1), 80, "MarkerFaceColor", mc_strainpeak, ...
        "MarkerEdgeColor", 'k', "Marker","square");

    % plot SGT
    for ii = 1:NSGB
        y_shift = SGB_x(ii)/1e3;
        p_SGB = plot(tmat_strain_event*1e3-onset_t_shift, Taumat2_removeoffset_smoothed(:, ii)/ampnorm + y_shift, "-", "Color", lc_strain, "LineWidth", lw_strain);
    end

    % plot SGT
    for ii = 1:NSGT
        y_shift = SGT_x(ii)/1e3;
        plot(tmat_strain_event*1e3-onset_t_shift, Taumat3_removeoffset_smoothed(:, ii)/ampnorm + y_shift, "-", "Color", lc_strain, "LineWidth", lw_strain);
    end

end

if ruptype == 2
    xlimit=[0, 120];
else
    xlimit=[0, 80];
end

% exceptions for x limit
if ismember(event_id, [11])
   xlimit=[0, 120];
end

if ismember(event_id, [14, 16, 18, 28, 36])
   xlimit=[0, 160];
end

% xlabel("Time [ms]");
% ylabel("Distance from western edge of fault [m]");
ylabel("Easting [m]");
xlim(xlimit);
ylim([-0.1, 4.1]);

% Plot scale of volt
scale_x = xlimit(2) * 0.1; %22.5;
scale_y = 0.1;
scale_len = 2.5; %1; %[V]
plot([scale_x, scale_x], [scale_y-(scale_len/ampnorm_AE)/2, scale_y+(scale_len/ampnorm_AE)/2], "-", "Color", lc_AE, "LineWidth", 2);
text(scale_x+0.1, scale_y, sprintf(" %.1fV", scale_len));

if ifPlotStrain
    % Plot scale of shear stress
    scale_x = xlimit(2) * 0.025; % 0.835; %23.5;
    scale_y = 0.1;
    scale_len = 0.1; %[MPa]
    plot([scale_x, scale_x], [scale_y-(scale_len/ampnorm)/2, scale_y+(scale_len/ampnorm)/2], "-", "Color", lc_strain, "LineWidth", 2);
    text(scale_x+0.1, scale_y, " 0.1MPa");
end


%---Plot timing of GP event---%
% select stick-slip experiment with ordinary events
Aevents = A((A.Var2==event_id) & (string(A.Var5)=={'Ordinary'}), :);

event_timing_before = 0;
for ievent = 1:size(Aevents, 1)
    Aevent = Aevents(ievent, :);
    event_timing = Aevent.Var4-Tstart;
    event_loc = Aevent.Var3;
    plot(event_timing*1e3-onset_t_shift, event_loc, "bo", "MarkerSize", 15, "LineWidth", 2);

    % abs(event_timing - event_timing_before)
    if ((abs(event_timing - event_timing_before)) < 2.0e-3) & (event_loc==event_loc_before) %1.0e-3
        % avoid overlap of text annotation
        text_shift = 3.6;
    else
        text_shift = 0;
    end

    ht = text(event_timing*1e3-onset_t_shift+text_shift, event_loc, string(Aevent.Var7)+"  ", "Color", gougeid_fontcolor, "FontSize", 18,...
         'HorizontalAlignment', 'right', 'FontWeight', "bold", 'VerticalAlignment', 'middle');
    % ht.Position(2) = ht.Position(2) - 1;
    event_timing_before = event_timing;
    event_loc_before = event_loc;
end

%% Plot gouge patch location
gougepatch_loc = 750 + (0:6)*500;
p_gouge = scatter(zeros(7, 1), gougepatch_loc/1e3, 80, "Marker", "o", "MarkerEdgeColor","k",...
    "MarkerFaceColor", mc_gougepatch, "LineWidth",1);
Child = get(p_gouge, 'Children');
set(Child, 'Clipping', 'on');

% annotate gouge patch name
GP_labels = "P" + string(1:7);
GP_label_x = 0.005 * (xlimit(2) - xlimit(1));
text(GP_label_x*ones(length(GP_labels), 1), gougepatch_loc/1e3-0.06, GP_labels, 'Units', 'data', 'HorizontalAlignment', 'left');

% %% Plot velocity reference
% vp = 6.2;
% vs = 3.6;
% vscale_len = 2.2;
% v_dt = 5; %[ms]
% v_t0= 6.0; %[ms]
% v_x0= 3.2; %[mm]
% 
% dx_vp = vp*v_dt;
% dx_vs = vs*v_dt;
% dx_vr = 0.92*vs*v_dt;
% 
% th_vp = atan2(dx_vp, v_dt);
% th_vs = atan2(dx_vs, v_dt);
% % th_vr = atan2(dx_vr, v_dt);
% 
% plot([v_t0, v_t0+vscale_len/tan(th_vp)], [v_x0, v_x0-vscale_len], "k-")
% plot([v_t0, v_t0+vscale_len/tan(th_vs)], [v_x0, v_x0-vscale_len], "k-")
% % plot([v_t0, v_t0+vscale_len/tan(th_vr)], [v_x0, -(v_x0+vscale_len)], "k-")
% 
% text(6, 0.92, "c_p", "FontSize",17);
% text(6.7, 0.92, "c_s", "FontSize",17);

%%

% titlestr = sprintf("FB03-%03d event %02d: T%4.4f-%4.4f[s]\nAE and Strain: rupvel=%.2f[m/s] ruptype=%d",...
%     expr_id, event_id, Tstart, Tstart+Tlength, T.nuc_rupturevel(event_id), T.eventtype(event_id));

% for master figure

titlestr = sprintf("Stick-slip event %02d: T%4.4f-%4.4f [s]", event_id, Tstart+onset_t_shift*1e-3, Tstart+xlimit(2)*1e-3+onset_t_shift*1e-3);

title(titlestr, "FontWeight", "normal");

if ifPlotStrain
    legend([p_AE, p_SGB], ["AE waveform", "Shear stress change"], "Location", "northwest");
else
    legend([p_AE], ["AE waveform"], "Location", "northwest");
end

% gca();
% pos = get(gca, 'Position');
% pos(1) = pos(1) + 0.03;
% set(gca, 'Position', pos)
set(gca,'XTickLabels',[])

%% plot macroscopic shear stress
subplot(12, 1, 11:12); hold on; box on;
plot(tmat_SS*1e3-onset_t_shift, SS_HFS_filtered, "k-", "LineWidth", 1);
xlabel("Time [ms]");
ylabel({"Macroscopic"; "Shear stress [MPa]"});
xlim(xlimit);
ylim([0.62, 0.7]);
% yticks([0.64, 0.65, 0.66]);

%% plot vline for the onset of main shock

% Update: read table from the GougeEventStats/Mainshocktiming

% init_SS_k = round(onset_t_shift*1e-3*Fs_strain)+1;
% init_SS = SS_HFS_filtered(init_SS_k);
% 
% grad_SS = gradient(SS_HFS_filtered(init_SS_k:end));
% 
% if ismember(event_id, [1])
%     onset_slope = 6e-6;
% elseif ismember(event_id, [35, 41, 48])
%     onset_slope = 5e-6;
% else
%     onset_slope = 3e-6;
% end
% onset_ind = init_SS_k+find(abs(grad_SS)>onset_slope, 1);

% read table
O = readtable(sprintf("../../../GougeEventStats/Mainshocktiming/data/MainshockTiming_FB03_%03d.csv", expr_id));

onset_t = O(event_id, :).mainshock_onset_all - Tstart;

xline(onset_t*1e3-onset_t_shift, "k--", "LineWidth",1.0);
% xline(onset_t*1e3, "b--", "LineWidth",1.0); % debug

%% Add figure annotation
% annotation('textbox', [0.065, 0.95, 0, 0], 'string', '(a)', "FontWeight","bold", "FontSize", 17);
% annotation('textbox', [0.065, 0.26, 0, 0], 'string', '(b)', "FontWeight","bold", "FontSize", 17);

%% Add east and west annotation
annotation('textbox', [0.095, 0.93, 0, 0], 'string', 'E', "FontWeight","normal", 'HorizontalAlignment', 'left');
annotation('textbox', [0.095, 0.275, 0, 0], 'string', 'W', "FontWeight","normal", 'HorizontalAlignment', 'left');

%%

figdir = "../figure";
if ~exist(figdir) mkdir(figdir); end

% figname = sprintf(figdir_jpg+"/plot_eventexample_master_fb03-%03d_event%02d_%d_wshearloading_filter%d.jpg", expr_id, event_id, ifPlotStrain, ifBandpass);
figname = sprintf(figdir_jpg+"/Supplementary_GPevents_FB03-%03d_stickslipevent%02d.jpg", expr_id, event_id);

if ifSavefig
    exportgraphics(fig, figname, 'Resolution',300);
end

% figname_eps = sprintf(figdir_eps+"/plot_eventexample_master_fb03-%03d_event%02d_%d_wshearloading_filter%d.eps", expr_id, event_id, ifPlotStrain, ifBandpass);
% if ifSavefig
%     exportgraphics(fig, figname_eps,'BackgroundColor','none','ContentType','vector');
% end
end
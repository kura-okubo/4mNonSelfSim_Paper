 % 01 plot static strain constraints for 1st revision rebuttal
clear all;
set(0,'DefaultTextFontsize',14, ...
    'DefaultTextFontname','Arial', ...
    'DefaultTextFontWeight','normal', ...
    'DefaultTextFontname','Arial', ...
    'DefaultAxesFontsize',14, ...
    'DefaultAxesFontname','Arial', ...
    'DefaultLineLineWidth', 1.0)

figdir="../figure";
if ~exist(figdir) mkdir(figdir); end

%% Read event data
% We focus on the largest target event: ID129 to show that
% S/N is not enough to apply the static strain inversion of GP event

load("../data/eventdata_FB03_087_event54.mat");


%% Read event onset table

T  = readtable("../../../../SourceInvFit//data/datacsv/AE_obs_location.csv");

expr_id = 87;
runID = sprintf('FB03_%03d', expr_id);

gougeevent_id = 129;

T_event = T(T.gougeevent_id == gougeevent_id, :);
T_origintime = T_event.origin_time;

%% Set constants
Fs_strain = 1e6; % sampling frequency of strain

NSGB=length(SGB_x);
NSGT=length(SGT_x);

%% Remove offset of strain
offsetind = 100;
Taumat2_removeoffset = Taumat2 - mean(Taumat2(1:offsetind, :));
Taumat3_removeoffset = Taumat3 - mean(Taumat3(1:offsetind, :));

%% Apply moving window average
SG_smooth_winlen = 100;
SG_smooth_winlen_t = SG_smooth_winlen/Fs_strain;
fpass = 1/SG_smooth_winlen_t;
dtpass = 1/fpass;
% Apply low-pass filter
[b, a] = butter(4, fpass/(Fs_strain/2), "low");
fprintf("Strain moving window smoothing with %4.4f[ms]\n" + ...
    "low-pass: %4.4f[kHz]\ndt: %4.1f[μs]\n", SG_smooth_winlen/Fs_strain*1e3, 1/SG_smooth_winlen_t/1e3, ...
    dtpass*1e6);

% Replace NaN to avoid the filter error
SGB_nanchan = isnan(Taumat2_removeoffset(1, :));
SGT_nanchan = isnan(Taumat3_removeoffset(1, :));
Taumat2_removeoffset(:, SGB_nanchan) = 0;
Taumat3_removeoffset(:, SGT_nanchan) = 0;

Taumat2_removeoffset_smoothed = filtfilt(b, a, Taumat2_removeoffset .* tukeywin(size(Taumat2_removeoffset, 1),0.01));
Taumat3_removeoffset_smoothed = filtfilt(b, a, Taumat3_removeoffset .* tukeywin(size(Taumat3_removeoffset, 1),0.01));

% refill NaN
Taumat2_removeoffset_smoothed(:, SGB_nanchan) = NaN;
Taumat3_removeoffset_smoothed(:, SGT_nanchan) = NaN;
%% Triming window
t_init = 0.0;
t_end = 0.115;
idx = find(tmat_strain_event >= t_init & tmat_strain_event <= t_end);

%% Plot strain time history (subplot version)
fig = figure(1);
fig.Units = 'point';
fig.Position = [0 500 1000 500];
clf(fig,'reset');

lc_SGB = "k";
lc_SGT = "k";

onset_t_shift = 15; % [ms]
plotamp = 5e-1;
plotind = 10; % span is needed, otherwise error in plotting due to too much data points.

% --- Subplot 1: full view ---
ax1 = subplot(1,2,1); hold on; box on;

% Plot SGB sensors
for i = 1:NSGB
    plot(ax1, ...
         tmat_strain_event(idx(1:plotind:end))*1e3 - onset_t_shift, ...
         Taumat2_removeoffset_smoothed(idx(1:plotind:end), i)/plotamp + SGB_x(i)/1e3, ...
         "Color", lc_SGB);
end

% Plot SGT sensors
for i = 1:NSGT
    plot(ax1, ...
         tmat_strain_event(idx(1:plotind:end))*1e3 - onset_t_shift, ...
         Taumat3_removeoffset_smoothed(idx(1:plotind:end), i)/plotamp + SGT_x(i)/1e3, ...
         "Color", lc_SGT);
end

ylabel(ax1, "Easting [m]");
xlabel(ax1, "Time [ms]");

% Set limits for full view
xlimit = [0 100];
xlim(ax1, xlimit);
ylim(ax1, [-0.1 4.1]);


% Plot scale of shear stress
scale_x = 10;
scale_y = 0.05;
scale_len = 0.1; %[MPa]
plot([scale_x, scale_x], [scale_y-(scale_len/plotamp)/2, scale_y+(scale_len/plotamp)/2], "-", "Color", lc_SGT, "LineWidth", 2);
text(scale_x+0.1, scale_y, " 0.1MPa");

% Plot AE timinig
T_GPevent = T_origintime-Tstart;
X_GPevent = T_event.X;
plot(T_GPevent*1e3-onset_t_shift, X_GPevent, "o", "MarkerSize", 10, "MarkerEdgeColor","r", "LineWidth",1);
plot([T_GPevent*1e3-onset_t_shift, T_GPevent*1e3-onset_t_shift], [-0.1, 4.2], "r--");

ylabel("Easting [m]");

title(ax1, "Full view", "FontWeight", "normal");

% --- Subplot 2: zoomed view (40–50 ms, 0–4 m) ---

zoom_init = 40;
zoom_end = 50;
idx_zoom = find(tmat_strain_event*1e3 - onset_t_shift >= zoom_init & tmat_strain_event*1e3 - onset_t_shift  <= zoom_end);


ax2 = subplot(1,2,2); hold on; box on;

plotind2 = 1;
plotamp2 = 2e-1;

% Plot SGB sensors
for i = 1:NSGB
    plot(ax2, ...
         tmat_strain_event(idx_zoom(1:plotind2:end))*1e3 - onset_t_shift, ...
         Taumat2_removeoffset_smoothed(idx_zoom(1:plotind2:end), i)/plotamp2 + SGB_x(i)/1e3, ...
         "Color", lc_SGB);
end

% Plot SGT sensors
for i = 1:NSGT
    plot(ax2, ...
         tmat_strain_event(idx_zoom(1:plotind2:end))*1e3 - onset_t_shift, ...
         Taumat3_removeoffset_smoothed(idx_zoom(1:plotind2:end), i)/plotamp2 + SGT_x(i)/1e3, ...
         "Color", lc_SGT);
end


% Plot AE timinig
% T_GPevent = T_origintime-Tstart;
% X_GPevent = T_event.X;
% plot(T_GPevent*1e3-onset_t_shift, X_GPevent, "o", "MarkerSize", 20, "MarkerEdgeColor","r", "LineWidth",1);
plot([T_GPevent*1e3-onset_t_shift, T_GPevent*1e3-onset_t_shift], [-0.1, 4.2], "r--");


ylabel(ax2, "Easting [m]");
xlabel(ax2, "Time [ms]");

% Set limits for full view
xlimit_zoom = [40, 50];
xlim(ax2, xlimit_zoom);
ylim(ax2, [-0.1 4.1]);


% Plot scale of shear stress
scale_x = 41;
scale_y = 0.2;
scale_len = 0.1; %[MPa]
plot([scale_x, scale_x], [scale_y-(scale_len/plotamp2)/2, scale_y+(scale_len/plotamp2)/2], "-", "Color", lc_SGT, "LineWidth", 2);
text(scale_x+0.1, scale_y, " 0.1MPa");

title(ax2, "Zoomed view (40–50 ms)", "FontWeight", "normal");

text(ax1, -0.15, 1.02, "(a)", ...
    "Units", "normalized", ...
    "FontWeight", "bold", ...
    "FontSize",16);
text(ax2, -0.15, 1.02, "(b)", ...
    "Units", "normalized", ...
    "FontWeight", "bold", ...
    "FontSize",16);

% --- Overall title (optional) ---
titlestr = sprintf("Stick-slip event %02d: Shear stress change\nT%4.4f-%4.4f [s] Low-pass filter at %.0f kHz", T_event.stickslip_id, Tstart+onset_t_shift*1e-3, Tstart+xlimit(2)*1e-3+onset_t_shift*1e-3, fpass/1e3);
sgtitle(titlestr);

figname = sprintf(figdir+"/StaticChange_Shearstress_FB03-%03d_stickslipevent%02d_ID%d.png", expr_id, T_event.stickslip_id, T_event.gougeevent_id);

exportgraphics(fig, figname, 'Resolution',300);
%%


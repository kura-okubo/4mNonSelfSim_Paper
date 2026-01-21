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

%% Remove offset of strain
offsetind = 100;
Dmat_event_removeoffset = Dmat_event - mean(Dmat_event(1:offsetind, :));

%% Triming window
t_init = 0.0;
t_end = 0.115;
idx = find(tmat_slip_event >= t_init & tmat_slip_event <= t_end);

%% Plot Slip
fig = figure(1);
fig.Units = 'point';
fig.Position = [0 500 1000 500];
clf(fig,'reset');

lc_D = "k";
NDisp = length(Disp_x);

onset_t_shift = 15; % [ms]
plotamp = 0.5e-1;
plotind = 1; % span is needed, otherwise error in plotting due to too much data points.

% --- Subplot 1: full view ---
ax1 = subplot(1,2,1); hold on; box on;

% Plot Gap sensors
for i = 1:NDisp
    plot(ax1, ...
         tmat_slip_event(idx(1:plotind:end))*1e3 - onset_t_shift, ...
         Dmat_event_removeoffset(idx(1:plotind:end), i)/plotamp + Disp_x(i)/1e3, ...
         "Color", lc_D);
end

ylabel(ax1, "Easting [m]");
xlabel(ax1, "Time [ms]");

% Set limits for full view
xlimit = [0 100];
xlim(ax1, xlimit);
ylim(ax1, [-0.1 4.1]);

%
% Plot scale of slip
scale_x = 10;
scale_y = 0.05;
scale_len = 0.01;
plot([scale_x, scale_x], [scale_y-(scale_len/plotamp)/2, scale_y+(scale_len/plotamp)/2], "-", "Color", lc_D, "LineWidth", 2);
text(scale_x+0.1, scale_y, sprintf(" %.0f μm", scale_len*1e3));

% Plot AE timinig
T_GPevent = T_origintime-Tstart;
X_GPevent = T_event.X;
plot(T_GPevent*1e3-onset_t_shift, X_GPevent, "o", "MarkerSize", 10, "MarkerEdgeColor","r", "LineWidth",1);
plot([T_GPevent*1e3-onset_t_shift, T_GPevent*1e3-onset_t_shift], [-0.1, 4.2], "r--");

ylabel("Easting [m]");

title(ax1, "Full view", "FontWeight", "normal");

%%
% --- Subplot 2: zoomed view (40–50 ms, 0–4 m) ---

zoom_init = 40;
zoom_end = 50;
idx_zoom = find(tmat_slip_event*1e3 - onset_t_shift >= zoom_init & tmat_slip_event*1e3 - onset_t_shift  <= zoom_end);

ax2 = subplot(1,2,2); hold on; box on;

plotind2 = 1;
plotamp2 = 0.1e-1;

% Plot Gap sensors
for i = 1:NDisp
    plot(ax2, ...
         tmat_slip_event(idx_zoom(1:plotind2:end))*1e3 - onset_t_shift, ...
         Dmat_event_removeoffset(idx_zoom(1:plotind2:end), i)/plotamp2 + Disp_x(i)/1e3, ...
         "Color", lc_D);
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


% Plot scale of slip
scale_x = 41;
scale_y = 0.5;
scale_len = 0.01;
plot([scale_x, scale_x], [scale_y-(scale_len/plotamp2)/2, scale_y+(scale_len/plotamp2)/2], "-", "Color", lc_D, "LineWidth", 2);
text(scale_x+0.1, scale_y, sprintf(" %.0f μm", scale_len*1e3));

title(ax2, "Zoomed view (40–50 ms)", "FontWeight", "normal");

% --- Overall title (optional) ---
titlestr = sprintf("Stick-slip event %02d: Slip \nT%4.4f-%4.4f [s] No filter", T_event.stickslip_id, Tstart+onset_t_shift*1e-3, Tstart+xlimit(2)*1e-3+onset_t_shift*1e-3);
sgtitle(titlestr);

figname = sprintf(figdir+"/StaticChange_Slip_FB03-%03d_stickslipevent%02d_ID%d.png", expr_id, T_event.stickslip_id, T_event.gougeevent_id);

exportgraphics(fig, figname, 'Resolution',300);



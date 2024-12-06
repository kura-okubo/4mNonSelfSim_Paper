% evaluate the local slip velocity at the gouge-mediated event
% 2024.3.11 Kurama Okubo

clear all; close all;
set(0,'DefaultTextFontsize',14, ...
    'DefaultTextFontname','Arial', ...
    'DefaultTextFontWeight','normal', ...
    'DefaultTextFontname','Arial', ...
    'DefaultAxesFontsize',14, ...
    'DefaultAxesFontname','Arial', ...
    'DefaultLineLineWidth', 1.0)

figdir="../figure/slipvel_M0_master";
if ~exist(figdir) mkdir(figdir); end

addpath("../../../utils/matlabcode_biax_v03");

%% Load the slip velocity data computed at 01_compute_slipvelocity_allevents.m

expr_id = 87;
runID = sprintf('FB03-%03d', expr_id);

ifSavefig = true;

%% Load slip and slipvelocity data
event_id = 29;
load(sprintf("../data/FB03-%03d_slipvelocity_event%02d.mat", expr_id, event_id));


% Time shift for the slip event 17
Tshift = 12e-3;%22e-3;
Tstart = Tstart + Tshift;

tmat_AE_event_decim = tmat_AE_event_decim - Tshift;
tmat_slip_event = tmat_slip_event - Tshift;
% Plot parameter
xlimit=[0, 80];

%% Load the AE event timing
M = readtable("../data/slipvelocity_and_M0_master.csv", 'VariableNamingRule', 'preserve');
df_event = M(M.stickslip_id==event_id, :);
gougeevent_ot = df_event.origin_time - Tstart; % Tstart is already shifted
gougeevent_x = 1.75; % we use the gouge patch location instead of the relocated hypocenter


%%
% lc_vel = sns_colorpalette(1, 1);

%%
fig = figure(1); clf; hold on; box on;
fig.Units = 'point';
fig.Position = [0 800 900 760];

N = length(Disp_x);

amp_norm_disp = 0.1;

% Plot slip
for i = 1:N
    plot(tmat_slip_event*1e3, Dmat_event_filtered(:, i)/amp_norm_disp + Disp_x(i)/1e3, "k--");
end

amp_norm_vel = 10;

% plot slip velocity
for i = 1:N
    plot(tmat_slip_event*1e3, SlipVelmat_filtered(:, i)/amp_norm_vel + Disp_x(i)/1e3, "-", "Color", "k", "LineWidth", 2.0);
end

% plot location of gap sensor
plot(zeros(length(Disp_x), 1), Disp_x/1e3, 's', 'Color', "k", 'MarkerSize', 10, "MarkerFaceColor", "w", "LineWidth", 1.0);

% plot the timing of gouge event
plot(gougeevent_ot*1e3, gougeevent_x, "ko", "MarkerSize", 10,...
    "MarkerFaceColor","w", "LineWidth", 1.5);


% Plot scale
scale_x_disp = 12;
scale_disp = 40e-3; %[mm]
plot([scale_x_disp, scale_x_disp], [Disp_x(N)/1e3, Disp_x(N)/1e3+scale_disp/amp_norm_disp], "k--", "Marker","_");
text(scale_x_disp, 4.1, sprintf(" %.0fμm", scale_disp*1e3));

scale_x_vel = 4;
scale_vel = 4.0; %[mm/s]
plot([scale_x_vel, scale_x_vel], [Disp_x(N)/1e3, Disp_x(N)/1e3+scale_vel/amp_norm_vel], "k-", "Marker","_");
text(scale_x_vel, 4.1, sprintf(" %.0fmm/s", scale_vel));


%% Plot velocity reference
vp = 6.2;
vs = 3.6;
vscale_len = 2.1;
v_dt = 5; %[ms]
v_t0= 15.0; %[ms]
v_x0= 3.2; %[mm]

dx_vp = vp*v_dt;
dx_vs = vs*v_dt;
dx_vs2 = 0.06*vs*v_dt;
dx_vr = 0.92*vs*v_dt;

th_vp = atan2(dx_vp, v_dt);
th_vs = atan2(dx_vs, v_dt);
th_vs2 = atan2(dx_vs2, v_dt);
% th_vr = atan2(dx_vr, v_dt);

% plot([v_t0, v_t0+vscale_len/tan(th_vp)], [v_x0, v_x0-vscale_len], "k-")
% plot([v_t0, v_t0+vscale_len/tan(th_vs)], [v_x0, v_x0-vscale_len], "k-")
plot([v_t0, v_t0+vscale_len/tan(th_vs2)], [v_x0, v_x0-vscale_len], "k-")
% plot([v_t0, v_t0+vscale_len/tan(th_vr)], [v_x0, -(v_x0+vscale_len)], "k-")

% text(6, 0.92, "c_p", "FontSize",17);
% text(6.7, 0.92, "c_s", "FontSize",17);
text(14, 1.82, "0.06c_s", "FontSize",17);



xlim(xlimit);
ylim([0, 4.5]);

xlabel("Time [ms]");
ylabel("Distance from western edge of fault [m]");

titlestr = sprintf("FB03-%03d stick-slip event %02d: T%4.4f-%4.4f[s]", expr_id, event_id, Tstart, Tstart+xlimit(2)*1e-3);
title(titlestr, "FontWeight", "normal");

%%

figname = sprintf(figdir+"/profile_slipvelocity_fb03-%03d_slipevent%02d.eps", expr_id, event_id);
exportgraphics(fig, figname, 'Resolution',80);

figname = sprintf(figdir+"/profile_slipvelocity_fb03-%03d_slipevent%02d.png", expr_id, event_id);
exportgraphics(fig, figname, 'Resolution',80);


% 04 plot event data
clear all;
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

expr_id = 87;
runID = sprintf('FB03_%03d', expr_id);
T = readtable(sprintf("../../Experiments/DetectEvent/data/p02_eventtype_fb03-%03d.csv", expr_id));

event_type = zeros(length(size(T, 1)), 3); %1. event type %2. rupture vel %3. event start time from strain
event_id = 29;

event_datdir="../data/eventdata_master";

load(event_datdir + sprintf("/eventdata_FB03_%03d_event%02d.mat", expr_id, event_id));

% %% Load gouge event timing
M = readtable("../../Others/LocalSlipVelocity/data/slipvelocity_and_M0_master.csv", 'VariableNamingRule', 'preserve');

df_event=M(M.stickslip_id==event_id, :);

%% Load event table from visually picked timing
% As the other events are not relocated yet, we use the visually picked
% time to annotate the timing
% The onset does not change in the order of 1ms, so it should be valid to
% investigate the incompleteness of the moment budget.

A = readtable("../../Experiments/DetectEvent/data/p06_visual_pick_gougeevents_merged.csv", 'NumHeaderLines',5);

Aevents = A((A.Var2==event_id) & (string(A.Var5)=={'Ordinary'}), :);

%% Add event at G5
% In the stick-slip event 29, a gouge event generated at G5 is not listed
% in the event table because the waveform is clipped due to the overload. 
% However, just for the sake of visuallization, we added this event in 
% the plot of slip accumulation.

T_ev29_G5 = 130.7299; % gouge event time for additional event on event 29.

C = {{'fb03-087'}, 29, 2.75, T_ev29_G5, {'Ordinary'}, 1, 502, 1, 501};
Aadditionalevents = cell2table(C, 'VariableNames', Aevents.Properties.VariableNames);
Aevents = [Aevents; Aadditionalevents];

%% Shift the start time
% Trim the time before the event nucleation
Tshift = 22e-3; % We syncronize this value to the master plot of AE and strain
xlimit=[0, 90]; %[ms]

Tstart = Tstart + Tshift;

Tshift_k = round(Tshift * Fs_slip);
Tend_k = round(xlimit(2) * 1e-3 * Fs_slip);

%%
tmat_slip_shifted = tmat_slip_event(Tshift_k+1:(Tshift_k+1) + Tend_k) - Tshift;
Dmat_event_shifted = Dmat_event(Tshift_k+1:(Tshift_k+1) + Tend_k, :);


%%
% Dmat_event_shifted = Dmat_event(:end)


% 
% lc_AE = 'k'; %sns_colorpalette(1, 1);
% lc_strain = [200/255, 0, 0]; %sns_colorpalette(1, 4);
% mc_strainpeak = sns_colorpalette(1, 1);
% lc_rupture = 'k';
mc_gougepatch = "w"; % sns_colorpalette(1, 9);

%--------------------------------%
%% 1. Plot accumulation of slip
%--------------------------------%

fig = figure(1); clf; hold on; box on;
fig.Units = 'point';
fig.Position = [0 800 350 300];

dt_event = 0.2e-3; %0.2e-3;
dt_event_k = round(dt_event*Fs_slip);
Tlength = xlimit(2) * 1e-3;
Tevent_k = round(Tlength*Fs_slip);

fprintf("Colored lines are plotted every %4.2f[ms]\n", dt_event*1e3);

% cd = colormap(flipud(viridis)); % take your pick (doc colormap)
cd = colormap(flipud(magma)); % take your pick (doc colormap)

lc_intp_event = @(x) interp1(linspace(0, Tlength+1e-8, length(cd)),cd, x);
%%

% remove offset in the begininng
Dmat_event_span = Dmat_event_shifted(1:dt_event_k:end, :);
Dmat_event_span = Dmat_event_span - mean(Dmat_event_span(1:5, :), 1);
Ndat_event = size(Dmat_event_span, 1);
tmat_event = tmat_slip_shifted(1:dt_event_k:end);


% for j = 1:Ndat_event
%     plot(Disp_x/1e3, Dmat_event_span(j, :), '-', 'Color', lc_intp_event(tmat_event(j)), "LineWidth", 1);
% end

for j = 1:Ndat_event

    Lw=1.0; 
    linecolor=lc_intp_event(tmat_event(j)); %0.5;
    
    plot(Disp_x/1e3, Dmat_event_span(j, :)*1e3, '-', 'Color', linecolor, "LineWidth", Lw);
end

% Plot bold line as the step
boldlinestep = 10; %[ms]
boldlinestep_k = round(boldlinestep*1e-3/dt_event);

for j = boldlinestep_k+1:boldlinestep_k:Ndat_event

    Lw=1.5;
    linecolor='k';  %3.5;

    plot(Disp_x/1e3, Dmat_event_span(j, :)*1e3, '-', 'Color', linecolor, "LineWidth", Lw);
    
    if tmat_event(j)*1e3<65
        text(3.88, Dmat_event_span(j, end)*1e3, sprintf("%.0f", tmat_event(j)*1e3), "FontSize", 9.5);
    end
end

text(3.88, 36.5, "ms", "FontSize", 9.);

%% plot gouge event timing
% mc_gougeevent = [255/255, 56/255, 16/255];
mc_gougeevent = '#FFD700'; %'#D55E00'; % foreshock_color

%1. relocated onset just for debug
gougeevent_ot = df_event.origin_time - Tstart; % Tstart is already shifted with the Tshift
gougeevent_x =  df_event.X; %for event 72 %1.75;
gougeevent_kt = find(tmat_event >= gougeevent_ot, 1); % use the t of the decimated data
gougeevent_accumulatedslip = Dmat_event_span(gougeevent_kt, :);
gougeevent_accumulatedslip_atpatch = interp1(Disp_x/1e3, gougeevent_accumulatedslip*1e3, gougeevent_x);
% plot(gougeevent_x, gougeevent_accumulatedslip_atpatch, "^", "MarkerEdgeColor","k",...
    % "MarkerFaceColor", "b", "MarkerSize", 18);
fprintf("Accumulated slip at gouge event: %.4f[um]\n", gougeevent_accumulatedslip_atpatch);


%% Plot gouge events from the visually picked onset
for ievent = 1:size(Aevents, 1)
    Aevent = Aevents(ievent, :);
    event_timing = Aevent.Var4-Tstart
    event_loc = Aevent.Var3;
    gougeevent_kt = find(tmat_event >= event_timing, 1); % use the t of the decimated data
    gougeevent_accumulatedslip = Dmat_event_span(gougeevent_kt, :);
    gougeevent_accumulatedslip_atpatch = interp1(Disp_x/1e3, gougeevent_accumulatedslip*1e3, event_loc)
    plot(event_loc, gougeevent_accumulatedslip_atpatch, "pentagram", "MarkerEdgeColor","k",...
        "MarkerFaceColor", mc_gougeevent, "MarkerSize", 9, "LineWidth", 0.75);
    text(event_loc+0.1, gougeevent_accumulatedslip_atpatch, sprintf("%d", Aevent.Var7));
end

%%

% plot sensor location
plot(Disp_x/1e3, zeros(length(Disp_x), 1), 's', 'Color', "k", 'MarkerSize', 6, "MarkerFaceColor",...
    [0.6, 0.6, 0.6], "LineWidth", 0.75);


% plot gougepatch
gougepatch_loc = 750 + (0:6)*500;
plot(gougepatch_loc/1e3, zeros(length(gougepatch_loc), 1)-1.25, 'o', 'Color', 'k',...
    'MarkerSize', 5, "MarkerFaceColor", mc_gougepatch, "LineWidth", 0.75);

xlabel("Easting [m]");
ylabel("Slip [μm]");

% title(sprintf("Slip: plot every %.2fms", dt_event*1e3));
% titlestr = sprintf("FB03-%03d stick-slip event %02d: T%4.4f-%4.4f[s]",...
    % expr_id, event_id, Tstart, Tstart+xlimit(2)*1e-3);
% title(titlestr, "FontWeight", "normal");

hcb = colorbar;
hcb.FontSize = 10;

ylabel(hcb, 'Time [ms]');
caxis([0, Tlength]*1e3);

% ref:https://www.mathworks.com/matlabcentral/answers/1449709-adjusting-width-of-horizontal-colorbar
ax = gca();
axPos = ax.Position;
pos = hcb.Position;
pos(1) = 0.80;
pos(3) = 0.03;
% pos(4) = 0.2;
hcb.Position = pos;
ax.Position = axPos;

xlim([0, 4.1]);
ylim([-2, 45]);
%     ylim([-0.001, 0.035]);

figdir = "../figure";
if ~exist(figdir) mkdir(figdir); end

% set(gcf, 'color', 'none');    
% set(gca, 'color', 'none');

figname = sprintf(figdir+"/plot_cumulativeslip_master_fb03-%03d_event%02d_forFig2.pdf", expr_id, event_id);

if ifSavefig
    print(fig, figname, '-dpdf', '-loose', '-vector');
    % exportgraphics(fig, figname);
end


figname = sprintf(figdir+"/plot_cumulativeslip_master_fb03-%03d_event%02d_forFig2.png", expr_id, event_id);

if ifSavefig
    exportgraphics(fig, figname, 'Resolution',80);
end

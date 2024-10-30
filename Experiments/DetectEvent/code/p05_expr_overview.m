% 05 plot overview of stick-slip experiments
clear all;
set(0,'DefaultTextFontsize',16, ...
    'DefaultTextFontname','Arial', ...
    'DefaultTextFontWeight','normal', ...
    'DefaultTextFontname','Arial', ...
    'DefaultAxesFontsize',16, ... 
    'DefaultAxesFontname','Arial', ...
    'DefaultLineLineWidth', 1.0)

figdir="../figure";
if ~exist(figdir) mkdir(figdir); end
addpath("../../../utils/matlabcode_biax_v03");

%% Load macro data
load("../../MacroData/data/MacroData_raw.mat");
load("../../MacroData/data/SlipData_raw.mat");

expr_id = 87;
ifSavefig = true;

Mex = M.(sprintf("FB03_%03d", expr_id));
Slip_ex = S.(sprintf("FB03_%03d", expr_id));
Fs_strain = 1/(Mex.tmat(2)-Mex.tmat(1));
Fs_slip = 1/(Slip_ex.tmat(2)-Slip_ex.tmat(1));

runID = sprintf('FB03_%03d', expr_id);
T = readtable(sprintf("../data/p02_eventtype_fb03-%03d.csv", expr_id));
Nevent = size(T,1);

eventtime = T.event_starttime;

%% 1. Macroscopic 
fig = figure(1); clf; hold on;
fig.Units = 'point';
fig.Position = [0 800 1000 800];

xlimit=[70, 130]; %[30, 260];
xline_width = 1.5;

t = tiledlayout(5, 1);
t.TileSpacing = 'compact';
t.Padding = 'loose';

nexttile; box on; hold on;

plot([eventtime eventtime], ylim(gca()), "k:");

plot(Mex.tmat, Mex.FCM, "k-");
xlim(xlimit);
ylim([0.395, 0.415]);
ylabel("Macroscopic τ/σ");

% Plot annotation for the detailed plots
if expr_id==87
    text(eventtime(14), 0.4105, " (i)", "FontWeight","bold");
    text(eventtime(17), 0.412, " (ii)", "FontWeight","bold");
    % text(eventtime(18), 0.411, " (ii)", "FontWeight","bold");
end

title(sprintf("FB03-%03d", expr_id));

%% 2. Drop of shear loading at event
% figure(2); clf; hold on;
% plot(Mex.tmat, Mex.SS, "k-");
% xlim(xlimit);

event_range = [0e-3, 170e-3]; % [s] window length from before to after the event
event_range_k = round(event_range * Fs_strain);
sheardrop = zeros(Nevent, 2);

for i=1:Nevent
% i = 12;
    st = T.event_starttime(i);
    ievent = find(Mex.tmat >=st, 1);
    SS_before = Mex.SS(ievent-event_range_k(1));
    SS_after = Mex.SS(ievent+event_range_k(2));

    SS_before_all(i, 1) = Mex.tmat(ievent-event_range_k(1));
    SS_before_all(i, 2) = SS_before;
    SS_before_all(i, 3) = Mex.FCM(ievent-event_range_k(1));

    SS_after_all(i, 1) = Mex.tmat(ievent+event_range_k(2));
    SS_after_all(i, 2) = SS_after;
    SS_after_all(i, 3) = Mex.FCM(ievent+event_range_k(2));

    sheardrop(i, 1) = st;
    sheardrop(i, 2) = SS_before - SS_after;
end

plot(SS_before_all(:, 1), SS_before_all(:, 3), "o", "MarkerSize", 8, "MarkerEdgeColor", "k", "MarkerFaceColor","w");
plot(SS_after_all(:, 1), SS_after_all(:, 3), "v", "MarkerSize", 8, "MarkerEdgeColor", "k", "MarkerFaceColor","w");


%% Plot rupture nucleation location
rupture_type = zeros(Nevent, 2);
rupture_type_char = strings(Nevent, 1);

for i=1:Nevent
    st = T.event_starttime(i);
    ruptype = T.eventtype(i);
    if ruptype == 1
        rup_y = 3;
        rup_char = "E";
    elseif ruptype == 2
        rup_y = 1;
        rup_char = "W";
    elseif ruptype == 3
        rup_y = 2;
        rup_char = "B";
    end

    rupture_type(i, 1) = st;
    rupture_type(i, 2) = rup_y;
    rupture_type_char(i, 1) = rup_char;
end

%%
nexttile; box on; hold on;

plot([eventtime eventtime], [0.5, 3.5], "k:");

ax = gca();
plot(eventtime, "k:", "LineWidth", xline_width);

plotind = (xlimit(1) < rupture_type(:, 1)) & (rupture_type(:, 1) <  xlimit(2));

ht = text(ax, rupture_type(plotind, 1), rupture_type(plotind, 2), rupture_type_char(plotind, 1),...
    "FontSize", 18, "FontWeight","bold", "HorizontalAlignment", "center", "BackgroundColor", "w");
xlim(xlimit);
ylim([0.5, 3.5]);
yticks([1, 2, 3]);
yticklabels([]);
ylabel("Nucleation location");


dim = [.87 .63 .06 .08];
str = {'East','Both','West'};
% annotation('textbox',dim,'String',str, "BackgroundColor", "w", ...
%     "HorizontalAlignment", "center", "VerticalAlignment","middle");

%% Plot rupture velocity of nucleation
col_E= sns_colorpalette(1, 4);
col_W= sns_colorpalette(1, 1);
col_B= sns_colorpalette(1, 8);

Trup = readtable(sprintf("../data/p02_eventtype_fb03-%03d.csv", expr_id)); % Need to update with p04: p04_renewedeventtype_fb03-%03d.csv

nuc_rupvel = zeros(Nevent, 3);

for i=1:Nevent
    st = Trup.event_starttime(i);
    nuc_rupvel(i, 1) = st;
    nuc_rupvel(i, 2) = sign(Trup.nuc_rupturevel(i));
    nuc_rupvel(i, 3) = abs(Trup.nuc_rupturevel(i));
end

%%
nexttile; box on; hold on;
plot([eventtime eventtime], [5, 500], "k:");

rupvel_E_ind = (nuc_rupvel(:, 2) == -1);
rupvel_W_ind = (nuc_rupvel(:, 2) == 1);
rupvel_B_ind = (isnan(nuc_rupvel(:, 2)));

% Plot rupture type 1 from east
h1 = plot(nuc_rupvel(rupvel_E_ind, 1), nuc_rupvel(rupvel_E_ind, 3), "s", "MarkerSize", 13, "MarkerEdgeColor", "k", "MarkerFaceColor", col_E);
h2 = plot(nuc_rupvel(rupvel_W_ind, 1), nuc_rupvel(rupvel_W_ind, 3), "^", "MarkerSize", 13, "MarkerEdgeColor", "k", "MarkerFaceColor", col_W);
set(gca, 'YScale', 'log');
% grid on;

% Plot average of east and west

mean_E_inds = (xlimit(1) < h1.XData) & (h1.XData < xlimit(2));
rupvel_mean_E = mean(h1.YData(mean_E_inds));
mean_W_inds = (xlimit(1) < h2.XData) & (h2.XData < xlimit(2));
rupvel_mean_W = mean(h2.YData(mean_W_inds));

yline(rupvel_mean_E, "-", sprintf("    %.1f m/s", rupvel_mean_E), "LabelHorizontalAlignment","left",...
    "FontSize", 14, "Color", col_E);
yline(rupvel_mean_W, "-", sprintf("    %.1f m/s", rupvel_mean_W), "LabelHorizontalAlignment","left",...
    "LabelVerticalAlignment", "top", "FontSize", 14, "Color", col_W);

xlim(xlimit);
ylim([5, 500]);
yticks([10, 100]);

ylabel({'Rupture velocity';'[m/s]'});

%%
nexttile; box on; hold on;
plot([eventtime eventtime], ylim(gca()), "k:");

h2 = bar(sheardrop(:, 1), sheardrop(:, 2), 1.2, "FaceColor", "flat", "EdgeColor","k", "LineWidth", 1.0);

for i = 1:Nevent
    if Trup.eventtype(i) == 1
        h2.CData(i, :) = col_E;
    elseif Trup.eventtype(i) == 2
        h2.CData(i, :) = col_W;
    else
        h2.CData(i, :) = col_B;
    end
end

xlim(xlimit);
ylim([0.01, 0.025]);
ylabel("Δτ [MPa]");

%% Slip
% We compute the average and standard deviation of slip from the gap
% sensors

event_range_slip = [0e-3, 300e-3]; % [s] window length from before to after the event
event_range_slip_k = round(event_range * Fs_slip);
average_slip = zeros(Nevent, 3);

for i=1:Nevent
% i = 12;
    st = T.event_starttime(i);
    ievent = find(Slip_ex.tmat >=st, 1);
    Slip_before = Slip_ex.Dmat(ievent-event_range_slip_k(1), :);
    Slip_after = Slip_ex.Dmat(ievent+event_range_slip_k(2), :);
    Slip_all = Slip_after - Slip_before;
    
    % Remove the broken sensor, which shows zero slip
    Slip_all(Slip_all<1e-3) = []; % ignore the slip less than 1μm
    average_slip(i, 1) = st;
    average_slip(i, 2) = mean(Slip_all);
    average_slip(i, 3) = std(Slip_all);
end


%%
nexttile; box on; hold on;
plot([eventtime eventtime], [25, 40], "k:");

% errorbar(average_slip(, 1), average_slip(:, 2)*1e3, average_slip(:, 3)*1e3, "ko",...
%     "MarkerSize", 8, "MarkerFaceColor", col, "MarkerEdgeColor","k", "LineWidth", 1.0,...
%     "CapSize", 8);

cols = [col_E; col_W; col_B];
inds = {rupvel_E_ind, rupvel_W_ind, rupvel_B_ind};
for i = 1:3
    plotinds = cell2mat(inds(i));
    errorbar(average_slip(plotinds, 1), average_slip(plotinds, 2)*1e3, average_slip(plotinds, 3)*1e3, "ko",...
    "MarkerSize", 12, "MarkerFaceColor", cols(i, :), "MarkerEdgeColor","k", "LineWidth", 1.0,...
    "CapSize", 8);
end

xlim(xlimit);
ylim([25, 40]);
ylabel("Slip [μm]");

%%
xlabel("Time [s]");

%% Reference: https://jp.mathworks.com/matlabcentral/fileexchange/41701-y-labels-alignment-in-subplots

allAxes = findall(fig,'type','axes');
xpos = [-38, -38, -25, -38, -38];
for k = 1:5
    allAxes(k).YLabel.Units='pixel';
    allAxes(k).YLabel.Position(1)=xpos(k);
end

%% annotation
annotation('textbox', [0.016, 0.95, 0, 0], 'string', '(a)', "FontWeight","bold");
annotation('textbox', [0.016, 0.77, 0, 0], 'string', '(b)', "FontWeight","bold");
annotation('textbox', [0.016, 0.60, 0, 0], 'string', '(c)', "FontWeight","bold");
annotation('textbox', [0.016, 0.43, 0, 0], 'string', '(d)', "FontWeight","bold");
annotation('textbox', [0.016, 0.26, 0, 0], 'string', '(e)', "FontWeight","bold");

%%
figname = sprintf("../figure/p05_overview_FB03-%03d.png", expr_id);

if ifSavefig
    exportgraphics(fig, figname, 'Resolution',300);
end

%% Save eps
figname = sprintf("../figure/p05_overview_FB03-%03d.eps", expr_id);

if ifSavefig
    exportgraphics(fig, figname);
end

%% save macroscopic stats
% save(sprintf("../data/Macroscopic_stats_%s.mat", runID));
Stickslip_eventID = (1:Nevent)';
Event_starttime = T.event_starttime;
Nucleation_loc = rupture_type_char;
Rupture_velocity = nuc_rupvel(:, 3);
Stressdrop = sheardrop(:, 2);
Slip_mean = average_slip(:, 2);
Slip_std = average_slip(:, 3);

Tstats = table(Stickslip_eventID,Event_starttime,Nucleation_loc,Rupture_velocity,Stressdrop,Slip_mean,Slip_std);

write(Tstats, sprintf("../data/Macroscopic_stats_%s.csv", runID));
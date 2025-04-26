% Script to plot the GP catalog with normal stress distribution on the
% fault
% 2025.03.18 updated for master plot; add foreshock/aftershock label


clear all;
set(0,'DefaultTextFontsize',10, ...
    'DefaultTextFontname','Arial', ...
    'DefaultTextFontWeight','normal', ...
    'DefaultTextFontname','Arial', ...
    'DefaultAxesFontsize',10, ... 
    'DefaultAxesFontname','Arial', ...
    'DefaultLineLineWidth', 1.0)

figdir="../figure";
if ~exist(figdir) mkdir(figdir); end
addpath("../../../utils/matlabcode_biax_v03");

set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'})


%% Load decimated data
expr_id = 87;
D_all = load(sprintf("../../../Experiments/MacroData/data/DecimatedData_raw_paper.mat"), "D");

D = D_all.D.(sprintf("FB03_%03d", expr_id)); % extract experiment data

%% Load catalog

% Read the picked csv file to superimposed the detected events.
A = readtable("../../../Experiments/DetectEvent/data/p06_visual_pick_gougeevents_merged.csv", 'NumHeaderLines',5);

%%
% filter by expr_id
idx_expr = string(A.Var1) == sprintf("fb03-%03d", expr_id);
A_expr = A(idx_expr,:);

%% Update 2025.3.19 add foreshock/aftershock label

% Aftershock list visually detected to double check the label
% Note that this reference includes the events on the all gouge patches
ref_aftershock_ids = [3,4,7,15,16,17,32,33,43,46,51,52,56,57,62,63,67,68,69,76,77,80,87,89,...
                  90,92,93,100,102,103,104,105,106,107,108,109,111,114,116,118,120,123,124,128,132]; % we visually categorized the aftershocks.

if ~ismember("eventlabel", A_expr.Properties.VariableNames)
    A_expr.eventlabel = strings(height(A_expr), 1); % Initialize as empty strings
end

for i = 1:height(A_expr)
    event_id = A_expr.Var7(i);  % Extract the event_id correctly
    event_type = A_expr.Var5(i);
    if strcmp(event_type{1}, 'LFE')
        A_expr.eventlabel(i) = "N";
    elseif ismember(event_id, ref_aftershock_ids)
        A_expr.eventlabel(i) = "A";  % Assign value correctly
    else
        A_expr.eventlabel(i) = "F";
    end
end


%% Read rupture type
B = readtable(sprintf("../../../Experiments/DetectEvent/data/p02_eventtype_fb03-%03d.csv", expr_id));

% check if the eventtype is consistent with the picked csv datasheet.
for i = 1:B.event_id(end)
    ruptype_event = A.Var6(A.Var2 == i);
    if ~isempty(ruptype_event)
        assert(all(ruptype_event == B(i, "eventtype").Variables));
    end
end

%% reorder table to plot stacked bar
event_c = zeros(7, 2);
gouge_loc = [0.75, 1.25, 1.75, 2.25, 2.75, 3.25, 3.75];
% ruptypes = ["Ordinary", "LFE"];
eventlabels = ["F", "A"];

for i = 1:length(gouge_loc)
    for j = 1:length(eventlabels)
        idx_loc = (A_expr.Var3 == gouge_loc(i)) & (string(A_expr.eventlabel) == eventlabels(j));
        A_loc = A_expr(idx_loc,:);
        count = size(A_loc, 1);
        event_c(i, j) = count;
    end
end

%% Plot normal stress
% Coordinate of strain gauges
% Load event data, including shifted strain gouges coordinates
C = load("../../../PlotEvent/data/eventdata_master/eventdata_FB03_087_event29.mat");

SGB_x = C.SGB_x; % biaxial
SGT_x = C.SGT_x; % triaxial

% %% Correct the relative location of the top rock specimen
% if expr_id == 87
%     xtop_shift = 82.46; %u_east (at t=0) = 13.57mm, u0 = -3.97
% else
%     error("the relative location of top block is not defined.");
% end
% 
% % adjust to the origin of the bottom rock specimen (west)
% SG3_x = xtop_shift + [170:240:4000, 230:240:4000];
% SG2_x = xtop_shift + [290:240:4000, 110:240:3900];

%% Remove the offset
offsetrange = 10;

D.Snmat_offsetremoved = D.Snmat - mean(D.Snmat(1:offsetrange, :)); 
D.Spmat_offsetremoved = D.Spmat - mean(D.Spmat(1:offsetrange, :)); 
D.Taumat2_offsetremoved = D.Taumat2 - mean(D.Taumat2(1:offsetrange, :)); 
D.Taumat3_offsetremoved = D.Taumat3 - mean(D.Taumat3(1:offsetrange, :)); 


%% Plot normal stress along fault
fig = figure(1);
fig.Units = 'point';
fig.Position = [0 500 500 450];
clf(fig,'reset'); cla(fig,'reset'); hold on; box on; 

% subplot(2, 1, 1); hold on; box on;
tiledlayout(2,1, 'Padding', 'loose', 'TileSpacing', 'compact'); 

mc_gougepatch = sns_colorpalette(1, 9);


%-----------------%
% Plot the catalog
%-----------------%
% subplot(2, 1, 2); hold on; box on; grid on;
nexttile; box on; hold on;

% grid on;
% b = bar(gouge_loc, event_c, "stacked", "EdgeColor", "k", "LineWidth",1);
% Plot only ordinary events
b = bar(gouge_loc, event_c, "stacked", "EdgeColor", "k", "LineWidth",0.75);

b(1).BarWidth = 0.31;
b(1).FaceColor = [0.7 0.7 0.7]; %sns_colorpalette(1, 1);
b(2).FaceColor = "w"; %sns_colorpalette(1, 2);
b(1).FaceAlpha = 1.0;
b(2).FaceAlpha = 0.75;

xlim([0, 4.1]);
ylim([0, 50]);

% Plot gouge patch location
p = plot(gouge_loc, 0.2 * ones(length(gouge_loc), 1), "o", "MarkerFaceColor", mc_gougepatch, ...
    "MarkerEdgeColor", "k", "MarkerSize", 7);

% Annotate patch ID
GPlabels = "P" + string(1:7);
gouge_loc = [0.75, 1.25, 1.75, 2.25, 2.75, 3.25, 3.75];
text(gouge_loc - 0.22, 1.5 * ones(size(gouge_loc)), GPlabels);

% legend([b(1), b(2), p], ["Ordinary event", "LFE", "gouge patch"], "Location", "Northwest");
legend([b(1),  b(2), p], ["Foreshock", "Aftershock", "GP Location"], "Location", "Northwest");

xlabel("Easting [m]");
ylabel("Event count");

% grid on;

%---------------------%
% plot normal stress
%---------------------%

nexttile;
ylimit = [0, 5];
sn_apply = 2; 
bar_width = 0.5;

hold on; box on; %grid on;

% nt = 97658;
Tplot = 150.0; % [s]
nt = find(D.tmat > Tplot, 1); % plot the time at 100

h1 = bar(SGT_x(17:end)/1e3, D.Snmat(nt, 17:end), bar_width, "FaceColor", sns_colorpalette(1, 4), "FaceAlpha", 0.3, "DisplayName", "SGT 17-32 (north)");
h2 = bar(SGT_x(1:16)/1e3, D.Snmat(nt, 1:16), bar_width, "FaceColor", sns_colorpalette(1, 1), "FaceAlpha", 0.3, "DisplayName", "SGT 1-16 (south)");

title(sprintf("fb03-087: Time at %.1f [s]", D.tmat(nt)));
xlabel("Easting [m]");
ylabel("Normal stress [MPa]");
xlim([0, 4.1]);
ylim(ylimit);

set(gcf, 'Color', 'w');
% exportgraphics(gcf, sprintf("../figure/normalstress_FB03_087.png"), "Resolution", 150);

% Compute average of north and south

Snmat_avg = zeros(16, 2);

for i = 1:16
    Snmat_avg(i, 1) = mean([SGT_x(i), SGT_x(i+16)]);
    Snmat_avg(i, 2) = mean([D.Snmat(nt, i), D.Snmat(nt, i+16)]);
end

h3 = plot(Snmat_avg(:, 1)/1e3, Snmat_avg(:, 2), "s-", "Color","k", "MarkerSize", 7,...
    "MarkerFaceColor", sns_colorpalette(1, 1), "LineWidth", 1, "DisplayName", "Mean");

% Plot gouge patch location
h4 = plot(gouge_loc, 0*ones(length(gouge_loc), 1),  "o", "MarkerFaceColor", mc_gougepatch, ...
    "MarkerEdgeColor", "k", "MarkerSize", 7, "DisplayName", "GP Location");

% annotate GP name
% text(gouge_loc - 0.16, 0.16 * ones(size(gouge_loc)), GPlabels);

legend([h2, h1, h3, h4], "Location", "northwest");


% Plot figure labels
text(-1, 3, 'a', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'right', 'Position', [-0.08, 2.37], 'Units', 'normalized');
text(-1, 3, 'b', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'right', 'Position', [-0.08, 1.07], 'Units', 'normalized');


set(gcf, 'Color', 'w');
exportgraphics(gcf, sprintf("../figure/FigS_GPactivity_fb03-%03d.png", expr_id), "Resolution", 70);
exportgraphics(gcf, sprintf("../figure/FigS_GPactivity_fb03-%03d.pdf", expr_id));

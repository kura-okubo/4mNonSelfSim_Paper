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

% Load the local slip velocity
E = readtable(sprintf("../data/localslipvel_fb03-%03d.csv", expr_id), 'VariableNamingRule', 'preserve');
Nevent = size(E, 1);

% add the column of seismic moment
E.("M0_best") = zeros(Nevent, 1);
E.("slip") = zeros(Nevent, 1);


%% Load the inverted seismic moment
% M = readtable(sprintf("../../SourceInversion/data/datacsv/gridsearch_bestparam_M0andTR_fb03-%03d.csv", expr_id));
% M = readtable(sprintf("../../EnergyBudget/data/gougeevent_sourceparameters_withLc.csv"));
waterlevel = 0.3;
denoisemethod = "detrend";
M_all = readtable(sprintf("../../../ComputeScaling/data/05_STFstats/SourceParam_meanstd_fb03-087_G3_wlv_%.2f_denoisemethod_%s.csv", waterlevel, denoisemethod));
M = M_all(M_all.Qinv_quart==50, :);


%% set label of foreshock and aftershock
E.("eventtiming") = zeros(Nevent, 1);

% aftershock_ids = [3, 35, 43, 52, 58, 71, 82, 90, 91, 97, 99]; % we visually categorized the aftershocks.
aftershock_ids = [3,4,7,15,16,17,32,33,43,46,51,52,56,57,62,63,67,68,69,76,77,80,87,89,90,92,93,100,102,103,104,105,106,107,108,109,111,114,116,118,120,123,124,128,132];

for i = 1:Nevent
    if ismember(E(i, :).gougeevent_id, aftershock_ids)
        E(i, :).("eventtiming") = 1;
    end
end

%% syncronize the tables
assert(Nevent == size(M, 1));
%%
for i = 1:Nevent
    gougeevent_id = M(i, :).gougeevent_id;
    E(E.gougeevent_id==gougeevent_id, :).("M0_best") = M(i, :).M0_mean;
end

%%
ctype = 1;
c_foreshock = sns_colorpalette(ctype, 4);
c_aftershock = sns_colorpalette(ctype, 1);

%%
M02Mw = @(x) (log10(x) - 9.105) * 2.0 / 3.0;

%%
fig = figure(1); clf; hold on; box on;
fig.Units = 'point';
fig.Position = [0 800 600 500];

E_foreshock = E(E.eventtiming==0, :);
E_aftershock = E(E.eventtiming==1, :);

ms = 10;

ylimit = [-2.0 0.5];

plot(E_foreshock.("slipvelocity[mm/s]"), log10(E_foreshock.("M0_best")),...
    "Marker", "v", "LineStyle","None", "MarkerSize", ms, "MarkerEdgeColor", "k", "MarkerFaceColor", c_foreshock);

plot(E_aftershock.("slipvelocity[mm/s]"), log10(E_aftershock.("M0_best")),...
    "Marker", "o", "LineStyle","None", "MarkerSize", ms, "MarkerEdgeColor", "k", "MarkerFaceColor", c_aftershock);

xlabel("Local slip velocity [mm/s]");
ylabel("log_{10}(M_0) [Nm]");
xlim([-0.5, 4.0]);
ylim(ylimit);
% yticks(-2.0:0.2:0.2)

% Annotate the 3 events for the case study
%%
% casestudy_events = [30, 91, 51];
% for ce_id = casestudy_events
%     E_ce = E(E.gougeevent_id==ce_id, :);
%     if E_ce.eventtiming==0
%         marker = "v";
%     else
%         marker = "o";
%     end
% 
%     % plot(E_ce.("slipvelocity[mm/s]"), log10(E_ce.("M0_best")),...
%     % "Marker", marker, "LineStyle","None",'HandleVisibility','off', "MarkerSize", ms, "MarkerEdgeColor", "k",...
%     % "LineWidth", 2, "MarkerFaceColor", "None");
% 
% end
% 
% annotation(fig,'textarrow',[0.471666666666667 0.45],[0.479 0.502],...
%     'String',{'ID:30'});
% annotation(fig,'textarrow',[0.4 0.411666666666667],[0.623 0.658],...
%     'String',{'ID:91'});
% annotation(fig,'textarrow',[0.586666666666667 0.565],...
%     [0.739 0.766],'String',{'ID:51'});


%%

yyaxis right;

ax = gca();
ax.YColor = 'k';

ms_2 = 12;

h1 = plot(E_foreshock.("slipvelocity[mm/s]"), M02Mw(E_foreshock.("M0_best")),...
    "Marker", "x", "LineStyle","None", "MarkerSize", ms_2, "MarkerEdgeColor", "y", "MarkerFaceColor", "k");

h2 = plot(E_aftershock.("slipvelocity[mm/s]"), M02Mw(E_aftershock.("M0_best")),...
    "Marker", "+", "LineStyle","None", "MarkerSize", ms_2, "MarkerEdgeColor", "y", "MarkerFaceColor", "k");

ax.XLim=[-0.5, 4.0];
ax.YLim=[M02Mw(10^ylimit(1)),M02Mw(10^ylimit(2))];
ax.YLabel.String="Mw";

h1.Visible = "off";
h2.Visible = "off";

%%
legend(["Foreshock", "Aftershock", "", ""] ,'Location', "southeast", "FontSize", 15);

figname = sprintf(figdir+"/slipvelocity_and_M0_fb03-%03d.eps", expr_id);
exportgraphics(fig, figname, 'Resolution',80);

figname = sprintf(figdir+"/slipvelocity_and_M0_fb03-%03d.png", expr_id);
exportgraphics(fig, figname, 'Resolution',80);

%% Dump the table
writetable(E, "../data/slipvelocity_and_M0_master.csv");



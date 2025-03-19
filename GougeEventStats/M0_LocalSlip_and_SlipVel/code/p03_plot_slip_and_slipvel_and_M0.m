% evaluate the local slip and slip velocity at the gouge-mediated event
% 2024.3.11 Kurama Okubo
% 2025.1.15 update computing the cumulative slip at the gouge event 

clear all; close all;
set(0,'DefaultTextFontsize',14, ...
    'DefaultTextFontname','Arial', ...
    'DefaultTextFontWeight','normal', ...
    'DefaultTextFontname','Arial', ...
    'DefaultAxesFontsize',14, ...
    'DefaultAxesFontname','Arial', ...
    'DefaultLineLineWidth', 1.0)

set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'})

figdir="../figure/p03_slip_slipvel_M0_master";
if ~exist(figdir) mkdir(figdir); end

addpath("../../../utils/matlabcode_biax_v03");

%% Load the slip velocity data computed at 01_compute_slipvelocity_allevents.m

expr_id = 87;
runID = sprintf('FB03-%03d', expr_id);

ifSavefig = true;

% Load the local slip velocity
E = readtable(sprintf("../data/localslip_and_slipvel_fb03-%03d.csv", expr_id), 'VariableNamingRule', 'preserve');
Nevent = size(E, 1);

%% set label of foreshock and aftershock
E.("eventtiming") = zeros(Nevent, 1);

for i = 1:Nevent
    if E(i, :).event_label == "A"
        E(i, :).("eventtiming") = 1;
    end
end


%%
ctype = 1;
c_foreshock = sns_colorpalette(ctype, 4);
c_aftershock = sns_colorpalette(ctype, 1);

%% Plot only the case with Nvalidsensors=4
Nvalidsensors_thresh = 4;

%% 1. Plot local slip vs M0
fig = figure(1); clf; hold on; box on;
fig.Units = 'point';
fig.Position = [0 800 600 500];

E_foreshock = E((E.eventtiming==0 & E.Nvalidsensors>=Nvalidsensors_thresh), :);
E_aftershock = E((E.eventtiming==1 & E.Nvalidsensors>=Nvalidsensors_thresh), :);

ms = 10;

% ylimit = [-1.2, 0.2];

% Plot only foreshock
% plot(E_foreshock.("cumulativelocalslip[mm]")*1e3, log10(E_foreshock.("M0")),...
% plot(log10(E_foreshock.("cumulativelocalslip[mm]")), log10(E_foreshock.("M0")),...


% plot(E_foreshock.("cumulativelocalslip[mm]")*1e3, E_foreshock.("M0"),...
%     "Marker", "v", "LineStyle","None", "MarkerSize", ms, "MarkerEdgeColor", "k", "MarkerFaceColor", c_foreshock);

plot(E_foreshock.("M0"), E_foreshock.("cumulativelocalslip[mm]")*1e3,...
    "Marker", "v", "LineStyle","None", "MarkerSize", ms, "MarkerEdgeColor", "k", "MarkerFaceColor", c_foreshock);

% annotate the event ID
% for i = 1:size(E_foreshock, 1)
%     text(E_foreshock.("cumulativelocalslip[mm]")(i)*1e3, E_foreshock.("M0")(i), sprintf("%d:%d", E_foreshock.("stickslip_id")(i), E_foreshock.("event_id")(i)));
% end

% plot(E_aftershock.("cumulativelocalslip[mm]"), log10(E_aftershock.("M0")),...
%     "Marker", "o", "LineStyle","None", "MarkerSize", ms, "MarkerEdgeColor", "k", "MarkerFaceColor", c_aftershock);

xscale('log');
yscale('log');

ylabel("Local cumulative slip [μm]");
xlabel("M_0 [Nm]");

% ylabel("log_{10}(M_0) [Nm]");
ylim([0.6, 10]);
xlim([1e-2, 5]);

% ylim(ylimit);
% yticks(-1.2:0.2:0.2)

figname = sprintf(figdir+"/slipandM0_master_fb03-%03d.png", expr_id);
if ifSavefig
    exportgraphics(fig, figname, 'Resolution',100);
end

%% Plot local slip velocity vs M0
fig = figure(2); clf; hold on; box on;
fig.Units = 'point';
fig.Position = [0 800 600 500];

E_foreshock = E(E.eventtiming==0, :);
E_aftershock = E(E.eventtiming==1, :);

ms = 10;

ylimit = [-1.2, 0.2];

plot(E_foreshock.("localslipvelocity[mm/s]"), log10(E_foreshock.("M0")),...
    "Marker", "v", "LineStyle","None", "MarkerSize", ms, "MarkerEdgeColor", "k", "MarkerFaceColor", c_foreshock);

plot(E_aftershock.("localslipvelocity[mm/s]"), log10(E_aftershock.("M0")),...
    "Marker", "o", "LineStyle","None", "MarkerSize", ms, "MarkerEdgeColor", "k", "MarkerFaceColor", c_aftershock);

xlabel("Local slip velocity [mm/s]");
ylabel("log_{10}(M_0) [Nm]");
xlim([-0.5, 4.0]);
ylim(ylimit);
yticks(-1.2:0.2:0.2)

% Annotate the 3 events for the case study
% %%
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

h1 = plot(E_foreshock.("localslipvelocity[mm/s]"), E_foreshock.("Mw"),...
    "Marker", "x", "LineStyle","None", "MarkerSize", ms_2, "MarkerEdgeColor", "y", "MarkerFaceColor", "k");

h2 = plot(E_aftershock.("localslipvelocity[mm/s]"), E_aftershock.("Mw"),...
    "Marker", "+", "LineStyle","None", "MarkerSize", ms_2, "MarkerEdgeColor", "y", "MarkerFaceColor", "k");

%%
M02Mw = @(x) (log10(x) - 9.105) * 2.0 / 3.0;

ax.XLim=[-0.5, 4.0];
ax.YLim=[M02Mw(10^ylimit(1)),M02Mw(10^ylimit(2))];
ax.YLabel.String="Mw";

h1.Visible = "off";
h2.Visible = "off";

%%
legend(["Foreshock", "Aftershock", "", ""] ,'Location', "southeast", "FontSize", 15);

% figname = sprintf(figdir+"/slipvelocity_and_M0_fb03-%03d.eps", expr_id);
% exportgraphics(fig, figname, 'Resolution',400);

figname = sprintf(figdir+"/slipvelocity_and_M0_fb03-%03d.png", expr_id);
exportgraphics(fig, figname, 'Resolution',100);




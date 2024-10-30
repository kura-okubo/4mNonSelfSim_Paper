% compute statistics of the experiments 
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


%% Load macro data
load("../../MacroData/data/MacroData_raw.mat");
% load("../../MacroData/data/SlipData_raw.mat");

expr_id = 87 %87;
ifSavefig = true;
%%
Mex = M.(sprintf("FB03_%03d", expr_id));
fs = 1/(Mex.tmat(2)-Mex.tmat(1));

%% convert from shear stress to shear load [N]
% The output of shear load cell is in [N]. When reading the data,
% we divide it by (4*0.1) to obtain the average shear stress on the fault.
% Here, we put back it from shear stress to shear loading as Mex.SF.

Mex.SF = Mex.SS * (4*0.1); %[MN]

%%
runID = sprintf('FB03_%03d', expr_id);
T = readtable(sprintf("../../DetectEvent/data/p02_eventtype_fb03-%03d.csv", expr_id));
Nevent = size(T,1);

eventtime = T.event_starttime;
%% 1. Macroscopic 

fig = figure(1); clf; hold on;
fig.Units = 'point';
fig.Position = [0 800 1000 600];

% compute x and y limits
xlimit = [0, round(eventtime(end)*0.1 + 1) * 10];
ylim_range = find((100<Mex.tmat) & (Mex.tmat<200)); 
ylim_shearload = [mean(Mex.SF(ylim_range))-0.02, mean(Mex.SF(ylim_range))+0.02];

xline_width = 1.5;

t = tiledlayout(3, 1);
t.TileSpacing = 'compact';
t.Padding = 'compact';

nexttile; box on; hold on;

plot([eventtime eventtime], ylim(gca()), "k:");

plot(Mex.tmat, Mex.SF, "k-");
xlim(xlimit);
ylim(ylim_shearload);
% ylim([0.245, 0.28]);
ylabel("Shear load [MN]");

%% 2. Drop of shear loading at event
% figure(2); clf; hold on;
% plot(Mex.tmat, Mex.SS, "k-");
% xlim(xlimit);

event_range = [0e-3, 150e-3]; % [s] window length from before to after the event
event_range_k = round(event_range * fs);
sheardrop = zeros(Nevent, 2);

for i=1:Nevent
% i = 12;
    st = T.event_starttime(i);
    ievent = find(Mex.tmat >=st, 1);
    SF_before = Mex.SF(ievent-event_range_k(1));
    SF_after = Mex.SF(ievent+event_range_k(2));

    SF_before_all(i, 1) = Mex.tmat(ievent-event_range_k(1));
    SF_before_all(i, 2) = SF_before;
    SF_before_all(i, 3) = ievent-event_range_k(1);

    SF_after_all(i, 1) = Mex.tmat(ievent+event_range_k(2));
    SF_after_all(i, 2) = SF_after;
    SF_after_all(i, 3) = ievent+event_range_k(2);

    sheardrop(i, 1) = st;
    sheardrop(i, 2) = SF_before - SF_after; %[MN]
end

plot(SF_before_all(:, 1), SF_before_all(:, 2), "o", "MarkerSize", 8, "MarkerEdgeColor", "k", "MarkerFaceColor","w");
plot(SF_after_all(:, 1), SF_after_all(:, 2), "v", "MarkerSize", 8, "MarkerEdgeColor", "k", "MarkerFaceColor","w");

%% Compute macroscopic shear loading rate
% Process flow
% 1. linear fiting between the event
% 2. compute the R^2 and threshold out with it
% 3. store the shear loading rate by its slope

rsq_thresh = 0.995;

loadrate_all = NaN * zeros(Nevent-1, 1);

for i=1:Nevent-1
% i = 16;
linfit_ind = find((SF_after_all(i, 1) < Mex.tmat) & (Mex.tmat <SF_before_all(i+1, 1)));  
p = polyfit(Mex.tmat(linfit_ind), Mex.SS(linfit_ind), 1);
yfit = p(1)*Mex.tmat(linfit_ind) + p(2);
yresid = Mex.SS(linfit_ind) - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(linfit_ind)-1) * var(Mex.SS(linfit_ind));
rsq = 1 - SSresid/SStotal;
 
% figure(2); clf; box on; hold on;
% plot(Mex.tmat(linfit_ind), Mex.SS(linfit_ind), "k-");
% plot(Mex.tmat(linfit_ind), yfit, "r-");
% title(sprintf("event %02d R^2 = %4.6f", i, rsq));

    if rsq>rsq_thresh
        loadrate_all(i) = p(1); %[MPa/s]

        % overplot the slope
        plot(Mex.tmat(linfit_ind), yfit*(4*0.1), "r-", "LineWidth", 1.5);
    end
end

Nloadrate_valid = nnz(~isnan(loadrate_all));
% disp(sprintf("Valid load rate span: %d", Nloadrate_valid));
title(sprintf("FB03-%03d Nevent:%d Valid load rate span: %d, Macro shear loading rate: %4.1f±%4.1f [kPa/s]",...
    expr_id, Nevent, Nloadrate_valid, mean(loadrate_all, 'omitnan')*1e3,  std(loadrate_all, 'omitnan')*1e3));


%% Compute displacement of west LDT at the event 
% The target of LDT after the main rupture is vibrated causing the
% oscillation in the measurement. Thus, we compute the temporal average of LDT.

nexttile; box on; hold on;

xlim(xlimit);
ylim([-0.1, 2.2]);
plot([eventtime eventtime], ylim(gca()), "k:");

h1 = plot(Mex.tmat, -Mex.DX(:, 1), "k-");
h2 = plot(Mex.tmat, -Mex.DX(:, 2), "r-");
%%
average_twin = 100e-3; %[s]
average_twin_k = round(average_twin * fs);

westdisp_all = zeros(Nevent, 1);
eastdisp_all = zeros(Nevent, 1);

for i=1:Nevent % Skip first event
% i = 24;
stind = SF_before_all(i, 3);
etind = SF_after_all(i, 3);
avg_st_DXw = mean(Mex.DX(stind-average_twin_k:stind, 1));
avg_et_DXw = mean(Mex.DX(etind:etind+average_twin_k, 1));

avg_st_DXe = mean(Mex.DX(stind-average_twin_k:stind, 2));
avg_et_DXe = mean(Mex.DX(etind:etind+average_twin_k, 2));

LDTdisp_west = -(avg_et_DXw - avg_st_DXw); %[mm]
LDTdisp_east = -(avg_et_DXe - avg_st_DXe); %[mm]

westdisp_all(i) = LDTdisp_west;
eastdisp_all(i) = LDTdisp_east;
% 
% plot(SS_before_all(i, 1), avg_st_DX, "ko")
% plot(SS_after_all(i, 1), avg_et_DX, "kv")
% text(SS_after_all(i, 1), avg_et_DX, sprintf("%4.1f", LDTdisp_west*1e3), "FontSize",12, "Color", "r");
end

slip_average = (westdisp_all + eastdisp_all)./2;
maxslip = -(avg_et_DXw + avg_et_DXe)/2;

% Drop the first event as the hold effect is large
slip_average_drop1st = slip_average(2:end);

title(sprintf("Average slip: %4.1f±%4.1f[μm], Total slip: %4.1f[mm]",...
    mean(slip_average_drop1st)*1e3, std(slip_average_drop1st)*1e3, maxslip));
ylabel("Displacement [mm]");

legend([h1 h2], ["West LDS", "East LDS"],'Location','northwest');

%% Compute stiffness
% we compute the stiffness by del sigma/del u at west

stiffness_all = zeros(Nevent, 1);

for i=1:Nevent % Skip first event
    stiffness_all(i) = sheardrop(i, 2)/(westdisp_all(i)*1e-3); %[MN/m]
end

%%
nexttile; box on; hold on;

title(sprintf("Stiffness %4.1f±%4.1f [MN/m]", mean(stiffness_all), std(stiffness_all)));
xlim(xlimit);
ylim([150, 400]);
plot(xlim(gca()), [mean(stiffness_all), mean(stiffness_all)], "b-");
plot([eventtime eventtime], ylim(gca()), "k:");
plot(SF_after_all(:, 1), stiffness_all, "ko", "MarkerFaceColor","w", "MarkerSize", 8);
ylabel("Stiffness [MN/m]");

%%
figname = sprintf("../figure/experiment_stats_FB03-%03d.png", expr_id);
exportgraphics(fig, figname, "Resolution", 80);
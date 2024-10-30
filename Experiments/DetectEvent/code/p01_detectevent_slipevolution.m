% 01 detect events from slip evolution
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

%%
movmean_winlen = 1e-3; %[s]
t_ave_transient = 2e-3; %[s] the time span to evaluate the offset before and after the time snap.

%% Load from mat file
% Run Experiments/MacroData/code/Slip_onfault.m to load the slip time
% history from raw data and save it in the mat format.

load("../../MacroData/data/SlipData_raw.mat");

%%
%Location of displacement sensors
x_disp = 270 + (0:15)*240;

Nsensor = length(x_disp);

Fs = 50e3; % sampling frequency

% for expr_id = [87, 94, 99] %[87, 94]
expr_id = 87;
runID = sprintf("FB03_%03d", expr_id);
%%
clear Dmat_trim_smooth tmat_trim;
tmat_trim = S.(sprintf("FB03_%03d", expr_id)).tmat;
Dmat_trim_smooth = movmean(S.(sprintf("FB03_%03d", expr_id)).Dmat, round(movmean_winlen*Fs));

%%
plot(tmat_trim(1:100:end), Dmat_trim_smooth(1:100:end, :), ".-");

%% Event detection
[Ndat,Nchan] = size(Dmat_trim_smooth);
t_ave_transient_k = round(t_ave_transient * Fs);

tmat_transient = zeros(Ndat-(t_ave_transient_k+1), 1);
Dmat_transient = zeros(Ndat-(t_ave_transient_k+1), Nchan);

for i = 1:Ndat-(t_ave_transient_k+1)
    tmat_transient(i) = tmat_trim(i);
    Dmat_transient(i, :) = Dmat_trim_smooth(i+t_ave_transient_k, :) - Dmat_trim_smooth(i, :); 
end

% %% plot
% step = 100;
% 
% figure(1);clf;
% [C,h] = contourf(x_disp, tmat_transient(1:step:end), Dmat_transient(1:step:end, :), 11);
% set(h,'LineColor','none')
% colorbar;
% caxis([0, 100e-6]);

%% plot surf at events
step=100;

[X,Y] = meshgrid(x_disp, tmat_transient(1:step:end));

figure(2);clf;
surf(X, Y, Dmat_transient(1:step:end, :));

% Plot time history of events
figure(3);clf;
plot(tmat_transient, max(abs(Dmat_transient), [], 2));

%% Search events
max_transient_offset=max(abs(Dmat_transient), [], 2); % get maximum transient slip over sensors
FW = 50; % window npts before event
BW = 500; % window npts after event

lc_long = "blue"; %[0/255, 75/255, 130/255];
lc_short = "black"; %[1, 0, 0];

plottime_index=[];
event_initid = [];
lc_index=[];

dt_long = 1.0;
dt_short = 1e-3;

dt_long_k = round(dt_long*Fs);
dt_short_k = round(dt_short*Fs);

event_slip_threshold = 0.002; %[mm]

it = 1; % current time stamp id
plottime_index(end+1) = it;
lc_index(end+1) = 1;

ievent = 0;
aftershock_flag = 0; % used to find successive events

assert ((FW + BW) > dt_short_k);

while it < Ndat-dt_long_k

    % trim the time series
    istart = it;
    iend = it + dt_long_k;
    tr = abs(max_transient_offset(istart:iend));

    if max(tr) < event_slip_threshold % no event
        plottime_index(end+1) = iend;
        lc_index(end+1) = 1;

        it = it + dt_long_k + 1;
        aftershock_flag = false;
        continue;

    else % with event
%         fprintf("%d\n", it)
        event_max = find(tr>event_slip_threshold,1) + it - 1; % find the first data point of event
        event_start = event_max - FW;
        event_end   = event_max + BW;
        %---plot with shorttime step---%
        if aftershock_flag
%             disp("debug")
            for k = it+dt_short_k:dt_short_k:event_end
                plottime_index(end+1) = k;
                lc_index(end+1) = 2;     
            end
        else
            event_initid(end+1) = event_start; % register new event timing
            for k = event_start:dt_short_k:event_end
%                 fprintf("k %d", k);
                plottime_index(end+1) = k;
                lc_index(end+1) = 2;
            end
        end
        it = plottime_index(end) + 1;
%         fprintf("%d, %d, %d\n", it, event_end, plottime_index(end));
        aftershock_flag = true;    
        %------------------------------%
    end

end

%% Plot slip accumulation
fig = figure(4);
fig.Units = 'point';
fig.Position = [0 500 1000 1200];
clf(fig,'reset'); cla(fig,'reset'); hold on; box on; 

for i = 1:length(plottime_index)
    pit = plottime_index(i);

    if i == 1
        plot(x_disp/1e3, Dmat_trim_smooth(pit, :), '-o', 'Color', lc_long, 'MarkerSize', 12);
        continue
    end

    if lc_index(i) == 2;
        plot(x_disp/1e3, Dmat_trim_smooth(pit, :), '-', 'Color', lc_short, 'LineWidth', 0.1);
    end

end

for i = 1:length(plottime_index)
    pit = plottime_index(i);
    if lc_index(i) == 1;
        plot(x_disp/1e3, Dmat_trim_smooth(pit, :), '-', 'Color', lc_long, 'LineWidth', 3);
    end
end

xlabel("Sensor location [m]");
ylabel("Displacement [mm]");

legend({sprintf("dt=%4.1f[s]", dt_long) , sprintf("dt=%4.1f[ms]", dt_short*1e3)});

title(sprintf("%s total %d events", runID, length(event_initid)));

% hcb = colorbar;
% ylabel(hcb, 'Time [s]');
% caxis([0, Tmax]);

xlim([0, 4]);
ylim([0, 2.5]);

% exportgraphics(gcf, sprintf("../figure/%s_displacement_evolution_adaptivetime_all.png", runID), "Resolution",300);

%% Save the event time to plot the AE waveform

tevent_start = tmat_trim(event_initid);
writematrix(tevent_start,sprintf('../data/%s_p01_eventstarttime.csv', runID)); 



%% Plot slip accumulation simple
fig = figure(10);
fig.Units = 'point';
fig.Position = [0 500 800 800];
clf(fig,'reset'); cla(fig,'reset'); hold on; box on; 

for i = 1:length(plottime_index)
    pit = plottime_index(i);

    if i == 1
        plot(x_disp/1e3, Dmat_trim_smooth(pit, :), '-o', 'Color', "k", 'MarkerSize', 12);
        continue
    end

    if lc_index(i) == 2;
%         plot(x_disp/1e3, Dmat_trim_smooth(pit, :), '-', 'Color', lc_short, 'LineWidth', 0.1);
    end

end

for i = 1:length(plottime_index)
    pit = plottime_index(i);
    if lc_index(i) == 1;
        plot(x_disp/1e3, Dmat_trim_smooth(pit, :), '-', 'Color', "k", 'LineWidth', 3);
    end
end

xlabel("Sensor location [m]");
ylabel("Displacement [mm]");

% legend({sprintf("dt=%4.1f[s]", dt_long) , sprintf("dt=%4.1f[ms]", dt_short*1e3)});

% title(sprintf("Plot every %.1fms", 1/Fs*span*1e3));

% hcb = colorbar;
% ylabel(hcb, 'Time [s]');
% caxis([0, Tmax]);

xlim([0, 4]);
ylim([0, 2.5]);

title(sprintf("fb03-%03d", expr_id));
% exportgraphics(gcf, sprintf("../figure/%s_displacement_evolution_adaptivetime_all_simple.png", runID), "Resolution",300);
exportgraphics(gcf, sprintf("../figure/%s_displacement_evolution_adaptivetime_all_simple.png", runID), "Resolution", 80);

% end

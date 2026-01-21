% Evaluate local normal stress on GP for revision
% 2025.12.07 Kurama Okubo
% This script is made with the support of ChatGPT.

clear all; close all;
set(0,'DefaultTextFontsize',14, ...
    'DefaultTextFontname','Arial', ...
    'DefaultTextFontWeight','normal', ...
    'DefaultTextFontname','Arial', ...
    'DefaultAxesFontsize',14, ...
    'DefaultAxesFontname','Arial', ...
    'DefaultLineLineWidth', 1.0)

set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'})

figdir = "../debug_figure/";
if ~exist(figdir, "dir"); mkdir(figdir); end
addpath("../../../../utils/matlabcode_biax_v03");

datadir = "../data/";
if ~exist(datadir, "dir"); mkdir(datadir); end

%% Load the master data
expr_id = 87;
runID   = sprintf('FB03_%03d', expr_id);

T = readtable(sprintf("../../../../Experiments/DetectEvent/data/p02_eventtype_fb03-%03d.csv", expr_id));

A = readtable("../../../../Experiments/DetectEvent/data/p06_visual_pick_gougeevents_merged.csv", ...
              'NumHeaderLines', 5);

%% GP (P3) configuration
gougepatch_x        = 1750;            % [mm] location of P3
SGT_search_range    = 100;             % [mm] search window for SGT
flatjack_x          = 250 + 500*(0:7); % [mm] flat-jack positions
flatjack_search_range = 100;           % [mm]

%% Time window for averaging
average_twinlen = 2e-3;   % 2 ms total
half_win        = average_twinlen / 2;

%% Directory containing event waveform data
event_datdir = sprintf("/Volumes/Okuboetal2025_masterHDD/4mBIAX_eventdata_master/p03_eventdata_FB03_%03d", expr_id);

%% Output table for GP events at P3
GP_normalstress = table();
% columns will be:
% expr_id, stickslip_id, gougeevent_id, event_time, event_loc_m, ...
% mean_SGT, mean_flatjack

%% Loop through all stick-slip events
for event_ind = 1:height(T)
% event_ind = 29;

    event_id = T.event_id(event_ind);
    fprintf("Processing stick-slip event %d (row %d)\n", event_id, event_ind);

    % Load event-specific .mat file
    matfile = fullfile(event_datdir, ...
        sprintf("eventdata_FB03_%03d_event%02d.mat", expr_id, event_id));

    if ~isfile(matfile)
        fprintf("  -> MAT file not found, skipping.\n");
        continue;
    end

    load(matfile);  
    % Assumes the .mat file contains:
    %   tmat_strain_event, Snmat, SGT_x,
    %   tmat_macro, NPmacro, Tstart, etc.

    % Find SGT and flat-jack indices near P3
    SGT_inds      = find(abs(SGT_x      - gougepatch_x) <= SGT_search_range);
    flatjack_inds = find(abs(flatjack_x - gougepatch_x) <= flatjack_search_range);

    if isempty(SGT_inds) || isempty(flatjack_inds)
        fprintf("  -> No sensors near P3, skipping.\n");
        continue;
    end

    % Extract GP events associated with this stick-slip event (Ordinary only)
    Aevents = A(A.Var2 == event_id & strcmp(string(A.Var5), "Ordinary"), :);

    if isempty(Aevents)
        fprintf("  -> No Ordinary GP event for this stick-slip.\n");
        continue;
    end

    % Containers for debug plotting (all P3 GPs within this stick-slip)
    gp_times           = [];
    mean_SGT_list      = [];
    mean_flatjack_list = [];
%%
    % Loop over GP events for this stick-slip
    for ievent = 1:height(Aevents)

        event_loc_m  = Aevents.Var3(ievent);          % [m]
        event_timing = Aevents.Var4(ievent) - Tstart; % [s]
        gouge_id     = Aevents.Var7(ievent);          % GP event id

        % Use only P3 (x = 1.75 m → 1750 mm)
        if abs(event_loc_m*1e3 - gougepatch_x) > 1e-6
            continue;
        end

        % Time window around the GP event
        tmin = event_timing - half_win;
        tmax = event_timing + half_win;

        idx_strain = find(tmat_strain_event >= tmin & tmat_strain_event <= tmax);
        idx_macro  = find(tmat_macro        >= tmin & tmat_macro        <= tmax);

        if isempty(idx_strain) || isempty(idx_macro)
            fprintf("  -> No sufficient data window for GP %d (time %.6f s).\n", gouge_id, event_timing);
            continue;
        end

        % Compute average normal stress in the window
        mean_SGT_this      = mean(Snmat(idx_strain, SGT_inds), 'all');
        mean_flatjack_this = mean(NPmacro(idx_macro, flatjack_inds) * 2/3, 'all');

        % Append one row to the output table for this GP event
        newRow = table( ...
            expr_id, ...
            event_id, ...
            gouge_id, ...
            event_timing, ...
            event_loc_m, ...
            mean_SGT_this, ...
            mean_flatjack_this, ...
            'VariableNames', { ...
                'expr_id', ...
                'stickslip_id', ...
                'gougeevent_id', ...
                'event_time', ...
                'event_loc_m', ...
                'mean_SGT', ...
                'mean_flatjack'});

        GP_normalstress = [GP_normalstress; newRow];

        % Store for debug plot
        gp_times           = [gp_times; event_timing];
        mean_SGT_list      = [mean_SGT_list; mean_SGT_this];
        mean_flatjack_list = [mean_flatjack_list; mean_flatjack_this];

    end % ievent loop

    % If no valid P3 GP for this stick-slip, skip plotting
    if isempty(gp_times)
        continue;
    end

    %% Debug plot for this stick-slip event
    %% Debug plot for this stick-slip event
    fig = figure(1); clf; 
    fig.Units    = 'point';
    fig.Position = [0 800 800 800];
    
    %% -------------------------
    %  Subplot (1): SGT & flat-jack & mean markers
    %% -------------------------
    subplot(2,1,1); hold on; box on;
    
    % --- Plot SGT traces near P3 ---
    h_SGT = gobjects(length(SGT_inds), 1);
    for k = 1:length(SGT_inds)
        h_SGT(k) = plot(tmat_strain_event, Snmat(:, SGT_inds(k)), '-', 'LineWidth', 1.0);
    end
    set(h_SGT(2:end), 'HandleVisibility', 'off');
    
    % --- Plot flat-jack traces near P3 ---
    h_FJ = gobjects(length(flatjack_inds), 1);
    for k = 1:length(flatjack_inds)
        h_FJ(k) = plot(tmat_macro, NPmacro(:, flatjack_inds(k))*2/3, '--', 'LineWidth', 1.0);
    end
    set(h_FJ(2:end), 'HandleVisibility', 'off');
    
    % Time window markers
    ylim([0.8, 2.0]);
    yL = ylim;
    for k = 1:length(gp_times)
        t_ev = gp_times(k);
        tmin = t_ev - half_win;
        tmax = t_ev + half_win;
    
        plot([tmin tmin], yL, ':k', 'LineWidth', 0.8, 'HandleVisibility','off');
        plot([tmax tmax], yL, ':k', 'LineWidth', 0.8, 'HandleVisibility','off');
    end
    
    % Mean markers
    h_meanSGT = plot(gp_times, mean_SGT_list, 'o', ...
        'MarkerSize', 7, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'k');
    
    h_meanFJ = plot(gp_times, mean_flatjack_list, 's', ...
        'MarkerSize', 7, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k');
    
    xlabel('Time [s]');
    ylabel('Normal stress [MPa]');
    title(sprintf('Event %d (P3 GP) – Local normal stress', event_id));
    
    legend([h_SGT(1), h_FJ(1), h_meanSGT, h_meanFJ], ...
        {'SGT near P3', 'Flat-jack near P3', ...
         'Mean SGT (P3 GPs)', 'Mean flat-jack (P3 GPs)'}, ...
        'Location', 'best', 'Box', 'off');
    
    
    %% -------------------------
    %  Subplot (2): Macroscopic shear stress history
    %% -------------------------
    subplot(2,1,2); hold on; box on;
    
    plot(tmat_macro, SSmacro, 'k-', 'LineWidth', 1.1);
    
    xlabel('Time [s]');
    ylabel('Macroscopic shear stress [MPa]');
    title('Macroscopic shear stress history during stick-slip');
    
    ylim([0.62 0.70]);   % Requested y-range

    % Save debug figure
    saveas(fig, fullfile(figdir, sprintf('fb03-%03d_normalstress_P3_event%02d.png', expr_id, event_id)));

end % event_ind loop

%% Save merged GP-event table
outfile = sprintf("%s/gp_normalstress_fb03-%03d_P3.csv", datadir, expr_id);
writetable(GP_normalstress, outfile);
fprintf("GP-normal-stress table written to: %s\n", outfile);

%%
SGT_inds
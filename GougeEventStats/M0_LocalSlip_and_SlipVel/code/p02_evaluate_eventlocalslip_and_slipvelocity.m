% evaluate the local slip velocity at the gouge event
% 2024.3.8 Kurama Okubo
% 2025.1.15 update computing the cumulative slip at the gouge event 
% 2025.1.28 update computing the time history of slip and slip velocity on
% the gouge

clear all; close all;
set(0,'DefaultTextFontsize',14, ...
    'DefaultTextFontname','Arial', ...
    'DefaultTextFontWeight','normal', ...
    'DefaultTextFontname','Arial', ...
    'DefaultAxesFontsize',14, ...
    'DefaultAxesFontname','Arial', ...
    'DefaultLineLineWidth', 1.0)

set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'})

datadir="../data/p02_slipinterpolation";
if ~exist(datadir) mkdir(datadir); end

figdir="../figure/p02_gougeevent_slipvel";
if ~exist(figdir) mkdir(figdir); end

addpath("../../../utils/matlabcode_biax_v03");

%% Load the slip velocity data computed at 01_compute_slipvelocity_allevents.m

expr_id = 87;
runID = sprintf('FB03-%03d', expr_id);

ifSavefig = true;


% The data below is for the gouge events
patch_id="G3";
Q_quart = "Q50";
fi_gougeevents = sprintf("../../../GougeEventCatalog/data/gougeeventcatalog__fb03-%03d__%s__%s.csv", expr_id, patch_id, Q_quart);
E = readtable(fi_gougeevents);
Nevent = size(E, 1);

% We use the origin time from the relocation of the gouge event to evaluate
% the local slip velocity.
% E = readtable("../../SourceInversion/data/datacsv/AE_obs_location.csv");
% Nevent = size(E, 1);

% add the column of slip velocity
E.("cumulativelocalslip[mm]") = zeros(Nevent, 1);
E.("localslipvelocity[mm/s]") = zeros(Nevent, 1);

stickslip_id_list = unique(E.stickslip_id);

for kk = 1:length(stickslip_id_list)
    slipevent_id = stickslip_id_list(kk);
    % slipevent_id = 29;

    fprintf("Start processing stick-slip event: %02d\n", slipevent_id);

    fname = sprintf("%s_slipvelocity_event%03d.mat", runID, slipevent_id);

    V = load("../data/p01_slipandslipvel/"+fname);
    E_slipevent = E(E.stickslip_id==slipevent_id, :);
    
    %%

    NT = length(V.tmat_slip_event);

    fig = figure(1); clf; hold on; box on;
    fig.Units = 'point';
    fig.Position = [0 800 800 600];

    xlimit=[0, 200]; %[ms]

    Tlength = xlimit(2)*1e-3;
    cd = colormap(flipud(viridis)); % line color by time
    lc_intp_event = @(x) interp1(linspace(0, Tlength+1e-8, length(cd)),cd, x);

    Lw = 1.5;

    plotspan = 10;

    for i = 1:plotspan:NT
        linecolor=lc_intp_event(V.tmat_slip_event(i)); %0.5;
        plot(V.Disp_x, V.SlipVelmat_filtered(i, :), "-", 'Color', linecolor, "LineWidth", Lw);
        % pause(0.1);
    end


    %% Evaluate the slip and slip velocity at the gouge event location and timing

    % We linear interpolate in both space and time of the slip and slip velocity
    % to evaluate the local slip velocity at the origin time and at the gouge
    % patch location
    % 2025.1.25 update to compute the slip

    [X, Y] = meshgrid(V.Disp_x*1e-3, V.Tstart + V.tmat_slip_event);

    % loop for each gouge event
    for ll = 1:size(E_slipevent, 1)

        gougeevent_id =  E_slipevent(ll, :).event_id;
        origin_time =  E_slipevent(ll, :).event_onset_time;

        % assert if the origin time is close to the visually picked time
        origin_time_vis = V.Aevents(V.Aevents.Var7==gougeevent_id, :).Var4;
        assert(abs(origin_time-origin_time_vis) < 1e-3);

        origin_ind = find(V.Tstart + V.tmat_slip_event >= origin_time, 1);
        fprintf("slip and slipvel at t=%4.4f, event origin time=%4.4f.\n", V.Tstart + V.tmat_slip_event(origin_ind), origin_time);

        % interpolation
        gougepatch_loc_x =  E_slipevent(ll, :).location;
        slip_interpolated = interp2(X, Y, V.Dmat_event_filtered, gougepatch_loc_x, origin_time,'linear'); % [mm]
        slipvel_interpolated = interp2(X, Y, V.SlipVelmat_filtered, gougepatch_loc_x, origin_time,'linear'); % [mm/s]

        % Update: compute time history of slip and slip velocity
        tr_tmat = V.Tstart + V.tmat_slip_event;
        Ndata = length(tr_tmat);
        tr_slip_gouge = zeros(Ndata, 1);
        tr_slipvel_gouge = zeros(Ndata, 1);
        
        
        for i = 1:Ndata
            tr_slip_gouge(i) = interp1(V.Disp_x*1e-3, V.Dmat_event_filtered(i, :), gougepatch_loc_x, 'linear');
            tr_slipvel_gouge(i) = interp1(V.Disp_x*1e-3, V.SlipVelmat_filtered(i, :), gougepatch_loc_x, 'linear');
            % debug
            % if mod(i, 25)==0
            %     figure(9); clf; hold on;
            %     plot(V.Disp_x*1e-3, V.Dmat_event_filtered(i, :), "k-");
            %     plot(gougepatch_loc_x, tr_slip_gouge(i), "ro");
            %     ylim([0, 40e-3]);
            %     ylabel("Slip [mm]")
            %     yyaxis right;
            %     plot(V.Disp_x*1e-3, V.SlipVelmat_filtered(i, :), "r-");
            %     plot(gougepatch_loc_x, tr_slipvel_gouge(i), "bv");
            %     ylim([0, 6]);
            %     ylabel("Slip velocity [mm/s]")
            % 
            %     title(sprintf("T=%12.6f",tr_tmat(i)));
            %     pause(0.01);
            % end
        end

        % plot figure to check if the interpolation is reasonable
        fig = figure(2); clf; hold on; box on;
        fig.Units = 'point';
        fig.Position = [0 800 800 500];

        subplot(2, 1, 1); hold on; box on;
        Nline = 10; % plot before and after the gouge event
        lc = turbo(Nline);
        for i = 1:Nline
            plot(V.Disp_x/1e3, V.SlipVelmat_filtered(origin_ind+i-round(Nline/2), :), "o-", "Color", lc(i, :));
        end

        % Plot interpolated velocity
        plot(gougepatch_loc_x, slipvel_interpolated, "k", "MarkerSize", 20, "MarkerFaceColor", "r", "Marker","pentagram");

        title(sprintf("Event%02d slipvel time %4.4f[s] origin %4.4f[s] event local slip vel %4.4f[mm/s]", ...
            gougeevent_id, V.Tstart + V.tmat_slip_event(origin_ind), origin_time, slipvel_interpolated));

        xlim([0, 4.1]);
        ylim([-0.5, 4]);
        xlabel("Distance from western edge of fault [m]");
        ylabel("Slip velocity [mm/s]");

        %----------------------%
        % Plot interpolated slip
        %----------------------%
        subplot(2, 1, 2); hold on; box on;
        for i = 1:Nline
            plot(V.Disp_x/1e3, V.Dmat_event_filtered(origin_ind+i-round(Nline/2), :)*1e3, "o-", "Color", lc(i, :));
        end

        plot(gougepatch_loc_x, slip_interpolated*1e3, "k", "MarkerSize", 20, "MarkerFaceColor", "r", "Marker","pentagram");

        title(sprintf("Event%02d slipvel time %4.6f[s] origin %4.6f[s]\nevent local slip vel %4.2f[mm/s] M0 %4.2f[Nm]", ...
            gougeevent_id, V.Tstart + V.tmat_slip_event(origin_ind), origin_time, slipvel_interpolated, E_slipevent(ll, :).M0));

        xlim([0, 4.1]);
        ylim([-2, 15]);
        xlabel("Distance from western edge of fault [m]");
        ylabel("Cumulative slip [μm]");

        figname = sprintf(figdir+"/gougeevent_slipvel_fb03-%03d_event%02d.png", expr_id, gougeevent_id);
        if ifSavefig
            exportgraphics(fig, figname, 'Resolution',80);
        end

        % add the slip velocity to the table
        E(E.event_id == gougeevent_id, :).("cumulativelocalslip[mm]") = slip_interpolated;
        E(E.event_id == gougeevent_id, :).("localslipvelocity[mm/s]") = slipvel_interpolated;

        %% Save data for the comparison to the seismic moment
        foname=sprintf(datadir+"/FB03-%03d_slip_and_slipvel_atgouge_gougeevent%03d.mat", expr_id, gougeevent_id);
        save(foname, "tr_tmat", "gougepatch_loc_x", "tr_slip_gouge", "tr_slipvel_gouge");

    end

end

writetable(E, sprintf('../data/localslip_and_slipvel_fb03-%03d.csv', expr_id),'Delimiter',',');
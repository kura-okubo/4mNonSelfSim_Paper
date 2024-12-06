% evaluate the local slip velocity at the gouge-mediated event
% 2024.3.8 Kurama Okubo

clear all; close all;
set(0,'DefaultTextFontsize',14, ...
    'DefaultTextFontname','Arial', ...
    'DefaultTextFontWeight','normal', ...
    'DefaultTextFontname','Arial', ...
    'DefaultAxesFontsize',14, ...
    'DefaultAxesFontname','Arial', ...
    'DefaultLineLineWidth', 1.0)

figdir="../figure/debug_gougeevent_slipvel";
if ~exist(figdir) mkdir(figdir); end

addpath("../../../utils/matlabcode_biax_v03");

%% Load the slip velocity data computed at 01_compute_slipvelocity_allevents.m

expr_id = 87;
runID = sprintf('FB03-%03d', expr_id);

ifSavefig = true;

% We use the origin time from the relocation of the gouge event to evaluate
% the local slip velocity.
E = readtable("../../../SourceInvFit/data/datacsv/AE_obs_location.csv");
Nevent = size(E, 1);

% add the column of slip velocity
E.("slipvelocity[mm/s]") = zeros(Nevent, 1);

for kk = 1:length(E.stickslip_id)
    % slipevent_id = 35;
    slipevent_id = E.stickslip_id(kk);
    fprintf("Start processing stick-slip event: %02d\n", slipevent_id);

    fname = sprintf("%s_slipvelocity_event%02d.mat", runID, slipevent_id);

    V = load("../data/"+fname);
    E_slipevent = E(E.stickslip_id==slipevent_id, :);
    
    %%

    NT = length(V.tmat_slip_event);

    % fig = figure(1); clf; hold on; box on;
    % fig.Units = 'point';
    % fig.Position = [0 800 800 600];
    % 
    % xlimit=[0, 200]; %[ms]
    % 
    % Tlength = xlimit(2)*1e-3;
    % cd = colormap(flipud(viridis)); % line color by time
    % lc_intp_event = @(x) interp1(linspace(0, Tlength+1e-8, length(cd)),cd, x);
    % 
    % Lw = 1.5;
    % 
    % plotspan = 100;
    % 
    % for i = 1:plotspan:NT
    %     linecolor=lc_intp_event(V.tmat_slip_event(i)); %0.5;
    %     plot(V.Disp_x, V.SlipVelmat_filtered(i, :), "-", 'Color', linecolor, "LineWidth", Lw);
    %     % pause(0.1);
    % end


    %% Evaluate the slip velocity at the gouge event timing

    % We linear interpolate in both space and time of the slip velocity
    % to evaluate the local slip velocity at the origin time and at the gouge
    % patch location


    [X, Y] = meshgrid(V.Disp_x*1e-3, V.Tstart + V.tmat_slip_event);

    for i = 1:size(E_slipevent, 1)

        gougeevent_id =  E_slipevent(i, :).gougeevent_id;
        origin_time =  E_slipevent(i, :).origin_time;

        % assert if the origin time is close to the visually picked time
        gougepatch_loc_x =  V.Aevents(V.Aevents.Var7==gougeevent_id, :).Var3;
        origin_time_vis = V.Aevents(V.Aevents.Var7==gougeevent_id, :).Var4;
        assert(abs(origin_time-origin_time_vis) < 1e-3);

        origin_ind = find(V.Tstart + V.tmat_slip_event >= origin_time, 1);
        fprintf("slipvel t=%4.4f, event origin time=%4.4f.\n", V.Tstart + V.tmat_slip_event(origin_ind), origin_time);

        slipvel_interpolated = interp2(X, Y, V.SlipVelmat_filtered, gougepatch_loc_x, origin_time,'linear'); % [mm/s]

        % plot figure to check if the interpolation is reasonable

        fig = figure(2); clf; hold on; box on;
        fig.Units = 'point';
        fig.Position = [0 800 800 500];

        Nline = 10;
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


        figname = sprintf(figdir+"/gougeevent_slipvel_fb03-%03d_event%02d.png", expr_id, gougeevent_id);
        if ifSavefig
            exportgraphics(fig, figname, 'Resolution',80);
        end

        % add the slip velocity to the table
        E(E.gougeevent_id == gougeevent_id, :).("slipvelocity[mm/s]") = slipvel_interpolated;
    end

end

writetable(E, sprintf('../data/localslipvel_fb03-%03d.csv', expr_id),'Delimiter',',');
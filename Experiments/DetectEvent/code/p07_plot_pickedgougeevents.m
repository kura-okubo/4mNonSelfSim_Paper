% 06 plot picked foreshocks from the csv datasheet
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

%% Load the compiled data

ifSavefig = true;
ifPlotAE = true; % true if plotting AE

% ifBandpass = 1; % Apply band-pass to find ordinary events.
% ifLFEpass = 0; % Set low-frequency range to detect LFEs.
% event_id = 1;
for event_id = 1:56

    expr_id = 87;
    runID = sprintf('FB03_%03d', expr_id);
    T = readtable(sprintf("../data/p02_eventtype_fb03-%03d.csv", expr_id));

    event_type = zeros(length(size(T, 1)), 3); %1. event type %2. rupture vel %3. event start time from strain

    % Read the picked csv file to superimposed the detected events.
    % A = readtable("../data/p06_visual_pick_gougeevents.csv", 'NumHeaderLines',5);
    A = readtable("../data/p06_visual_pick_gougeevents_merged.csv", 'NumHeaderLines',5);

    % filter by expr_id
    idx_expr = string(A.Var1) == sprintf("fb03-%03d", expr_id);
    A_expr = A(idx_expr,:);


    %%
    % for event_id = 1:size(T, 1)

    % filter by the event
    idx_ev = A.Var2 == event_id;
    if ~isempty(A_expr)
        A_event = A_expr(idx_ev,:);
    else
        A_event = [];
    end

    foreshock_pt = 0; %[ms] visually pick from the plot

    event_datdir=sprintf("/Volumes/4mGouge_WorkHDD/FB03data/4mBIAX_paper_tmp/p03_eventdata_FB03_%03d", expr_id);

    load(event_datdir + sprintf("/eventdata_FB03_%03d_event%02d.mat", expr_id, event_id));

    fig = figure(1); clf; hold on;
    fig.Units = 'point';
    fig.Position = [0 800 1600 1200];

    %--------------------------------------------%
    %% Plot AE and strain with rupture velocity
    %--------------------------------------------%

    % load event type
    assert(Tstart == T.event_starttime(event_id));

    % set id list to plot longer xlimit
    if expr_id == 87
        xlimit_long_list = [11, 23, 45];
    elseif expr_id == 94
        xlimit_long_list = [13];
    end

    if T.eventtype(event_id) == 2 || ismember(event_id, xlimit_long_list)% || T.eventtype(event_id) == 3
        AE_xlimit = [0, 200];
    else
        AE_xlimit = [0, 80];
    end

    %--------------%
    %% --Plot AE---%
    %--------------%

    bandpass_order = 4;
    fs_read = 1e7;

    % if ifLFEpass
    fmin_LFE = 2e4;
    fmax_LFE = 6e4;
    % else
    fmin_ORD = 1e5;
    fmax_ORD = 1e6; %6e5;
    % end

    ampnorm_AE = 5; %5 %20; %5e-5;

    % apply bandpass filter
    % if ifBandpass
    AEdatmat_raw = AEdatmat;

    f_nyq = fs_read/2;
    [b_ORD,a_ORD] = butter(bandpass_order, [fmin_ORD fmax_ORD]/f_nyq, "bandpass");
    AEdatmat_ORD = filtfilt(b_ORD,a_ORD,AEdatmat);

    [b_LFE,a_LFE] = butter(bandpass_order, [fmin_LFE fmax_LFE]/f_nyq, "bandpass");
    AEdatmat_LFE = filtfilt(b_LFE,a_LFE,AEdatmat);


    % clf(fig,'reset'); cla(fig,'reset'); hold on;
    % box on; grid on;
    Nsensors_AE = size(AEdatmat, 2);
    [~, AE_sortedinds] = sort(AEsensor_x);

    lc_AE = "k";

    % Plot raw data
    AEdata_plot = {AEdatmat_raw AEdatmat_ORD, AEdatmat_LFE};

    for l = 1:3
        subplot(3, 1, l); hold on; box on;

        dat = cell2mat(AEdata_plot(l));

        for i = 1:Nsensors_AE
            ii = AE_sortedinds(i);
            y_shift = AEsensor_x(ii)/1e3;
            p1 = plot(tmat_AE_event*1e3, dat(:, ii)/ampnorm_AE + y_shift, "-", "Color", lc_AE);
        end

        p1.DataTipTemplate.DataTipRows(1).Format = '%f'; % x

        % Plot scale of volt
        scale_x = AE_xlimit(2) * 0.92; %22.5;
        scale_y = 0.2;
        scale_len = 2; %[V]
        plot([scale_x, scale_x], [scale_y-(scale_len/ampnorm_AE)/2, scale_y+(scale_len/ampnorm_AE)/2], "k-", "LineWidth", 3);
        text(scale_x+1, scale_y, "2V");

        if l==1
            titlestr = sprintf("FB03-%03d event %02d Tstart:%.4f : AE with picked events: rupvel=%.2f[m/s] ruptype=%d",...
                expr_id, event_id, Tstart, T.nuc_rupturevel(event_id), T.eventtype(event_id));
        elseif l==2
            titlestr = sprintf("Ordinary event: bandpass %.0f-%.0fkHz", fmin_ORD/1e3, fmax_ORD/1e3);
        else
            titlestr = sprintf("LFE: bandpass %.0f-%.0fkHz", fmin_LFE/1e3, fmax_LFE/1e3);
        end

        title(titlestr);

        %---------------%
        %% --Plot strain--%
        %---------------%

        % NOTE: we recompute the polyfit to superimpose the linear regression of rupture
        % velocity to the plot.
        % Then, we renew the rupture velocity of nucleation.

        %%
        NSGB = size(Taumat2, 2);
        NSGT = size(Taumat3, 2);

        %% Plot along fault
        offsetind = 100;
        Taumat2_removeoffset = Taumat2 - mean(Taumat2(1:offsetind, :));
        Taumat3_removeoffset = Taumat3 - mean(Taumat3(1:offsetind, :));

        %% apply moving window average
        smooth_winlen = 100;
        Taumat2_removeoffset_smoothed = movmean(Taumat2_removeoffset, smooth_winlen, 1);
        Taumat3_removeoffset_smoothed = movmean(Taumat3_removeoffset, smooth_winlen, 1);


        %%
        ampnorm = 0.4; %5e-5;

        lc = [0.6, 0, 0, 0.5];

        for ii = 1:NSGB
            y_shift = SGB_x(ii)/1e3;
            %     plot(tmat_strain_event*1e3, Taumat2_removeoffset_smoothed(:, ii)/ampnorm + y_shift, "-", "Color", lc, "LineWidth", 1.2);
        end

        % plot SGT
        for ii = 1:NSGT
            y_shift = SGT_x(ii)/1e3;
            %     plot(tmat_strain_event*1e3, Taumat3_removeoffset_smoothed(:, ii)/ampnorm + y_shift, "-", "Color", lc, "LineWidth", 1.2);
        end

        % xline((AE_eventtime - Tstart) * 1e3, "--", "Color", "blue")
        % plot((AE_eventtime - Tstart) * 1e3, eventloc, "ro", "Markersize", 10)

        xlabel("Time [ms]");
        ylabel("x along fault [m]");
        % xlim([0, Tlength*1e3]);
        xlim(AE_xlimit);
        ylim([-0.1, 4.1]);

        % % Plot scale of shear stress
        % scale_x = AE_xlimit(2) * 0.85; %23.5;
        % scale_y = 0.2;
        % scale_len = 0.1; %[MPa]
        % plot([scale_x, scale_x], [scale_y-(scale_len/ampnorm)/2, scale_y+(scale_len/ampnorm)/2], "r-", "LineWidth", 3);
        % text(scale_x+0.2, scale_y, "0.1MPa");


        % Plot the x time

        % xline(foreshock_pt, "b--", "LineWidth", 1.5)
        %
        % sgtitlestr = sprintf("FB03-%03d event %02d: T%4.4f-%4.4f[s] Pick: %4.6f[s]", ...
        %     expr_id, event_id, Tstart, Tstart+AE_xlimit(2)/1e3, Tstart + foreshock_pt/1e3);
        %
        % % hSG = sgtitle(sgtitlestr, "FontWeight", "bold", "FontSize", 22);
        % text(0, 4.35, sgtitlestr, "FontWeight", "bold", "FontSize", 21);
        % pan on;

        % Plot picked event from csv file
        %  size(A_event)
        for i=1:size(A_event, 1)
            patchloc = A_event.Var3(i);
            tpick = A_event.Var4(i) - Tstart;
            AEtype = string(A_event.Var5(i));

            if AEtype == "Ordinary"
                mc = "r";
                mt = "o";
                lc = "r";
            elseif AEtype == "LFE"
                mc = "b";
                mt = "v";
                lc = "b";
            end

            xline(tpick*1e3, "--", "LineWidth", 1.0, "Color", lc);
            plot(tpick*1e3, patchloc, "Marker", mt, "MarkerEdgeColor",mc, "MarkerSize", 7);

        end

    end

    figdir = sprintf("../figure/p07_search_gougeevents_FB03_%03d_pickedevents", expr_id);
    if ~exist(figdir) mkdir(figdir); end


    figname = sprintf(figdir+"/p07_picked_gougeevents_FB03_%03d_event%02d.png", expr_id, event_id);

    if ifSavefig

        exportgraphics(fig, figname, 'Resolution',80); %150);
    end

    % To pick the event, get the data from figure and type the command below
    % Tstart + cursor_info.Position(1)/1e3

end

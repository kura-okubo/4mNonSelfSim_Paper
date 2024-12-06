% compute slip velocity profiles for all the event
% to investigate the relashionship with the seismic moment
% 2024.2.4 Kurama Okubo

clear all; close all;
set(0,'DefaultTextFontsize',14, ...
    'DefaultTextFontname','Arial', ...
    'DefaultTextFontWeight','normal', ...
    'DefaultTextFontname','Arial', ...
    'DefaultAxesFontsize',14, ...
    'DefaultAxesFontname','Arial', ...
    'DefaultLineLineWidth', 1.0)

figdir="../figure/plot_slipvel_master";
if ~exist(figdir) mkdir(figdir); end
figdir_debug="../figure/debug_slipvelhistory";
if ~exist(figdir_debug) mkdir(figdir_debug); end
addpath("../../../utils/matlabcode_biax_v03");

%% Load the master data

ifSavefig = true;

ifPlotContour = true; % this option is used to output figure in vector format.

expr_id = 87;
runID = sprintf('FB03_%03d', expr_id);
T = readtable(sprintf("../../../Experiments/DetectEvent/data/p02_eventtype_fb03-%03d.csv", expr_id));

event_type = zeros(length(size(T, 1)), 3); %1. event type %2. rupture vel %3. event start time from strain

%%
for event_ind = 1:length(T.event_id)
    % event_ind = 29;
    event_id = T.event_id(event_ind);
    fprintf("Start processing event %d\n", event_id);

    % event_datdir="../../../Experiments/DetectEvent/data/eventdata_master";
    event_datdir=sprintf("/Volumes/4mGouge_WorkHDD/FB03data/4mBIAX_paper_tmp/p03_eventdata_FB03_%03d", expr_id);

    load(event_datdir + sprintf("/eventdata_FB03_%03d_event%02d.mat", expr_id, event_id));

    %% Read event catalog
    A = readtable("../../../Experiments/DetectEvent/data/p06_visual_pick_gougeevents_merged.csv", 'NumHeaderLines',5);


    %% compute slip velocity using gradient
    dt_slip = 1/Fs_slip;

    % Apply low-pass filter to the slip history before compute the gradient
    % We apply the noraml moving-window average to avoid the artifact of
    % bandpass filter.

    lowpass_winlen = 2.0e-3; %[s]

    fprintf("Strain moving window smoothing with %4.4f[ms] low-pass: %4.4f[kHz]\n", lowpass_winlen*1e3, 1/lowpass_winlen/1e3);
    %%

    Dmat_event_offsetremoved = Dmat_event - mean(Dmat_event(1:5, :), 1);
    Dmat_event_filtered = zeros(size(Dmat_event_offsetremoved));


    % f_nyq = Fs_slip/2;
    % bandpass_order = 4;
    % [b,a] = butter(bandpass_order, fmax/f_nyq, "low");
    % freqz(b, a, [], Fs_slip);
    % ax = findall(gcf, 'Type', 'axes');
    % set(ax, 'XScale', 'log');
    % Dmat_event_filtered = filter(b,a,Dmat_event);

    for i=1:size(Dmat_event, 2)
        % Dmat_event_filtered(:, i) = movmean_causal(Dmat_event_offsetremoved(:, i), round(lowpass_winlen*Fs_slip));
        Dmat_event_filtered(:, i) = movmean(Dmat_event_offsetremoved(:, i), round(lowpass_winlen*Fs_slip));
    end

    fig = figure(1); clf; hold on; box on;
    fig.Units = 'point';
    fig.Position = [0 800 1200 1200];

    index = reshape(1:size(Dmat_event, 2), 2, 8).';

    for i = 1:size(Dmat_event, 2)

        subplot(size(Dmat_event, 2)/2, 2, index(i)); hold on; box on;

        % Compute slip velocity using MATLAB function gradient (central difference)
        % remove offset
        SlipVelmat = (gradient(Dmat_event_offsetremoved', dt_slip))'; % [mm/s]
        SlipVelmat_filtered = (gradient(Dmat_event_filtered', dt_slip))';

        yyaxis left
        plot(tmat_slip_event*1e3, Dmat_event_offsetremoved(:, i)*1e3, "b", "DisplayName","slip raw");
        plot(tmat_slip_event*1e3, Dmat_event_filtered(:, i)*1e3, "k-", "DisplayName","slip smoothed");
        % plot(tmat_slip_event, debug_singletrace, "b--");
        %
        if mod(i, 8)==0
            ylabel("Slip [μm]");
        end

        yyaxis right
        plot(tmat_slip_event*1e3, SlipVelmat(:, i), "-", "Color",[1, 0, 0, 0.3], "DisplayName","slipvel raw");
        plot(tmat_slip_event*1e3, SlipVelmat_filtered(:, i), "k-", "DisplayName","slipvel with smoothed slip");
        % plot(tmat_slip_event, SlipVelmat_filtered_2(:, disp_plotid), "b-");

        % ylabel("Slip velocity [mm/s]");
        % xlabel("Time [ms]");
        % legend("Location","southeast");
        title("gap sensor "+i);

        if mod(i, 8)==0
            ylabel("Slip velocity [mm/s]");
            xlabel("Time [ms]");
        end
    end

    sgtitle(sprintf("FB03\\_%03d stick-slip event %02d", expr_id, event_id));

    figname = sprintf(figdir_debug+"/slipandslipvel_fb03-%03d_event%02d.png", expr_id, event_id);
    if ifSavefig
        exportgraphics(fig, figname, 'Resolution',80);
    end

    %% Plot contour of slip velocity
    fig = figure(2); clf; hold on; box on;
    fig.Units = 'point';
    fig.Position = [0 800 800 800];

    t_plotshift = 0;%0.022; %0.025; % shift time to plot from t=0

    %
    map = viridis;
    % map(1, :) = [1, 1, 1];
    %
    colormap(map);

    Tmax_plot = 150; %[ms] max y limit

    if ifPlotContour

        h = imagesc(Disp_x/1e3, (tmat_slip_event-t_plotshift)*1e3, SlipVelmat_filtered);
        % set(h, 'EdgeColor', 'none');

        % [C,h] = contourf(Disp_x/1e3, (tmat_slip_event-t_plotshift)*1e3, SlipVelmat_filtered, 105);
        % set(h,'LineColor','none')

        % annotation('textbox',[.9 .78 .1 .2], 'String', 'x10^{-3}', 'EdgeColor','none');

        xlabel("Distance from western edge of fault [m]");
        ylabel("Time since nucleation of rupture [ms]");

        set(gca,'ColorScale','log');
        clim([0.2, 5]);

        ch = colorbar();
        ch.Ticks = [0.1, 0.2, 0.5, 1, 2, 5];

        %
        % ch.Ruler.Exponent = 3;
        % ch.Ruler.Scale = 'log';

        % ch.Ruler.Exponent=10;                   % set the desired exponent
        % ch.Ruler.TickLabelFormat='%0.5f';       % fix up ugly default %g formatting
        % ch.Ruler.TickValues=[0.1e-3, 0.2e-3, 0.5e-3, 1e-3, 2e-3, 5e-3];       % fix up ugly default %g formatting

        % set(ch, 'YTickLabel', cellstr(num2str(reshape(get(ch, 'YTick'),[],1),'%2.1e')) )

    else
        set(gca,'ColorScale','log');
        clim([0.2, 5]);

        ch = colorbar();
        ch.Ticks = [0.1, 0.2, 0.5, 1, 2, 5];


    end

    ylabel(ch,'Slip velocity [mm/s]', "FontSize", 16);

    %---Plot AE waveform---%

    bandpass_order = 4;
    fmin_AE = 0.6e5;
    fmax_AE = 6e5;
    ampnorm_AE = 5; %4 %5e-5;

    % apply bandpass filter
    f_nyq = Fs_AE/2;
    [b,a] = butter(bandpass_order, [fmin_AE fmax_AE]/f_nyq, "bandpass");
    %     freqz(b, a, [], fs_read);
    %     ax=gca();
    %     set(ax, 'XScale', 'log');
    AEdatmat_filtered = filtfilt(b,a,AEdatmat);

    lc_AE = [1, 1, 1, 0.1];

    % plot(AEdatmat_filtered/ampnorm_AE + AEsensor_x'/1e3, (tmat_AE_event-t_plotshift)*1e3, "-", "Color", lc_AE);

    %---Emphasize foreshock events---%
    % event_ids = [29,30,31];
    % emp_sensors = {[13, 28, 29]; [9, 10, 25, 26]; [7,8,22,23]};
    % emp_pretrigger = 5e-4;
    % event_winlen = 2e-3;
    %
    % for i = 1:length(event_ids)
    %     Aevent = A(event_ids(i), :);
    %     st = Aevent.Var4-Tstart-emp_pretrigger;
    %     et = st + emp_pretrigger + event_winlen;
    %     plot_sensors = cell2mat(emp_sensors(i, :));
    %
    %     % find time index
    %     emp_tinds = find((st < tmat_AE_event) & (tmat_AE_event < et));
    %     % plot(AEdatmat_filtered(emp_tinds, plot_sensors)/ampnorm_AE + AEsensor_x(plot_sensors)'/1e3, (tmat_AE_event(emp_tinds)-t_plotshift)*1e3, "-", "Color", "r");
    %
    % end

    %---Plot single AE trace---%
    plotstep = 10;

    for AE_trid = 1:size(AEdatmat_filtered, 2)%[23, 27, 13];
        AEtr_lc = [sns_colorpalette(1, 4), 0.02];
        ampnorm_AEtr = 4.2; %4 %5e-5;
        plot(AEdatmat_filtered(1:plotstep:end, AE_trid)/ampnorm_AEtr + AEsensor_x(AE_trid)/1e3, (tmat_AE_event(1:plotstep:end)-t_plotshift)*1e3,...
            "-", "Color", AEtr_lc);

    end


    %---Plot timing of event---%
    % select stick-slip experiment with ordinary events
    Aevents = A((A.Var2==event_id) & (string(A.Var5)=={'Ordinary'}), :);

    for ievent = 1:size(Aevents, 1)
        Aevent = Aevents(ievent, :);
        event_timing = Aevent.Var4-Tstart;
        event_loc = Aevent.Var3;
        plot(event_loc, (event_timing-t_plotshift)*1e3, "wo", "MarkerSize", 10, "LineWidth", 1);
        ht = text(event_loc, (event_timing-t_plotshift)*1e3, string(Aevent.Var7), "Color", "w", "FontSize", 14,...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
        ht.Position(2) = ht.Position(2) - 1;
    end

    % plot(AEdatmat_filtered/ampnorm_AE + AEsensor_x'/1e3, (tmat_AE_event-t_plotshift)*1e3, "-", "Color", "r");

    %---Plot Strain waveform---%

    fmin_strain = 1e4;
    % fmax_strain = ;

    % % apply bandpass filter
    f_nyq = Fs_strain/2;
    [b,a] = butter(4, fmin_strain/f_nyq, "low");
    % freqz(b, a, [], Fs_strain);
    % ax=gca();
    % set(ax, 'XScale', 'log');

    offsetind = 100;
    Taumat2_removeoffset = Taumat2 - mean(Taumat2(1:offsetind, :));
    Taumat3_removeoffset = Taumat3 - mean(Taumat3(1:offsetind, :));

    % Replace NaN to avoid the filter error
    SGB_nanchan = isnan(Taumat2_removeoffset(1, :));
    SGT_nanchan = isnan(Taumat3_removeoffset(1, :));
    Taumat2_removeoffset(:, SGB_nanchan) = 0;
    Taumat3_removeoffset(:, SGT_nanchan) = 0;

    Taumat2_removeoffset_smoothed = filtfilt(b, a, Taumat2_removeoffset .* tukeywin(size(Taumat2_removeoffset, 1),0.05));
    Taumat3_removeoffset_smoothed = filtfilt(b, a, Taumat3_removeoffset .* tukeywin(size(Taumat3_removeoffset, 1),0.05));

    % refill NaN
    Taumat2_removeoffset_smoothed(:, SGB_nanchan) = NaN;
    Taumat3_removeoffset_smoothed(:, SGT_nanchan) = NaN;


    ampnorm_strain = 1;

    % plot(Taumat2_removeoffset_smoothed/ampnorm_strain + SGB_x/1e3, (tmat_strain_event-t_plotshift)*1e3, "-", "Color", "w");
    % plot(Taumat3_removeoffset_smoothed/ampnorm_strain + SGT_x/1e3, (tmat_strain_event-t_plotshift)*1e3, "-", "Color", "w");

    %----------------------%

    %---plot location of gap sensors---%
    slipsensor_color="w"; %"#4DD9FF";
    plot(Disp_x/1e3, zeros(size(Disp_x, 2), 1), "s", "MarkerSize", 10, ...
        "MarkerEdgeColor", "k", "MarkerFaceColor", slipsensor_color);

    xlim([0, 4.1]);
    ylim([0.0, Tmax_plot]); %45]);

    titlestr = sprintf("FB03-%03d stick-slip event %02d: T%.4f-%.4f[s]", expr_id, event_id, Tstart, Tstart+Tmax_plot*1e-3);
    title(titlestr);
    
    figname = sprintf(figdir+"/plot_slipvelocity_master_fb03-%03d_event%02d_contour%d.jpg", expr_id, event_id, ifPlotContour);
    % figname_eps = sprintf(figdir+"/plot_slipvelocity_master_fb03-%03d_event%02d_contour%d.eps", expr_id, event_id, ifPlotContour);

    if ifSavefig
        exportgraphics(fig, figname, 'Resolution',80);
        % exportgraphics(fig, figname_eps, 'ContentType', 'vector');
    end

    %% Save data for the comparison to the seismic moment
    foname=sprintf("../data/FB03-%03d_slipvelocity_event%02d.mat", expr_id, event_id);
    tmat_AE_event_decim = tmat_AE_event(1:plotstep:end);
    AEdatmat_filtered_decim = AEdatmat_filtered(1:plotstep:end, :);

    save(foname, "Disp_x", "tmat_slip_event", "Dmat_event_filtered", "SlipVelmat_filtered", "Aevents", "Tstart", ...
        "tmat_AE_event_decim", "AEdatmat_filtered_decim", "AEsensor_x");

end
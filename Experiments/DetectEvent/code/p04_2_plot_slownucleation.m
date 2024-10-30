% 04-2 plot slow nucleation events
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

% Range to compute rupture velocity
x_bound_typecheck = [1.25, 2.0];
xrange_east = [2.0, 3.0];
xrange_west = [0.5, 2.0];

expr_id = 87;
runID = sprintf('FB03_%03d', expr_id);
T = readtable(sprintf("../data/p02_eventtype_fb03-%03d.csv", expr_id));

event_type = zeros(length(size(T, 1)), 3); %1. event type %2. rupture vel %3. event start time from strain

% for event_id = 1:size(T, 1)
for event_id = [11, 14, 16, 18, 23, 28, 36]
    % event_id = 14;

    event_datdir=sprintf("/Volumes/4mGouge_WorkHDD/FB03data/4mBIAX_paper_tmp/p03_eventdata_FB03_%03d", expr_id);

    load(event_datdir + sprintf("/eventdata_FB03_%03d_event%02d.mat", expr_id, event_id));

    %--------------------------------%
    %% 1. Plot accumulation of slip
    %--------------------------------%

    fig = figure(1); clf; hold on;
    fig.Units = 'point';
    fig.Position = [0 800 1900 850];

    dt_event = 0.1e-3;
    dt_event_k = round(dt_event*Fs_slip);
    Tlength = event_winlen;
    Tevent_k = round(Tlength*Fs_slip);

    cd = colormap('turbo'); % take your pick (doc colormap)

    lc_intp_event = @(x) interp1(linspace(0, Tlength+1e-8, length(cd)),cd, x);
    %%

    % remove offset in the begininng
    Dmat_event_span = Dmat_event(1:dt_event_k:end, :);
    Dmat_event_span = Dmat_event_span - mean(Dmat_event_span(1:5, :), 1);
    Ndat_event = size(Dmat_event_span, 1);
    tmat_event = tmat_slip_event(1:dt_event_k:end);

    subplot(1, 6, [1 2]); box on; hold on;

    for j = 1:Ndat_event
        if mod(j, 5) ==0
            %j
            Lw=1.5; %3.5;
        else
            Lw=1.5; %0.5;
        end
        plot(Disp_x/1e3, Dmat_event_span(j, :), '-', 'Color', lc_intp_event(tmat_event(j)), "LineWidth", Lw);
    end

    % plot sensor location
    plot(Disp_x/1e3, zeros(length(Disp_x), 1)-0.001, 'v', 'Color', "k", 'MarkerSize', 10, "MarkerFaceColor", "w");

    xlabel("x along fault [m]");
    ylabel("Slip [mm]");
    title(sprintf("Slip: plot every %.2fms", dt_event*1e3));

    hcb = colorbar;
    hcb.FontSize = 17;
    ylabel(hcb, 'Time [ms]');
    caxis([0, Tlength]*1e3);

    xlim([0, 4]);
    ylim([-0.002, 0.05]);
    %     ylim([-0.001, 0.035]);


    %--------------------------------------------%
    %% 2. Plot AE and strain with rupture velocity
    %--------------------------------------------%
    subplot(1, 6, [3, 4, 5, 6]); box on; hold on;

    % load event type
    assert(Tstart == T.event_starttime(event_id));

    % if T.eventtype(event_id) == 2
        AE_xlimit = [0, 200]; %[0, 100]
    % else
        % AE_xlimit = [0, 40];
    % end

    %--Plot AE---%
    %%
    ifBandpass = false;
    bandpass_order = 4;
    fmin = 0.6e5;
    fmax = 6e5;
    ampnorm_AE = 20; %5e-5;

    % apply bandpass filter
    if ifBandpass
        f_nyq = fs_read/2;
        [b,a] = butter(bandpass_order, [fmin fmax]/f_nyq, "bandpass");
        %     freqz(b, a, [], fs_read);
        %     ax=gca();
        %     set(ax, 'XScale', 'log');
        AEdatmat = filtfilt(b,a,AEdatmat);

    else
        AEdatmat = AEdatmat;
    end

    % clf(fig,'reset'); cla(fig,'reset'); hold on;
    % box on; grid on;
    Nsensors_AE = size(AEdatmat, 2);
    [~, AE_sortedinds] = sort(AEsensor_x);

    for i = 1:Nsensors_AE
        ii = AE_sortedinds(i);
        y_shift = AEsensor_x(ii)/1e3;
        if ii<16
            lc = [0, 0, 0.6, 0.2]; %[0, 0.4470, 0.7410]; %blue
        else
            lc = [0, 0, 0.6, 0.2]; %[0.8500, 0.3250, 0.0980]; %red
        end

        if ifPlotAE
            plot(tmat_AE_event*1e3, AEdatmat(:, ii)/ampnorm_AE + y_shift, "-", "Color", lc);
        end
    end

    % xline((AE_eventtime - Tstart) * 1e3, "--", "Color", "r", "LineWidth", 2);

    % Plot scale of volt
    scale_x = AE_xlimit(2) * 0.92; %22.5;
    scale_y = 0.2;
    scale_len = 5; %[V]
    plot([scale_x, scale_x], [scale_y-(scale_len/ampnorm_AE)/2, scale_y+(scale_len/ampnorm_AE)/2], "k-", "LineWidth", 3);
    text(scale_x+0.2, scale_y, "5V");


    titlestr = sprintf("AE and Strain: rupvel=%.2f[m/s] ruptype=%d", T.nuc_rupturevel(event_id), T.eventtype(event_id));
    title(titlestr);

    %--Plot strain--%
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

    %% Compute rupture pattern and rupture velocity
    % 1. compute timing of maximum shear strain at trace
    % 2. qualify if the S/N of max strain is larger than threshold
    % 3. linear regression to evaluate the rupture initiation and its velocity
    SNratio_threshold = 300;
    noisewin = 100;
    R = zeros(NSGB+NSGT, 3); %1. coordinate[mm], 2. timing [s] 3. max strain value

    % search SGB
    for ii = 1:NSGB
        [maxeps, imax] = max(Taumat2_removeoffset_smoothed(:, ii));
        noiselevel = mean(abs(Taumat2_removeoffset_smoothed(1:100, ii)));
        SNratio = maxeps/noiselevel;
        if SNratio > SNratio_threshold
            R(ii, :) = [SGB_x(ii)/1e3, tmat_strain_event(imax), maxeps];
        else
            R(ii, :) = [NaN, NaN, NaN];
        end

    end

    % search SGT
    for ii = 1:NSGT
        [maxeps, imax] = max(Taumat3_removeoffset_smoothed(:, ii));
        noiselevel = mean(abs(Taumat3_removeoffset_smoothed(1:100, ii)));
        SNratio = maxeps/noiselevel;

        if SNratio > SNratio_threshold
            R(ii+NSGB, :) = [SGT_x(ii)/1e3, tmat_strain_event(imax), maxeps];
        else
            R(ii+NSGB, :) = [NaN, NaN, NaN];
        end

    end

    %% Compute linear regression

    % To identify the rupture from both sides, we compare the slopes associated
    % with easetern and western of x = 1250[mm], and check the difference of
    % them.
    % We also compute the rupture velocity of nucleation with different range
    % of location for east and west nucleations.

    R(any(isnan(R), 2), :) = [];

    R_east = R(R(:, 1) > x_bound_typecheck(2), :);
    R_west = R(R(:, 1) <= x_bound_typecheck(1), :);

    pall = polyfit(R(:, 2), R(:, 1), 1);

    peast = polyfit(R_east(:, 2), R_east(:, 1), 1);
    pwest = polyfit(R_west(:, 2), R_west(:, 1), 1);

    rupvel = pall(1); % [m/s]

    % Evaluate the rupture type
    % Case2 rupture propagate from west; Note that some events from west are very slow, and no data point within the
    % range of polyfit. In this case, they are categorized in type 2.
    if (pwest(1) == 0) && (peast(1) > 0)
        ruptype = 2;

        R_rupvel_nuc = R((xrange_west(1) <= R(:, 1)) & (R(:, 1) <= xrange_west(2)) , :);
        pvel = polyfit(R_rupvel_nuc(:, 2), R_rupvel_nuc(:, 1), 1);
        rupvel = pvel(1);

        % case 1. rupture propagate from east
    elseif (pwest(1) == 0) && (peast(1) < 0)
        ruptype = 1;
        R_rupvel_nuc = R((xrange_east(1) <= R(:, 1)) & (R(:, 1) <= xrange_east(2)) , :);
        pvel = polyfit(R_rupvel_nuc(:, 2), R_rupvel_nuc(:, 1), 1);
        rupvel = pvel(1);

        % case 3. if the signs of slope in east and west are different, rupture nucleated
        % from both side
    elseif sign(peast(1)) ~= sign(pwest(1))
        ruptype = 3;
        rupvel = NaN;

        % case 1. rupture propagate from east
    elseif sign(rupvel) == -1
        ruptype = 1;
        R_rupvel_nuc = R((xrange_east(1) <= R(:, 1)) & (R(:, 1) <= xrange_east(2)) , :);
        pvel = polyfit(R_rupvel_nuc(:, 2), R_rupvel_nuc(:, 1), 1);
        rupvel = pvel(1);
        % case 2. rupture propagate from west

    elseif sign(rupvel) == 1
        ruptype = 2;

        R_rupvel_nuc = R((xrange_west(1) <= R(:, 1)) & (R(:, 1) <= xrange_west(2)) , :);
        pvel = polyfit(R_rupvel_nuc(:, 2), R_rupvel_nuc(:, 1), 1);
        rupvel = pvel(1);
    end

    event_type(event_id, 1) = ruptype;
    event_type(event_id, 2) = rupvel;
    event_type(event_id, 3) = Tstart;

    %% Check if the estimated rupture velocity is close to p02
    % NOTE: As we trimmed the strain traces and the moving window average
    % results in a slightly different smoothing causing the shift of imax, the
    % slope of polyfit is not identical to the p02.
    % However, it should be close each other.
%     assert(T.eventtype(event_id) == ruptype)
%     if ~isnan(rupvel)
%         assert(abs(T.nuc_rupturevel(event_id) - rupvel)/abs(T.nuc_rupturevel(event_id)) < 0.3); % assert if the difference is less than 5 %.
%     end
    %--------------------------------------------%
    %%
    ampnorm = 0.4; %5e-5;

    % plot rupture velocity
    for i = 1:size(R, 1)
        plot(R(i, 2)*1e3, R(i, 1), "s", "MarkerSize", 12, "MarkerFaceColor","w", "MarkerEdgeColor","r");
        %     plot(R(i, 2)*1e3, R(i, 1) + R(i, 3)/ampnorm, "bo");
    end

    % plot SGB
    % lc = [191/255, 134/255, 134/255, 0.9];
    % lc = [191/255, 0, 0, 1.0];
    lc = "k";

    for ii = 1:NSGB
        %     if ii<16
        %         lc = "b"; %[0, 0.4470, 0.7410]; %blue
        %     else
        %         lc = "b"; %[0.8500, 0.3250, 0.0980]; %red
        %     end
        y_shift = SGB_x(ii)/1e3;
        plot(tmat_strain_event*1e3, Taumat2_removeoffset_smoothed(:, ii)/ampnorm + y_shift, "-", "Color", lc, "LineWidth", 1.2);
    end

    % plot SGT
    for ii = 1:NSGT
        %     if ii<16
        %         lc = [0, 0.4470, 0.7410, 0.6]; %blue
        %     else
        %         lc = [0, 0.4470, 0.7410, 0.6]; %[0.8500, 0.3250, 0.0980]; %red
        %     end
        y_shift = SGT_x(ii)/1e3;
        plot(tmat_strain_event*1e3, Taumat3_removeoffset_smoothed(:, ii)/ampnorm + y_shift, "-", "Color", lc, "LineWidth", 1.2);
    end


    if ruptype <= 2 % skip ruptyre == 3
        xq = linspace(min(R(:, 2)), max(R(:, 2)), 101);
        yq = polyval(pvel, xq);
        plot(xq*1e3, yq, "--", "Color", "r", "LineWidth", 1.5);
    end

    % xline((AE_eventtime - Tstart) * 1e3, "--", "Color", "blue")
    % plot((AE_eventtime - Tstart) * 1e3, eventloc, "ro", "Markersize", 10)

    xlabel("Time [ms]");
    ylabel("x along fault [m]");
    % xlim([0, Tlength*1e3]);
    xlim(AE_xlimit);
    ylim([-0.1, 4.1]);

    % Plot scale of shear stress
    scale_x = AE_xlimit(2) * 0.85; %23.5;
    scale_y = 0.2;
    scale_len = 0.1; %[MPa]
    plot([scale_x, scale_x], [scale_y-(scale_len/ampnorm)/2, scale_y+(scale_len/ampnorm)/2], "k-", "LineWidth", 3);
    text(scale_x+0.2, scale_y, "0.1MPa");

    %%
    % xline(15.2);
    % xline(19.2);
    sgtitlestr = sprintf("FB03-%03d event %02d: T%4.4f-%4.4f[s]", ...
        expr_id, event_id, Tstart, Tstart+Tlength);
    % hSG = sgtitle(sgtitlestr, "FontWeight", "bold", "FontSize", 22);
    text(-5, 4.35, sgtitlestr, "FontWeight", "bold", "FontSize", 21);
    pan on;

    gca();
    pos = get(gca, 'Position');
    pos(1) = pos(1) + 0.03;
    set(gca, 'Position', pos)


    figdir = sprintf("../figure/p04_plotevent_FB03_%03d", expr_id);
    if ~exist(figdir) mkdir(figdir); end


    figname = sprintf(figdir+"/p04_plot_event_FB03_%03d_event%02d_slownucleation.jpg", expr_id, event_id);

    if ifSavefig
        exportgraphics(fig, figname, 'Resolution',60); %150
    end

end
% %% Renew event type
% 
% event_id = transpose(1:size(event_type, 1));
% event_starttime = event_type(:, 3);
% eventtype = event_type(:, 1);
% nuc_rupturevel = event_type(:, 2);
% 
% T = table(event_id, event_starttime, eventtype, nuc_rupturevel);
% writetable(T,sprintf("../data/p04_renewedeventtype_fb03-%03d.csv", expr_id));

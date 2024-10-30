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
ifSavefig = false;
Fs = 1e6;
cs = 3600; % [m/s] shear wave velocity of rock sample
eventst_pretrig = 35e-3; %15e-3; % [s] window length of pretrigger of event

tevent_FW = 100e-3; % event window before (forward) starttime 
tevent_BW = 50e-3; % event window after (back) starttime 

xscale_loc = 135; %[ms] location of scale bar


% Range to compute rupture velocity
x_bound_typecheck = [1250, 2000];
xrange_east = [2000, 3000];
xrange_west = [500, 2000];
%%
SGB_x = horzcat(290+(0:15)*240, 110+(0:15)*240);
SGT_x = horzcat(170+(0:15)*240, 230+(0:15)*240);

NSGB=length(SGB_x);
NSGT=length(SGT_x);

%%
% for expr_id = [87, 94]
expr_id = 87;

%% Correct the relative location of the top rock specimen
if expr_id == 87
    xtop_shift = 82.46; %u_east (at t=0) = 13.57mm, u0 = -3.97
elseif expr_id == 94
    xtop_shift = 90.23;
elseif expr_id == 99
    xtop_shift = 89.91;
else
    error("the relative location of top block is not defined.");
end

% adjust to the origin of the bottom rock specimen (west)
SGB_x = SGB_x + xtop_shift;
SGT_x = SGT_x + xtop_shift;


%% Load event init time
tevent_start=importdata(sprintf("../data/FB03_%03d_p01_eventstarttime.csv", expr_id));

pname = sprintf('/Volumes/4mGouge_WorkHDD/FB03data/fb03-%03d/biax/', expr_id);
fname = sprintf('fb03-%03d', expr_id);

event_type = zeros(length(tevent_start), 3); %1. event type %2. rupture vel %3. event start time from strain
for ievent = 1:length(tevent_start)
% ievent = 5;

st = tevent_start(ievent);

%% Load strain

rdecim = 0;
Tstart = st - tevent_FW;
Tlength = tevent_FW + tevent_BW;

[tmat,Snmat,Taumat3,Taumat2,Spmat]=fmbiaxlocalSTHFS3(pname,fname, Tstart, Tlength, rdecim);

fs_read = 1/(tmat(2)-tmat(1));

%% Plot along fault
offsetind = 100;
Taumat2_removeoffset = Taumat2 - mean(Taumat2(1:offsetind, :));
Taumat3_removeoffset = Taumat3 - mean(Taumat3(1:offsetind, :));

%% apply moving window average
smooth_winlen = 100;
Taumat2_removeoffset_smoothed = movmean(Taumat2_removeoffset, smooth_winlen, 1);
Taumat3_removeoffset_smoothed = movmean(Taumat3_removeoffset, smooth_winlen, 1);

%% ------%%
% Compute rupture pattern and velocity
%--------%
% 1. compute timing of maximum shear strain at trace
% 2. qualify if the S/N of max strain is larger than threshold
% 3. linear regression to evaluate the rupture initiation and its velocity
SNratio_threshold = 300;
noisewin = 100;
R = zeros(NSGB+NSGT, 3); %1. coordinate[mm], 2. timing [s] 3. max strain value

% search SGB
for ii = 1:NSGB
% ii = 10;
[maxeps, imax] = max(Taumat2_removeoffset_smoothed(:, ii));

% figure(1); clf; hold on;
% plot(tmat*1e3, Taumat2_removeoffset_smoothed(:, ii));
noiselevel = mean(abs(Taumat2_removeoffset_smoothed(1:100, ii)));
SNratio = maxeps/noiselevel;
if SNratio > SNratio_threshold
    R(ii, :) = [SGB_x(ii), tmat(imax), maxeps];
%     plot(tmat(imax)*1e3, Taumat2_removeoffset_smoothed(imax, ii), "ro");

else
    R(ii, :) = [NaN, NaN, NaN];
%     plot(tmat(imax)*1e3, Taumat2_removeoffset_smoothed(imax, ii), "ko");
%     pause(1)
end

end

% search SGT
for ii = 1:NSGT
% ii = 10;
[maxeps, imax] = max(Taumat3_removeoffset_smoothed(:, ii));

% figure(1); clf; hold on;
% plot(tmat*1e3, Taumat3_removeoffset_smoothed(:, ii));

noiselevel = mean(abs(Taumat3_removeoffset_smoothed(1:100, ii)));
SNratio = maxeps/noiselevel;

if SNratio > SNratio_threshold
    R(ii+NSGB, :) = [SGT_x(ii), tmat(imax), maxeps];
%     plot(tmat(imax)*1e3, Taumat3_removeoffset_smoothed(imax, ii), "ro");

else
    R(ii+NSGB, :) = [NaN, NaN, NaN];
%     plot(tmat(imax)*1e3, Taumat3_removeoffset_smoothed(imax, ii), "ko");
%     pause(1)
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

rupvel = pall(1) * 1e-3; % [mm/s] -> [m/s]

% Evaluate the rupture type
% Case2 rupture propagate from west; Note that some events from west are very slow, and no data point within the
% range of polyfit. In this case, they are categorized in type 2.
if (pwest(1) == 0) && (peast(1) > 0) 
    ruptype = 2;

    R_rupvel_nuc = R((xrange_west(1) <= R(:, 1)) & (R(:, 1) <= xrange_west(2)) , :);
    pvel = polyfit(R_rupvel_nuc(:, 2), R_rupvel_nuc(:, 1), 1);
    rupvel = pvel(1)*1e-3;

% case 1. rupture propagate from east  
elseif (pwest(1) == 0) && (peast(1) < 0) 
    ruptype = 1;
    R_rupvel_nuc = R((xrange_east(1) <= R(:, 1)) & (R(:, 1) <= xrange_east(2)) , :);
    pvel = polyfit(R_rupvel_nuc(:, 2), R_rupvel_nuc(:, 1), 1);
    rupvel = pvel(1)*1e-3;

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
    rupvel = pvel(1)*1e-3;
% case 2. rupture propagate from west

elseif sign(rupvel) == 1
    ruptype = 2;

    R_rupvel_nuc = R((xrange_west(1) <= R(:, 1)) & (R(:, 1) <= xrange_west(2)) , :);
    pvel = polyfit(R_rupvel_nuc(:, 2), R_rupvel_nuc(:, 1), 1);
    rupvel = pvel(1)*1e-3;
end

event_type(ievent, 1) = ruptype;
event_type(ievent, 2) = rupvel;
event_type(ievent, 3) = Tstart + min(R(:, 2)) - eventst_pretrig; % event start time estimated from strain

%%
% figure(3); clf;
% plot(R(:, 2), R(:, 1), "o");

%%
fig = figure(2);
fig.Units = 'point';
fig.Position = [0 800 1900 800];
clf; hold on; box on;

ampnorm = 0.4e-3; %5e-5;
plot_event_st = (min(R(:, 2)) - eventst_pretrig);

% plot SGB
% lc = [191/255, 134/255, 134/255, 0.9];
lc = [191/255, 0, 0, 1.0];

for ii = 1:NSGB
    y_shift = SGB_x(ii);
    plot(tmat*1e3, Taumat2_removeoffset_smoothed(:, ii)/ampnorm + y_shift, "-", "Color", lc, "LineWidth", 1.0);
end

% plot SGT
for ii = 1:NSGT
    y_shift = SGT_x(ii);
    plot(tmat*1e3, Taumat3_removeoffset_smoothed(:, ii)/ampnorm + y_shift, "-", "Color", lc, "LineWidth", 1.0);
end

% plot rupture velocity
for i = 1:size(R, 1)
    plot(R(i, 2)*1e3, R(i, 1), "bs");
%     plot(R(i, 2)*1e3, R(i, 1) + R(i, 3)/ampnorm, "bo");
end


if ruptype <= 2 % skip ruptyre == 3       
    xq = linspace(plot_event_st, plot_event_st+100e-3, 101);
    yq = polyval(pvel, xq);
    plot(xq*1e3, yq, "-", "Color", [0 0.5 0]);
end

% plot event start time
xline(plot_event_st*1e3, "k--", "LineWidth", 0.5);

% xline((AE_eventtime - Tstart) * 1e3, "--", "Color", "blue")
% plot((AE_eventtime - Tstart) * 1e3, eventloc, "ro", "Markersize", 10)

xlabel("Time [ms]");
ylabel("x along fault [mm]");
% xlim([0, Tlength*1e3]);
xlim([0, Tlength*1e3]);
ylim([-100, 4100]);

% Plot scale of shear stress
scale_x = xscale_loc;
scale_y = 200;
scale_len = 0.1; %[MPa]
plot([scale_x, scale_x], [scale_y-(scale_len/ampnorm)/2, scale_y+(scale_len/ampnorm)/2], "r-", "LineWidth", 3);
text(scale_x+0.2, scale_y, "0.1MPa");

titlestr = sprintf("%s  Fs=%.1f[MHz] event %d T%4.4f-%4.4f[s] rupvel=%.2f[m/s] ruptype=%d",...
                fname, fs_read/1e6, ievent, Tstart, Tstart+Tlength, rupvel, ruptype);
title(titlestr);

figdir_event_strain=sprintf("../figure/p02_strain_%s", fname);
if ~exist(figdir_event_strain) mkdir(figdir_event_strain); end

figname = sprintf(figdir_event_strain+"/%s_event%02d_p02_strain.png", fname, ievent);

% figname = sprintf(figdir_event_ae+"/%s_event%ievent_AEonly_forLFE.png", runID, ievent);
if ifSavefig
    exportgraphics(fig, figname, 'Resolution',300);
end

end

%% Save event type

event_id = transpose(1:size(event_type, 1));
event_starttime = event_type(:, 3);
eventtype = event_type(:, 1);
nuc_rupturevel = event_type(:, 2);

T = table(event_id, event_starttime, eventtype, nuc_rupturevel);
writetable(T,sprintf("../data/p02_eventtype_%s.csv", fname));

% Script to compute ARX model of sensor response
% 2022/07/12 Kurama Okubo
% update: for master script
clear all;
set(0,'DefaultTextFontsize',14, ...
    'DefaultTextFontname','Arial', ...
    'DefaultTextFontWeight','normal', ...
    'DefaultTextFontname','Arial', ...
    'DefaultAxesFontsize',14, ... 
    'DefaultAxesFontname','Arial', ...
    'DefaultLineLineWidth', 1.0)
set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'})

addpath("../../AEsensor_Calibration_ARX/src");
%%
figdir_arx = "../figure/TraceARX/";
if ~isfolder(figdir_arx); mkdir(figdir_arx); end

%% Read input (LDV) and output (AEsensor) data

casename_list = ["frontcenter", "fronttop", "sidecenter"];

ic = 2; % Choose fronttop
ifscale = true;

casename = casename_list(ic);
LDV = readtable(sprintf("../../AEsensor_Calibration_DataConversion/data/LDV_%s_S200V.csv", casename));
AE = readtable(sprintf("../../AEsensor_Calibration_DataConversion/data/AE_%s_S100V.csv", casename));
assert(all(LDV.t == AE.t));

%%
S.t = LDV.t;
S.LDV = LDV.LDV/1.79; %As the source is 200V on LDV while 100V on AE, devide LDV by 1.79 from linearity analysis.
S.AE = AE.AE;

%%
figure(1);
plot(S.t, S.LDV*3000, S.t, S.AE);
% xlim([0, 100e-6]);

%%
%---Parameters used to estimate ARX model---%
fs = 1e7; % sampling frequency
tlen = 100e-6; % [s] data length % master data with 100 us
NT = round(tlen * fs)+1; % npts of data

u_input_raw = S.LDV(1:NT); % exgeneous input data
y_output_raw = S.AE(1:NT); % output data
tvec = S.t(1:NT)*1e6;

%% apply bandpass

fmin = 0.2e5; %[Hz]

% switch fmax with cases
if ic <=2
    fmax = 2e5; %[Hz] master data with 100kHz
else
    fmax = 1e5; %[Hz] master data with 100kHz
end

% butterworth filter
[b,a]=butter(3, [fmin,fmax]/(fs/2),'bandpass'); % Master data: chose 3rd order.

% chebychef type 2
% [b,a] = cheby2(3,10,fmax/(fs/2));

% detrend
u_input_raw_demean = detrend(u_input_raw - mean(u_input_raw(1:100)));
y_output_raw_demean = detrend(y_output_raw - mean(y_output_raw(1:100)));

% tapered
u_input_raw_tapered = u_input_raw_demean.*tukeywin(NT,0.1);
y_output_raw_tapered = y_output_raw_demean.*tukeywin(NT,0.1);

% bandpass filtered
u_input_filtered = filter(b, a, u_input_raw_tapered);
y_output_filtered = filter(b, a, y_output_raw_tapered);

% tapered
u_input_filtered = u_input_filtered.*tukeywin(NT,0.1);
y_output_filtered = y_output_filtered.*tukeywin(NT,0.1);

u_input = u_input_filtered;
y_output = y_output_filtered;

% No pre-processing if removing comment out below.
% u_input = u_input_raw;
% y_output = y_output_raw;

%%
fig = figure(2);
fig.Units = 'point';
fig.Position = [0 800 800 240];
clf(fig,'reset'); cla(fig,'reset'); hold on;
box on; grid on;

yyaxis left;
plot(tvec, u_input*1e3, "k-", "DisplayName","LDV input");
ylim([-0.15, 0.1500001]);
ylabel("Velocity [mm/s]");

yyaxis right;
plot(tvec, y_output, "r-", "DisplayName","AEsensor output");
ylim([-1.5, 1.5]);
ylabel("Amplitude [V]");

xlim([0, 80]);%[0, tlen*1e6])
xlabel("Time [μs]");
legend("Location","northeast");

ax = gca;
ax.YAxis(1).Color="k";
ax.YAxis(2).Color="r";

title(sprintf("Bandpass %4.2f-%4.1fMHz", fmin/1e6, fmax/1e6))

set(gcf, 'Color', 'w');
figname = sprintf("../figure/TraceARX/LDV_and_AE_filtered_%s.pdf", casename);
exportgraphics(fig, figname);


%% Plot raw data
fig = figure(3);
fig.Units = 'point';
fig.Position = [0 800 800 240];
clf(fig,'reset'); cla(fig,'reset'); hold on;
box on; grid on;

yyaxis left;
plot(tvec, u_input_raw_tapered*1e3, "k-", "DisplayName","LDV input");
ylim([-0.8, 0.8]);
ylabel("Velocity [mm/s]");

yyaxis right;
plot(tvec, y_output_raw_tapered, "r-", "DisplayName","AEsensor output");
ylim([-5, 5]);
ylabel("Amplitude [V]");

xlim([0, 100]);%[0, tlen*1e6])
xlabel("Time [μs]");
legend("Location","northeast");

ax = gca;
ax.YAxis(1).Color="k";
ax.YAxis(2).Color="r";

title("Detrended and Tapered Raw data");

set(gcf, 'Color', 'w');
figname = sprintf("../figure/TraceARX/LDV_and_AE_raw_%s.eps", casename);
exportgraphics(fig, figname);


%%
% estimate_arx(u_input, y_output, Na, Nb);
AIC_best = 1e9; %initialize with large AIC
Na_best = 0;
Nb_best = 0;

NP = 30; % number of parameters

Ndeg = 0;
for na = 1:NP
    for nb = 1:NP
        Ndeg = Ndeg + 1;
    end
end

% scaling amplitude with std
scaling_u = std(u_input);
scaling_y = std(y_output);
u_input_scaled = u_input./scaling_u;
y_output_scaled = y_output./scaling_y;
u_normfact = scaling_u/scaling_y; % multiply this to the normalized output

all_AIC=zeros(Ndeg, 3);
AIC_mask=zeros(NP, NP);

i = 0;
for na = 1:NP
    for nb = 1:NP
        
        if ifscale
            [theta, AIC] = lsq_arx(u_input_scaled, y_output_scaled, na, nb);
        else
            [theta, AIC] = lsq_arx(u_input, y_output, na, nb);
        end
        
        i = i + 1;
        all_AIC(i, 1) = na;
        all_AIC(i, 2) = nb;
        all_AIC(i, 3) = AIC;

        
        tf_a = [1; theta(1:na)];
        tf_b = theta(na+1:end);
        [z,p,k] = tf2zpk(tf_b,tf_a);
        
        [warnmsg, msgid] = lastwarn;
        if strcmp(msgid,'MATLAB:nearlySingularMatrix') % inverse matrix is unstablly solved
%             na, nb
            warning('reset');
            fprintf("lsq unstable.");
            AIC_mask(na, nb) = 1.0;
            continue;
        end
        
        % skip if nb>=na
        if nb>na
            fprintf("mask nb>=na, improper.\n");
            AIC_mask(na, nb) = 1.0;
            continue;
        end
        
        if all(abs(z) <= 1) && all(abs(p) <= 1)
            fprintf("Na, Nb = (%d, %d) is stable\n", na, nb);
            if AIC < AIC_best
                AIC_best = AIC
                Na_best = na;
                Nb_best = nb;
            end
        else
            if AIC_mask(na, nb) == 0.0
                fprintf("%d %d; mask PAZ outside of unit circle.\n", na, nb);
            end
            AIC_mask(na, nb) = 1.0; % poles and zeros are outside of unit circle
        end
    end
end

%% Plot AIC
fig = figure(4); clf;
fig.Units = 'point';
fig.Position = [0 800 700 600];
% clf(fig,'reset'); cla(fig,'reset'); hold on;
box on; grid on;

x = all_AIC(:, 1);
y = all_AIC(:, 2);
z = all_AIC(:, 3);
N = NP; % Number Of Points Desired
xv = linspace(min(x), max(x), N);
yv = linspace(min(y), max(y), N);
[X,Y] = ndgrid(xv, yv);
Z = griddata(x, y, z, X, Y);

contourf(X, Y, Z, 100,'edgecolor','none');
axis('equal')

xlabel("Na");
ylabel("Nb");
zlabel("AIC");

colormap("jet");
% caxis([-5e4, -2e4]);
hc = colorbar;
hc.Label.String = 'AIC';

set(gcf, 'Color', 'w');
figname = sprintf("../figure/AIC_all_%s.pdf", casename);
exportgraphics(fig, figname);

%%

Na = Na_best; % 20; %3 % order of Auto-regressive model
Nb = Nb_best; %6; %2 % order of Moving-average associated with input

%NOTE: the preliminary result of MT inversion is donw with (Na, Nb) = (20, 6)

if ifscale
    [theta, AIC] = lsq_arx(u_input_scaled, y_output_scaled, Na, Nb);
else
    [theta, AIC] = lsq_arx(u_input, y_output, Na, Nb);
end
    
%% Deconvolution test to evaluate u_input
% Compute upred algebraically
upred = zeros(NT, 1);
b0 = theta(Na+1);
for k = 1:NT
    AY = 0;
    BU = 0;
    
    for l = 1:Na
        if k-l > 0
            AY = AY + theta(l) * y_output_raw_tapered(k-l);
        end
    end
    for l = Na+2:Na+Nb+1
        if k-l+Na+1 > 0
            BU = BU + theta(l) * upred(k-l+Na+1);
        end
    end

    upred(k) = (y_output_raw_tapered(k) + AY - BU)/b0;

end

% unscale upred
if ifscale
    upred = upred * u_normfact;
end

%%
fig = figure(5);
fig.Units = 'point';
fig.Position = [0 800 800 500];
clf(fig,'reset'); cla(fig,'reset'); hold on;
box on;
tvec = S.t(1:NT)*1e6;

% plot(tvec, u_input*2000, "k-");
plot(tvec, y_output_raw_tapered/8e3, "--", "Color", "b", "DisplayName", "AE output");
plot(tvec, u_input_raw_tapered, "k-", "DisplayName", "LDV");
plot(tvec, upred, "r-", "DisplayName","Predicted with AE sensor");

legend("Location","best");

% plot(tvec, real(upred_decon), "k-");
% plot(tvec, ypred1/k, "l-");
grid on;
xlim([10, 100]);
ylim([-0.6e-3, 0.6e-3]);
xlabel("Time [μs]");
ylabel("Velocity [mm/s]");

set(gcf, 'Color', 'w');
figname = "../figure/ARX_comparison_450_9.pdf";
exportgraphics(fig, figname);

%%
fig = figure(6);
fig.Units = 'point';
fig.Position = [0 800 800 500];
clf(fig,'reset'); cla(fig,'reset'); hold on;
box on; grid on;

fmin = 0.06e6;
fmax = 3e5;

% butterworth filter
[b,a]=butter(3, [fmin,fmax]/(fs/2),'bandpass');

% chebychef type 2
% % [b,a] = cheby2(6,40,fmax/(fs/2));

plot(tvec, filtfilt(b, a, u_input_raw_tapered.*hamming(NT))*1e3, "k-", "DisplayName", "LDV input");
plot(tvec, filtfilt(b, a, upred.*hamming(NT))*1e3, "r-", "DisplayName", "Modeled");

xlim([0, 100])
% ylim([-2e-1, 2e-1]);

xlabel("Time [μs]");
ylabel("Velocity [mm/s]");
legend()

set(gcf, 'Color', 'w');
figname = "../figure/ARX_comparison_model_450_9.pdf";
exportgraphics(fig, figname);

%%
tf_a = [1; theta(1:Na)]; % assuming a0 is 1
tf_b = theta(Na+1:end);
[z,p,k1] = tf2zpk(tf_b,tf_a);
k = k1/u_normfact;
assert(abs(real(prod(p)) - tf_a(end))<1e-6); % prod(p) should be tf_a[end]

%% Plot pole and zeros
fig = figure(7);
fig.Units = 'point';
fig.Position = [0 800 800 500];
clf(fig,'reset'); cla(fig,'reset'); hold on;
box on; grid on;

[hz1, hp1, ht1] = zplane(z, p);%with the unit circle
set(findobj(hz1, 'Type', 'line'), 'Color', 'r'); 
set(findobj(hp1, 'Type', 'line'), 'Color', 'b');
set(findobj(ht1, 'Type', 'line'), 'Color', 'k');
legend({"zeros", "poles", ""})
% H = tf(B,A); %set the transferfunction
% pzmap(H) % without the unit circle
title('Pole-zero map') 
set(gcf, 'Color', 'w');
figname = sprintf("../figure/pzplot_%s.pdf", casename);
exportgraphics(fig, figname);


%% dump data to use in MT inversion
foname = sprintf("../data/AE_resp_dataandcoef_%s.mat", casename);
save(foname, "u_input_raw", "y_output_raw", "tf_a", "tf_b", "z",...
    "p", "k", "Na_best", "Nb_best", "all_AIC", "AIC_mask", "AIC_best", "u_normfact");

% Script to find the main shock timing from macroscopic shear loading
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

%% load macro data
load("../data/MacroData_raw.mat", "M");

%% load event start time from gap sensor
st_gap = importdata("../../DetectEvent/data/FB03_087_p01_eventstarttime.csv");


%% Compute diff of shear loading
dt = M.FB03_087.tmat(2) - M.FB03_087.tmat(1);
del_shear = diff(M.FB03_087.SS)/dt;
del_tmat = M.FB03_087.tmat(1:end-1);

%% Find peaks of diff
[pks,locs] = findpeaks(abs(del_shear),'MinPeakHeight', 1, 'MinPeakDistance',100);
st_gap_tmp = del_tmat(locs); 
%% Compare to the start time from gap sensor
% search the event start times near the table from gap sensor
Nevent = length(st_gap); 
st_pks_thresh = zeros(Nevent,1);
st_locs_thresh = zeros(Nevent,1);
dt_thresh = 120e-3; % detection threshold near the start time 
for i = 1:Nevent
    st = st_gap(i);
    dt_ev = st_gap_tmp - st;
    stind = find(abs(dt_ev) < dt_thresh); 
    if length(stind)>1
        error("detection more than 2");
    end
    st_pks_thresh(i) = pks(stind);
    st_locs_thresh(i) = locs(stind);

end

%% Plot Main stick-slip experiments
fig = figure(1);
fig.Units = 'point';
fig.Position = [0 500 800 450];
clf(fig,'reset'); cla(fig,'reset'); hold on; box on; 


plot(M.FB03_087.tmat, M.FB03_087.SS, "k-", "DisplayName", "FB03-087");
xlabel("Time [s]");
ylabel("Shear stress [MPa]");
ylim([0.6, 0.7]);

yyaxis right;

plot(del_tmat, abs(del_shear), "r-", "DisplayName", "FB03-087 diff");
plot(del_tmat(st_locs_thresh), abs(st_pks_thresh), "ro", "DisplayName", "FB03-087 diff");

xlim([0, 280]);

ylabel("Diff. shear stress [MPa/s]");

exportgraphics(gcf, sprintf("../figure/Detect_mainshocks_FB03-087.png"), "Resolution", 80); %300);

%% Save mainshock timing

st_mainshock = del_tmat(st_locs_thresh);
writematrix(st_mainshock, '../data/FB03_087_mainshock_times.csv');



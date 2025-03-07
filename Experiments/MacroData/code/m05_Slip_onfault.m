% Script to plot the evolution of slip during the stick-slip events
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

%% Load slip

S = struct();

expr_ids = [87, 94]; %99]; %87, 94];
tends = [13e6, 11e6, 18e6]; % Trim the data to reduce the size
% for expr_id = [87]
for i = 1:length(expr_ids)
    clear Dmat
    expr_id = expr_ids(i);
    tend = tends(i);
    pname = sprintf('/Volumes/Okuboetal2025_masterHDD/FB03data/fb03-%03d/disp/', expr_id);
    runID = sprintf('fb03-%03d', expr_id);
    
    pname, runID
    [tmat,ddatmat]=fmbiaxdispdatread(pname,runID);
    Nsensor = size(ddatmat, 2);
    
     % convert from voltage to disp
    Gsens=load('../../../utils/matlabcode_biax_v03/GAPsensitivityS3.txt');
    offrange = 100;
    for m=1:Nsensor
        Dmat(:,m)=-(ddatmat(:,m)-mean(ddatmat(1:offrange,m)))/Gsens(m,2);
    end

    tmat = tmat(1:tend);
    Dmat = Dmat(1:tend, :);

    S.(sprintf("FB03_%03d", expr_id)).tmat=tmat;
    S.(sprintf("FB03_%03d", expr_id)).Dmat=Dmat;
     
end

save("../data/SlipData_raw.mat", "S", '-v7.3');

%% Load decimated data
load("../data/SlipData_raw.mat");

Nsensor = size(Dmat, 2);

%%
%Location of displacement sensors
x_disp = 270 + (0:15)*240;

Fs = 50e3; % sampling frequency
Tmax = 210; % maximum time to plot

%% Plot slip


% plot(S.FB03_094.tmat(:), S.FB03_094.Dmat(:, 8));
plot(S.FB03_087.tmat(:), S.FB03_087.Dmat(:, 8));

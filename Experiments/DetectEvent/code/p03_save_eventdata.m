% 03 compile and dump the event data

% We checked the reproducibility of the event data from raw continuous
% data by auxp03_check_reproducibility_eventdata.m

clear all;
set(0,'DefaultTextFontsize',16, ...
    'DefaultTextFontname','Arial', ...
    'DefaultTextFontWeight','normal', ...
    'DefaultTextFontname','Arial', ...
    'DefaultAxesFontsize',16, ...
    'DefaultAxesFontname','Arial', ...
    'DefaultLineLineWidth', 1.0)
set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'})

figdir="../figure";
if ~exist(figdir) mkdir(figdir); end
addpath("../../../utils/matlabcode_biax_v03");

%% Load event data

event_winlen = 200e-3; %[s] event window length;

Fs_slip = 50e3; % sampling frequency of slip
Fs_strain = 1e6; % sampling frequency of strain
Fs_AE = 1e7; % sampling frequency of AE

winlen_slip = round(event_winlen * Fs_slip);
winlen_strain = round(event_winlen * Fs_strain);
winlen_AE = round(event_winlen * Fs_AE);

Disp_x = 270 + (0:15)*240;
SGB_x = horzcat(290+(0:15)*240, 110+(0:15)*240);
SGT_x = horzcat(170+(0:15)*240, 230+(0:15)*240);

NSGB=length(SGB_x);
NSGT=length(SGT_x);

AEloc = readtable("../data/AEsensorlocation_onFB03_table.csv");
AEsensor_x = AEloc.North;
% [~, AE_sortedinds] = sort(AEsensor_x);

expr_id = 87;
runID = sprintf('FB03_%03d', expr_id);

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

% Load event start time estimated from strain

T = readtable(sprintf("../data/p02_eventtype_fb03-%03d.csv", expr_id));

% load slip data
load("../../MacroData/data/SlipData_raw.mat");

%% loop event

for i = 1:length(T.event_id) %1:3
    event_id = T.event_id(i);
    % event_id = 17;

    % Tstart_prewin = 20e-3; % include pretrigger window; update: We added
    % this in eventst_pretrig of p02.
    % Tstart = T.event_starttime(event_id) - Tstart_prewin;
    Tstart = T.event_starttime(event_id);

    % Trim slip
    st_ind = find(S.(runID).tmat >= Tstart, 1);
    tmat_slip_event = (0:winlen_slip-1)./Fs_slip;
    Dmat_event = S.(runID).Dmat(st_ind:st_ind+winlen_slip-1, :);

    % Trim strain
    pname_strain = sprintf('/Volumes/Okuboetal2025_masterHDD/FB03data/fb03-%03d/biax/', expr_id);
    fname = sprintf('fb03-%03d', expr_id);

    rdecim = 0;

    [tmat_strain_event,Snmat,Taumat3,Taumat2,Spmat]=fmbiaxlocalSTHFS3(pname_strain,fname, Tstart, event_winlen, rdecim);

    %% Trim AE
    fs_read = 1e7;
    pname_ae = sprintf('/Volumes/Okuboetal2025_masterHDD/FB03data/fb03-%03d/ae/', expr_id);
    [AEdatmat,tmat_AE_event,Ch_info]=SBENCHreader(Tstart,event_winlen,pname_ae,fname, 'fs_read', fs_read, 'pretrigger', false, 'rename', true);

    %% Load macroscopic data at Fs=1MHz
    [tmat_macro,Macro_datmat] = biaxrawdatreadneo(pname_strain,fname,1,Tstart,event_winlen);

    % convert the data to physical quantities
        
    % === Remove initial offset of displacement ===
    offrange=100;
    DX=Macro_datmat(:,10:11);
    DXwest=DX(:,1)-mean(DX(1:offrange,1));
    DXeast=-DX(:,2)+mean(DX(1:offrange,2));
    
    % === Mechanical data ===
    NPmacro=(Macro_datmat(:,2:9)-2.5e6)*1e-6;
    SSmacro=Macro_datmat(:,1)*1e-6/(4*0.1);

    % Save event data
    event_datdir=sprintf("/Volumes/Okuboetal2025_masterHDD/4mBIAX_eventdata_master/p03_eventdata_FB03_%03d", expr_id);
    if ~exist(event_datdir) mkdir(event_datdir); end

    save(event_datdir + sprintf("/eventdata_FB03_%03d_event%02d.mat", expr_id, event_id), ...
        "Tstart", "tmat_slip_event", "Dmat_event", "tmat_strain_event", "Snmat",...
        "Taumat3", "Taumat2", "Spmat", "tmat_AE_event", "AEdatmat", "tmat_macro", "DXwest", ...
        "DXeast", "NPmacro", "SSmacro",...
        "event_winlen", "Fs_slip", "Fs_strain","Fs_AE", "Disp_x", "SGB_x", "SGT_x",...
        "AEsensor_x");

end
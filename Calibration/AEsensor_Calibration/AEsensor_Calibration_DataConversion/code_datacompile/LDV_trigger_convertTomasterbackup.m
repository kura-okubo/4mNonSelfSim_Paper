% Script to convert pure binary to mat file of the LDV and AE data
% 2022/09/07 Kurama Okubo

clear all;
set(0,'DefaultTextFontsize',14, ...
    'DefaultTextFontname','Arial', ...
    'DefaultTextFontWeight','normal', ...
    'DefaultTextFontname','Arial', ...
    'DefaultAxesFontsize',14, ... 
    'DefaultAxesFontname','Arial', ...
    'DefaultLineLineWidth', 1.5)

%% LDV
fodir = "/Volumes/Okuboetal2025_masterHDD/AEsensorCalibration/AE_FreqResp_Mat_Masterdata";
if ~isfolder(fodir); mkdir(fodir); end

figdir = "../figure/debug_convertdata";
if ~isfolder(figdir); mkdir(figdir); end

sectionlength = 2^14; % length of trace set in Sbench configuration.
Nevent_LDV = 10000;

%%
casename_list=["master_S25V", "master_S100V", "master_S150V", "master_S200V",...
    "termination_1Mohm_S200V", "testdata_S200V", "testdata_side_S200V", "granite_S200V"];

for ic =1:length(casename_list)
    % ic = 4;
    casename = casename_list(ic);
    fprintf("start processing %s", casename);
    % pname_LDV = '/Volumes/HDD_Okubo_20220224_RAID1/BIAX/AE_FreqResp_MasterRawData/LDV';
    pname_LDV = '/Volumes/Okuboetal2025_masterHDD/AEsensorCalibration/AE_FreqResp_MasterRawData/LDV';
    fname_LDV = sprintf('FreqResp_LDV_%s', casename);
    [datmat,dattrig,tmat,Ch_info]=SBENCHreader_TrigStack(sectionlength,Nevent_LDV,pname_LDV,fname_LDV,'fs_read', 1e7);
    
    tr_LDV = mean(datmat(1:10000, :), 1);
    
    fig = figure(1);
    fig.Units = 'point';
    fig.Position = [0 500 800 500];
    clf(fig,'reset'); cla(fig,'reset'); hold on;
    box on;
    plot(tmat, tr_LDV, "k-");
    
    title(strrep(casename, "_", "\_"));
    figname = fullfile(figdir, casename+"_stackplot.png");
    exportgraphics(gcf, figname, "Resolution", 150);
    
    
    save(sprintf("%s/LDV_%s.mat", fodir, casename), "datmat", "tmat",'-v7');

end


%% save for linear prog

Nevent_linprog = 18000; % 2000 x 9 cases
input_ohm = 50; % input inpedence

V = [10, 25, 50, 75, 100, 125, 150, 175, 200]; %[V]

casename_list=["linearprog"];

for ic =1:length(casename_list)
    casename = casename_list(ic);
    fprintf("start processing %s", casename);
    % pname_LDV = '/Volumes/HDD_Okubo_20220224_RAID1/BIAX/AE_FreqResp_MasterRawData/LDV';
    pname_LDV = '/Volumes/Okuboetal2025_masterHDD/AEsensorCalibration/AE_FreqResp_MasterRawData/LDV';
    fname_LDV = sprintf('FreqResp_LDV_%s', casename);
    [datmat,dattrig,tmat,Ch_info]=SBENCHreader_TrigStack(sectionlength,Nevent_linprog,pname_LDV,fname_LDV,'fs_read', 1e7);
    
    fig = figure(1);
    fig.Units = 'point';
    fig.Position = [0 500 800 500];
    clf(fig,'reset'); cla(fig,'reset'); hold on;
    box on;

    for j = length(V):-1:1

        sind = 2000*(j-1)+1;
        eind = 2000*j; 

        tr_LDV = mean(datmat(sind:eind, :), 1);

        plot(tmat, tr_LDV, "-");
    end

    title(strrep(casename, "_", "\_"));
    figname = fullfile(figdir, casename+"_stackplot.png");
    exportgraphics(gcf, figname, "Resolution", 150);
    
    
    save(sprintf("%s/LDV_%s.mat", fodir, casename), "datmat", "tmat", 'V', '-v7.3');

end


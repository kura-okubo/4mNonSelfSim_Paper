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
Nevent_AE = 2000;

%%
casename_list=["master_S100V", "testdata_S100V", "testdata_side_S100V"];


for ic =1:length(casename_list)
    % ic = 4;
    casename = casename_list(ic);
    fprintf("start processing %s", casename);
    % pname_AE = '/Volumes/HDD_Okubo_20220224_RAID1/BIAX/AE_FreqResp_MasterRawData/AE';
    pname_AE = '/Volumes/Okuboetal2025_masterHDD/AEsensorCalibration/AE_FreqResp_MasterRawData/AE';
    fname_AE = sprintf('FreqResp_AE_%s', casename);
    [datmat,dattrig,tmat,Ch_info]=SBENCHreader_TrigStack(sectionlength,Nevent_AE,pname_AE,fname_AE,'fs_read', 1e7);
    
    tr_LDV = mean(datmat(1:Nevent_AE, :), 1);
    
    fig = figure(1);
    fig.Units = 'point';
    fig.Position = [0 500 800 500];
    clf(fig,'reset'); cla(fig,'reset'); hold on;
    box on;
    plot(tmat, tr_LDV, "k-");
    title(strrep(casename, "_", "\_"));
    figname = fullfile(figdir, casename+"_stackplot.png");
    exportgraphics(gcf, figname, "Resolution", 150);
    
    
    save(sprintf("%s/AE_%s.mat", fodir, casename), "datmat", "tmat",'-v7');

end

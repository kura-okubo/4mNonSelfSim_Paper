% Script to stack the signal of LDV from transducer on metal block.
% read data from the triggered traces
% 2022/07/01 Kurama Okubo
% 2022/08/10 Update file name for the validation
% 2022/09/05 Update for master data

clear all;
set(0,'DefaultTextFontsize',14, ...
    'DefaultTextFontname','Arial', ...
    'DefaultTextFontWeight','normal', ...
    'DefaultTextFontname','Arial', ...
    'DefaultAxesFontsize',14, ... 
    'DefaultAxesFontname','Arial', ...
    'DefaultLineLineWidth', 1.5)
%% Parameters
sectionlength = 2^14; % length of trace set in Sbench configuration.
Nevent_AE = 2000;
fs = 1e7;

if ~isfolder("../data"); mkdir("../data"); end

fodir = "/Volumes/Okuboetal2025_masterHDD/AEsensorCalibration/AE_FreqResp_Mat_Masterdata";
if ~isfolder(fodir); mkdir(fodir); end

figdir = "../figure/debug_convertdata";
if ~isfolder(figdir); mkdir(figdir); end

casename_list=["master_S100V", "testdata_S100V", "testdata_side_S100V"];
casename_list_out=["frontcenter_S100V", "fronttop_S100V", "sidecenter_S100V"];

for ic =1:length(casename_list)
% ic = 4;
casename = casename_list(ic);
finame = sprintf("%s/AE_%s.mat", fodir, casename);

A = load(finame);
tr_stacked = mean(A.datmat);

S.t = A.tmat;
S.AE =  mean(A.datmat)';
Stable = struct2table(S);
foname = sprintf("../data/AE_%s.csv", casename_list_out(ic));
writetable(Stable,foname);

end
